MODULE SCFs
  USE Parse
  USE InOut
  USE SCFKeys
  USE Overlay
  USE SCFKeys
  USE PunchHDF
  USE Thresholds
  USE OptionKeys
  USE Functionals
  USE ControlStructures
  IMPLICIT NONE 
  INTEGER HDFFileID,H5GroupID
CONTAINS
  !===============================================================================

  !===============================================================================
  SUBROUTINE SinglePoints(C)
    TYPE(Controls) :: C
    INTEGER        :: iBAS,iGEO
    !----------------------------------------------------------------------------!
    ! Just one geometry
    iGEO=1
    ! Init previous state
    C%Stat%Previous%I=(/0,1,1/)
    ! Archive MPIs and Geos 
    CALL MPIsArchive(C%Nams,C%Geos,C%MPIs)
    CALL InitClones(C%Nams,C%Geos)
    CALL GeomArchive(iGEO,C%Nams,C%Geos)    
    ! Loop over basis sets 
    DO iBAS=1,C%Sets%NBSets
       CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
       CALL SCF(iBAS,iGEO,C)
    ENDDO
  END SUBROUTINE SinglePoints
  !===============================================================================

  !===============================================================================
  SUBROUTINE SCF(cBAS,cGEO,C)
    TYPE(Controls)    :: C
    INTEGER,PARAMETER :: MaxSCFs=32
    INTEGER           :: cBAS,cGEO,iSCF
    !----------------------------------------------------------------------------!
    ! Compute one-electron matrices
    CALL OneEMats(cBAS,cGEO,C%Nams,C%Stat,C%Opts,C%MPIs)
    DO iSCF=0,MaxSCFs
       ! Do an SCF cycle
       IF(SCFCycle(iSCF,cBAS,cGEO,C%Nams,C%Stat,C%Opts,C%Geos,C%MPIs))RETURN
    ENDDO
    CALL MondoHalt(DRIV_ERROR,'Failed to converge SCF in ' &
         //TRIM(IntToChar(MaxSCFs))//' SCF iterations.')
  END SUBROUTINE SCF
  !===============================================================================
  !
  !===============================================================================
  FUNCTION SCFCycle(cSCF,cBAS,cGEO,N,S,O,G,M) 
    TYPE(FileNames)  :: N
    TYPE(State)      :: S
    TYPE(Options)    :: O
    TYPE(Geometries) :: G
    TYPE(Parallel)   :: M
    INTEGER          :: cSCF,cBAS,cGEO
    LOGICAL          :: SCFCycle
    !----------------------------------------------------------------------------!
    CALL DensityLogic(cSCF,cBAS,cGEO,S,O)
    CALL DensityBuild(N,S,M)
    IF(cSCF/=0) &
         S%Action=SCF_FOCK_BUILD
    CALL FockBuild(cSCF,cBAS,N,S,O,M)
    IF(cSCF/=0) &
         S%Action=SCF_SOLVE_SCF
    CALL SolveSCF(cBAS,N,S,O,M)
    IF(ConvergedQ(cSCF,cBAS,N,S,O,G))THEN
       SCFCycle=.TRUE.
    ELSE
       SCFCycle=.FALSE.
    ENDIF
    S%Previous%I=S%Current%I
!    CALL StateArchive(N,S)
  END FUNCTION SCFCycle
  !===============================================================================
  ! COMPUTE AN ENERGY GRADIENT
  !===============================================================================
  SUBROUTINE Force(cBAS,cGEO,N,O,S,G,M)
    TYPE(FileNames)  :: N
    TYPE(Options)    :: O
    TYPE(State)      :: S
    TYPE(Geometries) :: G
    TYPE(Parallel)   :: M
    TYPE(DBL_VECT)   :: GradE
    INTEGER          :: cBAS,cGEO,K,J,iATS,iCLONE
    CHARACTER(LEN=3) :: chGEO
    !----------------------------------------------------------------------------!
    S%Action='ForceEvaluation'
    ! The non-orthogonal response    
    CALL Invoke('SForce',N,S,M)
    ! Kinetic energy piece
    CALL Invoke('TForce',N,S,M)
    ! Build density with last DM
    CALL Invoke('MakeRho',N,S,M)
#ifdef PARALLEL
    CALL Invoke('ParaMakeRho',N,S,M)
#endif
    ! Coulomb part
    CALL Invoke('JForce',N,S,M)
    ! Exact Hartree-Fock exchange component
    IF(HasHF(O%Models(cBas)))  &
         CALL Invoke('XForce',N,S,M)
    ! DFT exchange corrleation term
    IF(HasDFT(O%Models(cBas))) &
         CALL Invoke('XCForce',N,S,M)
    ! Done, now retrieve forces 
    chGEO=IntToChar(cGEO)
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       CALL New(GradE,G%Clone(iCLONE)%NAtms*3)
       CALL Get(GradE,'GradE',Tag_O=chGEO)
       K=0
       G%Clone(iCLONE)%GradMax=Zero
       G%Clone(iCLONE)%GradRMS=Zero
       DO iATS=1,G%Clone(iCLONE)%NAtms
          IF(G%Clone(iCLONE)%CConstrain%I(iATS)==0)THEN
             DO J=1,3
                K=K+1
                G%Clone(iCLONE)%Vects%D(J,iATS)=GradE%D(K)
                G%Clone(iCLONE)%GradRMS=G%Clone(iCLONE)%GradRMS+GradE%D(K)**2
                G%Clone(iCLONE)%GradMax=MAX(G%Clone(iCLONE)%GradMax,ABS(GradE%D(K)))
             ENDDO
          ELSE
             ! Apply Cartesian constraints 
             DO J=1,3
                K=K+1
                GradE%D(K)=Zero
             ENDDO
          ENDIF
       ENDDO
       ! Put the zeroed forces back ...
       CALL Put(GradE,'GradE',Tag_O=chGEO)
       ! ... and close the group
       CALL CloseHDFGroup(HDF_CurrentID)
       G%Clone(iCLONE)%GradRMS=SQRT(G%Clone(iCLONE)%GradRMS)/DBLE(3*G%Clone(iCLONE)%NAtms)
       CALL Delete(GradE)
    ENDDO
    CALL CloseHDF(HDFFileID)
  END SUBROUTINE Force
  !===============================================================================
  ! BUILD A HGTF DENSITY BY HOOK OR BY CROOK
  !===============================================================================
  SUBROUTINE DensityBuild(N,S,M)
    TYPE(FileNames):: N
    TYPE(State)    :: S
    TYPE(Parallel) :: M    
    !----------------------------------------------------------------------------!
    IF(S%Action/=SCF_BASISSETSWITCH.AND. &
       S%Action/=SCF_DENSITY_NORMAL)     &
       CALL Invoke('P2Use',N,S,M)
    CALL Invoke('MakeRho',N,S,M)
#ifdef PARALLEL
    CALL Invoke('ParaMakeRho',N,S,M) 
#endif
  END SUBROUTINE DensityBuild
  !===============================================================================

  !===============================================================================
  SUBROUTINE DensityLogic(cSCF,cBAS,cGEO,S,O)
    TYPE(State)    :: S
    TYPE(Options)  :: O
    INTEGER        :: cSCF,cBAS,cGEO,pBAS
    !----------------------------------------------------------------------------!
    S%SubAction=""
    pBAS=S%Previous%I(2)
    S%Current%I=(/cSCF,cBAS,cGEO/)
    IF(O%Guess==GUESS_EQ_SUPR)THEN 
       O%Guess=0
       S%Previous%I=S%Current%I
       S%Action=SCF_SUPERPOSITION
    ELSEIF(O%Guess==GUESS_EQ_CORE)THEN
       O%Guess=0
       S%Previous%I=S%Current%I
       S%Action=SCF_GUESSEQCORE
    ELSEIF(O%Guess==GUESS_EQ_RESTART)THEN
       O%Guess=0
       S%Previous%I=S%Current%I
       S%Action=SCF_RESTART
    ELSEIF(cSCF==0.AND.cBAS==pBAS.AND.cGEO/=1)THEN
       S%Action=SCF_PROJECTION
    ELSEIF(cSCF==0.AND.cBAS/=pBAS)THEN
       S%Action=SCF_BASISSETSWITCH
       S%Previous%I(1)=S%Previous%I(1)+1
    ELSE
       S%Action=SCF_DENSITY_NORMAL
    ENDIF
  END SUBROUTINE DensityLogic
  !===============================================================================
  ! BUILD A FOCK MATRIX, POSSIBLY WITH DIIS EXTRAPOLATION
  !===============================================================================
  SUBROUTINE FockBuild(cSCF,cBAS,N,S,O,M)
    TYPE(FileNames):: N
    TYPE(State)    :: S
    TYPE(Options)  :: O
    TYPE(Parallel) :: M    
    INTEGER        :: cSCF,cBAS,Modl
    LOGICAL        :: DoDIIs
    !----------------------------------------------------------------------------!
    DoDIIS=cSCF>0
    Modl=O%Models(cBAS)
    CALL Invoke('QCTC',N,S,M)
    IF(S%SubAction/=SCF_GUESSEQCORE)THEN
       IF(HasHF(Modl)) CALL Invoke('ONX',N,S,M)
       IF(HasDFT(Modl))CALL Invoke('HiCu',N,S,M)
    ENDIF
    CALL Invoke('FBuild',N,S,M)
    IF(DoDIIS)CALL Invoke('DIIS',N,S,M)
  END SUBROUTINE FockBuild
  !===============================================================================
  !   Solve the SCF equations 
  !===============================================================================
  SUBROUTINE SolveSCF(cBAS,N,S,O,M)
    TYPE(FileNames):: N
    TYPE(State)    :: S
    TYPE(Options)  :: O
    TYPE(Parallel) :: M    
    INTEGER        :: cBAS
    !----------------------------------------------------------------------------!
    IF(O%Methods(cBAS)==RH_R_SCF)THEN
!       CALL LogSCF(S,'Solving SCF equations with Roothann-Hall')
       CALL Invoke('RHeqs',N,S,M)
    ELSEIF(O%Methods(cBAS)==SDMM_R_SCF) THEN
!       CALL LogSCF(S,'Solving SCF equations with SDMM')
       CALL Invoke('SDMM',N,S,M)
    ELSEIF(O%Methods(cBAS)==PM_R_SCF) THEN
!       CALL LogSCF(Current,'Solving SCF equations with PM')
       CALL Invoke('PM',N,S,M)
    ELSEIF(O%Methods(cBAS)==SP2_R_SCF) THEN
!       CALL LogSCF(Current,'Solving SCF equations with SP2')
       CALL Invoke('SP2',N,S,M)
    ELSEIF(O%Methods(cBAS)==SP4_R_SCF) THEN
!      CALL LogSCF(Current,'Solving SCF equations with SP4')
       CALL Invoke('SP4',N,S,M)
    ELSEIF(O%Methods(cBAS)==TS4_R_SCF) THEN
!       CALL LogSCF(Current,'Solving SCF equations with TS4')
       CALL Invoke('TS4',N,S,M)
    ELSE
       CALL MondoHalt(99,'No Method Chosen')
    ENDIF
    CALL Invoke('SCFstats',N,S,M)
  END SUBROUTINE SolveSCF
  !===============================================================================

  !===============================================================================
  SUBROUTINE OneEMats(cBAS,cGEO,N,S,O,M)
    TYPE(FileNames):: N
    TYPE(State)    :: S
    TYPE(Options)  :: O
    TYPE(Parallel) :: M
    INTEGER        :: cBAS,cGEO,pBAS
    LOGICAL, SAVE  :: DoPFFT=.TRUE.
    !----------------------------------------------------------------------------!
    pBAS=S%Previous%I(2)    
    S%Current%I=(/0,cBAS,cGEO/)
    S%Action='OneElectronMatrices'
    S%SubAction=""
#ifdef PERIODIC
    IF(pBAS/=cBAS)DoPFFT=.TRUE.
    IF(DoPFFT)THEN
       ! CALL LogSCF(Current,'Peridic Far-Field Tensor',.TRUE.)
       CALL Invoke('MakePFFT',N,S,M)
       DoPFFT=.FALSE.
    ENDIF
#endif
    CALL Invoke('MakeS',N,S,M)
    IF(O%Methods(cBAS)==RH_R_SCF)THEN
       CALL Invoke('LowdinO',N,S,M)
    ELSE
       CALL Invoke('AInv',N,S,M)
    ENDIF
    CALL Invoke('MakeT',N,S,M)
  END SUBROUTINE OneEMats
  !===============================================================================
  !
  !===============================================================================
  FUNCTION ConvergedQ(cSCF,cBAS,N,S,O,G)
    TYPE(FileNames)  :: N
    TYPE(State)      :: S
    TYPE(Options)    :: O
    TYPE(Geometries) :: G
    TYPE(Parallel)   :: M
    INTEGER          :: cSCF,cBAS,iGEO,iCLONE
    REAL(DOUBLE)     :: DIISA,DIISB,DDIIS,DIISQ,       &
                        DETOT,ETOTA,ETOTB,ETOTQ,ETEST, &
                        DDMAX,DMAXA,DMAXB,DMAXQ,DTEST
    LOGICAL          :: ConvergedQ
    CHARACTER(LEN=3) :: chGEO
    !----------------------------------------------------------------------------!
    ConvergedQ=.FALSE.
    IF(cSCF==0)RETURN
    ! Retrieve current statistics
    chGEO=IntToChar(iGEO)
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       ! Gather convergence parameters
       ! IN FUTURE CAN PROBABLY JUST GET AWAY WITH CURRENT AND PREVIOUS SCF CYCLE...
       CALL Get(EtotA,'Etot',StatsToChar(S%Previous%I))
       CALL Get(EtotB,'Etot',StatsToChar(S%Current%I))
       CALL Get(DMaxA,'DMax',StatsToChar(S%Previous%I))
       CALL Get(DMaxB,'DMax',StatsToChar(S%Current%I))
       IF(CSCF==0)THEN
          DIISA=1.D0
          DIISB=1.D0
       ELSEIF(CSCF==1)THEN
          CALL Get(DIISB,'diiserr',StatsToChar(S%Current%I))
          DIISA=DIISB
       ELSE
          CALL Get(DIISB,'diiserr',StatsToChar(S%Current%I))
          CALL Get(DIISA,'diiserr',StatsToChar(S%Previous%I))
       ENDIF
       CALL CloseHDFGroup(HDF_CurrentID)
       !  Absolute numbers
       dETot=ABS(ETotA-ETotB)
       dDMax=ABS(DMaxA-DMaxB)
       dDIIS=ABS(DIISA-DIISB)
       ! Quotients
       ETotQ=dETot/ABS(ETotB)
       DMaxQ=dDMax/ABS(DMaxB+1.D-50)
       DIISQ=dDIIS/ABS(DIISB+1.D-50)
       ! Load current energies
       G%Clone(iCLONE)%ETotal=ETotB
       ! Uncertainty in total energy
       !G%Clone(iCLONE)%DeltaE=ETotQ
    ENDDO
    CALL CloseHDF(HDFFileID)
    ! Convergence thresholds
    ETest=ETol(O%AccuracyLevels(cBAS))
    DTest=DTol(O%AccuracyLevels(cBAS))
    ! Convergence tests
    IF(((DMaxB<dTest.AND.ETotQ<ETest).OR.DMaxB<5D-1*dTest).AND.ETotB<ETotA)THEN
       Mssg='Normal SCF convergence.'
       ConvergedQ=.TRUE.
    ENDIF
    ! Accept convergence from wrong side if thresholds are tightend.
    IF(DMaxB<dTest*75D-2.AND.ETotQ<ETest*3D-1)THEN
       Mssg='Normal SCF convergence.'
       ConvergedQ=.TRUE.
    ENDIF
    ! Look for non-decreasing errors due to incomplete numerics
    IF(DIISQ<1.D-1.AND.DMaxQ<1.D-1.AND.CSCF>2)THEN
       IF(DIISB>DIISA.AND.DMaxB>DMaxA)THEN
          Mssg='SCF hit DIIS & DMax increase.'
          ConvergedQ=.TRUE.
       ENDIF
    ELSEIF(DIISQ<1.D-2.AND.DMaxQ<1.D-2.AND.CSCF>2)THEN
       IF(DIISB>DIISA)THEN
          Mssg='SCF hit DIIS increase.'
          ConvergedQ=.TRUE.
       ELSEIF(DMaxQ<1D-1.AND.DMaxB>DMaxA)THEN
          Mssg='SCF hit DMAX increase.'
          ConvergedQ=.TRUE.
       ENDIF
    ELSEIF((DIISQ<1D-3.OR.DMaxQ<1D-3).AND.CSCF>2)THEN
       Mssg='SCF convergence due to DIIS stagnation.'
       ConvergedQ=.TRUE.
    ENDIF
    ! Convergence announcement
    IF(ConvergedQ)THEN!.AND.PrintFlags%Key>DEBUG_NONE)THEN
       CALL OpenASCII(OutFile,Out)
       WRITE(Out,*)TRIM(Mssg)             
       CLOSE(Out)
    ENDIF
  END FUNCTION ConvergedQ
END MODULE SCFs




