MODULE SCFs
  USE Parse
  USE InOut

  USE LinAlg 
  USE GlobalObjects

  USE SCFKeys
  USE Overlay
  USE SCFKeys
  USE PunchHDF
  USE Numerics
  USE OptionKeys
  USE Functionals
  USE ControlStructures
  USE NEB
  IMPLICIT NONE 
  INTEGER HDFFileID,H5GroupID
CONTAINS
  !===============================================================================
  !
  !===============================================================================
  SUBROUTINE SinglePoints(C)
    TYPE(Controls) :: C
    INTEGER        :: iBAS,iGEO
    !----------------------------------------------------------------------------!
    ! Just one geometry
    iGEO=1
    ! Init previous state
    C%Stat%Previous%I=(/0,1,1/)
    ! Init groups
    CALL InitClones(C%Nams,C%MPIs,C%Sets,C%Geos)
    ! Loop over basis sets 
    DO iBAS=1,C%Sets%NBSets
       ! Archive 
       CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)
       CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
       ! Converge an SCF
       CALL SCF(iBAS,iGEO,C)
    ENDDO
  END SUBROUTINE SinglePoints
  !===============================================================================
  !
  !===============================================================================
  SUBROUTINE SCF(cBAS,cGEO,C)
    TYPE(Controls)    :: C
    TYPE(DBL_RNK2)    :: ETot,DMax,DIIS
    INTEGER,PARAMETER :: MaxSCFs=32
    INTEGER           :: cBAS,cGEO,iSCF
    !----------------------------------------------------------------------------!
    ! Compute one-electron matrices
    CALL OneEMats(cBAS,cGEO,C%Nams,C%Stat,C%Opts,C%MPIs)
    ! Allocate space for convergence statistics
    CALL New(ETot,(/MaxSCFs,C%Geos%Clones/),(/0,1/))
    CALL New(DMax,(/MaxSCFs,C%Geos%Clones/),(/0,1/))
    CALL New(DIIS,(/MaxSCFs,C%Geos%Clones/),(/0,1/))
    DO iSCF=0,MaxSCFs
       ! Do an SCF cycle
       IF(SCFCycle(iSCF,cBAS,cGEO, &
            C%Nams,C%Stat,C%Opts,C%Geos,C%MPIs,ETot,DMax,DIIS))THEN
          ! Free memory
          CALL Delete(ETot)
          CALL Delete(DMax)
          CALL Delete(DIIS)
          RETURN
       ENDIF
    ENDDO
    CALL MondoHalt(DRIV_ERROR,'Failed to converge SCF in ' &
         //TRIM(IntToChar(MaxSCFs))//' SCF iterations.')
  END SUBROUTINE SCF
  !===============================================================================
  !
  !===============================================================================
  FUNCTION SCFCycle(cSCF,cBAS,cGEO,N,S,O,G,M,ETot,DMax,DIIS) 
    TYPE(FileNames)  :: N
    TYPE(State)      :: S
    TYPE(Options)    :: O
    TYPE(Geometries) :: G
    TYPE(Parallel)   :: M
    TYPE(DBL_RNK2)   :: ETot,DMax,DIIS
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
    CALL StateArchive(N,S)
    SCFCycle=ConvergedQ(cSCF,cBAS,N,S,O,G,ETot,DMax,DIIS)
    S%Previous%I=S%Current%I
  END FUNCTION SCFCycle
  !===============================================================================
  ! COMPUTE AN ENERGY GRADIENT
  !===============================================================================
  SUBROUTINE Force(cBAS,cGEO,N,O,S,G,B,M)
    TYPE(FileNames)  :: N
    TYPE(Options)    :: O
    TYPE(State)      :: S
    TYPE(Geometries) :: G
    TYPE(Parallel)   :: M

    TYPE(BasisSets)  :: B

    TYPE(DBL_VECT)   :: GradE
    INTEGER          :: cBAS,cGEO,K,J,iATS,iCLONE
    CHARACTER(LEN=3) :: chGEO
    !----------------------------------------------------------------------------!
    ! Initialize the force vector in HDF, clone by clone
    chGEO=IntToChar(cGEO)
    CALL New(GradE,G%Clone(1)%NAtms*3)
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       GradE%D=BIG_DBL
       ! Put the initialized forces back ...
       CALL Put(GradE,'GradE',Tag_O=chGEO)
       ! ... and close the group
       CALL CloseHDFGroup(HDF_CurrentID)
       G%Clone(iCLONE)%GradRMS=SQRT(G%Clone(iCLONE)%GradRMS)/DBLE(3*G%Clone(iCLONE)%NAtms)
    ENDDO
    CALL CloseHDF(HDFFileID)
    ! Now evaluate the forces
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
    IF(HasHF(O%Models(cBas)))THEN
       CALL NXForce(cBAS,cGEO,N,G,B,S,M)
       ! CALL Invoke('XForce',N,S,M)
    ENDIF
    ! DFT exchange corrleation term
    IF(HasDFT(O%Models(cBas))) &
         CALL Invoke('XCForce',N,S,M)
    ! Load forces
    chGEO=IntToChar(cGEO)
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       CALL Get(GradE,'GradE',Tag_O=chGEO)
       K=0
       G%Clone(iCLONE)%Vects%D=Zero
       IF(O%Coordinates==GRAD_INTS_OPT) THEN
         DO iATS=1,G%Clone(iCLONE)%NAtms
           DO J=1,3;K=K+1
             G%Clone(iCLONE)%Vects%D(J,iATS)=GradE%D(K)
           ENDDO
         ENDDO
       ELSE
         DO iATS=1,G%Clone(iCLONE)%NAtms
            IF(G%Clone(iCLONE)%CConstrain%I(iATS)==0)THEN
               DO J=1,3;K=K+1
                  G%Clone(iCLONE)%Vects%D(J,iATS)=GradE%D(K)
               ENDDO
            ELSE
               K=K+3
            ENDIF
         ENDDO
       ENDIF
       ! Close the group
       CALL CloseHDFGroup(HDF_CurrentID)
       G%Clone(iCLONE)%GradRMS=SQRT(G%Clone(iCLONE)%GradRMS)/DBLE(3*G%Clone(iCLONE)%NAtms)
    ENDDO
    CALL CloseHDF(HDFFileID)

    ! NEB force projections
    IF(O%Grad==GRAD_TS_SEARCH_NEB) THEN
       CALL NEBForce(G,O)
    ENDIF
    ! Zero forces on contrained atoms and compute stats with projected forces
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       K=0
       GradE%D=Zero
       G%Clone(iCLONE)%GradMax=Zero
       G%Clone(iCLONE)%GradRMS=Zero
       DO iATS=1,G%Clone(iCLONE)%NAtms
          IF(G%Clone(iCLONE)%CConstrain%I(iATS)==0)THEN
             DO J=1,3;K=K+1
                GradE%D(K)=G%Clone(iCLONE)%Vects%D(J,iATS)
                G%Clone(iCLONE)%GradRMS=G%Clone(iCLONE)%GradRMS+GradE%D(K)**2
                G%Clone(iCLONE)%GradMax=MAX(G%Clone(iCLONE)%GradMax,ABS(GradE%D(K)))
             ENDDO
          ELSE
             K=K+3
          ENDIF
       ENDDO
       ! Put the zeroed forces back ...
       CALL Put(GradE,'GradE',Tag_O=chGEO)
       ! ... and close the group
       CALL CloseHDFGroup(HDF_CurrentID)
       G%Clone(iCLONE)%GradRMS=SQRT(G%Clone(iCLONE)%GradRMS)/DBLE(3*G%Clone(iCLONE)%NAtms)
    ENDDO
    ! Now close the HDF file ..
    CALL CloseHDF(HDFFileID)
    ! .. and clean up 
    CALL Delete(GradE)
  END SUBROUTINE Force
  !===============================================================================
  ! Numerically compute gradients of the exact HF exchange
  !===============================================================================
  SUBROUTINE NXForce(cBAS,cGEO,N,G,B,S,M)
    TYPE(FileNames)  :: N
    TYPE(Options)    :: O
    TYPE(State)      :: S
    TYPE(Geometries) :: G
    TYPE(BasisSets)  :: B
    TYPE(Parallel)   :: M
    TYPE(DBL_VECT)   :: GradE
    INTEGER          :: cBAS,cGEO,J,iATS,iCLONE
    CHARACTER(LEN=3) :: chGEO,chBAS,chSCF
    TYPE(BCSR),DIMENSION(G%Clones)   :: P
    TYPE(CRDS),DIMENSION(G%Clones)   :: GTmp
    TYPE(BCSR)                       :: K
    REAL(DOUBLE),DIMENSION(G%Clones,G%Clone(1)%NAtms*3) :: FX
    REAL(DOUBLE),DIMENSION(G%Clones,2)        :: EX
    INTEGER                          :: AtA,IX,II,IA,A1,A2,IS
    REAL(DOUBLE),PARAMETER           :: DDelta = 1.D-3
    CHARACTER(LEN=DCL)               :: TrixName
    !------------------------------------------------------------------------------
    chGEO=IntToChar(cGEO)
    chBAS=IntToChar(cBAS)
    chSCF=IntToChar(S%Current%I(1)+1)
    CALL New(BSiz,G%Clone(1)%NAtms)
    CALL New(OffS,G%Clone(1)%NAtms)
    CALL New(GradE,G%Clone(1)%NAtms*3)
    DO iCLONE=1,G%Clones
       ! Load globals 
       NAToms=G%Clone(1)%NAtms
       MaxAtms=B%MxAts(cBAS)
       MaxBlks=B%MxN0s(cBAS)
       MaxNon0=B%MxBlk(cBAS)
       NBasF=B%BSets(iCLONE,cBAS)%NBasF
       BSiz%I=B%BSiz(iCLONE,cBAS)%I
       OffS%I=B%OffS(iCLONE,cBAS)%I
       MaxBlkSize=0
       DO II=1,G%Clone(1)%NAtms; MaxBlkSize=MAX(MaxBlkSize,BSiz%I(II)); ENDDO
       ! Set temporary geometries
       GTmp(iCLONE)%NAtms=G%Clone(iCLONE)%NAtms
       CALL New_CRDS(GTmp(iCLONE))
       GTmp(iCLONE)%AbCarts%D=G%Clone(iCLONE)%AbCarts%D
       ! Get the density matrix for this clone
       TrixName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//'_Geom#'//TRIM(chGEO)//'_Base#'//TRIM(chBAS)//'_Cycl#'//TRIM(chSCF) &
            //'_Clone#'//TRIM(IntToChar(iCLONE))//'.D'
       CALL Get(P(iCLONE),TrixName)
    ENDDO
    chSCF=IntToChar(S%Current%I(1))
    FX=Zero
    DO AtA=1,G%Clone(1)%NAtms
       DO IX=1,3
          DO II=1,2
             DO iCLONE=1,G%Clones
                IF(II==1) THEN
                   G%Clone(iCLONE)%AbCarts%D(IX,AtA)=GTmp(iCLONE)%AbCarts%D(IX,AtA)+DDelta
                ELSEIF(II==2) THEN
                   G%Clone(iCLONE)%AbCarts%D(IX,AtA)=GTmp(iCLONE)%AbCarts%D(IX,AtA)-DDelta
                ENDIF
             ENDDO
             CALL GeomArchive(cBAS,cGEO,N,B,G)    
             CALL Invoke('ONX',N,S,M)
             DO iCLONE=1,G%Clones
                ! Load globals 
                NAToms=G%Clone(1)%NAtms
                MaxAtms=B%MxAts(cBAS)
                MaxBlks=B%MxN0s(cBAS)
                MaxNon0=B%MxBlk(cBAS)
                NBasF=B%BSets(iCLONE,cBAS)%NBasF
                BSiz%I=B%BSiz(iCLONE,cBAS)%I
                OffS%I=B%OffS(iCLONE,cBAS)%I
                ! Get the exact HF exchange matrix from disk
                TrixName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//'_Geom#'//TRIM(chGEO)//'_Base#'//TRIM(chBAS)//'_Cycl#'//TRIM(chSCF) &
                     //'_Clone#'//TRIM(IntToChar(iCLONE))//'.K'
                CALL Get(K,TrixName)
                EX(iCLONE,II)=Trace(P(iCLONE),K)
             ENDDO
          ENDDO
          IA=3*(AtA-1)+IX
          DO iCLONE=1,G%Clones
             FX(iCLONE,IA)=(EX(iCLONE,1)-EX(iCLONE,2))/(Two*DDelta)
          ENDDO
       ENDDO
    ENDDO
   ! Add in the forces to the global gradient and put back to HDF
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       CALL Get(GradE,'GradE',Tag_O=chGEO)
       GradE%D=GradE%D+FX(iCLONE,:)
       CALL Put(GradE,'GradE',Tag_O=chGEO)
       ! Close the group
       CALL CloseHDFGroup(HDF_CurrentID)
    ENDDO
    CALL CloseHDF(HDFFileID)
    DO iCLONE=1,G%Clones
       G%Clone(iCLONE)%AbCarts%D=GTmp(iCLONE)%AbCarts%D
       CALL Delete(GTmp(iCLONE))
       CALL Delete(P(iCLONE))
    ENDDO
    CALL Delete(GradE)
    CALL Delete(BSiz)
    CALL Delete(OffS)
    CALL Delete(K)
  END SUBROUTINE NXForce
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
       S%Previous%I=O%RestartState%I
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
       CALL MondoHalt(99,'Unknown method key = '//TRIM(IntToChar(O%Methods(cBAS))))
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
  FUNCTION ConvergedQ(cSCF,cBAS,N,S,O,G,ETot,DMax,DIIS)
    TYPE(FileNames)             :: N
    TYPE(State)                 :: S
    TYPE(Options)               :: O
    TYPE(Geometries)            :: G
    TYPE(Parallel)              :: M
    TYPE(DBL_RNK2)              :: ETot,DMax,DIIS
    INTEGER                     :: cSCF,cBAS,iGEO,iCLONE
    REAL(DOUBLE)                :: DIISA,DIISB,DDIIS,DIISQ,       &
                                   DETOT,ETOTA,ETOTB,ETOTQ,ETEST, &
                                   DDMAX,DMAXA,DMAXB,DMAXQ,DTEST
    LOGICAL,DIMENSION(G%Clones) :: Converged
    LOGICAL                     :: ConvergedQ
    CHARACTER(LEN=3)            :: chGEO
    !----------------------------------------------------------------------------!
    ! Convergence thresholds
    ETest=ETol(O%AccuracyLevels(cBAS))
    DTest=DTol(O%AccuracyLevels(cBAS))
    IF(cSCF==0)THEN
       ConvergedQ=.FALSE.
       RETURN
    ENDIF
    ! Accumulate current statistics
    chGEO=IntToChar(iGEO)
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       ! Gather convergence parameters
#ifdef PARALLEL_CLONES
       CALL Get(Etot%D(cSCF,iCLONE),'Etot')
       CALL Get(DMax%D(cSCF,iCLONE),'DMax')
       CALL Get(DIIS%D(cSCF,iCLONE),'DIISErr')
#else
       CALL Get(Etot%D(cSCF,iCLONE),'Etot',StatsToChar(S%Current%I))
       CALL Get(DMax%D(cSCF,iCLONE),'DMax',StatsToChar(S%Current%I))
       CALL Get(DIIS%D(cSCF,iCLONE),'DIISErr',StatsToChar(S%Current%I))
#endif
       CALL CloseHDFGroup(HDF_CurrentID)
       ! Load current energies
       G%Clone(iCLONE)%ETotal=ETot%D(cSCF,iCLONE)

       Converged(iCLONE)=.FALSE.
       IF(cSCF>1)THEN          
          ETotA=ETot%D(cSCF-1,iCLONE)
          ETotB=ETot%D(cSCF  ,iCLONE)
          DMaxA=DMax%D(cSCF-1,iCLONE)
          DMaxB=DMax%D(cSCF  ,iCLONE)
          DIISA=DIIS%D(cSCF-1,iCLONE)
          DIISB=DIIS%D(cSCF  ,iCLONE)
          ! Absolute numbers
          dETot=ABS(ETotA-ETotB)
          dDMax=ABS(DMaxA-DMaxB)
          dDIIS=ABS(DIISA-DIISB)
          ! Relative numbers (Quotients)
          ETotQ=dETot/ABS(ETotB)
          DMaxQ=dDMax/ABS(DMaxB+1.D-50)
          DIISQ=dDIIS/ABS(DIISB+1.D-50)
          !CALL OpenASCII(OutFile,Out)
          !WRITE(Out,*)'ETest = ',ETest
          !WRITE(Out,*)'DTest = ',DTest
          !WRITE(Out,*)'ETotQ = ',ETotQ
          !WRITE(Out,*)'ETotA = ',ETotA
          !WRITE(Out,*)'ETotB = ',ETotB
          !WRITE(Out,*)'DIISQ = ',DIISQ
          !WRITE(Out,*)'DMaxQ = ',DMaxQ
          !WRITE(Out,*)'DIISA = ',DIISA
          !WRITE(Out,*)'DIISB = ',DIISB
          !WRITE(Out,*)'DMaxA = ',DMaxA
          !WRITE(Out,*)'DMaxB = ',DMaxB
          !CLOSE(Out)
          ! Convergence tests
          IF(((DMaxB<dTest.AND.ETotQ<ETest).OR.DMaxB<5D-1*dTest).AND.ETotB<ETotA)THEN
             Converged(iCLONE)=.TRUE.
             Mssg='Normal SCF convergence.a'
          ENDIF
          ! Accept convergence from wrong side if DM thresholds are tightend.
          IF(DMaxB<dTest*75D-2.AND.ETotQ<ETest*3D-1)THEN
             !        IF(DMaxB<dTest*1D-1.AND.ETotQ<ETest*1D-1)THEN
             Converged(iCLONE)=.TRUE.
             Mssg='Normal SCF convergence.b'
          ENDIF
          ! Look for stall out if we have at least one digit in the DM
          IF(DMaxB<1.D-1)THEN
             ! Look for non-decreasing errors due to incomplete numerics
             IF(DIISQ<1.D-1.AND.DMaxQ<1.D-1.AND.cSCF>2)THEN
                IF(DIISB>DIISA.AND.DMaxB>DMaxA)THEN
                   Mssg='SCF hit DIIS & DMax increase.'
                   Converged(iCLONE)=.TRUE.
                ENDIF
             ELSEIF(DIISQ<1.D-2.AND.DMaxQ<1.D-2.AND.cSCF>2)THEN
                IF(DIISB>DIISA)THEN
                   Mssg='SCF hit DIIS increase. a'
                   Converged(iCLONE)=.TRUE.
                ELSEIF(DMaxQ<1D-1.AND.DMaxB>DMaxA)THEN
                   Mssg='SCF hit DIIS increase. b'
                   Converged(iCLONE)=.TRUE.
                ENDIF
             ELSEIF((DIISQ<1D-3.OR.DMaxQ<1D-3).AND.cSCF>2)THEN
                Mssg='SCF convergence due to DIIS stagnation.'
                Converged(iCLONE)=.TRUE.
             ENDIF
          ENDIF
       ENDIF
    ENDDO
    CALL CloseHDF(HDFFileID)
    IF(cSCF>1)ConvergedQ=.TRUE.
    DO iCLONE=1,G%Clones
       ConvergedQ=ConvergedQ.AND.Converged(iCLONE)
    ENDDO
    ! Convergence announcement
    IF(ConvergedQ)THEN!.AND.PrintFlags%Key>DEBUG_NONE)THEN
       CALL OpenASCII(OutFile,Out)
       WRITE(Out,*)TRIM(Mssg)
       WRITE(*,*)TRIM(Mssg)
!       WRITE(Out,*)'Normal SCF convergence.'
       CLOSE(Out)
    ENDIF
  END FUNCTION ConvergedQ



END MODULE SCFs




