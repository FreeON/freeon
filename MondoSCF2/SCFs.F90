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
  USE SetXYZ 
  IMPLICIT NONE 
  INTEGER HDFFileID,H5GroupID
  INTEGER,PARAMETER :: NOT_CONVERGE=345632
  INTEGER,PARAMETER :: SCF_STALLED =5345634
  INTEGER,PARAMETER :: DIIS_NOPATH =56778444
  INTEGER,PARAMETER :: DID_CONVERGE=293443454
CONTAINS
!===============================================================================
!
!===============================================================================
  SUBROUTINE SinglePoints(C)
    TYPE(Controls) :: C
    INTEGER        :: iBAS,iGEO,iBBegin
!----------------------------------------------------------------------------!
!   Loop over geometry
    DO iGEO = 1,1 !C%Geos%NGeom
!      Init previous state
       C%Stat%Previous%I=(/0,1,1/)
!      Init iBBegin
       iBBegin = 1
       IF(iGEO > 1) iBBegin = C%Sets%NBSets
!      Loop over basis sets 
       DO iBAS=iBBegin,C%Sets%NBSets
!         Archive 
          CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)
          CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
!         Converge an SCF
          CALL SCF(iBAS,iGEO,C)
       ENDDO
    ENDDO
  END SUBROUTINE SinglePoints
!===============================================================================
!
!===============================================================================
  SUBROUTINE SCF(cBAS,cGEO,C)
    TYPE(Controls)    :: C
    TYPE(DBL_RNK2)    :: ETot,DMax,DIIS
    INTEGER,PARAMETER :: MaxSCFs=256
    INTEGER           :: cBAS,cGEO,iSCF
    !----------------------------------------------------------------------------!
    CALL New(C%Stat%Action,1)
    ! Compute one-electron matrices
    CALL OneEMats(cBAS,cGEO,C%Nams,C%Sets,C%Stat,C%Opts,C%MPIs)
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
          CALL Delete(C%Stat%Action)
          RETURN
       ENDIF
    ENDDO
    CALL MondoHalt(DRIV_ERROR,'Failed to converge SCF in ' &
         //TRIM(IntToChar(MaxSCFs))//' SCF iterations.')
  END SUBROUTINE SCF
  !===============================================================================
  !
  !===============================================================================
  FUNCTION SCFCycle(cSCF,cBAS,cGEO,N,S,O,G,M,ETot,DMax,DIIS,CPSCF_O)
    TYPE(FileNames)  :: N
    TYPE(State)      :: S
    TYPE(Options)    :: O
    TYPE(Geometries) :: G
    TYPE(Parallel)   :: M
    TYPE(DBL_RNK2)   :: ETot,DMax,DIIS
    INTEGER          :: cSCF,cBAS,cGEO,iCLONE,Modl,WhatsUp
    LOGICAL,OPTIONAL :: CPSCF_O
    LOGICAL          :: DoDIIS,SCFCycle,DoCPSCF,DoODA,RebuildPostODA, &
                        DoMIX,RebuildPostMIX,Automatic1,Automatic2
    CHARACTER(LEN=128) :: Tmp
    !----------------------------------------------------------------------------!

    IF(cSCF==0)THEN
       DoODA=.TRUE.
       DoMIX=.FALSE.
       DoDIIS=.FALSE.
       Automatic1=.TRUE.
       Automatic2=.TRUE.
       ! Switching on/off rebuild doesn't work so well.  Lets be conservative for now
       RebuildPostODA=.TRUE.
    ENDIF
    ! Parse for strict ODA or DIIS invocation
    CALL OpenASCII(N%IFile,Inp)
    IF(cSCF>3.AND.OptKeyQ(Inp,CONVERGE_OPTION,CONVERGE_DODIIS))THEN
       DoODA=.FALSE.
       DoDIIS=.TRUE.
    ELSEIF(OptKeyQ(Inp,CONVERGE_OPTION,CONVERGE_DOODA))THEN
       DoODA=.TRUE.
       DoDIIS=.FALSE.
    ENDIF
    CLOSE(UNIT=Inp,STATUS='KEEP')
    ! Check semantics
    IF(DoDIIS.AND.(DoODA.OR.DoMIX))THEN
       CALL MondoHalt(DRIV_ERROR,'Logic failure 1 in SCFCycle')
    ENDIF
    ! Are we maybe solving CPSCF equations?
    IF(PRESENT(CPSCF_O))THEN
       DoCPSCF=CPSCF_O
    ELSE
       DoCPSCF=.FALSE.
    ENDIF
    ! The options...
    IF(DoCPSCF)THEN
       CALL DensityLogic(cSCF,cBAS,cGEO,S,O,CPSCF_O=.TRUE.)
       CALL DensityBuild(N,S,M)
       S%Action%C(1)=CPSCF_FOCK_BUILD
       CALL Invoke('QCTC',N,S,M)
       Modl=O%Models(cBAS)
       IF(HasHF(Modl)) CALL Invoke('ONX',N,S,M)
       IF(HasDFT(Modl))CALL Invoke('HiCu',N,S,M)
       CALL Invoke('FBuild',N,S,M)
       CALL Invoke('DDIIS',N,S,M)
    ELSEIF(DoDIIS)THEN
       CALL DensityLogic(cSCF,cBAS,cGEO,S,O)
       CALL DensityBuild(N,S,M)
       CALL FockBuild(cSCF,cBAS,N,S,O,M)
       IF(cSCF>0)CALL Invoke('DIIS',N,S,M)
       CALL SolveSCF(cBAS,N,S,O,M)
       CALL Invoke('SCFstats',N,S,M)
    ELSEIF(DoMIX)THEN
       RebuildPostMIX=.FALSE.
       IF((cSCF>1.AND.cGEO==1).OR.cSCF>0)THEN
          CALL DensityLogic(cSCF,cBAS,cGEO,S,O)
          CALL DensityBuild(N,S,M)
          CALL FockBuild(cSCF,cBAS,N,S,O,M)
          CALL SolveSCF(cBAS,N,S,O,M)
          CALL Invoke('SCFstats',N,S,M)
          S%Action%C(1)='Stanton MIX'
          CALL Invoke('MIX',N,S,M)
          CALL DensityLogic(cSCF,cBAS,cGEO,S,O)
          CALL DensityBuild(N,S,M)
          CALL FockBuild(cSCF,cBAS,N,S,O,M)
          CALL SolveSCF(cBAS,N,S,O,M)
       ELSE
          CALL DensityLogic(cSCF,cBAS,cGEO,S,O,CPSCF_O)
          CALL DensityBuild(N,S,M)
          CALL FockBuild(cSCF,cBAS,N,S,O,M)
          CALL SolveSCF(cBAS,N,S,O,M)
          CALL Invoke('SCFstats',N,S,M)
       ENDIF
    ELSEIF(DoODA)THEN
        IF(cSCF>1.OR.cSCF>0.AND.cGEO>1)THEN
          CALL DensityLogic(cSCF,cBAS,cGEO,S,O)
          CALL DensityBuild(N,S,M)
          CALL FockBuild(cSCF,cBAS,N,S,O,M)
          CALL SolveSCF(cBAS,N,S,O,M)
          Tmp=S%Action%C(1)
          S%Action%C(1)='Silent'
          CALL Invoke('SCFstats',N,S,M)
          S%Action%C(1)=Tmp
          CALL Invoke('ODA',N,S,M)
          IF(RebuildPostODA.AND.HasDFT(O%Models(cBAS)))THEN
             ! Rebuild non-linear KS matrix
             CALL DensityLogic(cSCF,cBAS,cGEO,S,O)
             CALL DensityBuild(N,S,M)
             CALL Invoke('HiCu',N,S,M)
             CALL Invoke('FBuild',N,S,M)
             !CALL FockBuild(cSCF,cBAS,N,S,O,M)
          ENDIF
          CALL SolveSCF(cBAS,N,S,O,M)
          CALL Invoke('SCFstats',N,S,M)
       ELSE
          CALL DensityLogic(cSCF,cBAS,cGEO,S,O,CPSCF_O)
          CALL DensityBuild(N,S,M)
          CALL FockBuild(cSCF,cBAS,N,S,O,M)
          CALL SolveSCF(cBAS,N,S,O,M)
          CALL Invoke('SCFstats',N,S,M)
       ENDIF
    ELSE
       CALL MondoHalt(DRIV_ERROR,'Logic failure 1 in SCFCycle')
    ENDIF
    !
    CALL StateArchive(N,S)
    WhatsUp=ConvergedQ(cSCF,cBAS,N,S,O,G,ETot,DMax,DIIS,DoDIIS,DoODA,RebuildPostODA,CPSCF_O)
    S%Previous%I=S%Current%I
    !
    IF(.NOT.DoCPSCF.AND.Automatic1.AND.DoDIIS.AND.WhatsUp==SCF_STALLED)THEN
       ! DIIS didnt work, dont try it again ...
       DoDIIS=.FALSE.
       DoODA=.TRUE.
       SCFCycle=.FALSE.
       Automatic1=.FALSE.
    ELSEIF(.NOT.DoCPSCF.AND.Automatic1.AND.DoODA.AND.WhatsUp==DIIS_NOPATH)THEN
       ! DIIS might be ok, try it just once
       DoDIIS=.TRUE.
       DoODA=.FALSE.
       SCFCycle=.FALSE.
    ELSEIF(.NOT.DoCPSCF.AND.Automatic2.AND.DoODA.AND.WhatsUp==SCF_STALLED)THEN
       ! Mixing the KS energy and matrix didnt work, must rebuild KS from scratch
       DoDIIS=.FALSE.
       DoODA=.TRUE.
       SCFCycle=.FALSE.
       Automatic2=.FALSE.
       RebuildPostODA=.TRUE.
    ELSEIF(WhatsUp==DID_CONVERGE)THEN
       SCFCycle=.TRUE.
    ELSE
       SCFCycle=.FALSE.
    ENDIF
  END FUNCTION SCFCycle
  !===============================================================================
  !
  !===============================================================================
  FUNCTION ConvergedQ(cSCF,cBAS,N,S,O,G,ETot,DMax,DIIS,DoDIIS,DoODA,RebuildPostODA,CPSCF_O)
    TYPE(FileNames)             :: N
    TYPE(State)                 :: S
    TYPE(Options)               :: O
    TYPE(Geometries)            :: G
    TYPE(Parallel)              :: M
    TYPE(DBL_RNK2)              :: ETot,DMax,DIIS
    LOGICAL,OPTIONAL            :: CPSCF_O
    LOGICAL                     :: CPSCF,DoDIIS,DoODA,RebuildPostODA
    LOGICAL                     :: ALogic,BLogic,CLogic,DLogic,ELogic,A2Logic, &
                                   GLogic,QLogic,ILogic,OLogic,RLogic,FLogic
    INTEGER                     :: cSCF,cBAS,iGEO,iCLONE
    REAL(DOUBLE)                :: DIISA,DIISB,DDIIS,DIISQ,       &
         DETOT,ETOTA,ETOTB,ETOTQ,ETEST, &
         DDMAX,DMAXA,DMAXB,DMAXQ,DTEST,ETOTO,ODAQ
    INTEGER,DIMENSION(G%Clones) :: Converged
    INTEGER                     :: ConvergedQ,iSCF
    CHARACTER(LEN=DCL)            :: chGEO
    !----------------------------------------------------------------------------!
    !
    IF(PRESENT(CPSCF_O)) THEN
       CPSCF=CPSCF_O
    ELSE
       CPSCF=.FALSE.
    ENDIF
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
       IF(CPSCF) THEN
          CALL Get(Etot%D(cSCF,iCLONE),'Prop'    )
          CALL Get(DMax%D(cSCF,iCLONE),'DPrimMax')
          CALL Get(DIIS%D(cSCF,iCLONE),'DDIISErr')
       ELSE
          CALL Get(Etot%D(cSCF,iCLONE),'Etot')
          CALL Get(DMax%D(cSCF,iCLONE),'DMax')
          CALL Get(DIIS%D(cSCF,iCLONE),'DIISErr' )
          IF(DoODA.AND.cSCF>1)THEN
             CALL Get(ETotO,'ODAEnergy')
          ELSE
             ETotO=1D10
          ENDIF
       ENDIF
       ! Load current energies
       G%Clone(iCLONE)%ETotal=ETot%D(cSCF,iCLONE)
       Converged(iCLONE)=NOT_CONVERGE
       IF(cSCF>1)THEN
          ETotA=ETot%D(cSCF-1,iCLONE)
          ETotB=ETot%D(cSCF  ,iCLONE)
          DMaxA=DMax%D(cSCF-1,iCLONE)
          DMaxB=DMax%D(cSCF  ,iCLONE)
          DIISA=DIIS%D(cSCF-1,iCLONE)
          DIISB=DIIS%D(cSCF  ,iCLONE)
          IF(cSCF>2)THEN
             ETotA=1D10
!             DMaxA=1D10
             DIISA=1D10
             DO iSCF=2,cSCF-1
                ETotA=MIN(ETotA,ETot%D(iSCF,iCLONE))
!                DMaxA=MIN(DMaxA,DMax%D(iSCF,iCLONE))
                DIISA=MIN(DIISA,DIIS%D(iSCF,iCLONE))
             ENDDO
          ENDIF
          dETot=ABS(ETotA-ETotB)
          dDMax=ABS(DMaxA-DMaxB)
          dDIIS=ABS(DIISA-DIISB)
          ! Relative numbers
          ETotQ=dETot/ABS(ETotB)
          DMaxQ=dDMax/ABS(DMaxB+1.D-50)
          DIISQ=dDIIS/ABS(DIISB+1.D-50)
          IF(DoODA)THEN
             ODAQ=ABS(ETotB-ETotO)/ABS(ETotB)
          ELSE
             ODAQ=Zero
          ENDIF
          CALL OpenASCII(OutFile,Out)
          WRITE(Out,*)'ODAQ = ',ODAQ
          WRITE(Out,*)'ETotQ = ',ETotQ
          WRITE(Out,*)'DIISQ = ',DIISQ
          WRITE(Out,*)'DMaxQ = ',DMaxQ
          WRITE(Out,*)'ETOTO = ',ETotO
          WRITE(Out,*)'ETOTA = ',ETOTA
          WRITE(Out,*)'ETOTB = ',ETOTB
          WRITE(Out,*)'DIISA = ',DIISA
          WRITE(Out,*)'DIISB = ',DIISB
          WRITE(Out,*)'DMaxA = ',DMaxA
          WRITE(Out,*)'DMaxB = ',DMaxB
          Converged(iCLONE)=NOT_CONVERGE
          ! Convergence from above +/- expected delta
          ! relative to historical minimum energy
          ALogic=ETotB*(One+ETest)<ETotA
          ! Convergence from above +/- expected delta
          ! simply relative to previous energy
          A2Logic=ETot%D(cSCF  ,iCLONE)*(One+ETest)<ETot%D(cSCF-1,iCLONE)
          ! Met all criteria
          CLogic=DMaxB<DTest.AND.ETotQ<ETest
          ! Exceeded density criteria
          DLogic=DMaxB<5D-2*DTest
          ! Exceeded energy criteria
          ELogic=ETotQ<3D-2*ETest
          ! Quasi convergence from below (bad)
          QLogic=(.NOT.ALogic).AND.DLogic.AND.ELogic
          ! Going to wrong state with DIIS
          ILogic=DoDIIS.AND.DLogic.AND.(.NOT.ELogic)
          ! DIIS is oscillating
          OLogic=DoDIIS.AND.(.NOT.ALogic).AND.( (ETotQ>1D-4.AND.(DMaxQ>1D0.AND.DIISQ>1D0)).OR. &
                                                (ETotQ>2D-3.AND.(DMaxQ>1D0.OR.DIISQ>1D0)) )
          ! Maybe DIIS would be a good idea
          GLogic=DoODA.AND.cSCF>4.AND.DIISB<75D-5.AND.ETotQ<5D-5.AND.DMaxB<1D-1
          ! Turn on KS recomputation if ODA is not strictly decreasing
          RLogic=DoODA.AND.(cSCF>2.AND.ODAQ>1D-4.OR..NOT.ALogic) !1D1*ETotQ
          ! If we are increasing with ODA and rebuild is on, we are well fucked.
          FLogic=DoODA.AND.RebuildPostODA.AND..NOT.ALogic
          ! Sort through logic hopefully in the conditionally correct order ...


          WRITE(Out,*)' ETest  = ',ETest
          WRITE(Out,*)' DTest  = ',DTest
          WRITE(Out,*)' ALogic = ',ALogic
          WRITE(Out,*)' A2Logic = ',A2Logic
          WRITE(Out,*)' ELogic = ',ELogic
          WRITE(Out,*)' CLogic = ',CLogic
          WRITE(Out,*)' DLogic = ',DLogic
          WRITE(Out,*)' QLogic = ',QLogic
          WRITE(Out,*)' ILogic = ',ILogic
          WRITE(Out,*)' GLogic = ',GLogic
          WRITE(Out,*)' RLogic = ',RLogic
          WRITE(Out,*)' FLogic = ',FLogic

          CLOSE(Out)

          IF(ALogic.AND.CLogic)THEN
             Converged(iCLONE)=DID_CONVERGE
             Mssg='Normal SCF convergence.'
          ELSEIF(A2Logic.AND.DLogic)THEN
             Converged(iCLONE)=DID_CONVERGE
             Mssg='Convergence of density only'
          ELSEIF(ALogic.AND.ELogic)THEN
             Converged(iCLONE)=DID_CONVERGE
             Mssg='Convergence of energy only'
          ELSEIF(QLogic)THEN
             Converged(iCLONE)=DID_CONVERGE
             Mssg='Quasi convergence from wrong side.'
          ELSEIF(FLogic)THEN
             Mssg='ODA with rebuild not strictly decreasing, as good as we can do for now ...'
             Converged(iCLONE)=DID_CONVERGE
          ELSEIF(RLogic)THEN
             Mssg='Rebuilding Fockian post ODA'
             Converged(iCLONE)=SCF_STALLED
             ! We may have a corrupted history here, rest...
             ETot%D(0:iSCF-1,iCLONE)=BIG_DBL
!          ELSEIF(ILogic)THEN
!             Converged(iCLONE)=SCF_STALLED
!             Mssg='DIIS converging to wrong state'
          ELSEIF(OLogic)THEN
             Converged(iCLONE)=SCF_STALLED
             Mssg='DIIS oscillation'
          ELSEIF(GLogic)THEN
             Mssg='Turning DIIS on'
             Converged(iCLONE)=DIIS_NOPATH
          ENDIF
       ENDIF
       IF(DoDIIS.AND.DIISB<DIISA)THEN
          ! If DIIS is making progress, then turn on archivation of the density
          CALL Put(.TRUE.,'ArchiveDensity')
       ELSE
          ! otherwise, dont archive a potential instability
          CALL Put(.FALSE.,'ArchiveDensity')
       ENDIF       
       CALL CloseHDFGroup(HDF_CurrentID)
    ENDDO
    CALL CloseHDF(HDFFileID)
    IF(cSCF>1)ConvergedQ=NOT_CONVERGE
    DO iCLONE=1,G%Clones
       ConvergedQ=MAX(ConvergedQ,Converged(iCLONE))
    ENDDO
    ! Convergence announcement
    IF(ConvergedQ.NE.NOT_CONVERGE.AND.cSCF>2)THEN
       CALL OpenASCII(OutFile,Out)
       WRITE(Out,*)TRIM(Mssg)
       WRITE(*,*)TRIM(Mssg)
       CLOSE(Out)
    ENDIF
  END FUNCTION ConvergedQ
  !===============================================================================
  ! BUILD A HGTF DENSITY BY HOOK OR BY CROOK
  !===============================================================================
  SUBROUTINE DensityBuild(N,S,M)
    TYPE(FileNames):: N
    TYPE(State)    :: S
    TYPE(Parallel) :: M    
    !----------------------------------------------------------------------------!

    IF(TRIM(S%Action%C(1))/=SCF_BASISSETSWITCH.AND.   &
       TRIM(S%Action%C(1))/=SCF_DENSITY_NORMAL.AND.   &
       TRIM(S%Action%C(1))/=CPSCF_START_RESPONSE.AND. &
       TRIM(S%Action%C(1))/=CPSCF_DENSITY_NORMAL) THEN
       CALL Invoke('P2Use',N,S,M)
    ENDIF
    ! Build some density ...
    ! Hold on, this is fucked up--we have got to unify MakeRho!!!!!!
    ! Chee Kwan and CJ, please, pretty please....
    IF(S%Action%C(1)/=CPSCF_START_RESPONSE) THEN
       CALL Invoke('MakeRho',N,S,M)
#ifdef PARALLEL
       CALL Invoke('ParaMakeRho',N,S,M) 
#endif
    ENDIF
  END SUBROUTINE DensityBuild
  !===============================================================================

  !===============================================================================
  SUBROUTINE DensityLogic(cSCF,cBAS,cGEO,S,O,CPSCF_O)
    TYPE(State)      :: S
    TYPE(Options)    :: O
    INTEGER          :: cSCF,cBAS,cGEO,pBAS
    LOGICAL          :: DoCPSCF
    LOGICAL,OPTIONAL :: CPSCF_O
    !----------------------------------------------------------------------------!
    pBAS=S%Previous%I(2)
    S%Current%I=(/cSCF,cBAS,cGEO/)
    IF(PRESENT(CPSCF_O))THEN
       DoCPSCF=CPSCF_O
    ELSE
       DoCPSCF=.FALSE.
    ENDIF
    IF(DoCPSCF)THEN
       IF(O%Guess==GUESS_EQ_DIPOLE.AND.cSCF==0)THEN
          O%Guess=0
          S%Previous%I=S%Current%I
          S%Action%C(1)=CPSCF_START_RESPONSE
       ELSE
          S%Action%C(1)=CPSCF_DENSITY_NORMAL
       ENDIF
    ELSE
       IF(O%Guess==GUESS_EQ_SUPR)THEN 
          O%Guess=0
          S%Previous%I=S%Current%I
          S%Action%C(1)=SCF_SUPERPOSITION
       ELSEIF(O%Guess==GUESS_EQ_CORE)THEN
          O%Guess=0
          S%Previous%I=S%Current%I
          S%Action%C(1)=SCF_GUESSEQCORE
       ELSEIF(O%Guess==GUESS_EQ_RESTART)THEN
          O%Guess=0
          S%Previous%I=O%RestartState%I
          S%Action%C(1)=SCF_RESTART
       ELSEIF(cSCF==0.AND.cBAS==pBAS.AND.cGEO/=1)THEN
!          S%Action%C(1)=SCF_PROJECTION
          S%Action%C(1)=SCF_EXTRAPOLATE
       ELSEIF(cSCF==0.AND.cBAS/=pBAS)THEN
          S%Action%C(1)=SCF_BASISSETSWITCH
          S%Previous%I(1)=S%Previous%I(1)+1
       ELSE
          S%Action%C(1)=SCF_DENSITY_NORMAL
       ENDIF
    ENDIF
  END SUBROUTINE DensityLogic
  !===============================================================================
  ! BUILD A FOCK MATRIX
  !===============================================================================
  SUBROUTINE FockBuild(cSCF,cBAS,N,S,O,M)
    TYPE(FileNames):: N
    TYPE(State)    :: S
    TYPE(Options)  :: O
    TYPE(Parallel) :: M    
    REAL(DOUBLE)   :: Lambda
    INTEGER        :: cSCF,cBAS,Modl
    LOGICAL        :: DoDIIs
    !----------------------------------------------------------------------------!
    Modl=O%Models(cBAS)

    IF(S%Action%C(1).NE.CPSCF_START_RESPONSE) CALL Invoke('QCTC',N,S,M)

    IF(S%Action%C(1)/=SCF_GUESSEQCORE.AND.S%Action%C(1).NE.CPSCF_START_RESPONSE)THEN
       IF(HasHF(Modl)) CALL Invoke('ONX',N,S,M)
       IF(HasDFT(Modl))CALL Invoke('HiCu',N,S,M)
    ENDIF

    CALL Invoke('FBuild',N,S,M)
  END SUBROUTINE FockBuild
  !===============================================================================
  ! EXTRAPOLATE (DIIS)
  !===============================================================================
  SUBROUTINE Xtra(cSCF,cBAS,N,S,O,M)
    TYPE(FileNames):: N
    TYPE(State)    :: S
    TYPE(Options)  :: O
    TYPE(Parallel) :: M    
    REAL(DOUBLE)   :: Lambda
    INTEGER        :: cSCF,cBAS,Modl
    LOGICAL        :: DoDIIs
    !----------------------------------------------------------------------------!
    DoDIIS=cSCF>0
    IF(DoDIIS)THEN
       SELECT CASE(S%Action%C(1))
       CASE(CPSCF_SOLVE_SCF,CPSCF_START_RESPONSE,CPSCF_DENSITY_NORMAL,CPSCF_FOCK_BUILD)
          CALL Invoke('DDIIS',N,S,M)
       CASE DEFAULT
          CALL Invoke('DIIS',N,S,M)
       END SELECT
    ENDIF
  END SUBROUTINE Xtra
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

    IF(    S%Action%C(1)==CPSCF_SOLVE_SCF.OR. &
         & S%Action%C(1)==CPSCF_START_RESPONSE)THEN
       CALL Invoke('TC2Response',N,S,M)    
    ELSEIF(O%Methods(cBAS)==RH_R_SCF)THEN
       CALL Invoke('RHeqs',N,S,M)
    ELSEIF(O%Methods(cBAS)==SDMM_R_SCF) THEN
       CALL Invoke('SDMM',N,S,M)
    ELSEIF(O%Methods(cBAS)==PM_R_SCF) THEN
       CALL Invoke('PM',N,S,M)
    ELSEIF(O%Methods(cBAS)==SP2_R_SCF) THEN
       CALL Invoke('SP2',N,S,M)
    ELSEIF(O%Methods(cBAS)==SP4_R_SCF) THEN
       CALL Invoke('SP4',N,S,M)
    ELSEIF(O%Methods(cBAS)==TS4_R_SCF) THEN
       CALL Invoke('TS4',N,S,M)
    ELSE
       CALL MondoHalt(99,'Unknown method key = '//TRIM(IntToChar(O%Methods(cBAS))))
    ENDIF
!!$<<<<<<< SCFs.F90
!!$    !
!!$    SELECT CASE(S%Action%C(1))
!!$    CASE(CPSCF_SOLVE_SCF,CPSCF_START_RESPONSE)
!!$       CALL Invoke('CPSCFStatus',N,S,M)
!!$    CASE DEFAULT
!!$       CALL Invoke('SCFstats',N,S,M)
!!$    END SELECT
!!$    !
!!$=======
!!$>>>>>>> 1.42
  END SUBROUTINE SolveSCF
  !===============================================================================

  !===============================================================================
  SUBROUTINE OneEMats(cBAS,cGEO,N,B,S,O,M)
    TYPE(FileNames):: N
    TYPE(BasisSets):: B
    TYPE(State)    :: S
    TYPE(Options)  :: O
    TYPE(Parallel) :: M
    INTEGER        :: cBAS,cGEO,pBAS
    LOGICAL, SAVE  :: DoPFFT=.TRUE.
    !----------------------------------------------------------------------------!
    pBAS=S%Previous%I(2)    
    S%Current%I=(/0,cBAS,cGEO/)
    S%Action%C(1)='OneElectronMatrices'
    IF(pBAS/=cBAS)DoPFFT=.TRUE.
    IF(DoPFFT)THEN
       CALL Invoke('MakePFFT',N,S,M)
       DoPFFT=.FALSE.
    ENDIF
    CALL Invoke('MakeS',N,S,M)
    IF(O%Methods(cBAS)==RH_R_SCF)THEN
       CALL Invoke('LowdinO',N,S,M)
    ELSE
       CALL Invoke('AInv',N,S,M)
    ENDIF
    ! Kinetic energy matrix T
    CALL Invoke('MakeT',N,S,M)
    IF(B%BSets(cBAS,1)%HasECPs)THEN
       ! Make the ECP matrix U 
       CALL Invoke('MakeU',N,S,M)
    ENDIF
  END SUBROUTINE OneEMats
  !===============================================================================
  ! COMPUTE AN ENERGY GRADIENT
  !===============================================================================
  SUBROUTINE Force(cBAS,cGEO,N,O,S,G,B,M)
    TYPE(FileNames)    :: N
    TYPE(Options)      :: O
    TYPE(State)        :: S
    TYPE(Geometries)   :: G
    TYPE(Parallel)     :: M
    TYPE(BasisSets)    :: B
    INTEGER            :: cBAS,cGEO,K,J,iATS,iCLONE
    CHARACTER(LEN=DCL) :: chGEO    
    REAL(DOUBLE)       :: GradVal
    !----------------------------------------------------------------------------!
    CALL New(S%Action,1)
    ! Initialize the force vector in HDF, clone by clone
    chGEO=IntToChar(cGEO)
    DO iCLONE=1,G%Clones
       G%Clone(iCLONE)%Gradients%D=BIG_DBL
       G%Clone(iCLONE)%GradRMS=&
            SQRT(G%Clone(iCLONE)%GradRMS)/DBLE(3*G%Clone(iCLONE)%NAtms)
    ENDDO
    CALL GeomArchive(cBAS,cGEO,N,B,G)
    ! Now evaluate the forces
    S%Action%C(1)='ForceEvaluation'
    ! The non-orthogonal response    
    CALL Invoke('SForce',N,S,M)
    ! Kinetic energy piece
    CALL Invoke('TForce',N,S,M)
    IF(B%BSets(cBAS,1)%HasECPs)THEN
       ! Compute ECP component of the force
       CALL Invoke('UForce',N,S,M)
    ENDIF
    ! Build density with last DM
    CALL Invoke('MakeRho',N,S,M)
#ifdef PARALLEL
    CALL Invoke('ParaMakeRho',N,S,M)
#endif
    ! Coulomb part
    CALL Invoke('JForce',N,S,M)
    ! Exact Hartree-Fock exchange component
    IF(HasHF(O%Models(cBas)))THEN
    !  CALL NXForce(cBAS,cGEO,N,G,B,S,M)
    !  CALL Invoke('XForce',N,S,M)
       CALL Invoke('GONX2',N,S,M)
    ENDIF
    ! DFT exchange corrleation term
    IF(HasDFT(O%Models(cBas))) THEN
       CALL Invoke('XCForce',N,S,M)
    ENDIF
!
!   Constraint the Gradients
!
    chGEO=IntToChar(cGEO)
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       CALL Get(G%Clone(iCLONE)%Gradients,'Gradients',Tag_O=chGEO)
       IF(O%Coordinates/=GRAD_INTS_OPT) THEN
         DO iATS=1,G%Clone(iCLONE)%NAtms
            IF(G%Clone(iCLONE)%CConstrain%I(iATS)/=0)THEN
               G%Clone(iCLONE)%Gradients%D(:,iATS)=Zero
            ENDIF
         ENDDO
       ENDIF
!      Close the group
       CALL CloseHDFGroup(HDF_CurrentID)
       G%Clone(iCLONE)%GradRMS=SQRT(G%Clone(iCLONE)%GradRMS)/DBLE(3*G%Clone(iCLONE)%NAtms)
!      Print the Forces and BoxForces if O%Grad==GRAD_ONE_FORCE
       IF(O%Grad==GRAD_ONE_FORCE) THEN
          CALL Print_Force(G%Clone(iCLONE))
       ENDIF
    ENDDO
    CALL CloseHDF(HDFFileID)
!   NEB force projections
    IF(O%Grad==GRAD_TS_SEARCH_NEB) THEN
       CALL NEBForce(G,O)
    ENDIF
!   Zero forces on constrained atoms and compute stats with projected forces
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       G%Clone(iCLONE)%GradMax=Zero
       G%Clone(iCLONE)%GradRMS=Zero
       DO iATS=1,G%Clone(iCLONE)%NAtms
          IF(G%Clone(iCLONE)%CConstrain%I(iATS)==0)THEN
             DO J=1,3
                GradVal=G%Clone(iCLONE)%Gradients%D(J,iATS)
                G%Clone(iCLONE)%GradRMS=G%Clone(iCLONE)%GradRMS+GradVal
                G%Clone(iCLONE)%GradMax=MAX(G%Clone(iCLONE)%GradMax,ABS(GradVal))
             ENDDO
          ENDIF
       ENDDO
       ! Put the zeroed forces back ...
       ! ... and close the group
       CALL CloseHDFGroup(HDF_CurrentID)
       G%Clone(iCLONE)%GradRMS=SQRT(G%Clone(iCLONE)%GradRMS)/DBLE(3*G%Clone(iCLONE)%NAtms)
    ENDDO
    ! Now close the HDF file ..
    CALL CloseHDF(HDFFileID)
    CALL Delete(S%Action)
  END SUBROUTINE Force
  !===============================================================================
!!$<<<<<<< SCFs.F90
!!$  !
!!$  !===============================================================================
!!$  FUNCTION ConvergedQ(cSCF,cBAS,N,S,O,G,ETot,DMax,DIIS,CPSCF_O)
!!$    TYPE(FileNames)             :: N
!!$    TYPE(State)                 :: S
!!$    TYPE(Options)               :: O
!!$    TYPE(Geometries)            :: G
!!$    TYPE(Parallel)              :: M
!!$    TYPE(DBL_RNK2)              :: ETot,DMax,DIIS
!!$    LOGICAL,OPTIONAL            :: CPSCF_O
!!$    LOGICAL                     :: CPSCF
!!$    INTEGER                     :: cSCF,cBAS,iGEO,iCLONE
!!$    REAL(DOUBLE)                :: DIISA,DIISB,DDIIS,DIISQ,       &
!!$                                   DETOT,ETOTA,ETOTB,ETOTQ,ETEST, &
!!$                                   DDMAX,DMAXA,DMAXB,DMAXQ,DTEST
!!$    LOGICAL,DIMENSION(G%Clones) :: Converged
!!$    LOGICAL                     :: ConvergedQ
!!$    CHARACTER(LEN=DCL)            :: chGEO
!!$    !----------------------------------------------------------------------------!
!!$    !
!!$    IF(PRESENT(CPSCF_O)) THEN
!!$       CPSCF=CPSCF_O
!!$    ELSE
!!$       CPSCF=.FALSE.
!!$    ENDIF
!!$    ! Convergence thresholds
!!$    IF(CPSCF) THEN
!!$       ETest=RTol(O%AccuracyLevels(cBAS))
!!$    ELSE
!!$       ETest=ETol(O%AccuracyLevels(cBAS))
!!$    ENDIF
!!$    DTest=DTol(O%AccuracyLevels(cBAS))
!!$    IF(cSCF==0)THEN
!!$       ConvergedQ=.FALSE.
!!$       RETURN
!!$    ENDIF
!!$    ! Accumulate current statistics
!!$    chGEO=IntToChar(iGEO)
!!$    HDFFileID=OpenHDF(N%HFile)
!!$    DO iCLONE=1,G%Clones
!!$       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
!!$       ! Gather convergence parameters
!!$#ifdef PARALLEL_CLONES
!!$       IF(CPSCF) THEN
!!$          CALL Get(Etot%D(cSCF,iCLONE),'Prop'    )
!!$          CALL Get(DMax%D(cSCF,iCLONE),'DPrimMax')
!!$          CALL Get(DIIS%D(cSCF,iCLONE),'DDIISErr')
!!$       ELSE
!!$          CALL Get(Etot%D(cSCF,iCLONE),'Etot')
!!$          CALL Get(DMax%D(cSCF,iCLONE),'DMax')
!!$          CALL Get(DIIS%D(cSCF,iCLONE),'DIISErr' )
!!$       ENDIF
!!$#else
!!$       IF(CPSCF) THEN
!!$          CALL Get(Etot%D(cSCF,iCLONE),'Prop'    ,StatsToChar(S%Current%I))
!!$          !Etot%D(cSCF,iCLONE)=0.0d0 !vw MUST BE CHANGED
!!$          CALL Get(DMax%D(cSCF,iCLONE),'DPrimMax',StatsToChar(S%Current%I))
!!$          CALL Get(DIIS%D(cSCF,iCLONE),'DDIISErr',StatsToChar(S%Current%I))
!!$       ELSE
!!$          CALL Get(Etot%D(cSCF,iCLONE),'Etot',StatsToChar(S%Current%I))
!!$          CALL Get(DMax%D(cSCF,iCLONE),'DMax',StatsToChar(S%Current%I))
!!$          CALL Get(DIIS%D(cSCF,iCLONE),'DIISErr',StatsToChar(S%Current%I))
!!$       ENDIF
!!$#endif
!!$       CALL CloseHDFGroup(HDF_CurrentID)
!!$       ! Load current energies
!!$       G%Clone(iCLONE)%ETotal=ETot%D(cSCF,iCLONE)
!!$
!!$       Converged(iCLONE)=.FALSE.
!!$       IF(cSCF>1)THEN          
!!$          ETotA=ETot%D(cSCF-1,iCLONE)
!!$          ETotB=ETot%D(cSCF  ,iCLONE)
!!$          DMaxA=DMax%D(cSCF-1,iCLONE)
!!$          DMaxB=DMax%D(cSCF  ,iCLONE)
!!$          DIISA=DIIS%D(cSCF-1,iCLONE)
!!$          DIISB=DIIS%D(cSCF  ,iCLONE)
!!$          ! Absolute numbers
!!$          dETot=ABS(ETotA-ETotB)
!!$          dDMax=ABS(DMaxA-DMaxB)
!!$          dDIIS=ABS(DIISA-DIISB)
!!$          ! Relative numbers (Quotients)
!!$          ETotQ=dETot/ABS(ETotB)
!!$!IF(CPSCF) write(*,*) 'dETot',dETot,' ETest',ETest
!!$!IF(CPSCF) write(*,*) 'DMaxB',DMaxB,' dTest',dTest
!!$          DMaxQ=dDMax/ABS(DMaxB+1.D-50)
!!$          DIISQ=dDIIS/ABS(DIISB+1.D-50)
!!$          CALL OpenASCII(OutFile,Out)
!!$!          WRITE(Out,*)'ETest = ',ETest
!!$!          WRITE(Out,*)'DTest = ',DTest
!!$!          WRITE(Out,*)'ETotQ = ',ETotQ
!!$!          WRITE(Out,*)'ETotA = ',ETotA
!!$!          WRITE(Out,*)'ETotB = ',ETotB
!!$!          WRITE(Out,*)'DIISQ = ',DIISQ
!!$!          WRITE(Out,*)'DMaxQ = ',DMaxQ
!!$!          WRITE(Out,*)'DIISA = ',DIISA
!!$!          WRITE(Out,*)'DIISB = ',DIISB
!!$!          WRITE(Out,*)'DMaxA = ',DMaxA
!!$!          WRITE(Out,*)'DMaxB = ',DMaxB
!!$          CLOSE(Out)
!!$          ! Convergence tests
!!$!         IF(((DMaxB<dTest.AND.ETotQ<ETest).OR.DMaxB<5D-1*dTest))THEN
!!$          IF(CPSCF) THEN
!!$             IF((DMaxB<dTest.AND.ETotQ<ETest).OR.DMaxB<dTest/10.0d0)THEN
!!$                Converged(iCLONE)=.TRUE.
!!$                Mssg='Normal SCF convergence.a'
!!$             ENDIF
!!$          ELSE
!!$             IF(((DMaxB<dTest.AND.ETotQ<ETest).OR.DMaxB<5D-1*dTest).AND.ETotB<ETotA)THEN
!!$                Converged(iCLONE)=.TRUE.
!!$                Mssg='Normal SCF convergence.a'
!!$             ENDIF
!!$             ! Accept convergence from wrong side if DM thresholds are tightend.
!!$             IF(DMaxB<dTest*75D-2.AND.ETotQ<ETest*3D-1)THEN
!!$                !        IF(DMaxB<dTest*1D-1.AND.ETotQ<ETest*1D-1)THEN
!!$                Converged(iCLONE)=.TRUE.
!!$                Mssg='Normal SCF convergence.b'
!!$             ENDIF
!!$          ENDIF
!!$          ! Look for stall out if we have at least one consecutive digit in the DM
!!$          IF(DMaxB<1.D-1.AND.DMaxA<1.D-1)THEN
!!$             ! Look for non-decreasing errors due to incomplete numerics
!!$             IF(DIISQ<1.D-1.AND.DMaxQ<1.D-1.AND.cSCF>6)THEN
!!$                IF(DIISB>DIISA.AND.DMaxB>DMaxA)THEN
!!$                   Mssg='SCF hit DIIS & DMax increase.'
!!$                   Converged(iCLONE)=.TRUE.
!!$                ENDIF
!!$             ELSEIF(DIISQ<1.D-2.AND.DMaxQ<1.D-2.AND.cSCF>6)THEN
!!$                IF(DIISB>DIISA)THEN
!!$                   Mssg='SCF hit DIIS increase'
!!$                   Converged(iCLONE)=.TRUE.
!!$                ELSEIF(DMaxQ<1D-1.AND.DMaxB>DMaxA)THEN
!!$                   Mssg='SCF hit DMAX increase'
!!$                   Converged(iCLONE)=.TRUE.
!!$                ENDIF
!!$             ELSEIF((DIISQ<1D-4.OR.DMaxQ<1D-4).AND.cSCF>6)THEN
!!$                Mssg='SCF convergence due to DIIS stagnation.'
!!$                Converged(iCLONE)=.TRUE.
!!$             ENDIF
!!$          ENDIF
!!$       ENDIF
!!$    ENDDO
!!$    CALL CloseHDF(HDFFileID)
!!$    IF(cSCF>1)ConvergedQ=.TRUE.
!!$    DO iCLONE=1,G%Clones
!!$       ConvergedQ=ConvergedQ.AND.Converged(iCLONE)
!!$    ENDDO
!!$    ! Convergence announcement
!!$    IF(ConvergedQ)THEN!.AND.PrintFlags%Key>DEBUG_NONE)THEN
!!$       CALL OpenASCII(OutFile,Out)
!!$       WRITE(Out,*)TRIM(Mssg)
!!$       WRITE(*,*)TRIM(Mssg)
!!$!       WRITE(Out,*)'Normal SCF convergence.'
!!$       CLOSE(Out)
!!$    ENDIF
!!$  END FUNCTION ConvergedQ
!!$  !===============================================================================
!!$=======
!!$>>>>>>> 1.42
  ! Numerically compute gradients of the exact HF exchange
  !===============================================================================
  SUBROUTINE NXForce(cBAS,cGEO,N,G,B,S,M)
    TYPE(FileNames)  :: N
    TYPE(Options)    :: O
    TYPE(State)      :: S
    TYPE(Geometries) :: G
    TYPE(BasisSets)  :: B
    TYPE(Parallel)   :: M
    TYPE(DBL_RNK2)   :: GradAux1,GradAux2
    INTEGER          :: cBAS,cGEO,J,iATS,iCLONE
    CHARACTER(LEN=DCL) :: chGEO,chBAS,chSCF
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
    CALL New(GradAux1,(/3,G%Clone(1)%NAtms/))
    CALL New(GradAux2,(/3,G%Clone(1)%NAtms/))
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
                !
                ! Move the atom.
                IF(II==1) THEN
                   G%Clone(iCLONE)%AbCarts%D(IX,AtA)=GTmp(iCLONE)%AbCarts%D(IX,AtA)+DDelta
                ELSEIF(II==2) THEN
                   G%Clone(iCLONE)%AbCarts%D(IX,AtA)=GTmp(iCLONE)%AbCarts%D(IX,AtA)-DDelta
                ENDIF
             ENDDO
             !
             CALL GeomArchive(cBAS,cGEO,N,B,G)    
             ! vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw>>>
             ! Move back the atom.
             DO iCLONE=1,G%Clones
                IF(II==1) THEN
                   G%Clone(iCLONE)%AbCarts%D(IX,AtA)=GTmp(iCLONE)%AbCarts%D(IX,AtA)
                ELSEIF(II==2) THEN
                   G%Clone(iCLONE)%AbCarts%D(IX,AtA)=GTmp(iCLONE)%AbCarts%D(IX,AtA)
                ENDIF
             ENDDO
             ! vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw<<<
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
#ifdef NGONX_INFO!vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw
             write(*,*) 'FX(',AtA,IX,')=',FX(iCLONE,IA)
#endif!vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw
          ENDDO
       ENDDO
    ENDDO
#ifdef NGONX_INFO!vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw
    DO iCLONE=1,G%Clones
       WRITE(*,'(A,I3,A,E20.12)') ' Total XForce(',iCLONE,') =',SUM(FX(iCLONE,:))
    ENDDO
#endif!vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw
   ! Add in the forces to the global gradient and put back to HDF
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       CALL Get(GradAux1,'gradients',Tag_O=chGEO)
       CALL CartRNK1ToCartRNK2(FX(iCLONE,:),GradAux2%D)
       GradAux1%D=GradAux1%D+GradAux2%D     
       CALL Put(GradAux1,'gradients',Tag_O=chGEO)
       ! Close the group
       CALL CloseHDFGroup(HDF_CurrentID)
    ENDDO
    CALL CloseHDF(HDFFileID)
    DO iCLONE=1,G%Clones
       G%Clone(iCLONE)%AbCarts%D=GTmp(iCLONE)%AbCarts%D
       CALL Delete(GTmp(iCLONE))
       CALL Delete(P(iCLONE))
    ENDDO
    CALL Delete(GradAux1)
    CALL Delete(GradAux2)
    CALL Delete(BSiz)
    CALL Delete(OffS)
    CALL Delete(K)
  END SUBROUTINE NXForce
END MODULE SCFs




