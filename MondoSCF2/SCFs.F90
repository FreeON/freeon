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
    TYPE(Controls)       :: C
    TYPE(DBL_RNK2),SAVE  :: ETot,DMax,DIIS
    INTEGER,PARAMETER    :: MaxSCFs=256
    INTEGER              :: cBAS,cGEO,iSCF
    !----------------------------------------------------------------------------!
    CALL New(C%Stat%Action,1)
    ! Determine if there was a geomety or Basis Set Change
    CALL SameBasisSameGeom(cBAS,cGEO,C%Nams,C%Opts,C%Stat)
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
    TYPE(FileNames)    :: N
    TYPE(State)        :: S
    TYPE(Options)      :: O
    TYPE(Geometries)   :: G
    TYPE(Parallel)     :: M
    TYPE(DBL_RNK2)     :: ETot,DMax,DIIS
    INTEGER            :: cSCF,cBAS,cGEO,iCLONE,Modl,IConAls
    INTEGER,SAVE       :: SCF_STATUS
    LOGICAL,OPTIONAL   :: CPSCF_O
    LOGICAL            :: SCFCycle,DoCPSCF
    LOGICAL,SAVE       :: DIIS_FAIL,ODA_DONE
    REAL(DOUBLE)       :: DIISErr
    CHARACTER(LEN=128) :: Tmp
!----------------------------------------------------------------------------!
!   Initailize
    SCFCycle=.FALSE.
    IF(cSCF==0) THEN
       SCF_STATUS = NOT_CONVERGE
       ODA_DONE   = .FALSE.
       DIIS_FAIL  = .FALSE.
    ENDIF
!   Are we maybe solving CPSCF equations?
    IF(PRESENT(CPSCF_O))THEN
       DoCPSCF=CPSCF_O
    ELSE
       DoCPSCF=.FALSE.
    ENDIF
!   The options...
    IF(DoCPSCF)THEN
       CALL DensityLogic(cSCF,cBAS,cGEO,N,S,O,CPSCF_O=.TRUE.)
       CALL DensityBuild(N,S,M)
       IF(cSCF.EQ.0)S%Action%C(1)=CPSCF_START_RESPONSE
       IF(cSCF.GT.0)S%Action%C(1)=CPSCF_FOCK_BUILD
       IF(cSCF.GT.0)THEN
          CALL Invoke('QCTC',N,S,M)
          Modl=O%Models(cBAS)
          IF(HasHF(Modl)) CALL Invoke('ONX',N,S,M)
          IF(HasDFT(Modl)) THEN
             CALL Halt('SCFs: DFT-Response not yet supported!')
             CALL Invoke('HiCu',N,S,M)
          ENDIF
       ENDIF
       CALL Invoke('FBuild',N,S,M)
       IF(cSCF.GT.0) CALL Invoke('DDIIS',N,S,M)
       S%Action%C(1)=CPSCF_SOLVE_SCF
       CALL SolveSCF(cBAS,N,S,O,M)
       CALL Invoke('CPSCFStatus',N,S,M)
    ELSE
!      Logic for Algorithm Choice
       IF(O%ConAls(cBAS)==ODMIX_CONALS) THEN
          IF(ODA_DONE) THEN
             IConAls = DIIS_CONALS
          ELSE
             IF(SCF_STATUS==DIIS_NOPATH) THEN
                IConAls   = DIIS_CONALS
                ODA_DONE  = .TRUE.
                Mssg = 'Turning on DIIS'
                CALL OpenASCII(OutFile,Out)
                WRITE(Out,*)TRIM(Mssg)
                WRITE(*,*)TRIM(Mssg)
                CLOSE(Out)
             ELSE
                IConAls = ODA_CONALS
             ENDIF
          ENDIF
       ELSEIF(O%ConAls(cBAS)==DOMIX_CONALS) THEN
          IF(DIIS_FAIL) THEN
             IConAls = ODA_CONALS
          ELSE
             IF(SCF_STATUS==SCF_STALLED) THEN
                IConAls     = ODA_CONALS
                DIIS_FAIL   = .TRUE.
                Mssg = 'DIIS Failed: Turning on ODA'
                CALL OpenASCII(OutFile,Out)
                WRITE(Out,*)TRIM(Mssg)
                WRITE(*,*)TRIM(Mssg)
                CLOSE(Out)
             ELSE
                IConAls = DIIS_CONALS
             ENDIF
          ENDIF
       ELSE
          IConAls = O%ConAls(cBAS)
       ENDIF
!      Parse for strict ODA or DIIS Over-Ride
       CALL OpenASCII(N%IFile,Inp)
       IF(OptKeyQ(Inp,CONALS_OVRIDE,CONALS_ODA))  IConAls = ODA_CONALS
       IF(OptKeyQ(Inp,CONALS_OVRIDE,CONALS_DIIS)) IConAls = DIIS_CONALS
       CLOSE(Inp)
!      Defaults
       IF(cSCF < 1)                                         IConAls = NO_CONALS
       IF(O%ConAls(cBAS)==SMIX_CONALS        .AND. cSCF<2)  IConAls = NO_CONALS
       IF(S%Action%C(1) ==SCF_GUESSEQCORE    .AND. cSCF<2)  IConAls = NO_CONALS 
       IF(S%Action%C(1) ==SCF_BASISSETSWITCH .AND. cSCF<2)  IConAls = NO_CONALS 
       IF(S%Action%C(1) ==SCF_RWBSS          .AND. cSCF<2)  IConAls = NO_CONALS 
!      Select the Case
       SELECT CASE (IConAls)
       CASE (DIIS_CONALS)
          CALL DensityLogic(cSCF,cBAS,cGEO,N,S,O)
          CALL DensityBuild(N,S,M)
          CALL FockBuild(cSCF,cBAS,N,S,O,M)
          CALL Invoke('DIIS',N,S,M)
          CALL SolveSCF(cBAS,N,S,O,M)
          CALL Invoke('SCFstats',N,S,M)
       CASE (ODA_CONALS)
          CALL DensityLogic(cSCF,cBAS,cGEO,N,S,O)
          CALL DensityBuild(N,S,M)
          CALL FockBuild(cSCF,cBAS,N,S,O,M)
          CALL SolveSCF(cBAS,N,S,O,M)
          CALL Invoke('ODA',N,S,M)
          IF(HasDFT(O%Models(cBAS)))THEN
!            Rebuild non-linear KS matrix
             CALL DensityLogic(cSCF,cBAS,cGEO,N,S,O)
             CALL DensityBuild(N,S,M)
             CALL Invoke('HiCu',N,S,M)
             CALL Invoke('FBuild',N,S,M)
          ENDIF
          CALL SolveSCF(cBAS,N,S,O,M)
          CALL Invoke('SCFstats',N,S,M)
       CASE (SMIX_CONALS)
          CALL DensityLogic(cSCF,cBAS,cGEO,N,S,O)
          CALL DensityBuild(N,S,M)
          CALL FockBuild(cSCF,cBAS,N,S,O,M)
          CALL SolveSCF(cBAS,N,S,O,M)
          CALL Invoke('SCFstats',N,S,M)
          S%Action%C(1)='Stanton MIX'
          CALL Invoke('MIX',N,S,M)
          CALL DensityLogic(cSCF,cBAS,cGEO,N,S,O)
          CALL DensityBuild(N,S,M)
          CALL FockBuild(cSCF,cBAS,N,S,O,M)
          CALL SolveSCF(cBAS,N,S,O,M)
       CASE (NO_CONALS)          
          CALL DensityLogic(cSCF,cBAS,cGEO,N,S,O)
          CALL DensityBuild(N,S,M)
          CALL FockBuild(cSCF,cBAS,N,S,O,M)
          CALL SolveSCF(cBAS,N,S,O,M)
          CALL Invoke('SCFstats',N,S,M)
       CASE(TEST_CONALS)
          CALL MondoHalt(DRIV_ERROR,'Test Algorithm Not Implimented')
       CASE DEFAULT
          CALL MondoHalt(DRIV_ERROR,'Logic failure in SCFCycle')
       END SELECT 
    ENDIF
!   Archive and Check Status
    CALL StateArchive(N,S)
    SCF_STATUS=ConvergedQ(cSCF,cBAS,N,S,O,G,ETot,DMax,DIIS,IConAls,CPSCF_O)
    S%Previous%I=S%Current%I
    IF(SCF_STATUS==DID_CONVERGE) SCFCycle=.TRUE.
!
  END FUNCTION SCFCycle
  !===============================================================================
  !
  !===============================================================================
  FUNCTION ConvergedQ(cSCF,cBAS,N,S,O,G,ETot,DMax,DIIS,IConAls,CPSCF_O)
    TYPE(FileNames)             :: N
    TYPE(State)                 :: S
    TYPE(Options)               :: O
    TYPE(Geometries)            :: G
    TYPE(Parallel)              :: M
    TYPE(DBL_RNK2)              :: ETot,DMax,DIIS
    LOGICAL,OPTIONAL            :: CPSCF_O
    LOGICAL                     :: DoCPSCF,DoDIIS,DoODA,RebuildPostODA
    LOGICAL                     :: ALogic,BLogic,CLogic,DLogic,ELogic,A2Logic, &
                                   GLogic,QLogic,ILogic,OLogic,RLogic,FLogic
    INTEGER                     :: cSCF,cBAS,iGEO,iCLONE
    REAL(DOUBLE)                :: DIISA,DIISB,DDIIS,DIISQ,       &
         DETOT,ETOTA,ETOTB,ETOTQ,ETEST, &
         DDMAX,DMAXA,DMAXB,DMAXQ,DTEST,ETOTO,ODAQ
    INTEGER,DIMENSION(G%Clones) :: Converged
    INTEGER                     :: ConvergedQ,iSCF,IConAls
    CHARACTER(LEN=DCL)            :: chGEO
    !----------------------------------------------------------------------------!
    !
    IF(PRESENT(CPSCF_O)) THEN
       DoCPSCF=CPSCF_O
    ELSE
       DoCPSCF=.FALSE.
    ENDIF
    ! Convergence thresholds
    IF(DoCPSCF) THEN
       ETest=RTol(O%AccuracyLevels(cBAS))
       DTest=DTol(O%AccuracyLevels(cBAS))
       IF(cSCF==0)THEN
          ConvergedQ=NOT_CONVERGE
          RETURN
       ENDIF
       ! Accumulate current statistics
       chGEO=IntToChar(iGEO)
       HDFFileID=OpenHDF(N%HFile)
       DO iCLONE=1,G%Clones
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          ! Gather convergence parameters
          CALL Get(Etot%D(cSCF,iCLONE),'Prop'    )
          CALL Get(DMax%D(cSCF,iCLONE),'DPrimMax')
          CALL Get(DIIS%D(cSCF,iCLONE),'DDIISErr')
          CALL CloseHDFGroup(HDF_CurrentID)
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
             ! Absolute numbers
             dETot=ABS(ETotA-ETotB)
             dDMax=ABS(DMaxA-DMaxB)
             dDIIS=ABS(DIISA-DIISB)
             ! Relative numbers (Quotients)
             ETotQ=dETot/ABS(ETotB)
             !IF(CPSCF) write(*,*) 'dETot',dETot,' ETest',ETest
             !IF(CPSCF) write(*,*) 'DMaxB',DMaxB,' dTest',dTest
             DMaxQ=dDMax/ABS(DMaxB+1.D-50)
             DIISQ=dDIIS/ABS(DIISB+1.D-50)
             CALL OpenASCII(OutFile,Out)
             !          WRITE(Out,*)'ETest = ',ETest
             !          WRITE(Out,*)'DTest = ',DTest
             !          WRITE(Out,*)'ETotQ = ',ETotQ
             !          WRITE(Out,*)'ETotA = ',ETotA
             !          WRITE(Out,*)'ETotB = ',ETotB
             !          WRITE(Out,*)'DIISQ = ',DIISQ
             !          WRITE(Out,*)'DMaxQ = ',DMaxQ
             !          WRITE(Out,*)'DIISA = ',DIISA
             !          WRITE(Out,*)'DIISB = ',DIISB
             !          WRITE(Out,*)'DMaxA = ',DMaxA
             !          WRITE(Out,*)'DMaxB = ',DMaxB
             CLOSE(Out)
             ! Convergence tests
             !         IF(((DMaxB<dTest.AND.ETotQ<ETest).OR.DMaxB<5D-1*dTest))THEN
             IF((DMaxB<dTest.AND.ETotQ<ETest).OR.DMaxB<dTest/10.0d0)THEN
                Converged(iCLONE)=DID_CONVERGE
                Mssg='Normal CPSCF convergence'
             ENDIF
             ! Look for stall out if we have at least one consecutive digit in the DM
             IF(DMaxB<1.D-1.AND.DMaxA<1.D-1)THEN
                ! Look for non-decreasing errors due to incomplete numerics
                IF(DIISQ<1.D-1.AND.DMaxQ<1.D-1.AND.cSCF>6)THEN
                   IF(DIISB>DIISA.AND.DMaxB>DMaxA)THEN
                      Mssg='CPSCF hit DDIIS & DMax increase.'
                      Converged(iCLONE)=DID_CONVERGE
                   ENDIF
                ELSEIF(DIISQ<1.D-2.AND.DMaxQ<1.D-2.AND.cSCF>6)THEN
                   IF(DIISB>DIISA)THEN
                      Mssg='CPSCF hit DDIIS increase'
                      Converged(iCLONE)=DID_CONVERGE
                   ELSEIF(DMaxQ<1D-1.AND.DMaxB>DMaxA)THEN
                      Mssg='CPSCF hit DMax increase'
                      Converged(iCLONE)=DID_CONVERGE
                   ENDIF
                ELSEIF((DIISQ<1D-4.OR.DMaxQ<1D-4).AND.cSCF>6)THEN
                   Mssg='CPSCF convergence due to DDIIS stagnation.'
                   Converged(iCLONE)=DID_CONVERGE
                ENDIF
             ENDIF
          ENDIF
       ENDDO
       CALL CloseHDF(HDFFileID)
       ! IF(cSCF>1)ConvergedQ=NOT_CONVERGE
       ConvergedQ = DID_CONVERGE
       DO iCLONE=1,G%Clones
          ConvergedQ=MIN(ConvergedQ,Converged(iCLONE))
       ENDDO
       ! Convergence announcement
       IF(ConvergedQ.NE.NOT_CONVERGE.AND.cSCF>2)THEN!.AND.PrintFlags%Key>DEBUG_NONE)THEN
          CALL OpenASCII(OutFile,Out)
          WRITE(Out,*)TRIM(Mssg)
          WRITE(*,*)TRIM(Mssg)
          !WRITE(Out,*)'Normal CPSCF convergence.'
          CLOSE(Out)
       ENDIF
       !
    ELSE
       ! NORMAL HUMANS CONVERGENCE CRITERIA 
       IF(IConAls==DIIS_CONALS) THEN
          DoDIIS=.TRUE.
          DoODA =.FALSE.
       ELSEIF(IConAls==ODA_CONALS) THEN
          DoDIIS=.FALSE.
          DoODA =.TRUE.
       ELSE
          DoDIIS=.FALSE.
          DoODA =.FALSE.
       ENDIF
!        
       ETest=ETol(O%AccuracyLevels(cBAS))
       DTest=DTol(O%AccuracyLevels(cBAS))
!*******
       IF(cSCF<10)THEN
          ConvergedQ=NOT_CONVERGE!.FALSE.
          RETURN
       ENDIF
!*******
       ! Accumulate current statistics
       chGEO=IntToChar(iGEO)
       HDFFileID=OpenHDF(N%HFile)
       DO iCLONE=1,G%Clones
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          ! Gather convergence parameters
          CALL Get(Etot%D(cSCF,iCLONE),'Etot')
          CALL Get(DMax%D(cSCF,iCLONE),'DMax')
          CALL Get(DIIS%D(cSCF,iCLONE),'DIISErr' )
          IF(DoODA.AND.cSCF>1)THEN
             CALL Get(ETotO,'ODAEnergy')
          ELSE
             ETotO=1D10
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
             IF(PrintFlags%Key==DEBUG_MAXIMUM) THEN
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
                CLOSE(Out)
             ENDIF
             Converged(iCLONE)=NOT_CONVERGE
             ! Convergence from above +/- expected delta
             ! relative to historical minimum energy
             ALogic=ETotB*(One+ETest)<ETotA
             ! Convergence from above +/- expected delta
             ! simply relative to previous energy
             A2Logic=ETot%D(cSCF  ,iCLONE)*(One+ETest)<ETot%D(cSCF-1,iCLONE)
             ! Met all criteria
             CLogic=DMaxB<DTest.AND.ETotQ<ETest.AND.DMaxB.NE.Zero
             ! Exceeded density criteria
             DLogic=DMaxB<5D-2*DTest.AND.DMaxB.NE.Zero
             ! Exceeded energy criteria
             ELogic=ETotQ<3D-2*ETest.AND.DMaxB<1D-2
             ! Quasi convergence from below (bad)
             QLogic=(.NOT.ALogic).AND.DLogic.AND.ELogic
             ! Going to wrong state with DIIS
             ILogic=DoDIIS.AND.DLogic.AND.(.NOT.ELogic)
             ! DIIS is oscillating
             OLogic=DoDIIS.AND.(.NOT.ALogic).AND.( (ETotQ>1D-4.AND.(DMaxQ>1D0.AND.DIISQ>1D0)).OR. &
                  (ETotQ>2D-3.AND.(DMaxQ>1D0.OR.DIISQ>1D0)) )
             ! Maybe DIIS would be a good idea
             GLogic=DoODA.AND.cSCF>4.AND.DIISB<75D-5.AND.ETotQ<5D-5.AND.DMaxB<1D-1
             ! If we are increasing with ODA and rebuild is on, we are well fucked.
             FLogic=DoODA.AND..NOT.ALogic.AND.cSCF>3
             ! Sort through logic hopefully in the conditionally correct order ...
             IF(PrintFlags%Key==DEBUG_MAXIMUM) THEN
                CALL OpenASCII(OutFile,Out)
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
             ENDIF
!
             Mssg=" "
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
                Converged(iCLONE)=DID_CONVERGE
                Mssg='Warning: ODA not strictly decreasing'
             ELSEIF(OLogic)THEN
                Converged(iCLONE)=SCF_STALLED
                Mssg='DIIS oscillation'
             ELSEIF(GLogic)THEN
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
      !IF(cSCF>1)ConvergedQ=NOT_CONVERGE
       ConvergedQ=DID_CONVERGE
       DO iCLONE=1,G%Clones
          ConvergedQ=MIN(ConvergedQ,Converged(iCLONE))
       ENDDO
       ! Convergence announcement
       IF(Mssg .NE. " " .AND. cSCF >2)THEN
          CALL OpenASCII(OutFile,Out)
          WRITE(Out,*)TRIM(Mssg)
          WRITE(*,*)TRIM(Mssg)
          CLOSE(Out)
       ENDIF
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
    IF(TRIM(S%Action%C(1))/=SCF_DENSITY_NORMAL   .AND. &
       TRIM(S%Action%C(1))/=SCF_BASISSETSWITCH   .AND. &
       TRIM(S%Action%C(1))/=CPSCF_START_RESPONSE .AND. &
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
  !
  !===============================================================================
  SUBROUTINE DensityLogic(cSCF,cBAS,cGEO,N,S,O,CPSCF_O)
    TYPE(FileNames)    :: N
    TYPE(State)        :: S
    TYPE(Options)      :: O
    INTEGER            :: cSCF,cBAS,cGEO,pBAS,I,J
    LOGICAL            :: DoCPSCF
    LOGICAL,OPTIONAL   :: CPSCF_O
    !----------------------------------------------------------------------------!
    pBAS=S%Previous%I(2)
    S%Current%I=(/cSCF,cBAS,cGEO/)
    IF(PRESENT(CPSCF_O))THEN
       DoCPSCF=CPSCF_O
    ELSE
       DoCPSCF=.FALSE.
    ENDIF
!   Determine the Action to be Taken 
    IF(DoCPSCF)THEN
       IF(O%Guess==GUESS_EQ_DIPOLE.AND.cSCF==0)THEN
          O%Guess=0
          S%Previous%I=S%Current%I
          S%Action%C(1)=CPSCF_START_RESPONSE
       ELSE
          S%Action%C(1)=CPSCF_DENSITY_NORMAL
       ENDIF
    ELSE
       IF(O%Guess==GUESS_EQ_CORE)THEN
          O%Guess=0
          S%Previous%I     = S%Current%I
          S%Action%C(1)    = SCF_GUESSEQCORE
       ELSEIF(O%Guess==GUESS_EQ_SUPR)THEN 
          O%Guess=0
          S%Previous%I     = S%Current%I
          S%Action%C(1)    = SCF_SUPERPOSITION
       ELSEIF(O%Guess==GUESS_EQ_RESTART)THEN
          IF(S%SameBasis .AND. .NOT. S%SameGeom) THEN
             O%Guess=0
             S%Previous%I  = O%RestartState%I
             S%Action%C(1) = SCF_EXTRAPOLATE
          ELSEIF( .NOT. S%SameBasis) THEN
             O%Guess=0
             S%Previous%I  = O%RestartState%I
             S%Action%C(1) = SCF_RWBSS
          ELSE
             O%Guess=0
             S%Previous%I  = O%RestartState%I
             S%Action%C(1) = SCF_RESTART
          ENDIF
       ELSEIF(S%SameBasis .AND. .NOT.S%SameGeom)THEN
          O%Guess=0
          S%Action%C(1)    = SCF_EXTRAPOLATE
       ELSEIF(.NOT. S%SameBasis .OR. pBAS /= cBAS)THEN
          O%Guess=0
          S%Action%C(1)    = SCF_BASISSETSWITCH
          S%Previous%I(1)  = S%Previous%I(1)+1
       ELSE
          O%Guess=0
          S%Action%C(1)    = SCF_DENSITY_NORMAL
       ENDIF
    ENDIF
    S%SameBasis=.TRUE.
    S%SameGeom =.TRUE.
    S%SameCrds =.TRUE.
    S%SameLatt =.TRUE.
!
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
  END SUBROUTINE SolveSCF
!---------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------
  SUBROUTINE SameBasisSameGeom(cBAS,cGEO,N,O,S)
    TYPE(FileNames)    :: N
    TYPE(Options)      :: O
    TYPE(State)        :: S
    TYPE(BSET)         :: BS,BS_rs
    TYPE(CRDS)         :: GM,GM_rs
    REAL(DOUBLE)       :: MaxDiff
    CHARACTER(LEN=DCL) :: chBAS,chGEO    
    INTEGER            :: I,J,cBAS,cGEO,pBAS,pGEO
!
    pBAS=S%Previous%I(2)
    pGEO=S%Previous%I(3)
    S%Current%I=(/0,cBAS,cGEO/)
!
    S%SameCrds  = .TRUE.
    S%SameLatt  = .TRUE.
    S%SameGeom  = .TRUE.
    S%SameBasis = .TRUE.
!

    IF(O%Guess==GUESS_EQ_RESTART) THEN
       HDFFileID=OpenHDF(N%HFile)
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
       chBAS = IntToChar(S%Current%I(2))
       chGEO = IntToChar(S%Current%I(3))
       CALL Get(BS,Tag_O=chBAS)
       CALL Get(GM,Tag_O=chGEO)
       CALL CloseHDFGroup(HDF_CurrentID)
       CALL CloseHDF(HDFFileID)
!
       HDFFileID=OpenHDF(N%RFile)
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
       chBAS = IntToChar(O%RestartState%I(2))
       chGEO = IntToChar(O%RestartState%I(3))
       CALL Get(BS_rs,Tag_O=chBAS)
       CALL Get(GM_rs,Tag_O=chGEO)
       CALL CloseHDFGroup(HDF_CurrentID)
       CALL CloseHDF(HDFFileID)
!
       IF(BS%BName /= BS_rs%BName) S%SameBasis=.FALSE.
       MaxDiff=Zero
       DO I=1,GM%Natms
          MaxDiff=MAX(MaxDiff,ABS(GM%Carts%D(1,I)-GM_rs%Carts%D(1,I)) + &
               ABS(GM%Carts%D(2,I)-GM_rs%Carts%D(2,I)) + &
               ABS(GM%Carts%D(3,I)-GM_rs%Carts%D(3,I))) 
       ENDDO
       IF(MaxDiff>1D-8)     S%SameCrds=.FALSE.
       MaxDiff=Zero
       DO I=1,3;DO J=1,3
          MaxDiff = MAX(MaxDiff,ABS(GM%PBC%BoxShape%D(I,J)-GM_rs%PBC%BoxShape%D(I,J)))
       ENDDO;ENDDO 
       IF(MaxDiff>1D-8)     S%SameLatt=.FALSE.
       IF(.NOT. S%SameCrds) S%SameGeom=.FALSE. 
       IF(.NOT. S%SameLatt) S%SameGeom=.FALSE. 
    ELSE
       HDFFileID=OpenHDF(N%HFile)
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
       chBAS = IntToChar(cBAS)
       chGEO = IntToChar(cGEO)
       CALL Get(GM,Tag_O=chGEO)
       CALL Get(BS,Tag_O=chBAS)
       chBAS = IntToChar(pBAS)
       chGEO = IntToChar(pGEO)
       CALL Get(BS_rs,Tag_O=chBAS)
       CALL Get(GM_rs,Tag_O=chGEO)
       CALL CloseHDFGroup(HDF_CurrentID)
       CALL CloseHDF(HDFFileID)
!
       IF(BS%BName /= BS_rs%BName) S%SameBasis=.FALSE.
       MaxDiff=Zero
       DO I=1,GM%Natms
          MaxDiff=MAX(MaxDiff,ABS(GM%Carts%D(1,I)-GM_rs%Carts%D(1,I)) + &
               ABS(GM%Carts%D(2,I)-GM_rs%Carts%D(2,I)) + &
               ABS(GM%Carts%D(3,I)-GM_rs%Carts%D(3,I))) 
       ENDDO
       IF(MaxDiff>1D-8)            S%SameCrds=.FALSE.
       MaxDiff=Zero
       DO I=1,3;DO J=1,3
          MaxDiff = MAX(MaxDiff,ABS(GM%PBC%BoxShape%D(I,J)-GM_rs%PBC%BoxShape%D(I,J)))
       ENDDO;ENDDO 
       IF(MaxDiff>1D-8)            S%SameLatt=.FALSE.
       IF(.NOT. S%SameCrds)        S%SameGeom=.FALSE. 
       IF(.NOT. S%SameLatt)        S%SameGeom=.FALSE. 
    ENDIF
!
  END SUBROUTINE SameBasisSameGeom
  !===============================================================================
  !
  !===============================================================================
  SUBROUTINE OneEMats(cBAS,cGEO,N,B,S,O,M)
    TYPE(FileNames):: N
    TYPE(BasisSets):: B
    TYPE(State)    :: S
    TYPE(Options)  :: O
    TYPE(Parallel) :: M
    INTEGER        :: cBAS,cGEO,pBAS
    LOGICAL        :: DoPFFT
    !----------------------------------------------------------------------------!
!
    S%Action%C(1)='OneElectronMatrices'
!
    DoPFFT = .FALSE.
    IF(O%Guess==GUESS_EQ_CORE)     DoPFFT=.TRUE.
    IF(O%Guess==GUESS_EQ_SUPR)     DoPFFT=.TRUE.
    IF(O%Guess==GUESS_EQ_RESTART)  DoPFFT=.TRUE.
    IF(.NOT. S%SameLatt)           DoPFFT=.TRUE.
    IF(.NOT. S%SameBasis)          DoPFFT=.TRUE.
!
    IF(DoPFFT) THEN
       CALL Invoke('MakePFFT',N,S,M)
    ENDIF
!
    IF(O%Guess==GUESS_EQ_RESTART .AND.  .NOT. S%SameGeom) THEN
       S%Action%C(1)='RestartGeomSwitch'
       CALL Invoke('MakeS',N,S,M)
       S%Action%C(1)='OneElectronMatrices'
       CALL Invoke('MakeS',N,S,M)
    ELSE
       CALL Invoke('MakeS',N,S,M)
    ENDIF
    IF(O%Methods(cBAS)==RH_R_SCF)THEN
       CALL Invoke('LowdinO',N,S,M)
!       CALL Invoke('IRInv',N,S,M)
    ELSE
       CALL Invoke('AInv',N,S,M)
    ENDIF
    ! Kinetic energy matrix T
    CALL Invoke('MakeT',N,S,M)
    IF(B%BSets(1,cBAS)%HasECPs)THEN
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
    !CALL NTHessian(cBAS,cGEO,N,G,B,S,M)
    IF(B%BSets(1,cBAS)%HasECPs)THEN
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
    !   CALL NXForce(cBAS,cGEO,N,G,B,S,M)
    !  CALL Invoke('XForce',N,S,M)
       CALL Invoke('GONX',N,S,M)
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
!      Get Forces
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
                !If I remember well rms=\sqrt{\frac{\sum_{i=1}^N x_i^2}{N}}
                G%Clone(iCLONE)%GradRMS=G%Clone(iCLONE)%GradRMS+GradVal**2
                !old G%Clone(iCLONE)%GradRMS=G%Clone(iCLONE)%GradRMS+GradVal
                G%Clone(iCLONE)%GradMax=MAX(G%Clone(iCLONE)%GradMax,ABS(GradVal))
             ENDDO
          ELSE
            G%Clone(iCLONE)%Gradients%D(1:3,iATS)=Zero
          ENDIF
       ENDDO
       ! Put the zeroed forces back ...
       ! ... and close the group
       CALL CloseHDFGroup(HDF_CurrentID)
       !
       !If I remember well rms=\sqrt{\frac{\sum_{i=1}^N x_i^2}{N}}
       G%Clone(iCLONE)%GradRMS=SQRT(G%Clone(iCLONE)%GradRMS/DBLE(3*G%Clone(iCLONE)%NAtms))
       !old G%Clone(iCLONE)%GradRMS=SQRT(G%Clone(iCLONE)%GradRMS)/DBLE(3*G%Clone(iCLONE)%NAtms)
       !
    ENDDO
    ! Now close the HDF file ..
    CALL CloseHDF(HDFFileID)
    CALL Delete(S%Action)
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
       !
       ! Load current gradients, we need that, cause we save the geo through
       ! GeoArchive and G%..%Gradients have been det to BIG_DBL previously.
       !
       HDFFileID=OpenHDF(N%HFile)
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
       CALL Get(G%Clone(iCLONE)%PBC%LatFrc,'latfrc',Tag_O=chGEO)
       CALL Get(G%Clone(iCLONE)%Gradients,'Gradients',Tag_O=chGEO)
       CALL CloseHDFGroup(HDF_CurrentID)
       CALL CloseHDF(HDFFileID)
       !
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
    !
    CALL GeomArchive(cBAS,cGEO,N,B,G)
    !
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
!===============================================================================
! Clean the Scratch
!===============================================================================
  SUBROUTINE CleanScratch(C,iGEO)
    TYPE(Controls)                 :: C
    INTEGER                        :: iGEO
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: RemoveFile,chGEO
!    
    chGEO = IntToChar(iGEO)
    RemoveFile=TRIM(C%Nams%M_SCRATCH)//'*_Geom#'//TRIM(chGEO)//"_*.*"         
    CALL SYSTEM('/bin/rm -f  '//RemoveFile)
!
  END SUBROUTINE CleanScratch
!===============================================================================
  SUBROUTINE NTHessian(cBAS,cGEO,N,G,B,S,M)
    TYPE(FileNames)  :: N
    TYPE(Options)    :: O
    TYPE(State)      :: S
    TYPE(Geometries) :: G
    TYPE(BasisSets)  :: B
    TYPE(Parallel)   :: M
    TYPE(DBL_RNK2)   :: HT
    INTEGER          :: cBAS,cGEO,J,iATS,iCLONE
    CHARACTER(LEN=DCL) :: chGEO,chBAS,chSCF
    TYPE(BCSR) :: P
    TYPE(CRDS) :: GTmp
    REAL(DOUBLE),DIMENSION(3,G%Clone(1)%NAtms*3,2)        :: FT
    INTEGER                          :: AtA,AtB,IX,II,IA,IB,A1,A2,IS
    REAL(DOUBLE),PARAMETER           :: DDelta = 1.D-3
    CHARACTER(LEN=DCL)               :: TrixName
    !------------------------------------------------------------------------------
    write(*,*) 'NTHessian NTHessian NTHessian NTHessian NTHessian NTHessian'
    chGEO=IntToChar(cGEO)
    chBAS=IntToChar(cBAS)
    chSCF=IntToChar(S%Current%I(1)+1)
    CALL New(BSiz,G%Clone(1)%NAtms)
    CALL New(OffS,G%Clone(1)%NAtms)
    CALL New(HT,(/G%Clone(1)%NAtms*3,G%Clone(1)%NAtms*3/))
    HT%D=0.0D0
    !
    ! Load current gradients, we need that, cause we save the geo through
    ! GeoArchive and G%..%Gradients have been set to BIG_DBL previously.
    !
    HDFFileID=OpenHDF(N%HFile)
    HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
    CALL Get(G%Clone(1)%PBC%LatFrc,'latfrc',Tag_O=chGEO)
    CALL Get(G%Clone(1)%Gradients,'Gradients',Tag_O=chGEO)
    CALL CloseHDFGroup(HDF_CurrentID)
    CALL CloseHDF(HDFFileID)
    !
    ! Load globals 
    NAToms=G%Clone(1)%NAtms
    MaxAtms=B%MxAts(cBAS)
    MaxBlks=B%MxN0s(cBAS)
    MaxNon0=B%MxBlk(cBAS)
    NBasF=B%BSets(1,cBAS)%NBasF
    BSiz%I=B%BSiz(1,cBAS)%I
    OffS%I=B%OffS(1,cBAS)%I
    MaxBlkSize=0
    DO II=1,NAtoms
       MaxBlkSize=MAX(MaxBlkSize,BSiz%I(II))
    ENDDO
    ! Set temporary geometries
    GTmp%NAtms=G%Clone(1)%NAtms
    CALL New_CRDS(GTmp)
    GTmp%AbCarts%D=G%Clone(1)%AbCarts%D
    GTmp%Gradients%D=G%Clone(1)%Gradients%D
    ! Get the density matrix for this clone
    TrixName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//'_Geom#'//TRIM(chGEO)//'_Base#'//TRIM(chBAS)//'_Cycl#'//TRIM(chSCF) &
         //'_Clone#'//TRIM(IntToChar(1))//'.D'
    CALL Get(P,TrixName)
    chSCF=IntToChar(S%Current%I(1))
    DO AtA=1,NAtoms
       DO IX=1,3
          DO II=1,2
             !
             ! Move the atom.
             IF(II==1) THEN
                G%Clone(1)%AbCarts%D(IX,AtA)=GTmp%AbCarts%D(IX,AtA)+DDelta
             ELSEIF(II==2) THEN
                G%Clone(1)%AbCarts%D(IX,AtA)=GTmp%AbCarts%D(IX,AtA)-DDelta
             ENDIF
             !
             G%Clone(1)%Gradients%D=0.0D0
             CALL GeomArchive(cBAS,cGEO,N,B,G)
             ! vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw>>>
             ! Move back the atom.
             IF(II==1) THEN
                G%Clone(1)%AbCarts%D(IX,AtA)=GTmp%AbCarts%D(IX,AtA)
             ELSEIF(II==2) THEN
                G%Clone(1)%AbCarts%D(IX,AtA)=GTmp%AbCarts%D(IX,AtA)
             ENDIF
             ! vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw<<<
             CALL Invoke('TForce',N,S,M)
             ! Get the T Gradient
             HDFFileID=OpenHDF(N%HFile)
             HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
             CALL Get(G%Clone(1)%Gradients,'Gradients',Tag_O=chGEO)
             CALL CloseHDFGroup(HDF_CurrentID)
             CALL CloseHDF(HDFFileID)
             !copy grad
             FT(:,:,II)=0.0D0
             FT(:,:,II)=G%Clone(1)%Gradients%D(:,:)
             !             DO IA=1,NAtoms
             !                write(*,91919) 'G2(',1,',',IA,')=',G%Clone(1)%Gradients%D(1,IA),';'
             !                write(*,91919) 'G2(',2,',',IA,')=',G%Clone(1)%Gradients%D(2,IA),';'
             !                write(*,91919) 'G2(',3,',',IA,')=',G%Clone(1)%Gradients%D(3,IA),';'
             !91919        format(A,I2,A,I2,A,E26.15,A)
             !             ENDDO
          ENDDO
          IB=0
          IA=3*(AtA-1)+IX
          DO AtB=1,NAtoms
             IB=IB+1
             HT%D(IA,IB)=(FT(1,AtB,1)-FT(1,AtB,2))/(Two*DDelta)
             IB=IB+1
             HT%D(IA,IB)=(FT(2,AtB,1)-FT(2,AtB,2))/(Two*DDelta)
             IB=IB+1
             HT%D(IA,IB)=(FT(3,AtB,1)-FT(3,AtB,2))/(Two*DDelta)
          ENDDO
       ENDDO
    ENDDO
    !
    !
    G%Clone(1)%AbCarts%D=GTmp%AbCarts%D
    G%Clone(1)%Gradients%D=GTmp%Gradients%D
    !
    CALL GeomArchive(cBAS,cGEO,N,B,G)
    !
    CALL Print_DBL_RNK2(HT,'THessian',Unit_O=6)

    CALL Delete(HT)
    CALL Delete(GTmp)
    CALL Delete(P)
    CALL Delete(BSiz)
    CALL Delete(OffS)
  END SUBROUTINE NTHessian
!======================================================================================
END MODULE SCFs




