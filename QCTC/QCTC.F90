!    FAST O(N lg N) COMPUTATION OF THE COULOMB MATRIX
!    Authors:  Matt Challacombe and CJ Tymczak
!===============================================================================
PROGRAM QCTC
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE ProcessControl
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE AtomPairs
  USE BraBloks
  USE QCTCThresholds
  USE PoleTree
  USE Globals
  USE PBCFarField
  USE JGen
  USE NuklarE
  USE ParallelQCTC
#ifdef PARALLEL
  USE FastMatrices
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(FastMat),POINTER          :: J
  REAL(DOUBLE)                   :: TmBegJ,TmEndJ,TmJ
  TYPE(DBL_VECT)                 :: TmJArr
  INTEGER                        :: IErr
  REAL(DOUBLE)                   :: LE_Nuc_Tot
#else
  TYPE(BCSR)                     :: J
#endif
  TYPE(BCSR)                     :: T1,T2

  REAL(DOUBLE)                   :: E_Nuc_Tot
  TYPE(TIME)                     :: TimeMakeJ,TimeMakeTree,TimeNukE
  CHARACTER(LEN=4),PARAMETER     :: Prog='QCTC'
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
  TYPE(CRDS)                     :: GM_MM
  REAL(DOUBLE)                   :: MM_COUL,E_C_EXCL,CONVF  
  INTEGER                        :: I,K,UOUT !!!!
!------------------------------------------------------------------------------- 
  ETimer(:) = Zero
! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
!
  IF(HasQM()) THEN
!    Get basis set and geometry
     CALL Get(BS,Tag_O=CurBase)
     CALL Get(GM,Tag_O=CurGeom)
!    Allocations 
     CALL NewBraBlok(BS)
  ELSEIF(HasMM()) THEN
    CALL Get(GM_MM,Tag_O='GM_MM'//CurGeom)
  ENDIF
  IF(SCFActn=='InkFok')THEN
     CALL Get(Rho,'DeltaRho',Args,0)
     CALL Get(RhoPoles,'Delta')
  ELSE IF(SCFActn=='ForceEvaluation')THEN
     IF(MMOnly()) THEN
       CALL Get(Rho,'Rho',Args,Current(1))
       CALL Get(RhoPoles)
     ELSE
       CALL Get(Rho,'Rho',Args,1)
       CALL Get(RhoPoles)
     ENDIF
  ELSE
     CALL Get(Rho,'Rho',Args,0,Bcast_O=.TRUE.)
     CALL Get(RhoPoles)
  ENDIF
! Set thresholds local to QCTC (for PAC and MAC)
  CALL SetLocalThresholds(Thresholds%TwoE)
! Potentially overide local QCTC thresholds
  IF(Args%NI==8)THEN
     TauPAC=1D1**(-Args%I%I(7))
     TauMAC=1D1**(-Args%I%I(8))
     CALL OpenASCII(OutFile,Out)         
     Mssg=TRIM(ProcessName('QCTC'))//' TauPAC = '//TRIM(DblToShrtChar(TauPAC))
     WRITE(Out,*)TRIM(Mssg)
     Mssg=TRIM(ProcessName('QCTC'))//' TauMAC = '//TRIM(DblToShrtChar(TauMAC))
     WRITE(Out,*)TRIM(Mssg)
     CLOSE(Out)
  ELSE
     CALL OpenASCII(InpFile,Inp)         
     IF(OptDblQ(Inp,'TauPAC',TauPAC))THEN
       Mssg=TRIM(ProcessName('QCTC'))//' TauPAC = '//TRIM(DblToShrtChar(TauPAC))
       CALL OpenASCII(OutFile,Out)         
       WRITE(Out,*)TRIM(Mssg)
       CLOSE(Out)
     ENDIF
     IF(OptDblQ(Inp,'TauMAC',TauMAC))THEN
       Mssg=TRIM(ProcessName('QCTC'))//' TauMAC = '//TRIM(DblToShrtChar(TauMAC))
       CALL OpenASCII(OutFile,Out)         
       WRITE(Out,*)TRIM(Mssg)
       CLOSE(Out)
     ENDIF
     CLOSE(Inp)
  ENDIF
! Initialize the auxiliary density arrays
  CALL InitRhoAux
! Setup global arrays for computation of multipole tensors
  CALL MultipoleSetUp()
! Build the global PoleTree representation of the total density
  CALL Elapsed_Time(TimeMakeTree,'Init')
#ifdef PARALLEL
  CALL ParaRhoToPoleTree
#else
  CALL RhoToPoleTree
#endif

  CALL RhoToPoleTree
  CALL Elapsed_TIME(TimeMakeTree,'Accum')
#ifdef PARALLEL
  CALL EqualTimeSetUp()
#endif
! Set the electrostatic background
  CALL PBCFarFieldSetUp(PoleRoot,GM)
! Delete the auxiliary density arrays
  CALL DeleteRhoAux
! Delete the Density
  CALL Delete(Rho)
!
  IF(HasQM()) THEN
     ! Allocate J
#ifdef PARALLEL
     CALL New_FASTMAT(J,0,(/0,0/))
#else
     CALL New(J)
#endif
     ! Compute the Coulomb matrix J in O(N Lg N)
     CALL Elapsed_Time(TimeMakeJ,'Init')
#ifdef PARALLEL
     TmBegJ = MPI_WTime()
#endif
     CALL MakeJ(J)
#ifdef PARALLEL
     TmEndJ = MPI_WTime()
     TmJ = TmEndJ - TmBegJ
#endif
     CALL Elapsed_TIME(TimeMakeJ,'Accum')
#ifdef PARALLEL
  IF(SCFActn=='InkFok') THEN
    STOP 'InkFok in PARALLEL QCTC is not supported.'
  ENDIF
  CALL Redistribute_FASTMAT(J)
  CALL ET_Part
  CALL Set_BCSR_EQ_DFASTMAT(T1,J) ! T1 is allocated in Set_BCSR...
#else
     IF(SCFActn=='InkFok')THEN
        !    Add in correction if incremental J build
        CALL New(T1)
        CALL New(T2)
        CALL Get(T1,TrixFile('J',Args,-1))
        CALL Add(T1,J,T2)
        CALL Filter(T1,T2)
        CALL Delete(T2)
     ELSE
        CALL Filter(T1,J)
     ENDIF
#endif
     ! Put J to disk
     IF(SCFActn=='FockPrimeBuild'.OR.SCFActn=='StartResponse')THEN
        CALL Put(T1,TrixFile('JPrime'//TRIM(Args%C%C(3)),Args,0))
     ELSE
        CALL Put(T1,TrixFile('J',Args,0))
     ENDIF
     ! Compute the nuclear-total electrostatic energy in O(N Lg N)
     IF(SCFActn=='InkFok')THEN
        CALL Get(E_Nuc_Tot,'E_NuclearTotal')
        CALL Elapsed_Time(TimeNukE,'Init')
        E_Nuc_Tot=E_Nuc_Tot+NukE(GM)
        CALL Elapsed_Time(TimeNukE,'Accum')
     ELSE     
        CALL Elapsed_Time(TimeNukE,'Init')
#ifdef PARALLEL
        CALL NukE_ENPart(GM)
        LE_Nuc_Tot=NukE(GM)
        E_Nuc_Tot = Reduce(LE_Nuc_Tot)
#else
        E_Nuc_Tot=NukE(GM)
#endif
        CALL Elapsed_Time(TimeNukE,'Accum')
     ENDIF

     CALL Put(E_Nuc_Tot,'E_NuclearTotal')
     CALL Put(E_Nuc_Tot,'E_NuclearTotal',StatsToChar(Current))
  ENDIF 
!-------------------------------------------------------------------------------
! QM calculations
  IF(HasMM()) THEN
     CALL Elapsed_Time(TimeNukE,'Init')
     MM_COUL = NukE(GM_MM)
     CALL Elapsed_Time(TimeNukE,'Accum')
     CALL Put(MM_COUL,'MM_COUL')
  ENDIF
!*******
!!$  OPEN(99,FILE='Timing_QCTC.dat',STATUS='UNKNOWN',POSITION='APPEND')
!!$  WRITE(99,7) TRIM(CurGeom),TRIM(CurBase),TRIM(SCFCycl),TimeNukE%CPUS,TimeNukE%WALL
!!$  WRITE(99,8) TRIM(CurGeom),TRIM(CurBase),TRIM(SCFCycl),TimeMakeJ%CPUS,TimeMakeJ%WALL
!!$  WRITE(99,9) TRIM(CurGeom),TRIM(CurBase),TRIM(SCFCycl),TimeMakeTree%CPUS,TimeMakeTree%WALL
!!$7 FORMAT(A2,'  ',A2,'  ',A2,'  QCTC.TimeNukE     = ',F12.4,1X,F12.4)
!!$8 FORMAT(A2,'  ',A2,'  ',A2,'  QCTC.TimeMakeJ    = ',F12.4,1X,F12.4)
!!$9 FORMAT(A2,'  ',A2,'  ',A2,'  QCTC.TimeMakeTree = ',F12.4,1X,F12.4)
!!$  CLOSE(99)
!*******
!-------------------------------------------------------------------------------
! Printing
!  CALL PChkSum(T1,'J['//TRIM(SCFCycl)//']',Prog,Unit_O=6)
  IF(SCFActn=='FockPrimeBuild'.OR.SCFActn=='StartResponse')THEN
     CALL PChkSum(T1,'J'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint( T1,'J'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']')
     CALL Plot(   T1,'J'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']')
  ELSE
     CALL PChkSum(T1,'J['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint( T1,'J['//TRIM(SCFCycl)//']')
     CALL Plot(   T1,'J['//TRIM(SCFCycl)//']')
  ENDIF
! Print Periodic Info
  CALL Print_Periodic(GM,Prog)
#ifdef PARALLEL
  CALL Delete_FastMat1(J)
#else
  CALL Delete(J)
#endif
  CALL Delete(T1)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(Args)
  CALL Delete(RhoPoles)
#ifdef PARALLEL
  CALL New(TmJArr,NPrc)
  CALL MPI_Gather(TmJ,1,MPI_DOUBLE_PRECISION,TmJArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  IF(MyID == ROOT) THEN
! Output needs a lot of work, and should not go to STDOUT
! Also, these statistics were written long ago by      
! Elapsed_TIME(T,Init_O,Proc_O) in PrettyPrint.  Why create
! another routine to do this????
!    CALL PImbalance(TmJArr,NPrc,Prog_O='MakeJ')
  ENDIF
  CALL Delete(TmJArr)
#endif
! didn't count flops, any accumulation is residual from matrix routines
  PerfMon%FLOP=Zero 
! Shutdown 
  CALL ShutDown(Prog)
END PROGRAM QCTC
