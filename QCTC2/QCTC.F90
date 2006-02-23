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
  USE NukularE
!  USE ParallelQCTC
  USE Clock
  USE TreeWalk
  IMPLICIT NONE
  TYPE(BCSR)                     :: J
  TYPE(BCSR)                     :: T1,T2
  REAL(DOUBLE)                   :: E_Nuc_Tot,SdvErrorJ,MaxErrorJ,JMExact
  REAL(DOUBLE)                   :: QCTC_TotalTime_Start
  TYPE(TIME)                     :: TimeMakeJ,TimeMakeTree,TimeNukE
  CHARACTER(LEN=4),PARAMETER     :: Prog='QCTC'
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
  TYPE(CRDS)                     :: GM_MM
  REAL(DOUBLE)                   :: MM_COUL,E_C_EXCL,CONVF  
  INTEGER                        :: I,K,UOUT !!!!
  !  REAL(DOUBLE),EXTERNAL          :: MondoTimer
  !------------------------------------------------------------------------------- 
  QCTC_TotalTime_Start=MTimer()

  ! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  ! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
  ! Allocations 
  CALL NewBraBlok(BS)
  ! What Action to Take
  IF(SCFActn=='InkFok')THEN
     CALL Get(Rho,'DeltaRho',Args,0)
     CALL Get(RhoPoles,'Delta')
  ELSE IF(SCFActn=='ForceEvaluation')THEN
     IF(MMOnly()) THEN
        CALL Get(Rho,'Rho',Args,Current(1))
        CALL Get(RhoPoles)
     ELSE
        CALL Get(Rho,'Rho',Args,1,Bcast_O=.TRUE.)
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
  ClusterSize=32
  MaxPoleEll=MIN(2*(BS%NASym+4),14)
  WRITE(*,*)' ClusterSize used in QCTC = ',ClusterSize
  WRITE(*,*)' Multipole Expansion used in QCTC = ',MaxPoleEll
  ! Initialize addressing for tensor contraction loops
  CALL TensorIndexingSetUp()
  ! Setup global arrays for computation of multipole tensors ...
  CALL MultipoleSetUp()
  ! Initialize some counters
  ! This preliminary density shit is very sloppy...
  MaxTier=0
  RhoLevel=0
  PoleNodes=0
  ! Initialize the root node   
  CALL NewPoleNode(PoleRoot,0)
  ! Initialize the auxiliary density arrays
  CALL InitRhoAux
  ! Here we set the max multipole expansion relative
  ! to the max angular symmetry
  NInts=0
  NPrim=0
  NFarAv=0
  NNearAv=0
  ! Build the global PoleTree representation of the total density
  CALL RhoToPoleTree
  ! Set the electrostatic background
  CALL PBCFarFieldSetUp(PoleRoot,GM)
  ! Delete the auxiliary density arrays
  CALL DeleteRhoAux
  ! Delete the Density
  CALL Delete(Rho)
  ! Allocate J
  CALL New(J)
  ! Compute the Coulomb matrix J in O(N Lg N)


  CALL Elapsed_Time(TimeMakeJ,'Init')

  JWalk_Time=0D0
  Integral_Time=0D0
  Multipole_Time=0D0

  ALLOCATE(NNearCount(1:CS_IN%NCells))
  NNearCount=0D0
  CALL MakeJ(J)

  K=0
  DO I=1,CS_IN%NCells
     WRITE(*,*)I,NNearCount(I)
     IF(NNearCount(I)==0D0)K=K+1
  ENDDO
  WRITE(*,*)' % of NoPAC = ',DBLE(K)/DBLE(CS_IN%NCells)

  WRITE(*,11)' Decompos_Time = ',Decompose_Time
  WRITE(*,11)' TreeMake_Time = ',TreeMake_Time
  WRITE(*,11)' JWalking_Time = ',JWalk_Time
  WRITE(*,11)' Integral_Time = ',Integral_Time
  WRITE(*,11)' Multipol_Time = ',Multipole_Time
  WRITE(*,11)' Total J Time  = ',Decompose_Time+TreeMake_Time+JWalk_Time+Multipole_Time+Integral_Time
  WRITE(*,11)' Total JWalks  = ',DBLE(NPrim)
  WRITE(*,11)' Av  Ints/Prim = ',DBLE(NInts)/DBLE(NPrim)
  WRITE(*,11)' Av  # NF/Prim = ',DBLE(NNearAv)/DBLE(NPrim)
  WRITE(*,11)' Av  # FF/Prim = ',DBLE(NFarAv)/DBLE(NPrim)
  WRITE(*,11)' Time per INode= ',Integral_Time/DBLE(NNearAv)
  WRITE(*,11)' Time per MNode= ',Multipole_Time/DBLE(NFarAv)
11 FORMAT(A20,D12.6)

  CALL Elapsed_TIME(TimeMakeJ,'Accum')
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
     E_Nuc_Tot=NukE(GM)
     CALL Elapsed_Time(TimeNukE,'Accum')
  ENDIF
  CALL Put(E_Nuc_Tot,'E_NuclearTotal',StatsToChar(Current))
  !-------------------------------------------------------------------------------
  ! Printing
  IF(SCFActn=='FockPrimeBuild'.OR.SCFActn=='StartResponse')THEN
     CALL PChkSum(T1,'J'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint( T1,'J'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']')
     CALL Plot(   T1,'J'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']')
  ELSE
     CALL PChkSum(T1,'J['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint( T1,'J['//TRIM(SCFCycl)//']')!,Unit_O=6)
     CALL Plot(   T1,'J['//TRIM(SCFCycl)//']')
  ENDIF
  ! Print Periodic Info
  CALL Print_Periodic(GM,Prog)
  CALL Delete(J)
  CALL Delete(T1)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(Args)
  CALL Delete(RhoPoles)
  ! didn't count flops, any accumulation is residual from matrix routines
  PerfMon%FLOP=Zero 
  ! Shutdown 


  WRITE(*,11)' QCTC Total Time = ',MTimer()-QCTC_TotalTime_Start

  CALL OpenASCII("Times.dat",111)
  WRITE(111,12)DBLE(NAtoms),DBLE(NInts)/DBLE(NPrim),DBLE(NNearAv)/DBLE(NPrim),DBLE(NFarAv)/DBLE(NPrim), &
               Decompose_Time,TreeMake_Time,JWalk_Time,Integral_Time,Multipole_Time,MTimer()-QCTC_TotalTime_Start
12 FORMAT(100(" ",D12.6))
  CLOSE(111)



  CALL ShutDown(Prog)
END PROGRAM QCTC
