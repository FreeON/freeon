! COMPUTE THE EXCHANGE CORRELATION MATRIX $K_{xc}$ IN O(N)
! USING ADAPTIVE, HIERARCHICAL CUBATURE  AND K-D BINARY TREE DATA STRUCTURES 
! FOR EFFICIENT RANGE QUERRIES OF THE DENSITY AND GRID
! Author: Matt Challacombe
!-------------------------------------------------------------------------------
PROGRAM HaiKu
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
  USE Thresholding
  USE HiCuThresholds
  USE RhoTree
  USE CubeTree
  USE Functionals
  USE KxcGen
  USE BoundingBox
#ifdef PARALLEL 
  USE ParallelHiCu
  USE FastMatrices
#endif
  IMPLICIT NONE
  TYPE(ARGMT)                    :: Args
#ifdef PARALLEL
  INTEGER::I
  REAL(DOUBLE)                   :: StartTm,EndTm,TotTm,CumuTm
#endif
#ifdef PARALLEL
  TYPE(FastMat),POINTER          :: Kxc
#else
  TYPE(BCSR)                     :: Kxc
#endif
  TYPE(BCSR)                     :: T1
  TYPE(TIME)                     :: TimeRhoToGrid,TimeGridToMat
  REAL(DOUBLE)                   :: Electrons
  CHARACTER(LEN=3)               :: SCFCycle
  CHARACTER(LEN=4),PARAMETER     :: Prog='HiCu'
  CHARACTER(LEN=12),PARAMETER    :: Sub1='HiCu.RhoTree' 
  CHARACTER(LEN=12),PARAMETER    :: Sub2='HiCu.GridGen' 
  CHARACTER(LEN=12),PARAMETER    :: Sub3='HiCu.MakeKxc' 
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg 
  TYPE(BBox)                     ::WBox
  REAL(DOUBLE)                   ::VolRho,VolExc
!-------------------------------------------------------------------------------
! Macro the start up
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
! Get basis set, geometry, thresholds and model type
  CALL Get(BS,CurBase)
  CALL Get(GM,CurGeom)
  NEl=GM%NElec
! Set local integration thresholds 
  CALL SetLocalThresholds(Thresholds%Cube)
#ifdef NewPAC
  CALL SetAACoef()
#endif
  ! Potentially overide local HiCu thresholds
  IF(Args%NI==8)THEN
     TauRel=1D1**(-Args%I%I(7))
     TauRho=1D1**(-Args%I%I(8))
#ifdef PARALLEL
     IF(MyId==0)THEN
#endif
     CALL OpenASCII(OutFile,Out)         
     Mssg=TRIM(ProcessName('HiCu'))//' TauRel = '//TRIM(DblToShrtChar(TauRel))
     WRITE(Out,*)TRIM(Mssg)
     Mssg=TRIM(ProcessName('HiCu'))//' TauRho = '//TRIM(DblToShrtChar(TauRho))
     WRITE(Out,*)TRIM(Mssg)
     CLOSE(Out)
#ifdef PARALLEL
     ENDIF
#endif
  ELSE
#ifdef PARALLEL
     IF(MyId==0)THEN
#endif
        CALL OpenASCII(InpFile,Inp)         
        IF(OptDblQ(Inp,'TauRel',TauRel))THEN
           Mssg=TRIM(ProcessName('HiCu'))//' TauRel = '//TRIM(DblToShrtChar(TauRel))
           CALL OpenASCII(OutFile,Out)         
           WRITE(Out,*)TRIM(Mssg)
           CLOSE(Out)
        ENDIF
        IF(OptDblQ(Inp,'TauRho',TauRho))THEN
           Mssg=TRIM(ProcessName('HiCu'))//' TauRho = '//TRIM(DblToShrtChar(TauRho))
           CALL OpenASCII(OutFile,Out)         
           WRITE(Out,*)TRIM(Mssg)
           CLOSE(Out)
        ENDIF
        CLOSE(Inp)
#ifdef PARALLEL
     ENDIF
#endif
  ENDIF
! Begin local performance accumulator for grid generation
  CALL Elapsed_Time(TimeRhoToGrid,'Init')
! Create the RhoTree (a 4-D BinTree)
#ifdef PARALLEL
  !! read in distributions, figure out the root bounding box
  CALL ParaInitRho(Args)
  CALL GetBBox()
  CALL AlignNodes()
  CALL SendBBox()
  CALL DistDist()
  CALL ParaRhoToTree()
#else
  CALL RhoToTree(Args)
#endif
! CALL New(Kxc)
#ifdef PARALLEL
  CALL New_FASTMAT(Kxc,0,(/0,0/))
#else
  CALL NewBCSR(Kxc)
#endif
  CALL NewBraBlok(BS)
#ifdef PARALLEL 
  CALL WorkBBOx(Kxc)
  CALL RepartitionVol()
#else
! Generate the CubeTree (a 3-D BinTree) 
  WBox%BndBox(1:3,1:2) = RhoRoot%Box%BndBox(1:3,1:2)
  CALL MakeBoxPeriodic(WBox)
  CALL CalCenterAndHalf(WBox)
  CALL GridGen(WBox,VolRho,VolExc)
  CALL Elapsed_TIME(TimeRhoToGrid,'Accum')
  IF(PrintFlags%Key>DEBUG_MEDIUM)THEN
     CALL PPrint(TimeRhoToGrid,Sub2)
     CALL PPrint(TimeRhoToGrid,Sub2,Unit_O=6)
  ENDIF
! Compute the exchange correlation matirix Kxc
! Begin local performance accumulator for matrix generation
  CALL Elapsed_Time(TimeGridToMat,'Init')
  CALL MakeKxc(Kxc,CubeRoot)
  CALL Elapsed_TIME(TimeGridToMat,'Accum')
  IF(PrintFlags%Key>DEBUG_MEDIUM)THEN
     CALL PPrint(TimeGridToMat,Sub3)
     CALL PPrint(TimeGridToMat,Sub3,Unit_O=6)
  ENDIF
#endif
! Delete the RhoTree
  CALL DeleteRhoTree(RhoRoot)
  CALL DeleteBraBlok()
! Put Kxc to disk
#ifdef PARALLEL
  CALL Redistribute_FASTMAT(Kxc)
  CALL AlignNodes()
  CALL Set_BCSR_EQ_DFASTMAT(T1,Kxc)
#else
  CALL Filter(T1,Kxc)
#endif

#ifdef PARALLEL
  CALL Put(T1,TrixFile('Kxc',Args,0))
#else
  CALL Put(Kxc,TrixFile('Kxc',Args,0))
#endif
! Put Exc to Info
#ifdef PARALLEL
    CALL Put(TotExc,'Exc',StatsToChar(Current))
#else
    CALL Put(Exc,'Exc',StatsToChar(Current))
#endif
! Printing
  CALL PChkSum(T1,'Kxc['//TRIM(SCFCycl)//']',Prog)
  CALL PPrint( T1,'Kxc['//TRIM(SCFCycl)//']')
  CALL Plot(   T1,'Kxc['//TRIM(SCFCycl)//']')

#ifdef PARALLEL
  CALL Delete_FastMat1(Kxc)
#else
  CALL Delete(Kxc)
#endif
!*******
!!$  OPEN(99,FILE='Timing_HiCu.dat',STATUS='UNKNOWN',POSITION='APPEND')
!!$  WRITE(99,7) TRIM(CurGeom),TRIM(CurBase),TRIM(SCFCycl),TimeRhoToGrid%CPUS,TimeRhoToGrid%WALL
!!$  WRITE(99,8) TRIM(CurGeom),TRIM(CurBase),TRIM(SCFCycl),TimeGridToMat%CPUS,TimeGridToMat%WALL
!!$  WRITE(99,*)
!!$7 FORMAT(A2,'  ',A2,'  ',A2,'  HiCu.RhoToGrid   = ',F12.4,1X,F12.4)
!!$8 FORMAT(A2,'  ',A2,'  ',A2,'  HiCu.GridToMat   = ',F12.4,1X,F12.4)
!!$  CLOSE(99)
!*******
  CALL Delete(T1)
  CALL Delete(BS)
  CALL Delete(GM)
! didn't count flops, any accumulation is residual from matrix routines
  PerfMon%FLOP=Zero 
  CALL ShutDown(Prog)
END PROGRAM HaiKu 
