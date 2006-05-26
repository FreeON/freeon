! COMPUTE THE EXCHANGE CORRELATION MATRIX $K_{xc}$ IN O(N)
! USING ADAPTIVE, HIERARCHICAL CUBATURE  AND K-D BINARY TREE DATA STRUCTURES 
! FOR EFFICIENT RANGE QUERRIES OF THE DENSITY AND GRID
! Author: Matt Challacombe
!-------------------------------------------------------------------------------
PROGRAM DIPMWavelet
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
  USE RhoTree
  USE Functionals
  USE BoundingBox 
  USE DIPMWThresholds
  USE DIPMWTree
  USE KxcGenDIPMW
!
  IMPLICIT NONE
  INTEGER                        :: I,J,K
  TYPE(ARGMT)                    :: Args
  TYPE(BCSR)                     :: Kxc
  TYPE(BCSR)                     :: T1
  TYPE(TIME)                     :: TimeRhoToTree,TimeDIPMWTree,TimeMakeKxc
  REAL(DOUBLE)                   :: Electrons,TmpDbl,Error
  CHARACTER(LEN=3)               :: SCFCycle
  CHARACTER(LEN=4),PARAMETER     :: Prog='DIPMW'
  CHARACTER(LEN=12),PARAMETER    :: Sub1='DIPMW.RhoTree' 
  CHARACTER(LEN=12),PARAMETER    :: Sub2='DIPMW.Tree' 
  CHARACTER(LEN=12),PARAMETER    :: Sub3='DIPMW.MakeKxc' 
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg 
!-------------------------------------------------------------------------------
! Macro the start up
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
! Get basis set, geometry, thresholds and model type
  CALL Get(BS,CurBase)
  CALL Get(GM,CurGeom)
! Set local integration thresholds 
  CALL SetLocalThresholdsDIPMW(Thresholds%Cube)
  TauDIPMW = 1.0D-3
  TauRho   = TauDIPMW*1.D-3
  WRITE(*,*) 'TauDIPMW = ',TauDIPMW
  WRITE(*,*) 'TauRho   = ',TauRho
  CALL SetAACoef()
! Begin local performance accumulator for grid generation
  CALL Elapsed_Time(TimeRhoToTree,'Init')
! Create the RhoTree (a 4-D BinTree)
  CALL RhoToTree(Args) 
  CALL Elapsed_Time(TimeRhoToTree,'Accum')
  WRITE(*,*) "RhoToTree = ",TimeRhoToTree%CPUS,TimeRhoToTree%WALL
! Genergate the Wavelet Representation of the XC potential
  CALL Elapsed_Time(TimeDIPMWTree,'Init')
  CALL DIPMWTree(20)
  CALL Elapsed_Time(TimeDIPMWTree,'Accum')
  WRITE(*,*) "RhoToPIGW = ",TimeDIPMWTree%CPUS,TimeDIPMWTree%WALL
! Allocate Kxc
  CALL NewBraBlok(BS)
  CALL New(Kxc)
! Make the Exchange Corelation Matrix
  CALL Elapsed_Time(TimeMakeKxc,'Init')
  CALL MakeKxcDIPMW(Kxc,PIGWRoot)
  CALL Elapsed_Time(TimeMakeKxc,'Accum')
  WRITE(*,*) " TimeMakeKxc = ",TimeMakeKxc%CPUS,TimeMakeKxc%WALL
!
  CALL PPrint(Kxc,'Kxc',Unit_O=6)
! Delete the RhoTree,PIGWRoot and BraBloks
  CALL DeleteRhoTree(RhoRoot)
  CALL DeleteDIPMWTree(PIGWRoot)
  CALL DeleteBraBlok()
! Put Kxc to disk
  CALL Filter(T1,Kxc)
  CALL Put(Kxc,TrixFile('Kxc',Args,0))
! Put Exc to Info
  CALL Put(XCEnergy,'Exc')
  CALL Put(XCEnergy,'Exc',StatsToChar(Current))
! Printing
  CALL PChkSum(T1,'Kxc['//TRIM(SCFCycl)//']',Prog)
  CALL PPrint( T1,'Kxc['//TRIM(SCFCycl)//']')
  CALL Plot(   T1,'Kxc['//TRIM(SCFCycl)//']')
  CALL Delete(Kxc)
  CALL Delete(T1)
  CALL Delete(BS)
  CALL Delete(GM)
! didn't count flops, any accumulation is residual from matrix routines
  PerfMon%FLOP=Zero 
  CALL ShutDown(Prog)
!
END PROGRAM DIPMWavelet
