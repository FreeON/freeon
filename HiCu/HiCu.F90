!    COMPUTE THE EXCHANGE CORRELATION MATRIX $K_{xc}$ IN O(N)
!    USING HIERARCHICAL CUBATURE BASED ON ADDAPTIVE CUBATURE 
!    AND K-D BINARY TREE DATA STRUCTURES FOR EFFICIENT RANGE QUERRIES
!
!    Author: Matt Challacombe
!------------------------------------------------------------------------------
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
  IMPLICIT NONE
  TYPE(ARGMT)                 :: Args
  TYPE(BCSR)                  :: Kxc,T1
  TYPE(TIME)                  :: TimeRhoToGrid,TimeGridToMat
  REAL(DOUBLE)                :: Electrons
  CHARACTER(LEN=3)            :: SCFCycle
  CHARACTER(LEN=4),PARAMETER  :: Prog='HiCu'
  CHARACTER(LEN=12),PARAMETER :: Sub1='HiCu.RhoTree' 
  CHARACTER(LEN=12),PARAMETER :: Sub2='HiCu.GridGen' 
  CHARACTER(LEN=12),PARAMETER :: Sub3='HiCu.MakeKxc' 
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg 
!---------------------------------------------------------------------------------------
! Macro the start up
  CALL StartUp(Args,Prog,Serial_O=.TRUE.)
! Get basis set, geometry, thresholds and model type
  CALL Get(BS,CurBase)
  CALL Get(GM,CurGeom)
  CALL Get(ModelChem,'ModelChemistry',CurBase)
  NEl=GM%NElec
#ifdef PERIODIC
! Calculate the Number of Cells
  CALL SetCellNumber(GM)
  CALL PPrint(CS_OUT,'CS_OUT',Prog)
#endif 
! Set local integration thresholds 
  CALL SetLocalThresholds(Thresholds%Cube)
! Begin local performance accumulator for grid generation
  CALL Elapsed_Time(TimeRhoToGrid,'Init')
! Create the RhoTree (a 4-D BinTree)
  CALL RhoToTree(Args)
! Generate the CubeTree (a 3-D BinTree) 
  CALL GridGen()
  CALL Elapsed_TIME(TimeRhoToGrid,'Accum')
  IF(PrintFlags%Key>DEBUG_MEDIUM)THEN
     CALL PPrint(TimeRhoToGrid,Sub2)
     CALL PPrint(TimeRhoToGrid,Sub2,Unit_O=6)
  ENDIF
! Delete the RhoTree
  CALL DeleteRhoTree(RhoRoot)
! Compute the exchange correlation matirix Kxc
  CALL New(Kxc)
  CALL NewBraBlok(BS)
! Begin local performance accumulator for matrix generation
  CALL Elapsed_Time(TimeGridToMat,'Init')
  CALL MakeKxc(Kxc,CubeRoot)
  CALL Elapsed_TIME(TimeGridToMat,'Accum')
  IF(PrintFlags%Key>DEBUG_MEDIUM)THEN
     CALL PPrint(TimeGridToMat,Sub3)
     CALL PPrint(TimeGridToMat,Sub3,Unit_O=6)
  ENDIF
  CALL DeleteBraBlok()
! Put Kxc to disk
  CALL Filter(T1,Kxc)
  CALL Put(Kxc,TrixFile('Kxc',Args,0))
! Put Exc to Info
  CALL Put(Exc,'Exc',Tag_O=SCFCycl)
! Printing
  CALL PChkSum(T1,'Kxc['//TRIM(SCFCycl)//']',Prog)
  CALL PPrint( T1,'Kxc['//TRIM(SCFCycl)//']')
  CALL Plot(   T1,'Kxc['//TRIM(SCFCycl)//']')
! Tidy up
  CALL Delete(Kxc)
  CALL Delete(T1)
  CALL Delete(BS)
  CALL Delete(GM)
! didn't count flops, any accumulation is residual
! from matrix routines
  PerfMon%FLOP=Zero 
! Shutdown 
  CALL ShutDown(Prog)
END PROGRAM HaiKu 
