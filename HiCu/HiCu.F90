!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    COMPUTE THE EXCHANGE CORRELATION MATRIX $K_{xc}$ USING HIERARCHICAL CUBATURE
!
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
  USE RhoTree
  USE CubeTree
  USE Functionals
  USE KxcGen
  IMPLICIT NONE
  TYPE(ARGMT)                 :: Args
  TYPE(CubeNode), POINTER     :: CubeRoot
  TYPE(BCSR)                  :: Kxc,T1
  TYPE(TIME)                  :: TimeRhoToTree,TimeGridGen, &
                                 TimeExcWalk,TimeKxcGen
  REAL(DOUBLE)                :: Electrons
  CHARACTER(LEN=3)            :: SCFCycle
  CHARACTER(LEN=15),PARAMETER :: Prog='HiCu           '
  CHARACTER(LEN=15),PARAMETER :: Sub1='HiCu.RhoTree   ' 
  CHARACTER(LEN=15),PARAMETER :: Sub2='HiCu.GridGen   ' 
  CHARACTER(LEN=15),PARAMETER :: Sub3='HiCu.MakeKxc   ' 
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg 
!---------------------------------------------------------------------------------------
! Macro the start up
  CALL StartUp(Args,Prog,Serial_O=.TRUE.)
! Get basis set, geometry, thresholds and model type
  CALL Get(BS,CurBase)
  CALL Get(GM,CurGeom)
  CALL NewBraKetBlok(BS)
  CALL Get(ModelChem,'ModelChemistry',CurBase)
  NEl=GM%NElec
#ifdef PERIODIC
! Calculate the Number of Cells
  CALL SetCellNumber(GM)
#endif
! Convert density to a 5-D BinTree
  CALL RhoToTree(Args)
! Generate the grid as a 3-D BinTree 
  CALL GridGen(CubeRoot)
! Delete the density
  CALL DeleteRhoTree(RhoRoot)
! Compute the exchange correlation matirix Kxc
  CALL MakeKxc(Kxc,CubeRoot)
! Put Kxc to disk
  CALL Filter(T1,Kxc)
  CALL Put(T1,TrixFile('Kxc',Args,0))
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
! Shutdown 
  CALL ShutDown(Prog)
END PROGRAM HaiKu 
