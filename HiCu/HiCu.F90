!------------------------------------------------------------------------------
!--  This code is part of the MondoSCF suite of programs for linear scaling 
!    electronic structure theory and ab initio molecular dynamics.
!
!--  Copyright (c) 2001, the Regents of the University of California.  
!    This SOFTWARE has been authored by an employee or employees of the 
!    University of California, operator of the Los Alamos National Laboratory 
!    under Contract No. W-7405-ENG-36 with the U.S. Department of Energy.  
!    The U.S. Government has rights to use, reproduce, and distribute this 
!    SOFTWARE.  The public may copy, distribute, prepare derivative works 
!    and publicly display this SOFTWARE without charge, provided that this 
!    Notice and any statement of authorship are reproduced on all copies.  
!    Neither the Government nor the University makes any warranty, express 
!    or implied, or assumes any liability or responsibility for the use of 
!    this SOFTWARE.  If SOFTWARE is modified to produce derivative works, 
!    such modified SOFTWARE should be clearly marked, so as not to confuse 
!    it with the version available from LANL.  The return of derivative works
!    to the primary author for integration and general release is encouraged. 
!    The first publication realized with the use of MondoSCF shall be
!    considered a joint work.  Publication of the results will appear
!    under the joint authorship of the researchers nominated by their
!    respective institutions. In future publications of work performed
!    with MondoSCF, the use of the software shall be properly acknowledged,
!    e.g. in the form "These calculations have been performed using MondoSCF, 
!    a suite of programs for linear scaling electronic structure theory and
!    ab initio molecular dynamics", and given appropriate citation.  
!------------------------------------------------------------------------------
!    Author: Matt Challacombe
!    COMPUTE THE EXCHANGE CORRELATION MATRIX $K_{xc}$ IN O(N)
!    USING HIERARCHICAL CUBATURE INVOLVING K-D BINARY TREES
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
  CALL Get(ModelChem,'ModelChemistry',CurBase)
  NEl=GM%NElec
#ifdef PERIODIC
! Calculate the Number of Cells
  CALL SetCellNumber(GM)
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
  CALL PPrint(TimeRhoToGrid,Sub2)
  CALL PPrint(TimeRhoToGrid,Sub2,Unit_O=6)
! Delete the RhoTree
  CALL DeleteRhoTree(RhoRoot)
! Compute the exchange correlation matirix Kxc
  CALL New(Kxc)
  CALL NewBraBlok(BS)
! Begin local performance accumulator for matrix generation
  CALL Elapsed_Time(TimeGridToMat,'Init')
  CALL MakeKxc(Kxc,CubeRoot)
  CALL Elapsed_TIME(TimeGridToMat,'Accum')
  CALL PPrint(TimeGridToMat,Sub3)
  CALL PPrint(TimeGridToMat,Sub3,Unit_O=6)
  CALL DeleteBraBlok()
! Put Kxc to disk
  CALL Filter(T1,Kxc)
  CALL Put(T1,TrixFile('Kxc',Args,0))
! Put Exc to Info
  CALL Put(Exc,'Exc',Tag_O=SCFCycl)
! Printing
!  CALL PChkSum(T1,'Kxc['//TRIM(SCFCycl)//']',Prog,Unit_O=6)
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
