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
!    COMPUTE DERIVATIVES OF THE EXCHANGE CORRELATION ENERGY $E_{xc}$ 
!    WRT TO NUCLEAR COORDINATES IN O(N) USING HIERARCHICAL CUBATURE
!------------------------------------------------------------------------------
PROGRAM XCForce
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
  USE dXCBlok
  IMPLICIT NONE
  TYPE(ARGMT)                 :: Args
  TYPE(BCSR)                  :: P
  TYPE(TIME)                  :: TimeRhoToTree,TimeGridGen
  TYPE(DBL_VECT)              :: XCFrc,Frc
  TYPE(AtomPair)              :: Pair
  REAL(DOUBLE)                :: Electrons,XCFrcChk
  INTEGER                     :: AtA,AtB,MA,A1,A2,JP,Q
  CHARACTER(LEN=3)            :: SCFCycle
  CHARACTER(LEN=15),PARAMETER :: Prog='XCForce        '
  CHARACTER(LEN=15),PARAMETER :: Sub1='XCForce.RhoTree' 
  CHARACTER(LEN=15),PARAMETER :: Sub2='XCForce.GridGen' 
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg 
#ifdef PERIODIC        
      INTEGER                  :: NCA,NCB
      REAL(DOUBLE)             :: Ax,Ay,Az,Bx,By,Bz
#endif     
!---------------------------------------------------------------------------------------
! Macro the start up
  CALL StartUp(Args,Prog,Serial_O=.TRUE.)
! Get basis set, geometry, thresholds and model type
  CALL Get(BS,CurBase)
  CALL Get(GM,CurGeom)
  CALL Get(ModelChem,'ModelChemistry',CurBase)
  NEl=GM%NElec
! For some reason need higher accuracy for 
! gradients.  This seems unessesary, needs fixing...
  Thresholds%Cube=Thresholds%Cube*1.D-1
#ifdef PERIODIC
! Calculate the Number of Cells
  CALL SetCellNumber(GM)
#endif
! Convert density to a 5-D BinTree
  CALL RhoToTree(Args)
! Generate the grid as a 3-D BinTree 
  CALL GridGen()
! Delete the density
  CALL DeleteRhoTree(RhoRoot)
! More allocations 
  CALL NewBraBlok(BS,Gradients_O=.TRUE.)
  CALL New(Frc,3*NAtoms)
  CALL New(XCFrc,3*NAtoms)
  CALL Get(P,TrixFile('D',Args,1))
!----------------------------------------------------------------------
! Compute the exchange-correlation contribution to the force in O(N)
!
  XCFrc%D=Zero
  DO AtA=1,NAtoms
     MA=BSiz%I(AtA)
     A1=3*(AtA-1)+1
     A2=3*AtA
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
        AtB=P%ColPt%I(JP)
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair))THEN
           Q=P%BlkPt%I(JP)
#ifdef PERIODIC
           Ax = Pair%A(1)
           Ay = Pair%A(2)
           Az = Pair%A(3)
           Bx = Pair%B(1)
           By = Pair%B(2)           
           Bz = Pair%B(3)
           DO NCA = 1,CS%NCells
              Pair%A(1) = Ax+CS%CellCarts%D(1,NCA)
              Pair%A(2) = Ay+CS%CellCarts%D(2,NCA)
              Pair%A(3) = Az+CS%CellCarts%D(3,NCA)
              DO NCB = 1,CS%NCells
                 Pair%B(1) = Bx+CS%CellCarts%D(1,NCB)
                 Pair%B(2) = By+CS%CellCarts%D(2,NCB)
                 Pair%B(3) = Bz+CS%CellCarts%D(3,NCB)
                 Pair%AB2  = (Pair%A(1)-Pair%B(1))**2+(Pair%A(2)-Pair%B(2))**2+(Pair%A(3)-Pair%B(3))**2
                 IF(TestAtomPair(Pair)) THEN
                    XCFrc%D(A1:A2)=XCFrc%D(A1:A2)+dXC(Pair,P%MTrix%D(Q:))
                 ENDIF
              ENDDO
           ENDDO
#else
           XCFrc%D(A1:A2)=XCFrc%D(A1:A2)+dXC(Pair,P%MTrix%D(Q:))
#endif
        ENDIF
     ENDDO
  ENDDO
!---------------------------------------------------------------
! Update forces
  CALL PPrint(XCFrc,'XCFrce')
  CALL Get(Frc,'GradE',Tag_O=CurGeom)
  Frc%D=Frc%D+XCFrc%D
  CALL Put(Frc,'GradE',Tag_O=CurGeom)
!  XCFrcChk=SQRT(DOT_PRODUCT(XCFrc%D,XCFrc%D))
!  WRITE(*,*)' XCFrcChk = ',XCFrcChk
!---------------------------------------------------------------
! Tidy up
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(P)
  CALL Delete(Frc)
  CALL Delete(XCFrc)
  CALL DeleteBraBlok(Gradients_O=.TRUE.)
! Shutdown 
  CALL ShutDown(Prog)
END PROGRAM XCForce
