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
!    FAST O(N Lg N) COMPUTATION OF GRADIENTS OF THE COULOMB ENERGY 
!    WRT TO NUCLEAR COORDINATES
!==============================================================================
PROGRAM JForce
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE Globals
  USE AtomPairs
  USE BraBloks
  USE PoleTree
#ifdef PERIODIC
  USE PBCFarField
#endif
  USE BlokTrPdJ
  USE BlokTrWdS
  USE NuklarE
#ifdef PARALLEL
  USE MondoMPI
  TYPE(DBCSR)                :: P
#else
  TYPE(BCSR)                 :: P
#endif
#ifdef PERIODIC        
  INTEGER                    :: NC
  REAL(DOUBLE)               :: Bx,By,Bz
  TYPE(DBL_VECT)             :: SFrc
#endif   
  INTEGER                    :: AtA,AtB,A1,A2,JP,Q
  TYPE(AtomPair)             :: Pair
  TYPE(DBL_VECT)             :: Frc,JFrc
  REAL(DOUBLE),DIMENSION(3)  :: FrcAB
  CHARACTER(LEN=6),PARAMETER :: Prog='JForce'
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)

  Thresholds%TwoE=1.D-10
  Thresholds%Dist=1.D-14


! Allocations 
  CALL NewBraBlok(BS,Gradients_O=.TRUE.)
  CALL New(Frc,3*NAtoms)
  CALL New(JFrc,3*NAtoms)
  CALL Get(P,TrixFile('D',Args,1))
! Get the Density for Poletree
  CALL Get(Rho,'Rho',Args,0)
! Setup global arrays for computation of multipole tensors
  CALL InitRhoAux
! Setup global arrays for computation of multipole tensors
  CALL MultipoleSetUp(FFEll2)
#ifdef PERIODIC
! Calculate the Number of Cells
  CALL SetCellNumber(GM)
! Set the electrostatic background 
  CALL PBCFarFieldSetUp(FFEll)
#endif
! Build the global PoleTree representation of the total density
  CALL RhoToPoleTree
! Delete the auxiliary density arrays
  CALL DeleteRhoAux
! Delete the Density
  CALL Delete(Rho)
!--------------------------------------------------------------
! Compute the Coulomb contribution to the force in O(N Lg N)
!
  JFrc%D=Zero
  DO AtA=1,NAtoms
     MA=BSiz%I(AtA)
     A1=3*(AtA-1)+1
     A2=3*AtA!
     JFrc%D(A1:A2)=Two*dNukE(AtA)
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
        AtB=P%ColPt%I(JP)
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair))THEN
           Q=P%BlkPt%I(JP)
#ifdef PERIODIC
           Bx = Pair%B(1)
           By = Pair%B(2)           
           Bz = Pair%B(3)
           DO NC=1,CS%NCells
              Pair%B(1) = Bx+CS%CellCarts%D(1,NC)
              Pair%B(2) = By+CS%CellCarts%D(2,NC)
              Pair%B(3) = Bz+CS%CellCarts%D(3,NC)
              Pair%AB2  = (Pair%A(1)-Pair%B(1))**2+(Pair%A(2)-Pair%B(2))**2+(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 JFrc%D(A1:A2)=JFrc%D(A1:A2)+Two*TrPdJ(Pair,P%MTrix%D(Q:))
              ENDIF
           ENDDO
#else
           JFrc%D(A1:A2)=JFrc%D(A1:A2)+Two*TrPdJ(Pair,P%MTrix%D(Q:))
#endif
        ENDIF
     ENDDO
  ENDDO
!---------------------------------------------------------------
! Update forces
  CALL PPrint(JFrc,'JFrce')
  CALL Get(Frc,'GradE',Tag_O=CurGeom)
  Frc%D=Frc%D+JFrc%D
  CALL Put(Frc,'GradE',Tag_O=CurGeom)
  JFrcChk=SQRT(DOT_PRODUCT(JFrc%D,JFrc%D))
!  WRITE(*,*)' JFrcChk = ',JFrcChk
!---------------------------------------------------------------
! Tidy up
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(P)
  CALL Delete(Frc)
  CALL Delete(JFrc)
  CALL DeleteBraBlok(Gradients_O=.TRUE.)
  CALL ShutDown(Prog)
END PROGRAM JForce
