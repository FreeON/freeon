!    FAST O(N Lg N) COMPUTATION OF GRADIENTS OF THE COULOMB ENERGY 
!    WRT TO NUCLEAR COORDINATES
!    Author: Matt Challacombe
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
