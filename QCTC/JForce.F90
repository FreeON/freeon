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
  TYPE(AtomPair)             :: Pair
  TYPE(DBL_VECT)             :: Frc,JFrc
#ifdef PERIODIC        
  REAL(DOUBLE),DIMENSION(3)  :: B
  INTEGER                    :: NC
#endif   
  INTEGER                    :: AtA,AtB,A1,A2,MA,NB,MN1,JP,Q
  REAL(DOUBLE)               :: JFrcChk
  CHARACTER(LEN=6),PARAMETER :: Prog='JForce'
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
! Allocations 
  CALL New(JFrc,3*NAtoms)
  CALL NewBraBlok(BS,Gradients_O=.TRUE.)
  CALL Get(P,TrixFile('D',Args,1))
! Get the Density for Poletree
  CALL Get(Rho,'Rho',Args,0)
! Set thresholds local to JForce (for PAC and MAC)
  CALL SetLocalThresholds(Thresholds%TwoE)
! Setup global arrays for computation of multipole tensors
  CALL InitRhoAux
! Setup global arrays for computation of multipole tensors
  CALL MultipoleSetUp(FFEll2)
! Build the global PoleTree representation of the total density
  CALL RhoToPoleTree
#ifdef PERIODIC
! Calculate the Number of Cells
  CALL SetCellNumber(GM)
! Set the electrostatic background 
  CALL PBCFarFieldSetUp(FFEll,PoleRoot)
#endif
! Delete the auxiliary density arrays
  CALL DeleteRhoAux
! Delete the Density
  CALL Delete(Rho)
!--------------------------------------------------------------------------------
! Compute the Coulomb contribution to the force in O(N Lg N)
!--------------------------------------------------------------------------------
  JFrc%D=Zero
  DO AtA=1,NAtoms
     MA=BSiz%I(AtA)
     A1=3*(AtA-1)+1
     A2=3*AtA
     JFrc%D(A1:A2)=dNukE(AtA)
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
        AtB=P%ColPt%I(JP)
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair))THEN
           Q=P%BlkPt%I(JP)
           NB=BSiz%I(AtB)
           MN1=MA*NB-1
#ifdef PERIODIC
           B=Pair%B
           DO NC=1,CS%NCells
              Pair%B=B+CS%CellCarts%D(:,NC)
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2 &
                      +(Pair%A(2)-Pair%B(2))**2 &
                      +(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair))THEN
                 JFrc%D(A1:A2)=JFrc%D(A1:A2)+TrPdJ(Pair,P%MTrix%D(Q:Q+MN1))
              ENDIF
           ENDDO
#else
           JFrc%D(A1:A2)=JFrc%D(A1:A2)+TrPdJ(Pair,P%MTrix%D(Q:Q+MN1))
#endif
        ENDIF
     ENDDO
  ENDDO
! Closed shell...
  JFrc%D=Two*JFrc%D
!--------------------------------------------------------------------------------
! Do some checksumming, resumming and IO 
!--------------------------------------------------------------------------------
  CALL PChkSum(JFrc,'JForce')  
! for temp debuging....
!  CALL PPrint(JFrc,'JForce',Unit_O=6)
!  CALL PPrint(JFrc,'JForce')
! Sum in contribution to total force
  CALL New(Frc,3*NAtoms)
  CALL Get(Frc,'GradE',Tag_O=CurGeom)
  Frc%D=Frc%D+JFrc%D
  CALL Put(Frc,'GradE',Tag_O=CurGeom)
!--------------------------------------------------------------------------------
! Tidy up
!--------------------------------------------------------------------------------
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(P)
  CALL Delete(Frc)
  CALL Delete(JFrc)
  CALL DeleteBraBlok(Gradients_O=.TRUE.)
  CALL ShutDown(Prog)
END PROGRAM JForce
