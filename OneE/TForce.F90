!    COMPUTE THE FORCE CORESPONDING TO THE DERIVATIVE OF THE 
!    KINETIC ENERGY MATRIX, TForce=2*Tr{P.dT}
!    Authors: Matt Challacombe and CJ Tymczak
!----------------------------------------------------------------------------------
PROGRAM TForce
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE AtomPairs
  USE BlokTrPdT
#ifdef PARALLEL
  USE MondoMPI
  TYPE(DBCSR)         :: P
#else
  TYPE(BCSR)          :: P
#endif
#ifdef PERIODIC 
  INTEGER             :: NC
  REAL(DOUBLE)        :: Bx,By,Bz
#endif
  TYPE(AtomPair)      :: Pair
  TYPE(BSET)          :: BS
  TYPE(CRDS)          :: GM
  TYPE(ARGMT)         :: Args
  INTEGER             :: Q,R,AtA,AtB,NN,iSwitch,IStrtP,IStopP,LP,JP,MB,MA,A1,A2
  TYPE(HGRho)         :: Rho
  TYPE(DBL_VECT)      :: TFrc,Frc
  CHARACTER(LEN=6),PARAMETER :: Prog='TForce'
!------------------------------------------------------------------------------------- 
! Start up macro
!
  CALL StartUp(Args,Prog)
!----------------------------------------------
! Get basis set and geometry
!
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
  CALL Get(P,TrixFile('D',Args,0))
  CALL New(Frc,3*NAtoms)
  CALL New(TFrc,3*NAtoms)
#ifdef PERIODIC
!-----------------------------------------------
! Calculate the Number of Cells
!
  CALL SetCellNumber(GM)
#endif
!---------------------------------------------------------------------------
! TForce=2*Tr{P.dT} (Extra 2 to account for symmetry of T in the trace)
  TFrc%D=Zero
  DO AtA=1,NAtoms
     A1=3*(AtA-1)+1
     A2=3*AtA
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
        AtB=P%ColPt%I(JP)
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair))THEN
           Q=P%BlkPt%I(JP)
#ifdef PERIODIC
           Bx = Pair%B(1)
           By = Pair%B(2)           
           Bz = Pair%B(3)
           DO NC=1,CS%NCells
              Pair%B(1)=Bx+CS%CellCarts%D(1,NC)
              Pair%B(2)=By+CS%CellCarts%D(2,NC)
              Pair%B(3)=Bz+CS%CellCarts%D(3,NC)
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2+(Pair%A(2)-Pair%B(2))**2+(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 TFrc%D(A1:A2)=TFrc%D(A1:A2)+Four*TrPdT(BS,Pair,P%MTrix%D(Q:))
              ENDIF
           ENDDO
#else
           TFrc%D(A1:A2)=TFrc%D(A1:A2)+Four*TrPdT(BS,Pair,P%MTrix%D(Q:))
#endif
        ENDIF
     ENDDO
  ENDDO
!---------------------------------------------------------------
! Update forces
  CALL PPrint(TFrc,'TForce')
  CALL Get(Frc,'GradE',Tag_O=CurGeom)
  Frc%D=Frc%D+TFrc%D
  CALL Put(Frc,'GradE',Tag_O=CurGeom)
  TFrcChk=SQRT(DOT_PRODUCT(TFrc%D,TFrc%D))
!------------------------------------------------------------------------------
! Tidy up 
  CALL Delete(P)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(TFrc)
  CALL ShutDown(Prog)
END PROGRAM TForce
