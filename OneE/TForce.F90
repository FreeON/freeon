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
  REAL(DOUBLE),DIMENSION(3) :: B
#endif
  TYPE(AtomPair)      :: Pair
  TYPE(BSET)          :: BS
  TYPE(CRDS)          :: GM
  TYPE(ARGMT)         :: Args
  INTEGER             :: Q,R,AtA,AtB,JP,MB,MA,NB,MN1,A1,A2
  TYPE(HGRho)         :: Rho
  TYPE(DBL_VECT)      :: TFrc,Frc
  REAL(DOUBLE)        :: TFrcChk
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
  CALL New(TFrc,3*NAtoms)
#ifdef PERIODIC
!-----------------------------------------------
! Calculate the Number of Cells
!
  CALL SetCellNumber(GM)
#endif
!--------------------------------------------------------------------------------
! TForce=2*Tr{P.dT} (Extra 2 to account for symmetry of T in the trace)
!--------------------------------------------------------------------------------
  TFrc%D=Zero
  DO AtA=1,NAtoms
     A1=3*(AtA-1)+1
     A2=3*AtA
     MA=BSiz%I(AtA)
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
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2  &
                      +(Pair%A(2)-Pair%B(2))**2  &
                      +(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 TFrc%D(A1:A2)=TFrc%D(A1:A2) &
                              +Four*TrPdT(BS,Pair,P%MTrix%D(Q:Q+MN1))
              ENDIF
           ENDDO
#else
           TFrc%D(A1:A2)=TFrc%D(A1:A2)  &
                        +Four*TrPdT(BS,Pair,P%MTrix%D(Q:Q+MN1))
#endif
        ENDIF
     ENDDO
  ENDDO
!--------------------------------------------------------------------------------
! Do some checksumming, resumming and IO 
!--------------------------------------------------------------------------------
  CALL PChkSum(TFrc,'TForce')  
! for tmp debuging ...
  CALL PPrint(TFrc,'TForce',Unit_O=6)
  CALL PPrint(TFrc,'TForce')
! Sum in contribution to total force
  CALL New(Frc,3*NAtoms)
  CALL Get(Frc,'GradE',Tag_O=CurGeom)
  Frc%D=Frc%D+TFrc%D
  CALL Put(Frc,'GradE',Tag_O=CurGeom)
!------------------------------------------------------------------------------
! Tidy up 
  CALL Delete(P)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(TFrc)
  CALL ShutDown(Prog)
END PROGRAM TForce
