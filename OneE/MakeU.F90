!    COMPUTE THE EFFECTIVE CORE POTENTIAL MATRIX U
!    Author: Matt Challacombe 
!------------------------------------------------------------
PROGRAM MakeU
  USE UBlok
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
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef LJDFLSJFLJSDF

#ifdef PARALLEL
  TYPE(DBCSR)                :: S,T1
#else
  TYPE(BCSR)                 :: S,T1
#endif
#ifdef PERIODIC 
  INTEGER                    :: NC
  REAL(DOUBLE),DIMENSION(3)  :: B
#endif
  TYPE(AtomPair)             :: Pair
  TYPE(BSET)                 :: BS
  TYPE(CRDS)                 :: GM
  TYPE(DBL_RNK4)             :: MD
  TYPE(ARGMT)                :: Args
  INTEGER                    :: P,R,AtA,AtB,NN                         
  CHARACTER(LEN=5),PARAMETER :: Prog='MakeS'
!--------------------------------------- 
! Start up macro

  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
!----------------------------------------------
! Get basis set and geometry

  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
#endif

END PROGRAM MakeU
