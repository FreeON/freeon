MODULE ONX2DataType
  !
  USE DerivedTypes
  USE ShellPairStruct
  !
  TYPE RNode
     INTEGER              :: iCell
     REAL(DOUBLE)         :: MaxInt
     TYPE(RNode), POINTER :: RNext
  END TYPE RNode
  !
  TYPE ANode
     INTEGER              :: Atom
     REAL(DOUBLE)         :: MaxInt
     TYPE(ANode), POINTER :: AtmNext
     TYPE(RNode), POINTER :: RNext
  END TYPE ANode
  !
  TYPE ANode2
     INTEGER :: Atom
     INTEGER :: NCell
     INTEGER     , DIMENSION(:), POINTER :: CellIdx
     REAL(DOUBLE), DIMENSION(:), POINTER :: SqrtInt
     TYPE(ANode2), POINTER :: AtmNext
  END TYPE ANode2
  !
  TYPE CList2
     TYPE(ANode2), POINTER :: GoList
  END TYPE CList2
  !
  TYPE CList
     TYPE(ANode), POINTER :: GoList
  END TYPE CList
  !
  TYPE ONX2OffSt
     INTEGER :: A,B,C,D
  END TYPE ONX2OffSt
  !
!!$  TYPE SmallAtomInfo
!!$     REAL(DOUBLE) :: Atm1X,Atm1Y,Atm1Z
!!$     REAL(DOUBLE) :: Atm2X,Atm2Y,Atm2Z
!!$  END TYPE SmallAtomInfo
!!$  !
!!$  TYPE ShellPair
!!$     INTEGER :: IntType
!!$     INTEGER :: L
!!$     TYPE(SmallAtomInfo) :: AtmInfo
!!$     REAL(DOUBLE), DIMENSION(5,3000) :: Cst
!!$     !REAL(DOUBLE), DIMENSION(:,:), POINTER :: Cst
!!$     !REAL(DOUBLE), DIMENSION(:,:), ALLOCATABLE :: Cst
!!$     !REAL(DOUBLE), DIMENSION(...) :: Cst
!!$  END TYPE ShellPair
  !
  TYPE AtomPr
     TYPE(ShellPair) :: SP
  END TYPE AtomPr
  !
!!$  TYPE AtomInfo
!!$     INTEGER      :: K1,K2,NCell
!!$     REAL(DOUBLE) :: Atm1X,Atm1Y,Atm1Z
!!$     REAL(DOUBLE) :: Atm2X,Atm2Y,Atm2Z
!!$     REAL(DOUBLE) :: R12
!!$     REAL(DOUBLE) :: Atm12X,Atm12Y,Atm12Z
!!$  END TYPE AtomInfo
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if 0
  TYPE ShellPair
     INTEGER :: N1
     INTEGER :: N2
     REAL(DOUBLE) :: At1x,At1y,At1z
     REAL(DOUBLE) :: At2x,At2y,At2z
     REAL(DOUBLE) :: R12
     REAL(DOUBLE), DIMENSION(10) :: Expt1
     REAL(DOUBLE), DIMENSION(10) :: Expt2
     REAL(DOUBLE), DIMENSION(10) :: Coef1
     REAL(DOUBLE), DIMENSION(10) :: Coef2
  END TYPE ShellPair
#endif

END MODULE ONX2DataType



