MODULE ShellPairStruct
  USE DerivedTypes
  INTEGER,PARAMETER :: PairLngth=5
  INTEGER,PARAMETER :: MaxSPairs=100
  !
  TYPE SmallAtomInfo
     REAL(DOUBLE) :: Atm1X,Atm1Y,Atm1Z
     REAL(DOUBLE) :: Atm2X,Atm2Y,Atm2Z
  END TYPE SmallAtomInfo
  !
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
  TYPE ShellPairG
     INTEGER :: IntType
     INTEGER :: L
     TYPE(SmallAtomInfo) :: AtmInfo
     REAL(DOUBLE), DIMENSION(7,MaxSPairs) :: Cst
  END TYPE ShellPairG
  !
  TYPE AtomInfo
     INTEGER      :: K1,K2,NFPair
     REAL(DOUBLE) :: Atm1X,Atm1Y,Atm1Z
     REAL(DOUBLE) :: Atm2X,Atm2Y,Atm2Z
     REAL(DOUBLE) :: R12
     REAL(DOUBLE) :: Atm12X,Atm12Y,Atm12Z
  END TYPE AtomInfo

  TYPE ShellPair
     INTEGER :: IntType
     INTEGER :: L
     TYPE(SmallAtomInfo) :: AtmInfo
     REAL(DOUBLE), DIMENSION(8,100) :: Cst
  END TYPE ShellPair

  TYPE AtomPr
     TYPE(ShellPair) :: SP
  END TYPE AtomPr

  TYPE AtomPrG
     TYPE(ShellPairG) :: SP
  END TYPE AtomPrG

  TYPE ONX2OffSt
     INTEGER :: A,B,C,D
  END TYPE ONX2OffSt
  !
END MODULE ShellPairStruct
