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
  TYPE ShellPairG
     INTEGER :: IntType
     INTEGER :: L
     LOGICAL :: Switch
     TYPE(SmallAtomInfo) :: AtmInfo
     REAL(DOUBLE), DIMENSION(10,MaxSPairs) :: Cst
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
     LOGICAL :: Switch
     TYPE(SmallAtomInfo) :: AtmInfo
     REAL(DOUBLE), DIMENSION(10,100) :: Cst
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
