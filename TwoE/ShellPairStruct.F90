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
  TYPE AtomInfo
     INTEGER      :: K1,K2,NCell
     REAL(DOUBLE) :: Atm1X,Atm1Y,Atm1Z
     REAL(DOUBLE) :: Atm2X,Atm2Y,Atm2Z
     REAL(DOUBLE) :: R12
     REAL(DOUBLE) :: Atm12X,Atm12Y,Atm12Z
  END TYPE AtomInfo

  TYPE ShellPair
     INTEGER :: IntType
     INTEGER :: L
     TYPE(SmallAtomInfo) :: AtmInfo
     REAL(DOUBLE), DIMENSION(7,3000) :: Cst
  END TYPE ShellPair

  TYPE AtomPr
     TYPE(ShellPair) :: SP
  END TYPE AtomPr

  TYPE ONX2OffSt
     INTEGER :: A,B,C,D
  END TYPE ONX2OffSt
  !
END MODULE ShellPairStruct
