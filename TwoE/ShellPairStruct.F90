MODULE ShellPairStruct
  USE DerivedTypes
  INTEGER,PARAMETER :: PairLngth=5
  INTEGER,PARAMETER :: MaxSPairs=100
  TYPE ShellPair
     INTEGER :: L
     REAL(DOUBLE), DIMENSION(PairLngth,MaxSPairs) :: SP
  END TYPE ShellPair
END MODULE ShellPairStruct
