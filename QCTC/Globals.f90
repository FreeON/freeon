MODULE Globals
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalObjects
!
  TYPE(BSET)        :: BS
  TYPE(CRDS)        :: GM
  TYPE(DBL_RNK4)    :: MD
  TYPE(ARGMT)       :: Args
  INTEGER,PARAMETER :: TimerSize = 40
  REAL(DOUBLE)      :: ETimer(TimerSize)
!
END MODULE Globals
