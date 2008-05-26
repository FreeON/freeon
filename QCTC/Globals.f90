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
  ! How many Gaussians per leaf node we can have

  INTEGER,PARAMETER ::   ClusterSize=2048
END MODULE Globals

