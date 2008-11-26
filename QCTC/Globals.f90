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
  ! How many Gaussians per leaf node we can have:
  INTEGER,PARAMETER ::   MaxCluster=2048
  ! At what do we stop recurring on the tree:
  INTEGER,PARAMETER ::   MinCluster=32
  ! Logical variable to control "electrons" only option
  ! Useful in case of RPA etc
  LOGICAL           :: NukesOn
END MODULE Globals

