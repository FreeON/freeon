MODULE Globals
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalObjects
!
  TYPE(BSET)      :: BS
  TYPE(CRDS)      :: GM
  TYPE(DBL_RNK4)  :: MD
  TYPE(ARGMT)     :: Args
  TYPE(HGRho)     :: Rho
END MODULE Globals
