MODULE ERIGlobals
  USE DerivedTypes
  USE GlobalScalars   
  IMPLICIT NONE
  TYPE(PrimPair)                  :: Prim
  TYPE(BBox)                      :: PBox
  REAL(DOUBLE)                    :: DP2
  REAL(DOUBLE)                    :: PoleSwitch
  INTEGER                         :: At
  REAL(DOUBLE),DIMENSION(1:HGLen) :: HGKet
  REAL(DOUBLE),DIMENSION(0:SPLen) :: SPKetC
  REAL(DOUBLE),DIMENSION(0:SPLen) :: SPKetS
  REAL(DOUBLE),DIMENSION(1000)     :: R
  REAL(DOUBLE),DIMENSION(500)      :: W
  REAL(DOUBLE),DIMENSION(0:50)     :: AuxR,G
CONTAINS
END MODULE ERIGlobals
