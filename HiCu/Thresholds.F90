
MODULE HiCuThresholds
   USE Derivedtypes
   USE GlobalScalars   
   USE GlobalObjects
   USE ProcessControl
   USE Indexing
   USE Parse
   USE InOut
   USE Macros
   USE Thresholding
   USE BoundingBox
   IMPLICIT NONE
   REAL(DOUBLE) :: TauRho
   REAL(DOUBLE) :: TauRel
   CONTAINS
      SUBROUTINE SetLocalThresholds(Tau)
         REAL(DOUBLE) :: Tau
         TauRel=Tau*1.D-0  ! Determines integration error
         TauRho=Tau*1.D-2  ! Deterines error of density on the grid
      END SUBROUTINE SetLocalThresholds
END MODULE
