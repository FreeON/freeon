MODULE DIPMWThresholds
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
   REAL(DOUBLE) :: TauDIPMW
!  Global Objects  
   TYPE(BSET)     :: BS  !  Global basis set
   TYPE(CRDS)     :: GM  !  Global molecular geometry
   CONTAINS
      SUBROUTINE SetLocalThresholdsDIPMW(Tau)
        REAL(DOUBLE) :: Tau
        TauRho  = 1.D-14     ! Reset rho overlap threshold
        TauDIPMW= Tau*1.D0   ! Deterines error of Wavelet representation
      END SUBROUTINE SetLocalThresholdsDIPMW
END MODULE
