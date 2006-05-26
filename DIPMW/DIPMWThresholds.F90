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
        TauRho  = Tau*1.D-3  ! Reset rho overlap threshold
        TauDIPMW= Tau*1.D0   ! Deterines error of Wavelet representation
      END SUBROUTINE SetLocalThresholdsDIPMW
END MODULE
