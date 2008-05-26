!    Local QCTC thresholds for tuning MAC and PAC
!    Author: Matt Challacombe
!------------------------------------------------------------------------------
MODULE QCTCThresholds
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
   USE PoleNodeType
   IMPLICIT NONE
   REAL(DOUBLE)            :: TauPAC
   REAL(DOUBLE)            :: TauMAC
   REAL(DOUBLE)            :: TauTwo
   CONTAINS
      SUBROUTINE SetLocalThresholds(Tau)
         REAL(DOUBLE) :: Tau
!        Penetration Acceptability Criterion (PAC) threshold
         TauPAC    = Tau!*1D2
!        Multipole Acceptability Criterion (MAC) threshold
         TauMAC    = Tau!*1D-1
         TauTwo    = Tau*1D-5
      END SUBROUTINE SetLocalThresholds
END MODULE
