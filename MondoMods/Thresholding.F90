!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe and  C. J. Tymczak
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    Compute and set intermediate thresholds
!
MODULE Thresholding
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
!  USE BoundingBox
!-------------------------------------------------  
!  Primary thresholds
!
   TYPE(TOLS), SAVE :: Thresholds
!-------------------------------------------------------------------------
! Intermediate thresholds
!
  REAL(DOUBLE), SAVE  :: MinExponent
  REAL(DOUBLE), SAVE  :: AtomPairDistanceThreshold  ! Atom pairs
  REAL(DOUBLE), SAVE  :: PrimPairDistanceThreshold  ! Prim pairs
  REAL(DOUBLE), SAVE  :: PenetratDistanceThreshold  ! Penetration threshold
  CONTAINS
!------------------------------------------------------------------------------
!    Set the Atom Pair Distance Threshhold: Zeta_Min*|A-B|^2 > -Log(Tau)
!
     SUBROUTINE SetAtomPairThresh(Tau)
        REAL(DOUBLE),INTENT(IN) :: Tau
        AtomPairDistanceThreshold=-LOG(Tau)/MinExponent
     END SUBROUTINE SetAtomPairThresh
!
     FUNCTION TestAtomPair(Pair)
        LOGICAL                   :: TestAtomPair
        TYPE(AtomPair)            :: Pair
        IF(Pair%AB2>AtomPairDistanceThreshold) THEN
           TestAtomPair = .FALSE.
        ELSE
           TestAtomPair = .TRUE.
        ENDIF
     END FUNCTION TestAtomPair
!--------------------------------------------------------------------------
!    Set the Primitive Pair Distance Threshhold: Xi_ab*Min*|A-B|^2 > -Log(Tau)
!
     SUBROUTINE SetPrimPairThresh(Tau)
        REAL(DOUBLE),INTENT(IN) :: Tau
        PrimPairDistanceThreshold=-LOG(Tau)
     END SUBROUTINE SetPrimPairThresh
!
     FUNCTION TestPrimPair(Xi,Dist)
        LOGICAL                :: TestPrimPair
        TYPE(AtomPair)         :: Pair
        REAL(DOUBLE)           :: Xi,Dist
        IF(Xi*Dist>PrimPairDistanceThreshold) THEN
           TestPrimPair = .FALSE.
        ELSE
           TestPrimPair = .TRUE.
        ENDIF
     END FUNCTION TestPrimPair
!--------------------------------------------------------------------------
!    Compute the extent of a Gaussian with exponent Zeta and amplitude Amp:
!
!    Amp*Exp[-Zeta*Extent^2] > Tau  
!
     SUBROUTINE SetPenetrationThresh(Tau)
        REAL(DOUBLE),INTENT(IN) :: Tau
        PenetratDistanceThreshold=-LOG(Tau)
     END SUBROUTINE SetPenetrationThresh
!
     FUNCTION GaussianExtent(Zeta,Amp)
        REAL(DOUBLE),INTENT(IN) :: Zeta,Amp
        REAL(DOUBLE)            :: GaussianExtent
        GaussianExtent=SQRT(MAX(1.D-10,PenetratDistanceThreshold+Amp)/Zeta)
     END FUNCTION GaussianExtent
!
END MODULE
