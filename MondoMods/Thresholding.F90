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
  USE InOut
!-------------------------------------------------  
!  Primary thresholds
!
   TYPE(TOLS), SAVE :: Thresholds
!-------------------------------------------------------------------------
! Intermediate thresholds
!
  REAL(DOUBLE), SAVE  :: MinRadialExponent
  REAL(DOUBLE), SAVE  :: MinDensityExponent
  REAL(DOUBLE), SAVE  :: AtomPairDistanceThreshold  ! Atom pairs
  REAL(DOUBLE), SAVE  :: PrimPairDistanceThreshold  ! Prim pairs
  REAL(DOUBLE), SAVE  :: PenetratDistanceThreshold  ! Penetration threshold
  CONTAINS
!------------------------------------------------------------------------------
!    Set and load global threholding parameters
!
     SUBROUTINE SetThresholds(CurBase)
         INTEGER          :: NExpt
         TYPE(DBL_VECT)   :: Expts
         CHARACTER(LEN=*) :: CurBase
!        Get the primary thresholds
         CALL Get(Thresholds,Tag_O=CurBase)
!        Get distribution exponents
         CALL Get(NExpt,'nexpt',Tag_O=CurBase)
         CALL New(Expts,NExpt)
         CALL Get(Expts,'dexpt',Tag_O=CurBase)
         MinDensityExponent=Half*Expts%D(1)
!        Xi_Min=Zeta*Zeta/(Zeta+Zeta)=Zeta_Min/2
         MinRadialExponent=Half*Expts%D(1)
         CALL Delete(Expts)
!        Set Atom-Atom thresholds
         CALL SetAtomPairThresh(Thresholds%Dist)
!        Set Prim-Prim thresholds
         CALL SetPrimPairThresh(Thresholds%Dist)
     END SUBROUTINE SetThresholds
!------------------------------------------------------------------------------
!    Set the Atom Pair Distance Threshhold: Zeta_Min*|A-B|^2 > -Log(Tau)
!
     SUBROUTINE SetAtomPairThresh(Tau)
        REAL(DOUBLE),INTENT(IN) :: Tau
        AtomPairDistanceThreshold=-LOG(Tau)/MinRadialExponent
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
     SUBROUTINE SetPenetrationThresh(Tau,ExpSwitch_O)
        REAL(DOUBLE),         INTENT(IN) :: Tau
        REAL(DOUBLE),OPTIONAL,INTENT(IN) :: ExpSwitch_O
        PenetratDistanceThreshold=-LOG(Tau)
        IF(PRESENT(ExpSwitch_O)) &
           PenetratDistanceThreshold=MIN(PenetratDistanceThreshold,ExpSwitch_O)
     END SUBROUTINE SetPenetrationThresh
!
     FUNCTION GaussianExtent(Zeta,Amp)
        REAL(DOUBLE),INTENT(IN) :: Zeta,Amp
        REAL(DOUBLE)            :: GaussianExtent
        GaussianExtent=SQRT(MAX(1.D-10,PenetratDistanceThreshold+Amp)/Zeta)
     END FUNCTION GaussianExtent
END MODULE
