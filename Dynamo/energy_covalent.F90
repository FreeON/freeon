!===============================================================================
!
!       Droits de reproduction et de diffusion réservés. © 2000 CEA/CNRS.
!              (Laboratoire de Dynamique Moléculaire/IBS/DSV) 2000
!
!===============================================================================
!
!                Copyright © 2000 CEA/CNRS. All Rights Reserved.
!              (Laboratoire de Dynamique Moléculaire/IBS/DSV) 2000
!
!===============================================================================
!                           The Covalent Energy Module
!===============================================================================
!
! . Subroutines (public):
!
!   ENERGY_ANGLE                   Calculate the angle energy.
!   ENERGY_BOND                    Calculate the bond energy.
!   ENERGY_DIHEDRAL                Calculate the dihedral energy.
!   ENERGY_IMPROPER                Calculate the improper energy.
!
! . Subroutines (private):
!
!   DIHEDRAL_ENERGY                Evaluate the dihedral energy.
!
!===============================================================================
MODULE ENERGY_COVALENT


! . Module declarations.
#ifdef MMech
USE DerivedTypes
USE InOut 
USE GlobalCharacters, Only: InfFile
#endif
USE DEFINITIONS, ONLY : DP
USE STATUS,      ONLY : ERROR

USE ATOMS,       ONLY : ATMCRD, ATMFIX, ATMIND, MM_NATOMS, NFREE
USE MM_TERMS

IMPLICIT NONE
PRIVATE
PUBLIC :: ENERGY_ANGLE, ENERGY_BOND, ENERGY_DIHEDRAL, ENERGY_IMPROPER

!===============================================================================
CONTAINS
!===============================================================================

#ifdef MMech
   !----------------------------------------------------
   SUBROUTINE ENERGY_ANGLE ( EANGLE, GRADIENT, HESSIAN, InfFile )
   !----------------------------------------------------
#else         
   !----------------------------------------------------
   SUBROUTINE ENERGY_ANGLE ( EANGLE, GRADIENT, HESSIAN )
   !----------------------------------------------------
#endif

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT) :: EANGLE

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:MM_NATOMS),              INTENT(INOUT), OPTIONAL :: GRADIENT
   REAL ( KIND = DP ), DIMENSION(1:(3*NFREE*(3*NFREE+1))/2), INTENT(INOUT), OPTIONAL :: HESSIAN

   ! . Local scalars.
   INTEGER            :: IANGLE, I, IFAC, INDEX, J, JFAC, K, KFAC
   LOGICAL            :: QGRADIENT, QHESSIAN
   REAL ( KIND = DP ) :: DDF, DF, DISP, DOTFAC, DTDX, D2TDX2, RIJ, RKJ, THETA

   ! . Local Hessian scalars.
   REAL ( KIND = DP ) :: D2IXX, D2IXY, D2IXZ, D2IYY, D2IYZ, D2IZZ, &
                         D2KXX, D2KXY, D2KXZ, D2KYY, D2KYZ, D2KZZ, &
                         D2IKXX, D2IKXY, D2IKXZ, D2IKYX, D2IKYY, D2IKYZ, D2IKZX, D2IKZY, D2IKZZ

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DRIJ, DRKJ, DTI, DTJ, DTK

   ! . Parameter declarations.
   REAL ( KIND = DP ), PARAMETER :: DOT_LIMIT = 1.0_DP - 1.0E-6_DP

#ifdef MMech
   TYPE(INT_VECT) :: ACTIVE_ANGLE
   CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
#endif

   ! . Initialization.
   EANGLE = 0.0_DP

   ! . Check the derivative options.
   QGRADIENT = PRESENT ( GRADIENT )
   QHESSIAN  = PRESENT ( HESSIAN  )

   ! . Check the consistency of the derivative options.
   IF ( QHESSIAN .AND. .NOT. QGRADIENT ) CALL ERROR ( "ENERGY_ANGLE", "First derivative argument missing" )

#ifdef MMech
!  IF(PRESENT(InfFile)) THEN
       IF(NANGLES/=0) THEN
       CALL NEW(ACTIVE_ANGLE,NANGLES)
       CALL Get(ACTIVE_ANGLE,'ACTIVE_ANGLE')
       ENDIF
!  ENDIF
#endif

   ! . Loop over the angles.
   DO IANGLE = 1,NANGLES

#ifdef MMech
      IF(ACTIVE_ANGLE%I(IANGLE)==1) THEN
#endif

      ! . Obtain the atoms for the angle.
      I = ANGLES(IANGLE)%I
      J = ANGLES(IANGLE)%J
      K = ANGLES(IANGLE)%K

      ! . Skip the term if all atoms are fixed.
      IF ( ATMFIX(I) .AND. ATMFIX(J) .AND. ATMFIX(K) ) CYCLE

      ! . Get the interatomic vectors.
      DRIJ = ATMCRD(1:3,I) - ATMCRD(1:3,J)
      DRKJ = ATMCRD(1:3,K) - ATMCRD(1:3,J)

      ! . Calculate their sizes.
      RIJ = SQRT ( DOT_PRODUCT ( DRIJ, DRIJ ) )
      RKJ = SQRT ( DOT_PRODUCT ( DRKJ, DRKJ ) )

      ! . Normalize the vectors.
      DRIJ = DRIJ / RIJ
      DRKJ = DRKJ / RKJ

      ! . Calculate the dot product between the two vectors.
      DOTFAC = DOT_PRODUCT ( DRIJ, DRKJ )

      ! . Ensure DOTFAC is between -1 and 1.
      DOTFAC = SIGN ( MIN ( ABS ( DOTFAC ), DOT_LIMIT ), DOTFAC )

      ! . Get the angle.
      THETA = ACOS ( DOTFAC )

      ! . Calculate some intermediate factors.
      DISP = THETA - ANGLES(IANGLE)%EQ
      DF   = ANGLES(IANGLE)%FC * DISP


      ! . Calculate the angle's contribution to the energy.
      EANGLE = EANGLE + ( DF * DISP )


      ! . Check for a gradient calculation.
      IF ( .NOT. QGRADIENT ) CYCLE

      ! . Calculate some factors for the first derivative calculation.
      DTDX = - 1.0_DP / SQRT ( 1.0_DP - DOTFAC * DOTFAC )
      DF   = 2.0_DP * DF * DTDX

      ! . Calculate the components of the derivatives.
      DTI  = ( DRKJ - DOTFAC * DRIJ ) / RIJ
      DTK  = ( DRIJ - DOTFAC * DRKJ ) / RKJ
      DTJ  = - ( DTI + DTK )

      ! . Calculate the gradients.
      GRADIENT(1:3,I) = GRADIENT(1:3,I) + DF * DTI
      GRADIENT(1:3,J) = GRADIENT(1:3,J) + DF * DTJ
      GRADIENT(1:3,K) = GRADIENT(1:3,K) + DF * DTK

      ! . Check for a Hessian calculation.
      IF ( .NOT. QHESSIAN ) CYCLE

      ! . Calculate some factors for the second derivative calculation.
      D2TDX2 = DOTFAC * DTDX * DTDX * DTDX
      DDF    = 2.0_DP * ANGLES(IANGLE)%FC * ( DTDX * DTDX + DISP * D2TDX2 )

      ! . Calculate some diagonal Hessian terms for I.
      D2IXX = - DF * ( 2.0_DP * DTI(1) * DRIJ(1) + DOTFAC * ( 1.0_DP - DRIJ(1) * DRIJ(1) ) / RIJ ) / RIJ
      D2IYY = - DF * ( 2.0_DP * DTI(2) * DRIJ(2) + DOTFAC * ( 1.0_DP - DRIJ(2) * DRIJ(2) ) / RIJ ) / RIJ
      D2IZZ = - DF * ( 2.0_DP * DTI(3) * DRIJ(3) + DOTFAC * ( 1.0_DP - DRIJ(3) * DRIJ(3) ) / RIJ ) / RIJ
      D2IXY = DF * ( DOTFAC * DRIJ(1) * DRIJ(2) / RIJ - DTI(1) * DRIJ(2) - DTI(2) * DRIJ(1) ) / RIJ
      D2IXZ = DF * ( DOTFAC * DRIJ(1) * DRIJ(3) / RIJ - DTI(1) * DRIJ(3) - DTI(3) * DRIJ(1) ) / RIJ
      D2IYZ = DF * ( DOTFAC * DRIJ(2) * DRIJ(3) / RIJ - DTI(2) * DRIJ(3) - DTI(3) * DRIJ(2) ) / RIJ

      ! . Calculate some diagonal Hessian terms for K.
      D2KXX = - DF * ( 2.0_DP * DTK(1) * DRKJ(1) + DOTFAC * ( 1.0_DP - DRKJ(1) * DRKJ(1) ) / RKJ ) / RKJ
      D2KYY = - DF * ( 2.0_DP * DTK(2) * DRKJ(2) + DOTFAC * ( 1.0_DP - DRKJ(2) * DRKJ(2) ) / RKJ ) / RKJ
      D2KZZ = - DF * ( 2.0_DP * DTK(3) * DRKJ(3) + DOTFAC * ( 1.0_DP - DRKJ(3) * DRKJ(3) ) / RKJ ) / RKJ
      D2KXY = DF * ( DOTFAC * DRKJ(1) * DRKJ(2) / RKJ - DTK(1) * DRKJ(2) - DTK(2) * DRKJ(1) ) / RKJ
      D2KXZ = DF * ( DOTFAC * DRKJ(1) * DRKJ(3) / RKJ - DTK(1) * DRKJ(3) - DTK(3) * DRKJ(1) ) / RKJ
      D2KYZ = DF * ( DOTFAC * DRKJ(2) * DRKJ(3) / RKJ - DTK(2) * DRKJ(3) - DTK(3) * DRKJ(2) ) / RKJ

      ! . Calculate the IK off-diagonal terms.
      D2IKXX = DF * ( DOTFAC * DRIJ(1) * DRKJ(1) - DRIJ(1) * DRIJ(1) - DRKJ(1) * DRKJ(1) + 1.0_DP ) / ( RIJ * RKJ )
      D2IKXY = DF * ( DOTFAC * DRIJ(1) * DRKJ(2) - DRIJ(1) * DRIJ(2) - DRKJ(1) * DRKJ(2)          ) / ( RIJ * RKJ )
      D2IKXZ = DF * ( DOTFAC * DRIJ(1) * DRKJ(3) - DRIJ(1) * DRIJ(3) - DRKJ(1) * DRKJ(3)          ) / ( RIJ * RKJ )
      D2IKYX = DF * ( DOTFAC * DRIJ(2) * DRKJ(1) - DRIJ(2) * DRIJ(1) - DRKJ(2) * DRKJ(1)          ) / ( RIJ * RKJ )
      D2IKYY = DF * ( DOTFAC * DRIJ(2) * DRKJ(2) - DRIJ(2) * DRIJ(2) - DRKJ(2) * DRKJ(2) + 1.0_DP ) / ( RIJ * RKJ )
      D2IKYZ = DF * ( DOTFAC * DRIJ(2) * DRKJ(3) - DRIJ(2) * DRIJ(3) - DRKJ(2) * DRKJ(3)          ) / ( RIJ * RKJ )
      D2IKZX = DF * ( DOTFAC * DRIJ(3) * DRKJ(1) - DRIJ(3) * DRIJ(1) - DRKJ(3) * DRKJ(1)          ) / ( RIJ * RKJ )
      D2IKZY = DF * ( DOTFAC * DRIJ(3) * DRKJ(2) - DRIJ(3) * DRIJ(2) - DRKJ(3) * DRKJ(2)          ) / ( RIJ * RKJ )
      D2IKZZ = DF * ( DOTFAC * DRIJ(3) * DRKJ(3) - DRIJ(3) * DRIJ(3) - DRKJ(3) * DRKJ(3) + 1.0_DP ) / ( RIJ * RKJ )

      ! . Calculate some index factors.
      IFAC = 3 * ( ATMIND(I) - 1 )
      JFAC = 3 * ( ATMIND(J) - 1 )
      KFAC = 3 * ( ATMIND(K) - 1 )

      ! . Calculate the II block of the hessian.
      IF ( ATMIND(I) > 0 ) THEN
         INDEX = ( IFAC * ( IFAC + 1 ) ) / 2 + IFAC + 1
         HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTI(1) + D2IXX
         INDEX = INDEX + IFAC + 1
         HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTI(2) + D2IXY
         HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTI(2) + D2IYY
         INDEX = INDEX + IFAC + 2
         HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTI(3) + D2IXZ
         HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTI(3) + D2IYZ
         HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(3) * DTI(3) + D2IZZ
      END IF

      ! . Calculate the JJ block of the hessian.
      IF ( ATMIND(J) > 0 ) THEN
         INDEX = ( JFAC * ( JFAC + 1 ) ) / 2 + JFAC + 1
         HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTJ(1) * DTJ(1) + ( D2IXX + D2KXX + 2.0_DP * D2IKXX )
         INDEX = INDEX + JFAC + 1
         HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTJ(1) * DTJ(2) + ( D2IXY + D2KXY + D2IKXY + D2IKYX )
         HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTJ(2) * DTJ(2) + ( D2IYY + D2KYY + 2.0_DP * D2IKYY )
         INDEX = INDEX + JFAC + 2
         HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTJ(1) * DTJ(3) + ( D2IXZ + D2KXZ + D2IKXZ + D2IKZX )
         HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTJ(2) * DTJ(3) + ( D2IYZ + D2KYZ + D2IKYZ + D2IKZY )
         HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTJ(3) * DTJ(3) + ( D2IZZ + D2KZZ + 2.0_DP * D2IKZZ )
      END IF

      ! . Calculate the KK block of the hessian.
      IF ( ATMIND(J) > 0 ) THEN
         INDEX = ( KFAC * ( KFAC + 1 ) ) / 2 + KFAC + 1
         HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTK(1) * DTK(1) + D2KXX
         INDEX = INDEX + KFAC + 1
         HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTK(1) * DTK(2) + D2KXY
         HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTK(2) * DTK(2) + D2KYY
         INDEX = INDEX + KFAC + 2
         HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTK(1) * DTK(3) + D2KXZ
         HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTK(2) * DTK(3) + D2KYZ
         HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTK(3) * DTK(3) + D2KZZ
      END IF

      ! . Calculate the IJ block of the hessian.
      IF ( ( ATMIND(I) > 0 ) .AND. ( ATMIND(J) > 0 ) ) THEN
         IF ( I > J ) THEN
            INDEX = ( IFAC * ( IFAC + 1 ) ) / 2 + JFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTJ(1) - D2IXX - D2IKXX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(1) * DTJ(2) - D2IXY - D2IKXY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(1) * DTJ(3) - D2IXZ - D2IKXZ
            INDEX = INDEX + IFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(2) * DTJ(1) - D2IXY - D2IKYX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTJ(2) - D2IYY - D2IKYY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(2) * DTJ(3) - D2IYZ - D2IKYZ
            INDEX = INDEX + IFAC + 2
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(3) * DTJ(1) - D2IXZ - D2IKZX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(3) * DTJ(2) - D2IYZ - D2IKZY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(3) * DTJ(3) - D2IZZ - D2IKZZ
         ELSE
            INDEX = ( JFAC * ( JFAC + 1 ) ) / 2 + IFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTJ(1) - D2IXX - D2IKXX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTJ(1) - D2IXY - D2IKYX
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(3) * DTJ(1) - D2IXZ - D2IKZX
            INDEX = INDEX + JFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTJ(2) - D2IXY - D2IKXY
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTJ(2) - D2IYY - D2IKYY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(3) * DTJ(2) - D2IYZ - D2IKZY
            INDEX = INDEX + JFAC + 2
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTJ(3) - D2IXZ - D2IKXZ
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTJ(3) - D2IYZ - D2IKYZ
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(3) * DTJ(3) - D2IZZ - D2IKZZ
         END IF
      END IF

      ! . Calculate the IK block of the hessian.
      IF ( ( ATMIND(I) > 0 ) .AND. ( ATMIND(K) > 0 ) ) THEN
         IF ( I > K ) THEN
            INDEX = ( IFAC * ( IFAC + 1 ) ) / 2 + KFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTK(1) + D2IKXX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(1) * DTK(2) + D2IKXY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(1) * DTK(3) + D2IKXZ
            INDEX = INDEX + IFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(2) * DTK(1) + D2IKYX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTK(2) + D2IKYY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(2) * DTK(3) + D2IKYZ
            INDEX = INDEX + IFAC + 2
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(3) * DTK(1) + D2IKZX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(3) * DTK(2) + D2IKZY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(3) * DTK(3) + D2IKZZ
         ELSE
            INDEX = ( KFAC * ( KFAC + 1 ) ) / 2 + IFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTK(1) + D2IKXX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTK(1) + D2IKYX
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(3) * DTK(1) + D2IKZX
            INDEX = INDEX + KFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTK(2) + D2IKXY
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTK(2) + D2IKYY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(3) * DTK(2) + D2IKZY
            INDEX = INDEX + KFAC + 2
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTK(3) + D2IKXZ
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTK(3) + D2IKYZ
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(3) * DTK(3) + D2IKZZ
         END IF
      END IF

      ! . Calculate the JK block of the hessian.
      IF ( ( ATMIND(J) > 0 ) .AND. ( ATMIND(K) > 0 ) ) THEN
         IF ( K > J ) THEN
            INDEX = ( KFAC * ( KFAC + 1 ) ) / 2 + JFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTK(1) * DTJ(1) - D2KXX - D2IKXX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTK(1) * DTJ(2) - D2KXY - D2IKYX
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTK(1) * DTJ(3) - D2KXZ - D2IKZX
            INDEX = INDEX + KFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTK(2) * DTJ(1) - D2KXY - D2IKXY
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTK(2) * DTJ(2) - D2KYY - D2IKYY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTK(2) * DTJ(3) - D2KYZ - D2IKZY
            INDEX = INDEX + KFAC + 2
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTK(3) * DTJ(1) - D2KXZ - D2IKXZ
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTK(3) * DTJ(2) - D2KYZ - D2IKYZ
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTK(3) * DTJ(3) - D2KZZ - D2IKZZ
         ELSE
            INDEX = ( JFAC * ( JFAC + 1 ) ) / 2 + KFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTK(1) * DTJ(1) - D2KXX - D2IKXX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTK(2) * DTJ(1) - D2KXY - D2IKXY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTK(3) * DTJ(1) - D2KXZ - D2IKXZ
            INDEX = INDEX + JFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTK(1) * DTJ(2) - D2KXY - D2IKYX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTK(2) * DTJ(2) - D2KYY - D2IKYY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTK(3) * DTJ(2) - D2KYZ - D2IKYZ
            INDEX = INDEX + JFAC + 2
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTK(1) * DTJ(3) - D2KXZ - D2IKZX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTK(2) * DTJ(3) - D2KYZ - D2IKZY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTK(3) * DTJ(3) - D2KZZ - D2IKZZ
         END IF
      END IF

#ifdef MMech
      ENDIF
#endif
   END DO
#ifdef MMech
     IF(NANGLES/=0) CALL DELETE(ACTIVE_ANGLE)
#endif

   END SUBROUTINE ENERGY_ANGLE

#ifdef MMech
   !----------------------------------------------------------
   SUBROUTINE ENERGY_BOND ( EBOND, VIRIAL, GRADIENT, HESSIAN, InfFile )
   !----------------------------------------------------------
#else
   !----------------------------------------------------------
   SUBROUTINE ENERGY_BOND ( EBOND, VIRIAL, GRADIENT, HESSIAN )
   !----------------------------------------------------------
#endif

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT)   :: EBOND
   REAL ( KIND = DP ), INTENT(INOUT) :: VIRIAL

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:MM_NATOMS),              INTENT(INOUT), OPTIONAL :: GRADIENT
   REAL ( KIND = DP ), DIMENSION(1:(3*NFREE*(3*NFREE+1))/2), INTENT(INOUT), OPTIONAL :: HESSIAN

   ! . Local scalars.
   INTEGER            :: I, IBOND, IFAC, INDEX, J, JFAC, SWAP
   LOGICAL            :: QGRADIENT, QHESSIAN
   REAL ( KIND = DP ) :: DF, DISP, HXX, HXY, HXZ, HYY, HYZ, HZZ, KDF, RIJ

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DRIJ

#ifdef MMech
   TYPE(INT_VECT) :: ACTIVE_BOND
   CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
!  IF(PRESENT(InfFile)) THEN
       IF(NBONDS/=0) THEN
       CALL NEW(ACTIVE_BOND,NBONDS)
       CALL Get(ACTIVE_BOND,'ACTIVE_BOND')
       ENDIF
!  ENDIF
#endif

   ! . Initialization.
   EBOND = 0.0_DP

   ! . Check the derivative options.
   QGRADIENT = PRESENT ( GRADIENT )
   QHESSIAN  = PRESENT ( HESSIAN  )

   ! . Check the consistency of the derivative options.
   IF ( QHESSIAN .AND. .NOT. QGRADIENT ) CALL ERROR ( "ENERGY_BOND", "First derivative argument missing" )

   ! . Loop over the bonds.
   DO IBOND = 1,NBONDS

#ifdef MMech
   IF(ACTIVE_BOND%I(IBOND)==1) THEN
#endif

      ! . Obtain the atoms for the bond.
      I = BONDS(IBOND)%I
      J = BONDS(IBOND)%J

      ! . Skip the term if all atoms are fixed.
      IF ( ATMFIX(I) .AND. ATMFIX(J) ) CYCLE

      ! . Calculate some intermediate factors.
      DRIJ = ATMCRD(1:3,I) - ATMCRD(1:3,J)
      RIJ  = SQRT ( DOT_PRODUCT ( DRIJ, DRIJ ) )
      DISP = RIJ - BONDS(IBOND)%EQ
      DF   = BONDS(IBOND)%FC * DISP

      ! . Calculate the bond's contribution to the bond energy.
      EBOND = EBOND + ( DF * DISP )

      ! . Check for a gradient calculation.
      IF ( .NOT. QGRADIENT ) CYCLE

      ! . Calculate the contribution to the virial.
      VIRIAL = VIRIAL + ( 2.0_DP * DF * RIJ )

      ! . Calculate the gradients.
      DF = ( 2.0_DP * DF ) / RIJ
      GRADIENT(1:3,I) = GRADIENT(1:3,I) + ( DF * DRIJ )
      GRADIENT(1:3,J) = GRADIENT(1:3,J) - ( DF * DRIJ )

      ! . Check for a Hessian calculation.
      IF ( .NOT. QHESSIAN ) CYCLE

      ! . Normalize the displacement vector.
      DRIJ = DRIJ / RIJ

      ! . Calculate a hessian factor.
      KDF = 2.0_DP * BONDS(IBOND)%FC - DF

      ! . Calculate the elements of the hessian block.
      HXX = DRIJ(1) * DRIJ(1) * KDF + DF
      HXY = DRIJ(1) * DRIJ(2) * KDF
      HXZ = DRIJ(1) * DRIJ(3) * KDF
      HYY = DRIJ(2) * DRIJ(2) * KDF + DF
      HYZ = DRIJ(2) * DRIJ(3) * KDF
      HZZ = DRIJ(3) * DRIJ(3) * KDF + DF

      ! . Calculate some index factors.
      IFAC = 3 * ( ATMIND(I) - 1 )
      JFAC = 3 * ( ATMIND(J) - 1 )

      ! . Calculate the II block of the hessian.
      IF ( ATMIND(I) > 0 ) THEN
         INDEX = ( IFAC * ( IFAC + 1 ) ) / 2 + IFAC + 1
         HESSIAN(INDEX)   = HESSIAN(INDEX)   + HXX
         INDEX = INDEX + IFAC + 1
         HESSIAN(INDEX)   = HESSIAN(INDEX)   + HXY
         HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + HYY
         INDEX = INDEX + IFAC + 2
         HESSIAN(INDEX)   = HESSIAN(INDEX)   + HXZ
         HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + HYZ
         HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + HZZ
      END IF

      ! . Calculate the JJ block of the hessian.
      IF ( ATMIND(J) > 0 ) THEN
	 INDEX = ( JFAC * ( JFAC + 1 ) ) / 2 + JFAC + 1
	 HESSIAN(INDEX)   = HESSIAN(INDEX)   + HXX
	 INDEX = INDEX + JFAC + 1
	 HESSIAN(INDEX)   = HESSIAN(INDEX)   + HXY
	 HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + HYY
	 INDEX = INDEX + JFAC + 2
	 HESSIAN(INDEX)   = HESSIAN(INDEX)   + HXZ
	 HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + HYZ
	 HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + HZZ
      END IF
         
      ! . Swap the I and J indices if I is less than J.
      IF ( I < J ) THEN
         SWAP = IFAC ; IFAC = JFAC ; JFAC = SWAP
      END IF

      ! . Calculate the IJ block of the hessian.
      IF ( ( ATMIND(I) > 0 ) .AND. ( ATMIND(J) > 0 ) ) THEN
	 INDEX = ( IFAC * ( IFAC + 1 ) ) / 2 + JFAC + 1
	 HESSIAN(INDEX)   = HESSIAN(INDEX)   - HXX
	 HESSIAN(INDEX+1) = HESSIAN(INDEX+1) - HXY
	 HESSIAN(INDEX+2) = HESSIAN(INDEX+2) - HXZ
	 INDEX = INDEX + IFAC + 1
	 HESSIAN(INDEX)   = HESSIAN(INDEX)   - HXY
	 HESSIAN(INDEX+1) = HESSIAN(INDEX+1) - HYY
	 HESSIAN(INDEX+2) = HESSIAN(INDEX+2) - HYZ
	 INDEX = INDEX + IFAC + 2
	 HESSIAN(INDEX)   = HESSIAN(INDEX)   - HXZ
	 HESSIAN(INDEX+1) = HESSIAN(INDEX+1) - HYZ
	 HESSIAN(INDEX+2) = HESSIAN(INDEX+2) - HZZ
      END IF

#ifdef MMech
   ENDIF
#endif

   END DO

#ifdef MMech
   IF(NBONDS/=0) CALL DELETE(ACTIVE_BOND)
#endif

   END SUBROUTINE ENERGY_BOND

#ifdef MMech
   !-----------------------------------------------------------
   SUBROUTINE ENERGY_DIHEDRAL ( EDIHEDRAL, GRADIENT, HESSIAN, InfFile )
   !-----------------------------------------------------------
#else
   !-----------------------------------------------------------
   SUBROUTINE ENERGY_DIHEDRAL ( EDIHEDRAL, GRADIENT, HESSIAN )
   !-----------------------------------------------------------
#endif

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT) :: EDIHEDRAL

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:MM_NATOMS),              INTENT(INOUT), OPTIONAL :: GRADIENT
   REAL ( KIND = DP ), DIMENSION(1:(3*NFREE*(3*NFREE+1))/2), INTENT(INOUT), OPTIONAL :: HESSIAN
#ifdef MMech
   TYPE(INT_VECT) :: ACTIVE_DIHEDRAL
   CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
!  IF(PRESENT(InfFile)) THEN
       IF(NDIHEDRALS/=0) THEN
       CALL NEW(ACTIVE_DIHEDRAL,NDIHEDRALS)
       CALL Get(ACTIVE_DIHEDRAL,'ACTIVE_DIHEDRAL')
       ENDIF
!  ENDIF
#endif

   ! . Call the energy evaluation routine with the appropriate lists.

#ifdef MMech
   CALL DIHEDRAL_ENERGY ( NDIHEDRALS, DIHEDRALS, EDIHEDRAL, GRADIENT, HESSIAN, ACTIVE_DIHEDRAL )
   IF(NDIHEDRALS/=0) CALL DELETE(ACTIVE_DIHEDRAL)
#else
   CALL DIHEDRAL_ENERGY ( NDIHEDRALS, DIHEDRALS, EDIHEDRAL, GRADIENT, HESSIAN )
#endif

   END SUBROUTINE ENERGY_DIHEDRAL

#ifdef MMech
   !-----------------------------------------------------------
   SUBROUTINE ENERGY_IMPROPER ( EIMPROPER, GRADIENT, HESSIAN, InfFile )
   !-----------------------------------------------------------
#else
   !-----------------------------------------------------------
   SUBROUTINE ENERGY_IMPROPER ( EIMPROPER, GRADIENT, HESSIAN )
   !-----------------------------------------------------------
#endif

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT) :: EIMPROPER

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:MM_NATOMS),              INTENT(INOUT), OPTIONAL :: GRADIENT
   REAL ( KIND = DP ), DIMENSION(1:(3*NFREE*(3*NFREE+1))/2), INTENT(INOUT), OPTIONAL :: HESSIAN

#ifdef MMech
   TYPE(INT_VECT) :: ACTIVE_IMPROPER
   CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
!  IF(PRESENT(InfFile)) THEN
       IF(NIMPROPERS/=0) THEN
       CALL NEW(ACTIVE_IMPROPER,NIMPROPERS)
       CALL Get(ACTIVE_IMPROPER,'ACTIVE_IMPROPER')
       ENDIF
!  ENDIF
#endif

   ! . Call the energy evaluation routine with the appropriate lists.
#ifdef MMech
   CALL DIHEDRAL_ENERGY ( NIMPROPERS, IMPROPERS, EIMPROPER, GRADIENT, HESSIAN, ACTIVE_IMPROPER )
   IF(NIMPROPERS/=0) CALL DELETE(ACTIVE_IMPROPER)
#else
   CALL DIHEDRAL_ENERGY ( NIMPROPERS, IMPROPERS, EIMPROPER, GRADIENT, HESSIAN )
#endif

   END SUBROUTINE ENERGY_IMPROPER

   !============================================================================
   ! . Private Subroutines.
   !============================================================================

   !---------------------------------------------------------------------------------
#ifdef MMech
   SUBROUTINE DIHEDRAL_ENERGY ( NDIHEDRALS, DIHEDRALS, EDIHEDRAL, GRADIENT, HESSIAN, ACTIVE )
#else
   SUBROUTINE DIHEDRAL_ENERGY ( NDIHEDRALS, DIHEDRALS, EDIHEDRAL, GRADIENT, HESSIAN )
#endif
   !---------------------------------------------------------------------------------

   ! . Dihedral data.
   INTEGER,                           INTENT(IN) :: NDIHEDRALS
   TYPE(DIHEDRAL_TERM), DIMENSION(:), INTENT(IN) :: DIHEDRALS

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT) :: EDIHEDRAL

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:MM_NATOMS),              INTENT(INOUT), OPTIONAL :: GRADIENT
   REAL ( KIND = DP ), DIMENSION(1:(3*NFREE*(3*NFREE+1))/2), INTENT(INOUT), OPTIONAL :: HESSIAN

   ! . Local scalars.
   INTEGER            :: IDIHE, I, IFAC, J, JFAC, K, KFAC, L, LFAC
   LOGICAL            :: QGRADIENT, QHESSIAN
   REAL ( KIND = DP ) :: COSNPHI, COSPHI, COSPHI2, D, DCOS, DDF, DF, DOTIJ, DOTLK, D2COS, E, R, RKJ, S

   ! . Local Hessian scalars.
   REAL ( KIND = DP ) :: D2IXX, D2IXY, D2IXZ, D2IYY, D2IYZ, D2IZZ, &
                         D2JXX, D2JXY, D2JXZ, D2JYY, D2JYZ, D2JZZ, &
                         D2KXX, D2KXY, D2KXZ, D2KYY, D2KYZ, D2KZZ, &
                         D2LXX, D2LXY, D2LXZ, D2LYY, D2LYZ, D2LZZ, &
                         D2IJXX, D2IJXY, D2IJXZ, D2IJYX, D2IJYY, D2IJYZ, D2IJZX, D2IJZY, D2IJZZ, &
                         D2IKXX, D2IKXY, D2IKXZ, D2IKYX, D2IKYY, D2IKYZ, D2IKZX, D2IKZY, D2IKZZ, &
                         D2ILXX, D2ILXY, D2ILXZ, D2ILYX, D2ILYY, D2ILYZ, D2ILZX, D2ILZY, D2ILZZ, &
                         D2JKXX, D2JKXY, D2JKXZ, D2JKYX, D2JKYY, D2JKYZ, D2JKZX, D2JKZY, D2JKZZ, &
                         D2JLXX, D2JLXY, D2JLXZ, D2JLYX, D2JLYY, D2JLYZ, D2JLZX, D2JLZY, D2JLZZ, &
                         D2KLXX, D2KLXY, D2KLXZ, D2KLYX, D2KLYY, D2KLYZ, D2KLZX, D2KLZY, D2KLZZ

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DDJ, DDK, DEJ, DEK, DIL, DR, DRIJ, DRKJ, DRLK, DS, DTI, DTJ, DTK, DTL

   ! . Parameter declarations.
   REAL ( KIND = DP ), PARAMETER :: DOT_LIMIT = 1.0_DP

#ifdef MMech
   TYPE(INT_VECT) :: ACTIVE
#endif

   ! . Initialization.
   EDIHEDRAL = 0.0_DP

   ! . Check the derivative options.
   QGRADIENT = PRESENT ( GRADIENT )
   QHESSIAN  = PRESENT ( HESSIAN  )

   ! . Check the consistency of the derivative options.
   IF ( QHESSIAN .AND. .NOT. QGRADIENT ) CALL ERROR ( "DIHEDRAL_ENERGY", "First derivative argument missing" )

   ! . Loop over the dihedrals.
   DO IDIHE = 1,NDIHEDRALS
#ifdef MMech
   IF(ACTIVE%I(IDIHE)==1) THEN
#endif

      ! . Obtain the atoms for the dihedral.
      I = DIHEDRALS(IDIHE)%I
      J = DIHEDRALS(IDIHE)%J
      K = DIHEDRALS(IDIHE)%K
      L = DIHEDRALS(IDIHE)%L

      ! . Skip the term if all atoms are fixed.
      IF ( ATMFIX(I) .AND. ATMFIX(J) .AND. ATMFIX(K) .AND. ATMFIX(L) ) CYCLE

      ! . Calculate some displacement vectors.
      DRIJ = ATMCRD(1:3,I) - ATMCRD(1:3,J)
      DRKJ = ATMCRD(1:3,K) - ATMCRD(1:3,J)
      DRLK = ATMCRD(1:3,L) - ATMCRD(1:3,K)

      ! . Calculate the size of the RKJ vector.
      RKJ = SQRT ( DOT_PRODUCT ( DRKJ, DRKJ ) )

      ! . Normalize the vector.
      DRKJ = DRKJ / RKJ

      ! . Calculate some intermediate dot products.
      DOTIJ = DOT_PRODUCT ( DRIJ, DRKJ )
      DOTLK = DOT_PRODUCT ( DRLK, DRKJ )

      ! . Calculate the DR and DS vectors.
      DR = DRIJ - DOTIJ * DRKJ
      DS = DRLK - DOTLK * DRKJ

      ! . Calculate the magnitudes of DR and DS.
      R = SQRT ( DOT_PRODUCT ( DR, DR ) )
      S = SQRT ( DOT_PRODUCT ( DS, DS ) )

      ! . Normalize DR and DS.
      DR = DR / R
      DS = DS / S

      ! . Calculate the dot product of the vectors.
      COSPHI = DOT_PRODUCT ( DR, DS )

      ! . Ensure DOTFAC is between -1 and 1.
      COSPHI = SIGN ( MIN ( ABS ( COSPHI ), DOT_LIMIT ), COSPHI )

      ! . Calculate the square.
      COSPHI2 = COSPHI * COSPHI

      ! . Calculate Cos(n phi) and its derivative for the appropriate periodicity.
      SELECT CASE ( DIHEDRALS(IDIHE)%PERIOD )
      CASE ( 0 ) ; COSNPHI = 1.0_DP
                   DCOS    = 0.0_DP
                   D2COS   = 0.0_DP
      CASE ( 1 ) ; COSNPHI = COSPHI
                   DCOS    = 1.0_DP
                   D2COS   = 0.0_DP
      CASE ( 2 ) ; COSNPHI = 2.0_DP * COSPHI2 - 1.0_DP
                   DCOS    = 4.0_DP * COSPHI
                   D2COS   = 4.0_DP
      CASE ( 3 ) ; COSNPHI = ( 4.0_DP * COSPHI2 - 3.0_DP ) * COSPHI
                   DCOS    = 12.0_DP * COSPHI2 - 3.0_DP
                   D2COS   = 24.0_DP * COSPHI
      CASE ( 4 ) ; COSNPHI = 8.0_DP * ( COSPHI2 - 1.0_DP ) * COSPHI2 + 1.0_DP
                   DCOS    = 16.0_DP * ( 2.0_DP * COSPHI2 - 1.0_DP ) * COSPHI
                   D2COS   = 16.0_DP * ( 6.0_DP * COSPHI2 - 1.0_DP )
      CASE ( 5 ) ; COSNPHI = ( ( 16.0_DP * COSPHI2 - 20.0_DP ) * COSPHI2 + 5.0_DP ) * COSPHI
                   DCOS    = 20.0_DP * ( 4.0_DP * COSPHI2 - 3.0_DP ) * COSPHI2 + 5.0_DP
                   D2COS   = 40.0_DP * ( 8.0_DP * COSPHI2 - 3.0_DP ) * COSPHI
      CASE ( 6 ) ; COSNPHI = ( ( 32.0_DP * COSPHI2 - 48.0_DP ) * COSPHI2 + 18.0_DP ) * COSPHI2 - 1.0_DP
                   DCOS    = ( 192.0_DP * ( COSPHI2 - 1.0_DP ) * COSPHI2 + 36.0_DP ) * COSPHI
                   D2COS   =   192.0_DP * ( 5.0_DP * COSPHI2 - 3.0_DP ) * COSPHI2 + 36.0_DP
      END SELECT

      ! . Calculate the dihedral's contribution to the energy.
      EDIHEDRAL = EDIHEDRAL + DIHEDRALS(IDIHE)%FC * ( 1.0_DP + DIHEDRALS(IDIHE)%PHASE * COSNPHI )

      ! . Check for a gradient calculation.
      IF ( .NOT. QGRADIENT ) CYCLE

      ! . Calculate an intermediate factor for the gradient calculation.
      DF  = DIHEDRALS(IDIHE)%FC * DIHEDRALS(IDIHE)%PHASE * DCOS

      ! . Calculate the first derivative vectors.
      D = DOTIJ / RKJ
      E = DOTLK / RKJ
      DTI   = ( DS - COSPHI * DR ) / R
      DTL   = ( DR - COSPHI * DS ) / S
      DTJ   =   DTI * ( D - 1.0_DP ) + E * DTL
      DTK   = - DTL * ( E + 1.0_DP ) - D * DTI

      ! . Calculate the gradient.
      GRADIENT(1:3,I) = GRADIENT(1:3,I) + DF * DTI
      GRADIENT(1:3,J) = GRADIENT(1:3,J) + DF * DTJ
      GRADIENT(1:3,K) = GRADIENT(1:3,K) + DF * DTK
      GRADIENT(1:3,L) = GRADIENT(1:3,L) + DF * DTL

      ! . Check for a Hessian calculation.
      IF ( .NOT. QHESSIAN ) CYCLE

      ! . Calculate an intermediate factor for the Hessian calculation.
      DDF = DIHEDRALS(IDIHE)%FC * DIHEDRALS(IDIHE)%PHASE * D2COS

      ! . Calculate the derivatives of DOTIJ and DOTLK.
      DIL =   DF * DRKJ / RKJ
      DDJ =   DF * ( ( 2.0_DP * D - 1.0_DP ) * DRKJ - DRIJ / RKJ ) / RKJ
      DDK = - DF * (   2.0_DP * D            * DRKJ - DRIJ / RKJ ) / RKJ
      DEJ =   DF * (   2.0_DP * E            * DRKJ - DRLK / RKJ ) / RKJ
      DEK = - DF * ( ( 2.0_DP * E + 1.0_DP ) * DRKJ - DRLK / RKJ ) / RKJ

      ! . Calculate some diagonal Hessian terms for I.
      D2IXX = DF * ( COSPHI * ( DR(1) * DR(1) + DRKJ(1) * DRKJ(1) - 1.0_DP ) / R - 2.0_DP * DR(1) * DTI(1) ) / R
      D2IYY = DF * ( COSPHI * ( DR(2) * DR(2) + DRKJ(2) * DRKJ(2) - 1.0_DP ) / R - 2.0_DP * DR(2) * DTI(2) ) / R
      D2IZZ = DF * ( COSPHI * ( DR(3) * DR(3) + DRKJ(3) * DRKJ(3) - 1.0_DP ) / R - 2.0_DP * DR(3) * DTI(3) ) / R
      D2IXY = DF * ( COSPHI * ( DR(1) * DR(2) + DRKJ(1) * DRKJ(2) ) / R - DR(1) * DTI(2) - DR(2) * DTI(1) ) / R
      D2IXZ = DF * ( COSPHI * ( DR(1) * DR(3) + DRKJ(1) * DRKJ(3) ) / R - DR(1) * DTI(3) - DR(3) * DTI(1) ) / R
      D2IYZ = DF * ( COSPHI * ( DR(2) * DR(3) + DRKJ(2) * DRKJ(3) ) / R - DR(2) * DTI(3) - DR(3) * DTI(2) ) / R

      ! . Calculate some diagonal Hessian terms for L.
      D2LXX = DF * ( COSPHI * ( DS(1) * DS(1) + DRKJ(1) * DRKJ(1) - 1.0_DP ) / S - 2.0_DP * DS(1) * DTL(1) ) / S
      D2LYY = DF * ( COSPHI * ( DS(2) * DS(2) + DRKJ(2) * DRKJ(2) - 1.0_DP ) / S - 2.0_DP * DS(2) * DTL(2) ) / S
      D2LZZ = DF * ( COSPHI * ( DS(3) * DS(3) + DRKJ(3) * DRKJ(3) - 1.0_DP ) / S - 2.0_DP * DS(3) * DTL(3) ) / S
      D2LXY = DF * ( COSPHI * ( DS(1) * DS(2) + DRKJ(1) * DRKJ(2) ) / S - DS(1) * DTL(2) - DS(2) * DTL(1) ) / S
      D2LXZ = DF * ( COSPHI * ( DS(1) * DS(3) + DRKJ(1) * DRKJ(3) ) / S - DS(1) * DTL(3) - DS(3) * DTL(1) ) / S
      D2LYZ = DF * ( COSPHI * ( DS(2) * DS(3) + DRKJ(2) * DRKJ(3) ) / S - DS(2) * DTL(3) - DS(3) * DTL(2) ) / S

      ! . Calculate the IL off-diagonal terms.
      D2ILXX = DF * ( COSPHI * DR(1) * DS(1) - DR(1) * DR(1) - DS(1) * DS(1) - DRKJ(1) * DRKJ(1) + 1.0_DP ) / ( R * S )
      D2ILXY = DF * ( COSPHI * DR(1) * DS(2) - DR(1) * DR(2) - DS(1) * DS(2) - DRKJ(1) * DRKJ(2)          ) / ( R * S )
      D2ILXZ = DF * ( COSPHI * DR(1) * DS(3) - DR(1) * DR(3) - DS(1) * DS(3) - DRKJ(1) * DRKJ(3)          ) / ( R * S )
      D2ILYX = DF * ( COSPHI * DR(2) * DS(1) - DR(2) * DR(1) - DS(2) * DS(1) - DRKJ(2) * DRKJ(1)          ) / ( R * S )
      D2ILYY = DF * ( COSPHI * DR(2) * DS(2) - DR(2) * DR(2) - DS(2) * DS(2) - DRKJ(2) * DRKJ(2) + 1.0_DP ) / ( R * S )
      D2ILYZ = DF * ( COSPHI * DR(2) * DS(3) - DR(2) * DR(3) - DS(2) * DS(3) - DRKJ(2) * DRKJ(3)          ) / ( R * S )
      D2ILZX = DF * ( COSPHI * DR(3) * DS(1) - DR(3) * DR(1) - DS(3) * DS(1) - DRKJ(3) * DRKJ(1)          ) / ( R * S )
      D2ILZY = DF * ( COSPHI * DR(3) * DS(2) - DR(3) * DR(2) - DS(3) * DS(2) - DRKJ(3) * DRKJ(2)          ) / ( R * S )
      D2ILZZ = DF * ( COSPHI * DR(3) * DS(3) - DR(3) * DR(3) - DS(3) * DS(3) - DRKJ(3) * DRKJ(3) + 1.0_DP ) / ( R * S )

      ! . Calculate the IJ off-diagonal terms.
      D2IJXX = D2IXX * ( D - 1.0_DP ) + D2ILXX * E + DIL(1) * DTI(1)
      D2IJXY = D2IXY * ( D - 1.0_DP ) + D2ILXY * E + DIL(1) * DTI(2)
      D2IJXZ = D2IXZ * ( D - 1.0_DP ) + D2ILXZ * E + DIL(1) * DTI(3)
      D2IJYX = D2IXY * ( D - 1.0_DP ) + D2ILYX * E + DIL(2) * DTI(1)
      D2IJYY = D2IYY * ( D - 1.0_DP ) + D2ILYY * E + DIL(2) * DTI(2)
      D2IJYZ = D2IYZ * ( D - 1.0_DP ) + D2ILYZ * E + DIL(2) * DTI(3)
      D2IJZX = D2IXZ * ( D - 1.0_DP ) + D2ILZX * E + DIL(3) * DTI(1)
      D2IJZY = D2IYZ * ( D - 1.0_DP ) + D2ILZY * E + DIL(3) * DTI(2)
      D2IJZZ = D2IZZ * ( D - 1.0_DP ) + D2ILZZ * E + DIL(3) * DTI(3)

      ! . Calculate the IK off-diagonal terms.
      D2IKXX = - D2IXX * D - D2ILXX * ( E + 1.0_DP ) - DIL(1) * DTI(1)
      D2IKXY = - D2IXY * D - D2ILXY * ( E + 1.0_DP ) - DIL(1) * DTI(2)
      D2IKXZ = - D2IXZ * D - D2ILXZ * ( E + 1.0_DP ) - DIL(1) * DTI(3)
      D2IKYX = - D2IXY * D - D2ILYX * ( E + 1.0_DP ) - DIL(2) * DTI(1)
      D2IKYY = - D2IYY * D - D2ILYY * ( E + 1.0_DP ) - DIL(2) * DTI(2)
      D2IKYZ = - D2IYZ * D - D2ILYZ * ( E + 1.0_DP ) - DIL(2) * DTI(3)
      D2IKZX = - D2IXZ * D - D2ILZX * ( E + 1.0_DP ) - DIL(3) * DTI(1)
      D2IKZY = - D2IYZ * D - D2ILZY * ( E + 1.0_DP ) - DIL(3) * DTI(2)
      D2IKZZ = - D2IZZ * D - D2ILZZ * ( E + 1.0_DP ) - DIL(3) * DTI(3)

      ! . Calculate the JL off-diagonal terms.
      D2JLXX = D2ILXX * ( D - 1.0_DP ) + D2LXX * E + DIL(1) * DTL(1)
      D2JLXY = D2ILXY * ( D - 1.0_DP ) + D2LXY * E + DIL(2) * DTL(1)
      D2JLXZ = D2ILXZ * ( D - 1.0_DP ) + D2LXZ * E + DIL(3) * DTL(1)
      D2JLYX = D2ILYX * ( D - 1.0_DP ) + D2LXY * E + DIL(1) * DTL(2)
      D2JLYY = D2ILYY * ( D - 1.0_DP ) + D2LYY * E + DIL(2) * DTL(2)
      D2JLYZ = D2ILYZ * ( D - 1.0_DP ) + D2LYZ * E + DIL(3) * DTL(2)
      D2JLZX = D2ILZX * ( D - 1.0_DP ) + D2LXZ * E + DIL(1) * DTL(3)
      D2JLZY = D2ILZY * ( D - 1.0_DP ) + D2LYZ * E + DIL(2) * DTL(3)
      D2JLZZ = D2ILZZ * ( D - 1.0_DP ) + D2LZZ * E + DIL(3) * DTL(3)

      ! . Calculate the KL off-diagonal terms.
      D2KLXX = - D2ILXX * D - D2LXX * ( E + 1.0_DP ) - DIL(1) * DTL(1)
      D2KLXY = - D2ILXY * D - D2LXY * ( E + 1.0_DP ) - DIL(2) * DTL(1)
      D2KLXZ = - D2ILXZ * D - D2LXZ * ( E + 1.0_DP ) - DIL(3) * DTL(1)
      D2KLYX = - D2ILYX * D - D2LXY * ( E + 1.0_DP ) - DIL(1) * DTL(2)
      D2KLYY = - D2ILYY * D - D2LYY * ( E + 1.0_DP ) - DIL(2) * DTL(2)
      D2KLYZ = - D2ILYZ * D - D2LYZ * ( E + 1.0_DP ) - DIL(3) * DTL(2)
      D2KLZX = - D2ILZX * D - D2LXZ * ( E + 1.0_DP ) - DIL(1) * DTL(3)
      D2KLZY = - D2ILZY * D - D2LYZ * ( E + 1.0_DP ) - DIL(2) * DTL(3)
      D2KLZZ = - D2ILZZ * D - D2LZZ * ( E + 1.0_DP ) - DIL(3) * DTL(3)

      ! . Calculate the JJ diagonal terms.
      D2JXX = D2IJXX * ( D - 1.0_DP ) + D2JLXX * E + DDJ(1) * DTI(1) + DEJ(1) * DTL(1)
      D2JYY = D2IJYY * ( D - 1.0_DP ) + D2JLYY * E + DDJ(2) * DTI(2) + DEJ(2) * DTL(2)
      D2JZZ = D2IJZZ * ( D - 1.0_DP ) + D2JLZZ * E + DDJ(3) * DTI(3) + DEJ(3) * DTL(3)
      D2JXY = D2IJXY * ( D - 1.0_DP ) + D2JLYX * E + DDJ(2) * DTI(1) + DEJ(2) * DTL(1)
      D2JXZ = D2IJXZ * ( D - 1.0_DP ) + D2JLZX * E + DDJ(3) * DTI(1) + DEJ(3) * DTL(1)
      D2JYZ = D2IJYZ * ( D - 1.0_DP ) + D2JLZY * E + DDJ(3) * DTI(2) + DEJ(3) * DTL(2)

      ! . Calculate the KK diagonal terms.
      D2KXX = - D2KLXX * ( E + 1.0_DP ) - D2IKXX * D - DDK(1) * DTI(1) - DEK(1) * DTL(1)
      D2KYY = - D2KLYY * ( E + 1.0_DP ) - D2IKYY * D - DDK(2) * DTI(2) - DEK(2) * DTL(2)
      D2KZZ = - D2KLZZ * ( E + 1.0_DP ) - D2IKZZ * D - DDK(3) * DTI(3) - DEK(3) * DTL(3)
      D2KXY = - D2KLXY * ( E + 1.0_DP ) - D2IKYX * D - DDK(1) * DTI(2) - DEK(1) * DTL(2)
      D2KXZ = - D2KLXZ * ( E + 1.0_DP ) - D2IKZX * D - DDK(1) * DTI(3) - DEK(1) * DTL(3)
      D2KYZ = - D2KLYZ * ( E + 1.0_DP ) - D2IKZY * D - DDK(2) * DTI(3) - DEK(2) * DTL(3)
         
      ! . Calculate the JK off-diagonal terms.
      D2JKXX = - D2IJXX * D - D2JLXX * ( E + 1.0_DP ) - DDJ(1) * DTI(1) - DEJ(1) * DTL(1)
      D2JKXY = - D2IJYX * D - D2JLXY * ( E + 1.0_DP ) - DDJ(1) * DTI(2) - DEJ(1) * DTL(2)
      D2JKXZ = - D2IJZX * D - D2JLXZ * ( E + 1.0_DP ) - DDJ(1) * DTI(3) - DEJ(1) * DTL(3)
      D2JKYX = - D2IJXY * D - D2JLYX * ( E + 1.0_DP ) - DDJ(2) * DTI(1) - DEJ(2) * DTL(1)
      D2JKYY = - D2IJYY * D - D2JLYY * ( E + 1.0_DP ) - DDJ(2) * DTI(2) - DEJ(2) * DTL(2)
      D2JKYZ = - D2IJZY * D - D2JLYZ * ( E + 1.0_DP ) - DDJ(2) * DTI(3) - DEJ(2) * DTL(3)
      D2JKZX = - D2IJXZ * D - D2JLZX * ( E + 1.0_DP ) - DDJ(3) * DTI(1) - DEJ(3) * DTL(1)
      D2JKZY = - D2IJYZ * D - D2JLZY * ( E + 1.0_DP ) - DDJ(3) * DTI(2) - DEJ(3) * DTL(2)
      D2JKZZ = - D2IJZZ * D - D2JLZZ * ( E + 1.0_DP ) - DDJ(3) * DTI(3) - DEJ(3) * DTL(3)

      ! . Calculate some index factors.
      IFAC = 3 * ( ATMIND(I) - 1 )
      JFAC = 3 * ( ATMIND(J) - 1 )
      KFAC = 3 * ( ATMIND(K) - 1 )
      LFAC = 3 * ( ATMIND(L) - 1 )

      ! . Calculate the diagonal blocks of the hessian.
      CALL DIAGONAL_BLOCK ( IFAC, D2IXX, D2IXY, D2IXZ, D2IYY, D2IYZ, D2IZZ, DTI )
      CALL DIAGONAL_BLOCK ( JFAC, D2JXX, D2JXY, D2JXZ, D2JYY, D2JYZ, D2JZZ, DTJ )
      CALL DIAGONAL_BLOCK ( KFAC, D2KXX, D2KXY, D2KXZ, D2KYY, D2KYZ, D2KZZ, DTK )
      CALL DIAGONAL_BLOCK ( LFAC, D2LXX, D2LXY, D2LXZ, D2LYY, D2LYZ, D2LZZ, DTL )

      ! . Calculate the off-diagonal blocks of the hessian.
      CALL OFF_DIAGONAL_BLOCK ( IFAC, JFAC, D2IJXX, D2IJXY, D2IJXZ, D2IJYX, D2IJYY, D2IJYZ, D2IJZX, D2IJZY, D2IJZZ, DTI, DTJ )
      CALL OFF_DIAGONAL_BLOCK ( IFAC, KFAC, D2IKXX, D2IKXY, D2IKXZ, D2IKYX, D2IKYY, D2IKYZ, D2IKZX, D2IKZY, D2IKZZ, DTI, DTK )
      CALL OFF_DIAGONAL_BLOCK ( IFAC, LFAC, D2ILXX, D2ILXY, D2ILXZ, D2ILYX, D2ILYY, D2ILYZ, D2ILZX, D2ILZY, D2ILZZ, DTI, DTL )
      CALL OFF_DIAGONAL_BLOCK ( JFAC, KFAC, D2JKXX, D2JKXY, D2JKXZ, D2JKYX, D2JKYY, D2JKYZ, D2JKZX, D2JKZY, D2JKZZ, DTJ, DTK )
      CALL OFF_DIAGONAL_BLOCK ( JFAC, LFAC, D2JLXX, D2JLXY, D2JLXZ, D2JLYX, D2JLYY, D2JLYZ, D2JLZX, D2JLZY, D2JLZZ, DTJ, DTL )
      CALL OFF_DIAGONAL_BLOCK ( KFAC, LFAC, D2KLXX, D2KLXY, D2KLXZ, D2KLYX, D2KLYY, D2KLYZ, D2KLZX, D2KLZY, D2KLZZ, DTK, DTL )

#ifdef MMech
   ENDIF
#endif

   END DO

   !=====================================================================
   CONTAINS
   !=====================================================================

      !--------------------------------------------------------------------------------
      SUBROUTINE DIAGONAL_BLOCK ( IFAC, D2IXX, D2IXY, D2IXZ, D2IYY, D2IYZ, D2IZZ, DTI )
      !--------------------------------------------------------------------------------

      ! . Scalar arguments.
      INTEGER,            INTENT(IN) :: IFAC
      REAL ( KIND = DP ), INTENT(IN) :: D2IXX, D2IXY, D2IXZ, D2IYY, D2IYZ, D2IZZ

      ! . Array arguments.
      REAL ( KIND = DP ), DIMENSION(1:3), INTENT(IN) :: DTI

      ! . Local scalars.
      INTEGER :: INDEX

      ! . Add in the elements to the diagonal block.
      IF ( IFAC >= 0 ) THEN
	 INDEX = ( IFAC * ( IFAC + 1 ) ) / 2 + IFAC + 1
	 HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTI(1) + D2IXX
	 INDEX = INDEX + IFAC + 1
	 HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTI(2) + D2IXY
	 HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTI(2) + D2IYY
	 INDEX = INDEX + IFAC + 2
	 HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTI(3) + D2IXZ
	 HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTI(3) + D2IYZ
	 HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(3) * DTI(3) + D2IZZ
      END IF

      END SUBROUTINE DIAGONAL_BLOCK

      !-----------------------------------------------------------------------------
      SUBROUTINE OFF_DIAGONAL_BLOCK ( IFAC, JFAC, D2IJXX, D2IJXY, D2IJXZ, &
                                                  D2IJYX, D2IJYY, D2IJYZ, &
                                                  D2IJZX, D2IJZY, D2IJZZ, DTI, DTJ )
      !-----------------------------------------------------------------------------

      ! . Scalar arguments.
      INTEGER,            INTENT(IN) :: IFAC, JFAC
      REAL ( KIND = DP ), INTENT(IN) :: D2IJXX, D2IJXY, D2IJXZ, &
                                        D2IJYX, D2IJYY, D2IJYZ, &
                                        D2IJZX, D2IJZY, D2IJZZ

      ! . Array arguments.
      REAL ( KIND = DP ), DIMENSION(1:3), INTENT(IN) :: DTI, DTJ

      ! . Local scalars.
      INTEGER :: INDEX

      ! . Add in the elements to the off-diagonal block.
      IF ( ( IFAC >= 0 ) .AND. ( JFAC >= 0 ) ) THEN
	 IF ( IFAC > JFAC ) THEN
            INDEX = ( IFAC * ( IFAC + 1 ) ) / 2 + JFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTJ(1) + D2IJXX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(1) * DTJ(2) + D2IJXY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(1) * DTJ(3) + D2IJXZ
            INDEX = INDEX + IFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(2) * DTJ(1) + D2IJYX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTJ(2) + D2IJYY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(2) * DTJ(3) + D2IJYZ
            INDEX = INDEX + IFAC + 2
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(3) * DTJ(1) + D2IJZX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(3) * DTJ(2) + D2IJZY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(3) * DTJ(3) + D2IJZZ
	 ELSE
            INDEX = ( JFAC * ( JFAC + 1 ) ) / 2 + IFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTJ(1) + D2IJXX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTJ(1) + D2IJYX
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(3) * DTJ(1) + D2IJZX
            INDEX = INDEX + JFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTJ(2) + D2IJXY
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTJ(2) + D2IJYY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(3) * DTJ(2) + D2IJZY
            INDEX = INDEX + JFAC + 2
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + DDF * DTI(1) * DTJ(3) + D2IJXZ
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + DDF * DTI(2) * DTJ(3) + D2IJYZ
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + DDF * DTI(3) * DTJ(3) + D2IJZZ
	 END IF
      END IF

      END SUBROUTINE OFF_DIAGONAL_BLOCK

   END SUBROUTINE DIHEDRAL_ENERGY

END MODULE ENERGY_COVALENT
