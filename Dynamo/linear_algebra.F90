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
!                           The Linear Algebra Module
!===============================================================================
!
! . Functions:
!
!   CROSS_PRODUCT                Calculate the cross-product of two 3-vectors.
!   DETERMINANT                  Calculate the determinant of a 3x3 matrix.
!   NORM                         The norm of a vector.
!   NORMALIZE                    Normalize a vector.
!
! . Subroutines:
!
!   PROJECT_OUT                  Make a vector orthogonal to a set of other ones.
!   SCHMIDT_ORTHOGONALIZE        Orthogonalize a set of vectors.  
!
!===============================================================================
MODULE LINEAR_ALGEBRA

! . Module declarations.
USE DEFINITIONS, ONLY : DP

IMPLICIT NONE
PUBLIC
PRIVATE :: NORMALIZATION_TOLERANCE
SAVE

! . Module parameters.
REAL ( KIND = DP ), PARAMETER :: NORMALIZATION_TOLERANCE = 1.0E-16_DP

!==============================================================================
CONTAINS
!==============================================================================

   !------------------------------
   FUNCTION CROSS_PRODUCT ( A, B )
   !------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:3) :: CROSS_PRODUCT

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(IN) :: A, B

   ! . Calculate the cross product A x B.
   CROSS_PRODUCT(1) = A(2) * B(3) - A(3) * B(2)
   CROSS_PRODUCT(2) = A(3) * B(1) - A(1) * B(3)
   CROSS_PRODUCT(3) = A(1) * B(2) - A(2) * B(1)

   END FUNCTION CROSS_PRODUCT

   !-------------------------
   FUNCTION DETERMINANT ( A )
   !-------------------------

   ! . Function declarations.
   REAL ( KIND = DP ) :: DETERMINANT

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:3), INTENT(IN) :: A

   ! . Calculate the determinant of A.
   DETERMINANT = A(1,1) * ( A(2,2) * A(3,3) - A(2,3) * A(3,2) ) - &
                 A(1,2) * ( A(2,1) * A(3,3) - A(3,1) * A(2,3) ) + &
                 A(1,3) * ( A(2,1) * A(3,2) - A(3,1) * A(2,2) )

   END FUNCTION DETERMINANT

   !-----------------------
   FUNCTION NORM ( VECTOR )
   !-----------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: VECTOR

   ! . Function declarations.
   REAL ( KIND = DP ) :: NORM

   ! . Local scalars.
   INTEGER :: N

   ! . Initialization.
   NORM = 0.0_DP

   ! . Determine and check the size of the vector.
   N = SIZE ( VECTOR )
   IF ( N <= 0 ) RETURN

   ! . Get the size of the vector.
   NORM = SQRT ( DOT_PRODUCT ( VECTOR(1:N), VECTOR(1:N) ) )

   ! . If NORM is too small make it zero.
   IF ( NORM < NORMALIZATION_TOLERANCE ) NORM = 0.0_DP

   END FUNCTION NORM

   !----------------------------
   FUNCTION NORMALIZE ( VECTOR )
   !----------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: VECTOR

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:SIZE(VECTOR)) :: NORMALIZE

   ! . Local scalars.
   INTEGER            :: N
   REAL ( KIND = DP ) :: FACT

   ! . Determine and check the size of the vector.
   N = SIZE ( VECTOR )
   IF ( N <= 0 ) RETURN

   ! . Get the NORM of the vector.
   FACT = NORM ( VECTOR )

   ! . Normalize the vector.
   IF ( FACT > 0.0_DP ) NORMALIZE = VECTOR / FACT

   END FUNCTION NORMALIZE

   !-----------------------------------------
   SUBROUTINE PROJECT_OUT ( VECTOR, VECTORS )
   !-----------------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:),   INTENT(INOUT) :: VECTOR
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN)    :: VECTORS

   ! . Scalar arguments.
   INTEGER            :: IVEC, N, NVEC
   REAL ( KIND = DP ) :: DOTFAC

   ! . Get the dimension of the vectors.
   N    = SIZE ( VECTOR,  1 )
   NVEC = SIZE ( VECTORS, 2 )

   ! . Return if there is dimension mismatch or the second set of vectors is null.
   IF ( ( N /= SIZE ( VECTORS, 1 ) ) .OR. ( NVEC <= 0 ) ) RETURN

   ! . Project out VECTORS from VECTOR.
   DO IVEC = 1,NVEC
      DOTFAC = DOT_PRODUCT ( VECTOR(1:N), VECTORS(1:N,IVEC) )
      VECTOR(1:N) = VECTOR(1:N) - DOTFAC * VECTORS(1:N,IVEC)
   END DO

   END SUBROUTINE PROJECT_OUT

   !---------------------------------------------------------
   SUBROUTINE SCHMIDT_ORTHOGONALIZE ( VECTORS, NINDEPENDENT )
   !---------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(OUT) :: NINDEPENDENT

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT) :: VECTORS

   ! . Scalar arguments.
   INTEGER            :: IVEC, JVEC, N, NVEC
   REAL ( KIND = DP ) :: DOTFAC

   ! . Initialization.
   NINDEPENDENT = 0

   ! . Get the dimensions of the vectors.
   N    = SIZE ( VECTORS, 1 )
   NVEC = SIZE ( VECTORS, 2 )

   ! . Return if there is any null dimension.
   IF ( ( N <= 0 ) .OR. ( NVEC <= 0 ) ) RETURN

   ! . Loop over the vectors.
   DO IVEC = 1,NVEC

      ! . Take out the contributions from vectors of lower index.
      DO JVEC = 1,NINDEPENDENT
         DOTFAC = DOT_PRODUCT ( VECTORS(1:N,IVEC), VECTORS(1:N,JVEC) )
         VECTORS(1:N,IVEC) = VECTORS(1:N,IVEC) - DOTFAC * VECTORS(1:N,JVEC)
      END DO

      ! . Get the norm of the vector.
      DOTFAC = NORM ( VECTORS(1:N,IVEC) )

      ! . Add the vector to the set.
      IF ( DOTFAC > 0.0_DP ) THEN

         ! . Increment the number of independent vectors.
         NINDEPENDENT = NINDEPENDENT + 1

         ! . Store the vector.
         VECTORS(1:N,NINDEPENDENT) = VECTORS(1:N,IVEC) / DOTFAC

      END IF

   END DO

   ! . Zero out the remaining vector components.
   VECTORS(1:N,NINDEPENDENT+1:NVEC) = 0.0_DP

   END SUBROUTINE SCHMIDT_ORTHOGONALIZE

END MODULE LINEAR_ALGEBRA

