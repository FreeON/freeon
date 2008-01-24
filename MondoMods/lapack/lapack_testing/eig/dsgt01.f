      SUBROUTINE DSGT01( ITYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D,
     $                   WORK, RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     modified August 1997, a new parameter M is added to the calling
*     sequence.
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            ITYPE, LDA, LDB, LDZ, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), D( * ), RESULT( * ),
     $                   WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DDGT01 checks a decomposition of the form
*
*     A Z   =  B Z D or
*     A B Z =  Z D or
*     B A Z =  Z D
*
*  where A is a symmetric matrix, B is
*  symmetric positive definite, Z is orthogonal, and D is diagonal.
*
*  One of the following test ratios is computed:
*
*  ITYPE = 1:  RESULT(1) = | A Z - B Z D | / ( |A| |Z| n ulp )
*
*  ITYPE = 2:  RESULT(1) = | A B Z - Z D | / ( |A| |Z| n ulp )
*
*  ITYPE = 3:  RESULT(1) = | B A Z - Z D | / ( |A| |Z| n ulp )
*
*  Arguments
*  =========
*
*  ITYPE   (input) INTEGER
*          The form of the symmetric generalized eigenproblem.
*          = 1:  A*z = (lambda)*B*z
*          = 2:  A*B*z = (lambda)*z
*          = 3:  B*A*z = (lambda)*z
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrices A and B is stored.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  M       (input) INTEGER
*          The number of eigenvalues found.  0 <= M <= N.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA, N)
*          The original symmetric matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB, N)
*          The original symmetric positive definite matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  Z       (input) DOUBLE PRECISION array, dimension (LDZ, M)
*          The computed eigenvectors of the generalized eigenproblem.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= max(1,N).
*
*  D       (input) DOUBLE PRECISION array, dimension (M)
*          The computed eigenvalues of the generalized eigenproblem.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N*N)
*
*  RESULT  (output) DOUBLE PRECISION array, dimension (1)
*          The test ratio as described above.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   ANORM, ULP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANSY
      EXTERNAL           DLAMCH, DLANGE, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSYMM
*     ..
*     .. Executable Statements ..
*
      RESULT( 1 ) = ZERO
      IF( N.LE.0 )
     $   RETURN
*
      ULP = DLAMCH( 'Epsilon' )
*
*     Compute product of 1-norms of A and Z.
*
      ANORM = DLANSY( '1', UPLO, N, A, LDA, WORK )*
     $        DLANGE( '1', N, M, Z, LDZ, WORK )
      IF( ANORM.EQ.ZERO )
     $   ANORM = ONE
*
      IF( ITYPE.EQ.1 ) THEN
*
*        Norm of AZ - BZD
*
         CALL DSYMM( 'Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO,
     $               WORK, N )
         DO 10 I = 1, M
            CALL DSCAL( N, D( I ), Z( 1, I ), 1 )
   10    CONTINUE
         CALL DSYMM( 'Left', UPLO, N, M, ONE, B, LDB, Z, LDZ, -ONE,
     $               WORK, N )
*
         RESULT( 1 ) = ( DLANGE( '1', N, M, WORK, N, WORK ) / ANORM ) /
     $                 ( N*ULP )
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        Norm of ABZ - ZD
*
         CALL DSYMM( 'Left', UPLO, N, M, ONE, B, LDB, Z, LDZ, ZERO,
     $               WORK, N )
         DO 20 I = 1, M
            CALL DSCAL( N, D( I ), Z( 1, I ), 1 )
   20    CONTINUE
         CALL DSYMM( 'Left', UPLO, N, M, ONE, A, LDA, WORK, N, -ONE, Z,
     $               LDZ )
*
         RESULT( 1 ) = ( DLANGE( '1', N, M, Z, LDZ, WORK ) / ANORM ) /
     $                 ( N*ULP )
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        Norm of BAZ - ZD
*
         CALL DSYMM( 'Left', UPLO, N, M, ONE, A, LDA, Z, LDZ, ZERO,
     $               WORK, N )
         DO 30 I = 1, M
            CALL DSCAL( N, D( I ), Z( 1, I ), 1 )
   30    CONTINUE
         CALL DSYMM( 'Left', UPLO, N, M, ONE, B, LDB, WORK, N, -ONE, Z,
     $               LDZ )
*
         RESULT( 1 ) = ( DLANGE( '1', N, M, Z, LDZ, WORK ) / ANORM ) /
     $                 ( N*ULP )
      END IF
*
      RETURN
*
*     End of DDGT01
*
      END
