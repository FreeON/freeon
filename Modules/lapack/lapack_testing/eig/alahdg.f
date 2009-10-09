      SUBROUTINE ALAHDG( IOUNIT, PATH )
*
*  -- LAPACK test routine (version 3.1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER*3       PATH
      INTEGER           IOUNIT
*     ..
*
*  Purpose
*  =======
*
*  ALAHDG prints header information for the different test paths.
*
*  Arguments
*  =========
*
*  IOUNIT  (input) INTEGER
*          The unit number to which the header information should be
*          printed.
*
*  PATH    (input) CHARACTER*3
*          The name of the path for which the header information is to
*          be printed.  Current paths are
*             GQR:  GQR (general matrices)
*             GRQ:  GRQ (general matrices)
*             LSE:  LSE Problem
*             GLM:  GLM Problem
*             GSV:  Generalized Singular Value Decomposition
*
*  =====================================================================
*
*     .. Local Scalars ..
      CHARACTER*3       C2
      INTEGER           ITYPE
*     ..
*     .. External Functions ..
      LOGICAL           LSAMEN
      EXTERNAL          LSAMEN
*     ..
*     .. Executable Statements ..
*
      IF( IOUNIT.LE.0 )
     $   RETURN
      C2 = PATH( 1: 3 )
*
*     First line describing matrices in this path
*
      IF( LSAMEN( 3, C2, 'GQR' ) ) THEN
         ITYPE = 1
         WRITE( IOUNIT, FMT = 9991 )PATH
      ELSE IF( LSAMEN( 3, C2, 'GRQ' ) ) THEN
         ITYPE = 2
         WRITE( IOUNIT, FMT = 9992 )PATH
      ELSE IF( LSAMEN( 3, C2, 'LSE' ) ) THEN
         ITYPE = 3
         WRITE( IOUNIT, FMT = 9993 )PATH
      ELSE IF( LSAMEN( 3, C2, 'GLM' ) ) THEN
         ITYPE = 4
         WRITE( IOUNIT, FMT = 9994 )PATH
      ELSE IF( LSAMEN( 3, C2, 'GSV' ) ) THEN
         ITYPE = 5
         WRITE( IOUNIT, FMT = 9995 )PATH
      END IF
*
*     Matrix types
*
      WRITE( IOUNIT, FMT = 9999 )'Matrix types: '
*
      IF( ITYPE.EQ.1 )THEN
         WRITE( IOUNIT, FMT = 9950 )1
         WRITE( IOUNIT, FMT = 9952 )2
         WRITE( IOUNIT, FMT = 9954 )3
         WRITE( IOUNIT, FMT = 9955 )4
         WRITE( IOUNIT, FMT = 9956 )5
         WRITE( IOUNIT, FMT = 9957 )6
         WRITE( IOUNIT, FMT = 9961 )7
         WRITE( IOUNIT, FMT = 9962 )8
      ELSE IF( ITYPE.EQ.2 )THEN
         WRITE( IOUNIT, FMT = 9951 )1
         WRITE( IOUNIT, FMT = 9953 )2
         WRITE( IOUNIT, FMT = 9954 )3
         WRITE( IOUNIT, FMT = 9955 )4
         WRITE( IOUNIT, FMT = 9956 )5
         WRITE( IOUNIT, FMT = 9957 )6
         WRITE( IOUNIT, FMT = 9961 )7
         WRITE( IOUNIT, FMT = 9962 )8
      ELSE IF( ITYPE.EQ.3 )THEN
         WRITE( IOUNIT, FMT = 9950 )1
         WRITE( IOUNIT, FMT = 9952 )2
         WRITE( IOUNIT, FMT = 9954 )3
         WRITE( IOUNIT, FMT = 9955 )4
         WRITE( IOUNIT, FMT = 9955 )5
         WRITE( IOUNIT, FMT = 9955 )6
         WRITE( IOUNIT, FMT = 9955 )7
         WRITE( IOUNIT, FMT = 9955 )8
      ELSE IF( ITYPE.EQ.4 )THEN
         WRITE( IOUNIT, FMT = 9951 )1
         WRITE( IOUNIT, FMT = 9953 )2
         WRITE( IOUNIT, FMT = 9954 )3
         WRITE( IOUNIT, FMT = 9955 )4
         WRITE( IOUNIT, FMT = 9955 )5
         WRITE( IOUNIT, FMT = 9955 )6
         WRITE( IOUNIT, FMT = 9955 )7
         WRITE( IOUNIT, FMT = 9955 )8
      ELSE IF( ITYPE.EQ.5 )THEN
         WRITE( IOUNIT, FMT = 9950 )1
         WRITE( IOUNIT, FMT = 9952 )2
         WRITE( IOUNIT, FMT = 9954 )3
         WRITE( IOUNIT, FMT = 9955 )4
         WRITE( IOUNIT, FMT = 9956 )5
         WRITE( IOUNIT, FMT = 9957 )6
         WRITE( IOUNIT, FMT = 9959 )7
         WRITE( IOUNIT, FMT = 9960 )8
      END IF
*
*     Tests performed
*
      WRITE( IOUNIT, FMT = 9999 )'Test ratios: '
*
      IF( ITYPE.EQ.1 ) THEN
*
*        GQR decomposition of rectangular matrices
*
         WRITE( IOUNIT, FMT = 9930 )1
         WRITE( IOUNIT, FMT = 9931 )2
         WRITE( IOUNIT, FMT = 9932 )3
         WRITE( IOUNIT, FMT = 9933 )4
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        GRQ decomposition of rectangular matrices
*
         WRITE( IOUNIT, FMT = 9934 )1
         WRITE( IOUNIT, FMT = 9935 )2
         WRITE( IOUNIT, FMT = 9932 )3
         WRITE( IOUNIT, FMT = 9933 )4
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        LSE Problem
*
         WRITE( IOUNIT, FMT = 9937 )1
         WRITE( IOUNIT, FMT = 9938 )2
      ELSE IF( ITYPE.EQ.4 ) THEN
*
*        GLM Problem
*
         WRITE( IOUNIT, FMT = 9939 )1
      ELSE IF( ITYPE.EQ.5 ) THEN
*
*        GSVD
*
         WRITE( IOUNIT, FMT = 9940 )1
         WRITE( IOUNIT, FMT = 9941 )2
         WRITE( IOUNIT, FMT = 9942 )3
         WRITE( IOUNIT, FMT = 9943 )4
         WRITE( IOUNIT, FMT = 9944 )5
      END IF
*
 9999 FORMAT( 1X, A )
 9991 FORMAT( / 1X, A3, ': GQR factorization of general matrices' )
 9992 FORMAT( / 1X, A3, ': GRQ factorization of general matrices' )
 9993 FORMAT( / 1X, A3, ': LSE Problem' )
 9994 FORMAT( / 1X, A3, ': GLM Problem' )
 9995 FORMAT( / 1X, A3, ': Generalized Singular Value Decomposition' )
*
 9950 FORMAT( 3X, I2, ': A-diagonal matrix  B-upper triangular' )
 9951 FORMAT( 3X, I2, ': A-diagonal matrix  B-lower triangular' )
 9952 FORMAT( 3X, I2, ': A-upper triangular B-upper triangular' )
 9953 FORMAT( 3X, I2, ': A-lower triangular B-diagonal triangular' )
 9954 FORMAT( 3X, I2, ': A-lower triangular B-upper triangular' )
*
 9955 FORMAT( 3X, I2, ': Random matrices cond(A)=100, cond(B)=10,' )
*
 9956 FORMAT( 3X, I2, ': Random matrices cond(A)= sqrt( 0.1/EPS ) ',
     $      'cond(B)= sqrt( 0.1/EPS )' )
 9957 FORMAT( 3X, I2, ': Random matrices cond(A)= 0.1/EPS ',
     $      'cond(B)= 0.1/EPS' )
 9959 FORMAT( 3X, I2, ': Random matrices cond(A)= sqrt( 0.1/EPS ) ',
     $      'cond(B)=  0.1/EPS ' )
 9960 FORMAT( 3X, I2, ': Random matrices cond(A)= 0.1/EPS ',
     $      'cond(B)=  sqrt( 0.1/EPS )' )
*
 9961 FORMAT( 3X, I2, ': Matrix scaled near underflow limit' )
 9962 FORMAT( 3X, I2, ': Matrix scaled near overflow limit' )
*
*
*     GQR test ratio
*
 9930 FORMAT( 3X, I2, ': norm( R - Q'' * A ) / ( min( N, M )*norm( A )',
     $       '* EPS )' )
 9931 FORMAT( 3X, I2, ': norm( T * Z - Q'' * B )  / ( min(P,N)*norm(B)',
     $       '* EPS )' )
 9932 FORMAT( 3X, I2, ': norm( I - Q''*Q )   / ( N * EPS )' )
 9933 FORMAT( 3X, I2, ': norm( I - Z''*Z )   / ( P * EPS )' )
*
*     GRQ test ratio
*
 9934 FORMAT( 3X, I2, ': norm( R - A * Q'' ) / ( min( N,M )*norm(A) * ',
     $       'EPS )' )
 9935 FORMAT( 3X, I2, ': norm( T * Q - Z'' * B )  / ( min( P,N ) * nor',
     $       'm(B)*EPS )' )
*
*     LSE test ratio
*
 9937 FORMAT( 3X, I2, ': norm( A*x - c )  / ( norm(A)*norm(x) * EPS )' )
 9938 FORMAT( 3X, I2, ': norm( B*x - d )  / ( norm(B)*norm(x) * EPS )' )
*
*     GLM test ratio
*
 9939 FORMAT( 3X, I2, ': norm( d - A*x - B*y ) / ( (norm(A)+norm(B) )*',
     $       '(norm(x)+norm(y))*EPS )' )
*
*     GSVD test ratio
*
 9940 FORMAT( 3X, I2, ': norm( U'' * A * Q - D1 * R ) / ( min( M, N )*',
     $       'norm( A ) * EPS )' )
 9941 FORMAT( 3X, I2, ': norm( V'' * B * Q - D2 * R ) / ( min( P, N )*',
     $       'norm( B ) * EPS )' )
 9942 FORMAT( 3X, I2, ': norm( I - U''*U )   / ( M * EPS )' )
 9943 FORMAT( 3X, I2, ': norm( I - V''*V )   / ( P * EPS )' )
 9944 FORMAT( 3X, I2, ': norm( I - Q''*Q )   / ( N * EPS )' )
      RETURN
*
*     End of ALAHDG
*
      END
