      SUBROUTINE GetBigInt(N,X,A)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      REAL*8 A(N)
      X=0.0D0
      DO I=1,N
         X=MAX(X,ABS(A(I)))
      ENDDO
      X=DSQRT(X)
      RETURN
      END
