      SUBROUTINE RLoader(A,N,X)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      REAL*8 A(N),X
      INTEGER N,I
      DO I=1,N
         A(I)=X
      ENDDO
      RETURN
      END
