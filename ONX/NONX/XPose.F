      SUBROUTINE XPose1C(M,N,A,B)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER I,J,M,N
      REAL*8 A(M,N)
      REAL*8 B(N,M)
      DO I=1,N
         DO J=I+1,M
            B(I,J)=A(J,I)
         ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE XPose2C(M,N,A,B)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER I,J,M,N
      REAL*8 A(M,N)
      REAL*8 B(N,M)
      DO I=1,N
         DO J=1,M
            B(I,J)=A(J,I)
         ENDDO
      ENDDO
      RETURN
      END

