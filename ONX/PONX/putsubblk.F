      SUBROUTINE PutSubBlk(Ind,N,NAR,NAC,NBR,NBC,RS,CS,A,B)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      REAL*8 A(NAR,NAC)
      REAL*8 B(N,NBR,NBC)
      INTEGER NAR,NAC
      INTEGER NBR,NBC,RS,CS
      INTEGER I,J,RS1,CS1,Ind,N
      RS1=RS-1
      CS1=CS-1
      DO J=1,NBC
         DO I=1,NBR
            A(I+RS1,J+CS1)=A(I+RS1,J+CS1)+B(Ind,I,J)
         ENDDO
      ENDDO
      RETURN
      END
