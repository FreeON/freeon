      SUBROUTINE INITDBL(N,A)
      IMPLICIT NONE
      INTEGER N
      REAL*8  A(N)
      REAL*8  B
      INTEGER M,MM,I
      MM=16
      B=0D0
      M=MOD(N,MM)
      DO I=1,M
         A(I)=B
      ENDDO
      IF(N.LT.MM)RETURN
      DO I=M+1,N,MM
         A(I )=B
         A(I+1)=B
         A(I+2)=B
         A(I+3)=B
         A(I+4)=B
         A(I+5)=B
         A(I+6)=B
         A(I+7)=B
         A(I+8)=B
         A(I+9)=B
         A(I+10)=B
         A(I+11)=B
         A(I+12)=B
         A(I+13)=B
         A(I+14)=B
         A(I+15)=B
      ENDDO
      RETURN
      END

      SUBROUTINE DBLAXPY(N,A,B)
      IMPLICIT NONE
      INTEGER  I,MM,M,N
      REAL*8   A(N),B(N)
      MM=16
      M=MOD(N,MM)
      DO I=1,M
         A(I)=A(I)+B(I)
      ENDDO
      IF(N.LT.MM)RETURN
      DO I=M+1,N,MM
         A(I   )=A(I   )+B(I   )
         A(I+1 )=A(I+1 )+B(I+1 )
         A(I+2 )=A(I+2 )+B(I+2 )
         A(I+3 )=A(I+3 )+B(I+3 )
         A(I+4 )=A(I+4 )+B(I+4 )
         A(I+5 )=A(I+5 )+B(I+5 )
         A(I+6 )=A(I+6 )+B(I+6 )
         A(I+7 )=A(I+7 )+B(I+7 )
         A(I+8 )=A(I+8 )+B(I+8 )
         A(I+9 )=A(I+9 )+B(I+9 )
         A(I+10)=A(I+10)+B(I+10)
         A(I+11)=A(I+11)+B(I+11)
         A(I+12)=A(I+12)+B(I+12)
         A(I+13)=A(I+13)+B(I+13)
         A(I+14)=A(I+14)+B(I+14)
         A(I+15)=A(I+15)+B(I+15)
      ENDDO
      RETURN
      END 

      SUBROUTINE DBLAXPZY(N,A,Z,B)
      IMPLICIT NONE
      INTEGER  I,MM,M,N
      REAL*8   Z,A(N),B(N)
      MM=16
      M=MOD(N,MM)
      DO I=1,M
         A(I)=A(I)+Z*B(I)
      ENDDO
      IF(N.LT.MM)RETURN
      DO I=M+1,N,MM
         A(I   )=A(I   )+Z*B(I   )
         A(I+1 )=A(I+1 )+Z*B(I+1 )
         A(I+2 )=A(I+2 )+Z*B(I+2 )
         A(I+3 )=A(I+3 )+Z*B(I+3 )
         A(I+4 )=A(I+4 )+Z*B(I+4 )
         A(I+5 )=A(I+5 )+Z*B(I+5 )
         A(I+6 )=A(I+6 )+Z*B(I+6 )
         A(I+7 )=A(I+7 )+Z*B(I+7 )
         A(I+8 )=A(I+8 )+Z*B(I+8 )
         A(I+9 )=A(I+9 )+Z*B(I+9 )
         A(I+10)=A(I+10)+Z*B(I+10)
         A(I+11)=A(I+11)+Z*B(I+11)
         A(I+12)=A(I+12)+Z*B(I+12)
         A(I+13)=A(I+13)+Z*B(I+13)
         A(I+14)=A(I+14)+Z*B(I+14)
         A(I+15)=A(I+15)+Z*B(I+15)
      ENDDO
      RETURN
      END 

