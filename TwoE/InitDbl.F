      SUBROUTINE INITDBL(N,A)
      IMPLICIT NONE
      INTEGER N
      REAL*8  A(N)
      REAL*8  B
      INTEGER M,I,I1,I2,I3,I4,I5,I6
      B=0D0
      M=MOD(N,7)
      DO I=1,M
         A(I)=B
      ENDDO
      IF(N.LT.7)RETURN
      DO I=M+1,N,7
         I1=I+1
         I2=I+2
         I3=I+3
         I4=I+4
         I5=I+5
         I6=I+6
         A(I )=B
         A(I1)=B
         A(I2)=B
         A(I3)=B
         A(I4)=B
         A(I5)=B
         A(I6)=B
      ENDDO
      RETURN
      END
      SUBROUTINE DBLAXPY(N,A,B)
      IMPLICIT NONE
      INTEGER  I,I1,I2,I3,I4,I5,I6,M,N
      REAL*8   A(N),B(N)
      M=MOD(N,7)
      DO I=1,M
         A(I)=A(I)+B(I)
      ENDDO
      IF(N.LT.7)RETURN
      DO I=M+1,N,7
         I1=I+1
         I2=I+2
         I3=I+3
         I4=I+4
         I5=I+5
         I6=I+6
         A(I )=A(I )+B(I)
         A(I1)=A(I1)+B(I)
         A(I2)=A(I2)+B(I)
         A(I3)=A(I3)+B(I)
         A(I4)=A(I4)+B(I)
         A(I5)=A(I5)+B(I)
         A(I6)=A(I6)+B(I)
      ENDDO
      RETURN
      END 

