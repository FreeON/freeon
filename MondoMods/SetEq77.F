C===================================================================
C     Transform a BCSR matrix to a dense matrix
C===================================================================
      SUBROUTINE BCSR_TO_DENS(NBasF,NAtoms,NBlks,NNon0,MSiz,OffS,
     >                        A,MTrix,RowPt,ColPt,BlkPt)
      IMPLICIT NONE
      INTEGER NBasF,NAtoms,NBlks,NNon0
      REAL*8  A(NBasF,NBasF),MTrix (NNon0)
      INTEGER RowPt(NAtoms+1),ColPt(NBlks),BlkPt(NBlks), 
     >        MSiz(NAtoms),OffS(NAtoms),
     >        I,J,K,L,P,S,II,JJ,JP,MA,NA
      DO 10 I=1,NBasF
         CALL DBL_VECT_EQ_DBL_SCLR(NBasF,A(1,I),0.0D0)
 10   CONTINUE
      DO 100 I=1,NAtoms
         MA=MSiz(I)
         II=OffS(I)
         DO 200 JP=RowPt(I),RowPt(I+1)-1
            J=ColPt(JP)
            P=BlkPt(JP)
            NA=MSiz(J)
            JJ=OffS(J)
            DO 300 K=0,NA-1
               DO 400 L=0,MA-1
                  A(II+L,JJ+K)=MTrix(P+K*MA+L)
 400           CONTINUE
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C==========================================================
C     Add a double vector to a double vector
C==========================================================
      SUBROUTINE DBL_VECT_PLS_DBL_VECT(N,A,B)
      IMPLICIT NONE
      INTEGER N
      REAL*8  A(N),B(N)
      INTEGER M,I,I1,I2,I3,I4,I5,I6
C----------------------------------------------
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
         A(I )=A(I )+B(I )
         A(I1)=A(I1)+B(I1)
         A(I2)=A(I2)+B(I2)
         A(I3)=A(I3)+B(I3)
         A(I4)=A(I4)+B(I4)
         A(I5)=A(I5)+B(I5)
         A(I6)=A(I6)+B(I6)
      ENDDO
      RETURN
      END
C==========================================================
C     Copy a double vector to a double vector
C==========================================================
      SUBROUTINE DBL_VECT_EQ_DBL_VECT(N,A,B)
      IMPLICIT NONE
      INTEGER N
      REAL*8  A(N),B(N)
      INTEGER M,I,I1,I2,I3,I4,I5,I6
C----------------------------------------------
      M=MOD(N,7)
      DO I=1,M
         A(I)=B(I)
      ENDDO
      IF(N.LT.7)RETURN
      DO I=M+1,N,7
         I1=I+1
         I2=I+2
         I3=I+3
         I4=I+4
         I5=I+5
         I6=I+6
         A(I )=B(I )
         A(I1)=B(I1)
         A(I2)=B(I2)
         A(I3)=B(I3)
         A(I4)=B(I4)
         A(I5)=B(I5)
         A(I6)=B(I6)
      ENDDO
      RETURN
      END
C==========================================================
C     Copy an integer vector to an integer vector
C==========================================================
      SUBROUTINE INT_VECT_EQ_INT_VECT(N,A,B)
      IMPLICIT NONE
      INTEGER N
      INTEGER A(N),B(N)
      INTEGER M,I,I1,I2,I3,I4,I5,I6
C----------------------------------------------
      M=MOD(N,7)
      DO I=1,M
         A(I)=B(I)
      ENDDO
      IF(N.LT.7)RETURN
      DO I=M+1,N,7
         I1=I+1
         I2=I+2
         I3=I+3
         I4=I+4
         I5=I+5
         I6=I+6
         A(I )=B(I )
         A(I1)=B(I1)
         A(I2)=B(I2)
         A(I3)=B(I3)
         A(I4)=B(I4)
         A(I5)=B(I5)
         A(I6)=B(I6)
      ENDDO
      RETURN
      END
C==========================================================
C     Initialize a double vector with a double scalar
C==========================================================
      SUBROUTINE DBL_VECT_EQ_DBL_SCLR(N,A,B)
      IMPLICIT NONE
      INTEGER N
      REAL*8  A(N),B
      INTEGER M,I,I1,I2,I3,I4,I5,I6
C----------------------------------------------
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
C==========================================================
C     Initialize an integer vector with an integer scalar
C==========================================================
      SUBROUTINE INT_VECT_EQ_INT_SCLR(N,A,B)
      IMPLICIT NONE
      INTEGER N
      INTEGER A(N),B
      INTEGER M,I,I1,I2,I3,I4,I5,I6
C----------------------------------------------
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


