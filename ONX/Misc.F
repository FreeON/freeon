      INTEGER FUNCTION LHGTF(L)
      INTEGER L
      LHGTF=(L+1)*(L+2)*(L+3)/6
      RETURN
      END
C
      INTEGER FUNCTION IStride(I,J)
      INTEGER I,J
      IStride=LHGTF(J)-LHGTF(I-1)
      RETURN
      END
C
      INTEGER FUNCTION LBegin(L)
      INTEGER L
      LBegin=(L*(L+1)*(L+2))/6+1
      RETURN
      END
C
      INTEGER FUNCTION LEnd(L)
      INTEGER L
      LEnd=LBegin(L)+L*(L+1)/2+L
      RETURN
      END
c
      INTEGER FUNCTION LMNDex(L,M,N)
      INTEGER L,M,N
      LMNDex=LBegin(L+M+N)+N*(2*(L+M+N)-N+3)/2+M ! Q-Chem ordering
      RETURN
      END

