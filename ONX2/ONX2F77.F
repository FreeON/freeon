C
      FUNCTION DGetAbsMax(N,A)
      INTEGER N
      REAL*8 A(N)
      REAL*8 DGetAbsMax
      INTEGER I
      DGetAbsMax=0.0D0
      DO I=1,N
         DGetAbsMax=MAX(DGetAbsMax,ABS(A(I)))
      ENDDO
      RETURN
      END
C     
