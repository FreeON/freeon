      SUBROUTINE BuffTbl(MaxK,BLength)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      INTEGER BLength(MaxK,MaxK)
      INTEGER MaxK
      DO I=1,MaxK
         DO J=1,I
            MX1=MXBUF1/(I*J)
            MX2=MXBUF2/I
            MX3=MXBATCH
            MXBT=MIN(MX1,MX2)
            MXBT=MIN(MX3,MXBT)
            BLength(I,J)=MXBT
            BLength(J,I)=MXBT
         ENDDO
      ENDDO
      RETURN
      END 
