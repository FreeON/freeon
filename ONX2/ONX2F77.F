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
      SUBROUTINE GetAbsDenBlk2(D,NR,NC,DM,NFR,NFC,IRB,IRE,ICB,ICE)
      IMPLICIT NONE
      INTEGER NR,NC,NFR,NFC
      INTEGER IRB(NFR),IRE(NFR),ICB(NFC),ICE(NFC)
      REAL*8 D(NR,NC),DM(NFR,NFC)
      INTEGER INFC,INFR,I,J,IRB1,ICB1,STDR,STDC
      REAL*8 TMP
      ICB1=0
      DO INFC=1,NFC
         IRB1=0
         STDC=ICE(INFC)-ICB(INFC)+1
         DO INFR=1,NFR
            STDR=IRE(INFR)-IRB(INFR)+1
            TMP=0.0D0
            DO J=ICB1+1,ICB1+STDC
               DO I=IRB1+1,IRB1+STDR
                  TMP=MAX(TMP,ABS(D(I,J)))
               ENDDO
            ENDDO
            DM(INFR,INFC)=TMP
            IRB1=IRB1+STDR
         ENDDO
         ICB1=ICB1+STDC
      ENDDO
      RETURN
      END
C     
      SUBROUTINE GetAbsDenBlk(D,NR,NC,DM,NFR,NFC,IRB,IRE,ICB,ICE)
      IMPLICIT NONE
      INTEGER NR,NC,NFR,NFC
      INTEGER IRB(NFR),IRE(NFR),ICB(NFC),ICE(NFC)
      REAL*8 D(NR,NC),DM(NFC,NFR)
      INTEGER INFC,INFR,I,J,IRB1,ICB1,STDR,STDC
      REAL*8 TMP
      ICB1=0
      DO INFC=1,NFC
         STDC=ICE(INFC)-ICB(INFC)+1
         IRB1=0
         DO INFR=1,NFR
            STDR=IRE(INFR)-IRB(INFR)+1
            TMP=0.0D0
            DO J=ICB1+1,ICB1+STDC
               DO I=IRB1+1,IRB1+STDR
                  TMP=MAX(TMP,ABS(D(I,J)))
               ENDDO
            ENDDO
            DM(INFC,INFR)=TMP
            IRB1=IRB1+STDR
         ENDDO
         ICB1=ICB1+STDC
      ENDDO
      RETURN
      END
C     
      INTEGER FUNCTION IBinSrch(IVec,IVal,NDim)
      IMPLICIT NONE
      INTEGER IVal,NDim
      INTEGER IVec(NDim)
      INTEGER IMin,IMax,IMid
      IMin=0
      IMax=NDim+1
      IBinSrch=-1
      DO WHILE(IMax-IMin.GT.1)
         IMid=(IMax+IMin)/2
         IF(IVal.GT.IVec(IMid)) THEN
            IMin=IMid
         ELSEIF(IVal.LT.IVec(IMid)) THEN
            IMax=IMid
         ELSE
            IBinSrch=IMid
            RETURN
         ENDIF
      ENDDO
      RETURN
      END FUNCTION IBinSrch
C     
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
C     
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
C
