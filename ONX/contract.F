      SUBROUTINE Contract(N,KBra,KKet,LngDrv,CDrv,CB,CK,C,U)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 C(N,*)
      REAL*8 U(N,KBra,KKet,*)
      REAL*8 CB(N,KBra,*)
      REAL*8 CK(N,KKet,*)
      INTEGER CDrv(4,*)

      IF(KBra.EQ.1.AND.KKet.EQ.1) THEN
      DO IDrv=1,LngDrv
         iC  = CDrv(1,IDrv)
         iP  = CDrv(2,IDrv)
         iQ  = CDrv(3,IDrv)
         iU  = CDrv(4,IDrv)
         IF(iP.EQ.0.AND.iQ.EQ.0) THEN
            DO I=1,N
               C(I,iC) = U(I,1,1,iU)
            ENDDO
         ELSEIF(iP.GT.0.AND.iQ.EQ.0) THEN
            DO I=1,N
               C(I,iC) = U(I,1,1,iU)*CB(I,1,iP)
            ENDDO
         ELSEIF(iP.EQ.0.AND.iQ.GT.0) THEN
            DO I=1,N
               C(I,iC) = U(I,1,1,iU)*CK(I,1,iQ)
            ENDDO
         ELSE
            DO I=1,N
               C(I,iC) = U(I,1,1,iU)*CB(I,1,iP)*CK(I,1,iQ)
            ENDDO
         ENDIF
      ENDDO
      ELSE
      DO IDrv=1,LngDrv
         iC  = CDrv(1,IDrv)
         iP  = CDrv(2,IDrv)
         iQ  = CDrv(3,IDrv)
         iU  = CDrv(4,IDrv)
         DO I=1,N
            C(I,iC)=0.0D0
         ENDDO
         IF(iP.EQ.0.AND.iQ.EQ.0) THEN
            DO K=1,KKet
            DO J=1,KBra
            DO I=1,N
               C(I,iC) = C(I,iC) + U(I,J,K,iU)
            ENDDO
            ENDDO
            ENDDO
         ELSEIF(iP.GT.0.AND.iQ.EQ.0) THEN
            DO K=1,KKet
            DO J=1,KBra
            DO I=1,N
               C(I,iC) = C(I,iC) + U(I,J,K,iU)*CB(I,J,iP)
            ENDDO
            ENDDO
            ENDDO
         ELSEIF(iP.EQ.0.AND.iQ.GT.0) THEN
            DO K=1,KKet
            DO J=1,KBra
            DO I=1,N
               C(I,iC) = C(I,iC) + U(I,J,K,iU)*CK(I,K,iQ)
            ENDDO
            ENDDO
            ENDDO
         ELSE
            DO K=1,KKet
            DO J=1,KBra
            DO I=1,N
               C(I,iC) = C(I,iC) + U(I,J,K,iU)*
     &                   CB(I,J,iP)*CK(I,K,iQ)
            ENDDO
            ENDDO
            ENDDO
         ENDIF
      ENDDO      
      ENDIF
      RETURN
      END
