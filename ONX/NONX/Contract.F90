SUBROUTINE Contract(N,KBra,KKet,NVRR,LngDrv,CDrv,CB,CK,C,U)
  USE DerivedTypes
  USE GlobalScalars
  USE ONXParameters
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: N,KBra,KKet,NVRR,LngDrv
  INTEGER, INTENT(IN)        :: CDrv(4,LngDrv)
  REAL(DOUBLE), INTENT(IN)   :: CB(KBra,3)
  REAL(DOUBLE), INTENT(IN)   :: CK(N,KKet,3)
  REAL(DOUBLE), INTENT(OUT)  :: C(N,NVRR)
  REAL(DOUBLE), INTENT(IN)   :: U(N,KBra,KKet,NVRR)
  INTEGER                    :: I,J,K,iC,iP,iQ,iU,IDrv

  IF (KBra==1.AND.KKet==1) THEN
    DO IDrv=1,LngDrv
      iC  = CDrv(1,IDrv)
      iP  = CDrv(2,IDrv)
      iQ  = CDrv(3,IDrv)
      iU  = CDrv(4,IDrv)
      IF (iP==0.AND.iQ==0) THEN
        DO I=1,N
          C(I,iC) = U(I,1,1,iU)
        END DO
      ELSEIF(iP.GT.0.AND.iQ.EQ.0) THEN
        DO I=1,N
          C(I,iC) = U(I,1,1,iU)*CB(1,iP)
        END DO
      ELSEIF(iP.EQ.0.AND.iQ.GT.0) THEN
        DO I=1,N
          C(I,iC) = U(I,1,1,iU)*CK(I,1,iQ)
        END DO
      ELSE
        DO I=1,N
          C(I,iC) = U(I,1,1,iU)*CB(1,iP)*CK(I,1,iQ)
        END DO
      END IF
    END DO ! IDrv
  ELSE
  DO IDrv=1,LngDrv
    iC  = CDrv(1,IDrv)
    iP  = CDrv(2,IDrv)
    iQ  = CDrv(3,IDrv)
    iU  = CDrv(4,IDrv)
    DO I=1,N
      C(I,iC)=0.0D0
    END DO
    IF(iP.EQ.0.AND.iQ.EQ.0) THEN
      DO K=1,KKet
        DO J=1,KBra
          DO I=1,N
            C(I,iC) = C(I,iC) + U(I,J,K,iU)
          END DO
        END DO
      END DO
    ELSEIF(iP.GT.0.AND.iQ.EQ.0) THEN
      DO K=1,KKet
        DO J=1,KBra
          DO I=1,N
            C(I,iC) = C(I,iC) + U(I,J,K,iU)*CB(J,iP)
          END DO
        END DO
      END DO
    ELSEIF(iP.EQ.0.AND.iQ.GT.0) THEN
      DO K=1,KKet
        DO J=1,KBra
          DO I=1,N
            C(I,iC) = C(I,iC) + U(I,J,K,iU)*CK(I,K,iQ)
          END DO
        END DO
      END DO
    ELSE
      DO K=1,KKet
        DO J=1,KBra
          DO I=1,N
            C(I,iC) = C(I,iC) + U(I,J,K,iU)*CB(J,iP)*CK(I,K,iQ)
          END DO
        END DO 
      END DO
    END IF
  END DO
  END IF

END SUBROUTINE Contract


