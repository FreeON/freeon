SUBROUTINE ContractG(N,KBra,KKet,NVRR,CB,CK,C,U,PrmBufB,DB,SB,GD)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE ONXParameters
  IMPLICIT NONE
  TYPE(DBuf), INTENT(IN)     :: DB
  TYPE(DSL), INTENT(IN)      :: SB
  TYPE(GradD), INTENT(IN)    :: GD
  INTEGER, INTENT(IN)        :: N,KBra,KKet,NVRR
  REAL(DOUBLE), INTENT(IN)   :: CB(KBra,3)
  REAL(DOUBLE), INTENT(IN)   :: CK(N,KKet,3)
  REAL(DOUBLE), INTENT(OUT)  :: C(N,NVRR*4)
  REAL(DOUBLE), INTENT(IN)   :: U(N,KBra,KKet,NVRR)
  REAL(DOUBLE), INTENT(IN)   :: PrmBufB(DB%MAXP,KBra+DB%MInfo)
  INTEGER                    :: I,J,K,L,iC,iP,iQ,iU,Ind,Ioff
  REAL(DOUBLE)               :: Za2,Zb2,Zc2,Zd2

  IF (N.EQ.0) RETURN
!----------------------------------------!
! Total contraction length of one        !
!----------------------------------------!
  IF (KBra==1.AND.KKet==1) THEN
    Za2=PrmBufB( 9,1)
    Zc2=PrmBufB(10,1)
!----------------------------------------!
! Contract integrals and scale by 2 * Za |
!----------------------------------------!
    DO L=1,GD%LenDa
      iC  = GD%GDrv1%I(1,L) 
      iP  = GD%GDrv1%I(2,L)
      iQ  = GD%GDrv1%I(3,L)
      iU  = GD%GDrv1%I(4,L)

      IF (iP==0.AND.iQ==0) THEN
        DO I=1,N
          C(I,iC) = U(I,1,1,iU)*Za2
        END DO
      ELSEIF(iP.GT.0.AND.iQ.EQ.0) THEN
        DO I=1,N
          C(I,iC) = U(I,1,1,iU)*CB(1,iP)*Za2
        END DO
      ELSEIF(iP.EQ.0.AND.iQ.GT.0) THEN
        DO I=1,N
          C(I,iC) = U(I,1,1,iU)*CK(I,1,iQ)*Za2
        END DO
      ELSE
        DO I=1,N
          C(I,iC) = U(I,1,1,iU)*CB(1,iP)*CK(I,1,iQ)*Za2
        END DO
      END IF
    END DO ! L
!----------------------------------------!
! Contract integrals and scale by 2 * Zc |
!----------------------------------------!
    DO L=GD%LenDa+1,GD%LenDc
      iC  = GD%GDrv1%I(1,L)
      iP  = GD%GDrv1%I(2,L)
      iQ  = GD%GDrv1%I(3,L)
      iU  = GD%GDrv1%I(4,L)

      IF (iP==0.AND.iQ==0) THEN
        DO I=1,N
          C(I,iC) = U(I,1,1,iU)*Zc2
        END DO
      ELSEIF(iP.GT.0.AND.iQ.EQ.0) THEN
        DO I=1,N
          C(I,iC) = U(I,1,1,iU)*CB(1,iP)*Zc2
        END DO
      ELSEIF(iP.EQ.0.AND.iQ.GT.0) THEN
        DO I=1,N
          C(I,iC) = U(I,1,1,iU)*CK(I,1,iQ)*Zc2
        END DO
      ELSE
        DO I=1,N
          C(I,iC) = U(I,1,1,iU)*CB(1,iP)*CK(I,1,iQ)*Zc2
        END DO
      END IF
    END DO ! L
!----------------------------------------!
! Contract integrals and scale by 2 * Zb |
!----------------------------------------!
    DO L=GD%LenDc+1,GD%LenDb
      iC  = GD%GDrv1%I(1,L)
      iP  = GD%GDrv1%I(2,L)
      iQ  = GD%GDrv1%I(3,L)
      iU  = GD%GDrv1%I(4,L)
      IF (iP==0.AND.iQ==0) THEN
        DO I=1,N
          Ind = SB%SLPrm%I(I)+8
          Zb2 = DB%PrmBuf%D(Ind)
          C(I,iC) = U(I,1,1,iU)*Zb2
        END DO
      ELSEIF(iP.GT.0.AND.iQ.EQ.0) THEN
        DO I=1,N
          Ind = SB%SLPrm%I(I)+8
          Zb2 = DB%PrmBuf%D(Ind)
          C(I,iC) = U(I,1,1,iU)*CB(1,iP)*Zb2
        END DO
      ELSEIF(iP.EQ.0.AND.iQ.GT.0) THEN
        DO I=1,N
          Ind = SB%SLPrm%I(I)+8
          Zb2 = DB%PrmBuf%D(Ind)
          C(I,iC) = U(I,1,1,iU)*CK(I,1,iQ)*Zb2
        END DO
      ELSE
        DO I=1,N
          Ind = SB%SLPrm%I(I)+8
          Zb2 = DB%PrmBuf%D(Ind)
          C(I,iC) = U(I,1,1,iU)*CB(1,iP)*CK(I,1,iQ)*Zb2
        END DO
      END IF
    END DO ! L
!----------------------------------------!
! Contract integrals (no scale factors)  |
!----------------------------------------!
    DO L=GD%LenDb+1,GD%LenDn
      iC  = GD%GDrv1%I(1,L)
      iP  = GD%GDrv1%I(2,L)
      iQ  = GD%GDrv1%I(3,L)
      iU  = GD%GDrv1%I(4,L)
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
    END DO ! L
  ELSE
!----------------------------------------!
! Total contraction length > 1           !
!----------------------------------------!
    C=0.0D0
!----------------------------------------!
! Contract integrals and scale by 2 * Za |
!----------------------------------------!
    DO L=1,GD%LenDa
      iC  = GD%GDrv1%I(1,L)
      iP  = GD%GDrv1%I(2,L)
      iQ  = GD%GDrv1%I(3,L)
      iU  = GD%GDrv1%I(4,L)

      IF (iP==0.AND.iQ==0) THEN
        DO K=1,KKet
          DO J=1,KBra
            Za2=PrmBufB(9,J)
            DO I=1,N
              C(I,iC) = C(I,iC)+U(I,J,K,iU)*Za2
            END DO
          END DO
        END DO
      ELSEIF(iP.GT.0.AND.iQ.EQ.0) THEN
        DO K=1,KKet
          DO J=1,KBra
            Za2=PrmBufB(9,J)
            DO I=1,N
              C(I,iC) = C(I,iC)+U(I,J,K,iU)*CB(J,iP)*Za2
            END DO
          END DO
        END DO
      ELSEIF(iP.EQ.0.AND.iQ.GT.0) THEN
        DO K=1,KKet
          DO J=1,KBra
            Za2=PrmBufB(9,J)
            DO I=1,N
              C(I,iC) = C(I,iC)+U(I,J,K,iU)*CK(I,K,iQ)*Za2
            END DO
          END DO
        END DO
      ELSE
        DO K=1,KKet
          DO J=1,KBra
            Za2=PrmBufB(9,J)
            DO I=1,N
              C(I,iC) = C(I,iC)+U(I,J,K,iU)*CB(J,iP)*CK(I,K,iQ)*Za2
            END DO
          END DO
        END DO
      END IF
    END DO ! L
!----------------------------------------!
! Contract integrals and scale by 2 * Zc |
!----------------------------------------!
    DO L=GD%LenDa+1,GD%LenDc
      iC  = GD%GDrv1%I(1,L)
      iP  = GD%GDrv1%I(2,L)
      iQ  = GD%GDrv1%I(3,L)
      iU  = GD%GDrv1%I(4,L)

      write(*,*) "C:iU=",iU

      IF (iP==0.AND.iQ==0) THEN
        DO K=1,KKet
          DO J=1,KBra
            Zc2=PrmBufB(10,J)
            DO I=1,N
              C(I,iC) = C(I,iC)+U(I,J,K,iU)*Zc2
            END DO
          END DO
        END DO
      ELSEIF(iP.GT.0.AND.iQ.EQ.0) THEN
        DO K=1,KKet
          DO J=1,KBra
            Zc2=PrmBufB(10,J)
            DO I=1,N
              C(I,iC) = C(I,iC)+U(I,J,K,iU)*CB(J,iP)*Zc2
            END DO
          END DO
        END DO
      ELSEIF(iP.EQ.0.AND.iQ.GT.0) THEN
        DO K=1,KKet
          DO J=1,KBra
            Zc2=PrmBufB(10,J)
            DO I=1,N
              C(I,iC) = C(I,iC)+U(I,J,K,iU)*CK(I,K,iQ)*Zc2
            END DO
          END DO
        END DO
      ELSE
        DO K=1,KKet
          DO J=1,KBra
            Zc2=PrmBufB(10,J)
            DO I=1,N
              C(I,iC) = C(I,iC)+U(I,J,K,iU)*CB(J,iP)*CK(I,K,iQ)*Zc2
            END DO
          END DO
        END DO
      END IF
    END DO ! L
!----------------------------------------!
! Contract integrals and scale by 2 * Zb |
!----------------------------------------!
    DO L=GD%LenDc+1,GD%LenDb
      iC  = GD%GDrv1%I(1,L)
      iP  = GD%GDrv1%I(2,L)
      iQ  = GD%GDrv1%I(3,L)
      iU  = GD%GDrv1%I(4,L)

      write(*,*) "B:iU=",iU

      IF (iP==0.AND.iQ==0) THEN
        DO K=1,KKet
          Ioff=DB%MAXP*(K-1)+8
          DO J=1,KBra
            DO I=1,N
              Ind=SB%SLPrm%I(I)+Ioff
              Zb2=DB%PrmBuf%D(Ind)
              C(I,iC) = C(I,iC)+U(I,J,K,iU)*Zb2
            END DO
          END DO
        END DO
      ELSEIF(iP.GT.0.AND.iQ.EQ.0) THEN
        DO K=1,KKet
          Ioff=DB%MAXP*(K-1)+8
          DO J=1,KBra
            DO I=1,N
              Ind=SB%SLPrm%I(I)+Ioff
              Zb2=DB%PrmBuf%D(Ind)
              C(I,iC) = C(I,iC)+U(I,J,K,iU)*CB(J,iP)*Zb2
            END DO
          END DO
        END DO
      ELSEIF(iP.EQ.0.AND.iQ.GT.0) THEN
        DO K=1,KKet
          Ioff=DB%MAXP*(K-1)+8
          DO J=1,KBra
            DO I=1,N
              Ind=SB%SLPrm%I(I)+Ioff
              Zb2=DB%PrmBuf%D(Ind)
              C(I,iC) = C(I,iC)+U(I,J,K,iU)*CK(I,K,iQ)*Zb2
            END DO
          END DO
        END DO
      ELSE
        DO K=1,KKet
          Ioff=DB%MAXP*(K-1)+8
          DO J=1,KBra
            DO I=1,N
              Ind=SB%SLPrm%I(I)+Ioff
              Zb2=DB%PrmBuf%D(Ind)
              C(I,iC) = C(I,iC)+U(I,J,K,iU)*CB(J,iP)*CK(I,K,iQ)*Zb2
            END DO
          END DO
        END DO
      END IF
    END DO ! L
!----------------------------------------!
! Contract integrals (no scale factors)  |
!----------------------------------------!
    DO L=GD%LenDb+1,GD%LenDn
      iC  = GD%GDrv1%I(1,L)
      iP  = GD%GDrv1%I(2,L)
      iQ  = GD%GDrv1%I(3,L)
      iU  = GD%GDrv1%I(4,L)
      IF (iP==0.AND.iQ==0) THEN
        DO K=1,KKet
          DO J=1,KBra
            DO I=1,N
              C(I,iC) = C(I,iC)+U(I,J,K,iU)
            END DO
          END DO
        END DO
      ELSEIF(iP.GT.0.AND.iQ.EQ.0) THEN
        DO K=1,KKet
          DO J=1,KBra
            DO I=1,N
              C(I,iC) = C(I,iC)+U(I,J,K,iU)*CB(J,iP)
            END DO
          END DO
        END DO
      ELSEIF(iP.EQ.0.AND.iQ.GT.0) THEN
        DO K=1,KKet
          DO J=1,KBra
            DO I=1,N
              C(I,iC) = C(I,iC)+U(I,J,K,iU)*CK(I,K,iQ)
            END DO
          END DO
        END DO
      ELSE
        DO K=1,KKet
          DO J=1,KBra
            DO I=1,N
              C(I,iC) = C(I,iC)+U(I,J,K,iU)*CB(J,iP)*CK(I,K,iQ)
            END DO
          END DO
        END DO
      END IF
    END DO ! L
  END IF
END SUBROUTINE ContractG
