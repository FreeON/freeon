SUBROUTINE Int1122(N,IntCode,CBra,CKet,DisBufB,PrmBufB,DB,IB,SB,C,U)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE GammaF2
  USE InvExp
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Distribution buffer stuff
!--------------------------------------------------------------------------------
  TYPE(DBuf),INTENT(IN)    :: DB            ! ONX distribution buffers
  TYPE(IBuf),INTENT(INOUT) :: IB            ! ONX 2-e eval buffers
  TYPE(DSL),INTENT(IN)     :: SB            ! ONX distribution pointers
  INTEGER       :: N,IntCode,CBra,CKet
  REAL(DOUBLE)  :: DisBufB(DB%MAXC)
  REAL(DOUBLE)  :: PrmBufB(DB%MAXP,CBra+DB%MInfo)
!--------------------------------------------------------------------------------
! Temporary space for computing 2-e integrals
!--------------------------------------------------------------------------------
  REAL(DOUBLE)  :: U(10),C(N,17)
  REAL(DOUBLE)  :: Cx,Cy,Cz
  REAL(DOUBLE)  :: Zeta,Up,Px,Py,Pz,cp1,cp2,cp3
  REAL(DOUBLE)  :: Eta, Uq,Qx,Qy,Qz,cq1,cq2,cq3
  REAL(DOUBLE)  :: EtaQx,EtaQy,EtaQz
  REAL(DOUBLE)  :: r1xZpE,ZpE,Rkk
  REAL(DOUBLE)  :: PAx,PAy,PAz,PQx,PQy,PQz
  REAL(DOUBLE)  :: QBx,QBy,QBz,WQx,WQy,WQz
  REAL(DOUBLE)  :: Wx,Wy,Wz,WPx,WPy,WPz
  REAL(DOUBLE)  :: T1,T2,T3,T4,TwoT
  REAL(DOUBLE)  :: R1,R2,R3,G0,G1,G2,ET
  REAL(DOUBLE)  :: Ts0,r1xE,r1x2E,rpxE,r1x2ZE
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  INTEGER       :: I,J,K,MG,ME
  INTEGER       :: I0,I1,I2

  IF(N.EQ.0) RETURN

  C=0.0D0

  IF(IntCode.EQ.1010202) THEN

  DO I=1,N
    I0 = SB%SLPrm%I(I)
    I2 = SB%SLDis%I(I)-4
    DO J=1,CKet
      I1=I0+(J-1)*DB%MAXP
      Eta = DB%PrmBuf%D(I1)
      Qx  = DB%PrmBuf%D(I1+1)
      Qy  = DB%PrmBuf%D(I1+2)
      Qz  = DB%PrmBuf%D(I1+3)
      Uq  = DB%PrmBuf%D(I1+4)
      cq1 = DB%PrmBuf%D(I1+5)
      cq2 = DB%PrmBuf%D(I1+6)
      cq3 = DB%PrmBuf%D(I1+7)
      r1xE  = 1.0D0/Eta
      r1x2E = 0.5D0*r1xE
      EtaQx = Eta*Qx
      EtaQy = Eta*Qy
      EtaQz = Eta*Qz
      QBx = Qx-DB%DisBuf%D(I2+7)
      QBy = Qy-DB%DisBuf%D(I2+8)
      QBz = Qz-DB%DisBuf%D(I2+9)
      DO K=1,CBra
        Zeta = PrmBufB(1,K)
        Px   = PrmBufB(2,K)
        Py   = PrmBufB(3,K)
        Pz   = PrmBufB(4,K)
        Up   = PrmBufB(5,K)
        ZpE  = Zeta+Eta
        r1xZpE = 1.0D0/ZpE
        r1x2ZE = 0.5D0*r1xZpE
        rpxE   = Zeta*r1xZpE
        Rkk  = Up*Uq*DSQRT(r1xZpE)
        Wx  = (Zeta*Px+EtaQx)*r1xZpE
        Wy  = (Zeta*Py+EtaQy)*r1xZpE
        Wz  = (Zeta*Pz+EtaQz)*r1xZpE
        WQx = Wx - Qx
        WQy = Wy - Qy
        WQz = Wz - Qz
        PQx = Px - Qx
        PQy = Py - Qy
        PQz = Pz - Qz
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE

        IF (T1<Gamma_Switch)THEN
          MG=AINT(T1*Gamma_Grid)
          ME=AINT(T1*Exp_Grid)
          T2=T1*T1
          T3=T2*T1
          T4=T2*T2
          TwoT=2.0D0*T1
          ET=Exp_0(ME)+T1*Exp_1(ME)+T2*Exp_2(ME)+T3*Exp_3(ME)+T4*Exp_4(ME)
          G2=F2_0(MG) +T1*F2_1(MG) +T2*F2_2(MG) +T3*F2_3(MG) +T4*F2_4(MG)
          G1=.333333333333333D+00*(TwoT*G2+ET)
          G0=TwoT*G1+ET
          R1=Rkk*G0
          R2=Rkk*G1
          R3=Rkk*G2
        ELSE
          T2=1.0D0/T1
          T3=DSQRT(T2)
          R1=Rkk*IB%GammaA%D(1)*T3
          T3=T3*T2
          R2=Rkk*IB%GammaA%D(2)*T3
          T3=T3*T2
          R3=Rkk*IB%GammaA%D(3)*T3
        ENDIF
        U(1)=R1
        U(6)=R2
        U(5)=R3
        U(2)=QBx*U(1)+WQx*U(6)
        U(3)=QBy*U(1)+WQy*U(6)
        U(4)=QBz*U(1)+WQz*U(6)
        U(7)=QBx*U(6)+WQx*U(5)
        U(8)=QBy*U(6)+WQy*U(5)
        U(9)=QBz*U(6)+WQz*U(5)
        Ts0=r1x2E*(U(1)-rpxE*U(6))
        U(5)=QBx*U(2)+WQx*U(7)+Ts0
        U(7)=QBy*U(3)+WQy*U(8)+Ts0
        U(10)=QBz*U(4)+WQz*U(9)+Ts0
        U(6)=QBx*U(3)+WQx*U(8)
        U(8)=QBx*U(4)+WQx*U(9)
        U(9)=QBy*U(4)+WQy*U(9)
        C(I,1)=C(I,1)+U(1)*cq1
        C(I,2)=C(I,2)+U(2)*cq3
        C(I,3)=C(I,3)+U(3)*cq3
        C(I,4)=C(I,4)+U(4)*cq3
        C(I,5)=C(I,5)+U(2)*cq2
        C(I,6)=C(I,6)+U(5)
        C(I,7)=C(I,7)+U(1)*cq3
        C(I,8)=C(I,8)+U(2)
        C(I,9)=C(I,9)+U(3)*cq2
        C(I,10)=C(I,10)+U(6)
        C(I,11)=C(I,11)+U(7)
        C(I,12)=C(I,12)+U(3)
        C(I,13)=C(I,13)+U(4)*cq2
        C(I,14)=C(I,14)+U(8)
        C(I,15)=C(I,15)+U(9)
        C(I,16)=C(I,16)+U(10)
        C(I,17)=C(I,17)+U(4)
      END DO ! K
    END DO ! J
  END DO ! I

  ELSEIF(IntCode.EQ.1010302) THEN

  DO I=1,N
    I0 = SB%SLPrm%I(I)
    I2 = SB%SLDis%I(I)-4
    DO J=1,CKet
      I1=I0+(J-1)*DB%MAXP
      Eta = DB%PrmBuf%D(I1)
      Qx  = DB%PrmBuf%D(I1+1)
      Qy  = DB%PrmBuf%D(I1+2)
      Qz  = DB%PrmBuf%D(I1+3)
      Uq  = DB%PrmBuf%D(I1+4)
      cq1 = DB%PrmBuf%D(I1+5)
      cq2 = DB%PrmBuf%D(I1+6)
      cq3 = DB%PrmBuf%D(I1+7)
      r1xE  = 1.0D0/Eta
      r1x2E = 0.5D0*r1xE
      EtaQx = Eta*Qx
      EtaQy = Eta*Qy
      EtaQz = Eta*Qz
      QBx = Qx-DB%DisBuf%D(I2+7)
      QBy = Qy-DB%DisBuf%D(I2+8)
      QBz = Qz-DB%DisBuf%D(I2+9)
      DO K=1,CBra
        Zeta = PrmBufB(1,K)
        Px   = PrmBufB(2,K)
        Py   = PrmBufB(3,K)
        Pz   = PrmBufB(4,K)
        Up   = PrmBufB(5,K)
        ZpE  = Zeta+Eta
        r1xZpE = 1.0D0/ZpE
        r1x2ZE = 0.5D0*r1xZpE
        rpxE   = Zeta*r1xZpE
        Rkk  = Up*Uq*DSQRT(r1xZpE)
        Wx  = (Zeta*Px+EtaQx)*r1xZpE
        Wy  = (Zeta*Py+EtaQy)*r1xZpE
        Wz  = (Zeta*Pz+EtaQz)*r1xZpE
        WQx = Wx - Qx
        WQy = Wy - Qy
        WQz = Wz - Qz
        PQx = Px - Qx
        PQy = Py - Qy
        PQz = Pz - Qz
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
        IF (T1<Gamma_Switch)THEN
          MG=AINT(T1*Gamma_Grid)
          ME=AINT(T1*Exp_Grid)
          T2=T1*T1
          T3=T2*T1
          T4=T2*T2
          TwoT=2.0D0*T1
          ET=Exp_0(ME)+T1*Exp_1(ME)+T2*Exp_2(ME)+T3*Exp_3(ME)+T4*Exp_4(ME)
          G2=F2_0(MG) +T1*F2_1(MG) +T2*F2_2(MG) +T3*F2_3(MG) +T4*F2_4(MG)
          G1=.333333333333333D+00*(TwoT*G2+ET)
          G0=TwoT*G1+ET
          R1=Rkk*G0
          R2=Rkk*G1
          R3=Rkk*G2
        ELSE
          T2=1.0D0/T1
          T3=DSQRT(T2)
          R1=Rkk*IB%GammaA%D(1)*T3
          T3=T3*T2
          R2=Rkk*IB%GammaA%D(2)*T3
          T3=T3*T2
          R3=Rkk*IB%GammaA%D(3)*T3
        ENDIF
        U(1)=R1
        U(6)=R2
        U(5)=R3
        U(2)=QBx*U(1)+WQx*U(6)
        U(3)=QBy*U(1)+WQy*U(6)
        U(4)=QBz*U(1)+WQz*U(6)
        U(7)=QBx*U(6)+WQx*U(5)
        U(8)=QBy*U(6)+WQy*U(5)
        U(9)=QBz*U(6)+WQz*U(5)
        Ts0=r1x2E*(U(1)-rpxE*U(6))
        U(5)=QBx*U(2)+WQx*U(7)+Ts0
        U(7)=QBy*U(3)+WQy*U(8)+Ts0
        U(10)=QBz*U(4)+WQz*U(9)+Ts0
        U(6)=QBx*U(3)+WQx*U(8)
        U(8)=QBx*U(4)+WQx*U(9)
        U(9)=QBy*U(4)+WQy*U(9)
        C(I,1)=C(I,1)+U(2)*cq2
        C(I,2)=C(I,2)+U(2)
        C(I,3)=C(I,3)+U(3)
        C(I,4)=C(I,4)+U(4)
        C(I,5)=C(I,5)+U(3)*cq2
        C(I,6)=C(I,6)+U(5)
        C(I,7)=C(I,7)+U(6)
        C(I,8)=C(I,8)+U(7)
        C(I,9)=C(I,9)+U(4)*cq2
        C(I,10)=C(I,10)+U(8)
        C(I,11)=C(I,11)+U(9)
        C(I,12)=C(I,12)+U(10)
      END DO ! K
    END DO ! J
  END DO ! I

  ELSEIF(IntCode.EQ.1010303) THEN

  DO I=1,N
    I0 = SB%SLPrm%I(I)
    I2 = SB%SLDis%I(I)-4
    DO J=1,CKet
      I1=I0+(J-1)*DB%MAXP
      Eta = DB%PrmBuf%D(I1)
      Qx  = DB%PrmBuf%D(I1+1)
      Qy  = DB%PrmBuf%D(I1+2)
      Qz  = DB%PrmBuf%D(I1+3)
      Uq  = DB%PrmBuf%D(I1+4)
      cq1 = DB%PrmBuf%D(I1+5)
      cq2 = DB%PrmBuf%D(I1+6)
      cq3 = DB%PrmBuf%D(I1+7)
      r1xE  = 1.0D0/Eta
      r1x2E = 0.5D0*r1xE
      EtaQx = Eta*Qx
      EtaQy = Eta*Qy
      EtaQz = Eta*Qz
      QBx = Qx-DB%DisBuf%D(I2+7)
      QBy = Qy-DB%DisBuf%D(I2+8)
      QBz = Qz-DB%DisBuf%D(I2+9)
      DO K=1,CBra
        Zeta = PrmBufB(1,K)
        Px   = PrmBufB(2,K)
        Py   = PrmBufB(3,K)
        Pz   = PrmBufB(4,K)
        Up   = PrmBufB(5,K)
        ZpE  = Zeta+Eta
        r1xZpE = 1.0D0/ZpE
        r1x2ZE = 0.5D0*r1xZpE
        rpxE   = Zeta*r1xZpE
        Rkk  = Up*Uq*DSQRT(r1xZpE)
        Wx  = (Zeta*Px+EtaQx)*r1xZpE
        Wy  = (Zeta*Py+EtaQy)*r1xZpE
        Wz  = (Zeta*Pz+EtaQz)*r1xZpE
        WQx = Wx - Qx
        WQy = Wy - Qy
        WQz = Wz - Qz
        PQx = Px - Qx
        PQy = Py - Qy
        PQz = Pz - Qz
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
        IF (T1<Gamma_Switch)THEN
          MG=AINT(T1*Gamma_Grid)
          ME=AINT(T1*Exp_Grid)
          T2=T1*T1
          T3=T2*T1
          T4=T2*T2
          TwoT=2.0D0*T1
          ET=Exp_0(ME)+T1*Exp_1(ME)+T2*Exp_2(ME)+T3*Exp_3(ME)+T4*Exp_4(ME)
          G2=F2_0(MG) +T1*F2_1(MG) +T2*F2_2(MG) +T3*F2_3(MG) +T4*F2_4(MG)
          G1=.333333333333333D+00*(TwoT*G2+ET)
          G0=TwoT*G1+ET
          R1=Rkk*G0
          R2=Rkk*G1
          R3=Rkk*G2
        ELSE
          T2=1.0D0/T1
          T3=DSQRT(T2)
          R1=Rkk*IB%GammaA%D(1)*T3
          T3=T3*T2
          R2=Rkk*IB%GammaA%D(2)*T3
          T3=T3*T2
          R3=Rkk*IB%GammaA%D(3)*T3
        ENDIF
        U(1)=R1
        U(6)=R2
        U(5)=R3
        U(2)=QBx*U(1)+WQx*U(6)
        U(3)=QBy*U(1)+WQy*U(6)
        U(4)=QBz*U(1)+WQz*U(6)
        U(7)=QBx*U(6)+WQx*U(5)
        U(8)=QBy*U(6)+WQy*U(5)
        U(9)=QBz*U(6)+WQz*U(5)
        Ts0=r1x2E*(U(1)-rpxE*U(6))
        U(5)=QBx*U(2)+WQx*U(7)+Ts0
        U(7)=QBy*U(3)+WQy*U(8)+Ts0
        U(10)=QBz*U(4)+WQz*U(9)+Ts0
        U(6)=QBx*U(3)+WQx*U(8)
        U(8)=QBx*U(4)+WQx*U(9)
        U(9)=QBy*U(4)+WQy*U(9)
        C(I,1)=C(I,1)+U(2)
        C(I,2)=C(I,2)+U(3)
        C(I,3)=C(I,3)+U(4)
        C(I,4)=C(I,4)+U(5)
        C(I,5)=C(I,5)+U(6)
        C(I,6)=C(I,6)+U(7)
        C(I,7)=C(I,7)+U(8)
        C(I,8)=C(I,8)+U(9)
        C(I,9)=C(I,9)+U(10)
      END DO ! K
    END DO ! J
  END DO ! I

  ELSE
    CALL Halt('Illegal IntCode in Int1122')
  ENDIF

END SUBROUTINE Int1122
