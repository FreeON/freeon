SUBROUTINE Int2121(N,IntCode,CBra,CKet,DisBufB,PrmBufB,DB,IB,SB,C,U)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
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
  REAL(DOUBLE)  :: U(16),C(N,16)
  REAL(DOUBLE)  :: Cx,Cy,Cz
  REAL(DOUBLE)  :: Zeta,Up,Px,Py,Pz,cp1,cp2,cp3
  REAL(DOUBLE)  :: Eta, Uq,Qx,Qy,Qz,cq1,cq2,cq3
  REAL(DOUBLE)  :: EtaQx,EtaQy,EtaQz
  REAL(DOUBLE)  :: r1xZpE,ZpE,Rkk,r1xZ,r1x2Z,rpxZ
  REAL(DOUBLE)  :: PAx,PAy,PAz,PQx,PQy,PQz
  REAL(DOUBLE)  :: QBx,QBy,QBz,WQx,WQy,WQz
  REAL(DOUBLE)  :: Wx,Wy,Wz,WPx,WPy,WPz
  REAL(DOUBLE)  :: T1,T2,T3,TwoT
  REAL(DOUBLE)  :: R1,R2,R3,G0,G1,G2,ET
  REAL(DOUBLE)  :: Ts0,r1xE,r1x2E,rpxE,r1x2ZE
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  INTEGER       :: I,J,K,M,l
  INTEGER       :: I0,I1,I2

  IF(N.EQ.0) RETURN

  C=0.0D0

  Cx=DisBufB( 8)
  Cy=DisBufB( 9)
  Cz=DisBufB(10)

  IF(IntCode.EQ.02010201) THEN

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
        cp1  = PrmBufB(6,K)
        cp2  = PrmBufB(7,K)
        cp3  = PrmBufB(8,K)
        ZpE  = Zeta+Eta
        r1xZpE = 1.0D0/ZpE
        r1x2ZE = 0.5D0*r1xZpE
        r1xZ   = 1.0D0/Zeta
        r1x2Z  = 0.5D0*r1xZ
        rpxZ   = Eta*r1xZpE
        rpxE   = Zeta*r1xZpE
        Rkk = Up*Uq*DSQRT(r1xZpE)
        PAx = Px-Cx
        PAy = Py-Cy
        PAz = Pz-Cz
        Wx  = (Zeta*Px+EtaQx)*r1xZpE
        Wy  = (Zeta*Py+EtaQy)*r1xZpE
        Wz  = (Zeta*Pz+EtaQz)*r1xZpE
        WPx = Wx-Px
        WPy = Wy-Py
        WPz = Wz-Pz
        WQx = Wx-Qx
        WQy = Wy-Qy
        WQz = Wz-Qz
        PQx = Px-Qx
        PQy = Py-Qy
        PQz = Pz-Qz
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE

        IF (T1<IB%Switch) THEN
          M=INT(T1*IB%Grid)
          T2=T1*T1
          T3=T2*T1
          TwoT=2.0D0*T1
          ET=IB%ET%D(0,M)+T1*IB%ET%D(1,M)+T2*IB%ET%D(2,M)+T3*IB%ET%D(3,M)
          G2=IB%GT%D(0,M)+T1*IB%GT%D(1,M)+T2*IB%GT%D(2,M)+T3*IB%GT%D(3,M)
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
        U(4)=R2
        U(2)=R3
        U(5)=QBx*U(1)+WQx*U(4)
        U(9)=QBy*U(1)+WQy*U(4)
        U(13)=QBz*U(1)+WQz*U(4)
        U(8)=QBx*U(4)+WQx*U(2)
        U(12)=QBy*U(4)+WQy*U(2)
        U(15)=QBz*U(4)+WQz*U(2)
        Ts0 = r1x2ZE*U(4)
        U(6)=PAx*U(5)+WPx*U(8)+ Ts0
        U(11)=PAy*U(9)+WPy*U(12)+ Ts0
        U(16)=PAz*U(13)+WPz*U(15)+ Ts0
        U(2)=PAx*U(1)+WPx*U(4)
        U(3)=PAy*U(1)+WPy*U(4)
        U(4)=PAz*U(1)+WPz*U(4)
        U(7)=PAy*U(5)+WPy*U(8)
        U(8)=PAz*U(5)+WPz*U(8)
        U(10)=PAx*U(9)+WPx*U(12)
        U(12)=PAz*U(9)+WPz*U(12)
        U(14)=PAx*U(13)+WPx*U(15)
        U(15)=PAy*U(13)+WPy*U(15)

        C(I,1)=C(I,1)+U(1)*cp1*cq1
        C(I,2)=C(I,2)+U(2)*cq1
        C(I,3)=C(I,3)+U(3)*cq1
        C(I,4)=C(I,4)+U(4)*cq1
        C(I,5)=C(I,5)+U(5)*cp1
        C(I,6)=C(I,6)+U(6)
        C(I,7)=C(I,7)+U(7)
        C(I,8)=C(I,8)+U(8)
        C(I,9)=C(I,9)+U(9)*cp1
        C(I,10)=C(I,10)+U(10)
        C(I,11)=C(I,11)+U(11)
        C(I,12)=C(I,12)+U(12)
        C(I,13)=C(I,13)+U(13)*cp1
        C(I,14)=C(I,14)+U(14)
        C(I,15)=C(I,15)+U(15)
        C(I,16)=C(I,16)+U(16)

      END DO ! K
    END DO ! J
  END DO ! I

  ELSEIF(IntCode.EQ.03010201) THEN

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
        cp1  = PrmBufB(6,K)
        cp2  = PrmBufB(7,K)
        cp3  = PrmBufB(8,K)
        ZpE  = Zeta+Eta
        r1xZpE = 1.0D0/ZpE
        r1x2ZE = 0.5D0*r1xZpE
        r1xZ   = 1.0D0/Zeta
        r1x2Z  = 0.5D0*r1xZ
        rpxZ   = Eta*r1xZpE
        rpxE   = Zeta*r1xZpE
        Rkk = Up*Uq*DSQRT(r1xZpE)
        PAx = Px-Cx
        PAy = Py-Cy
        PAz = Pz-Cz
        Wx  = (Zeta*Px+EtaQx)*r1xZpE
        Wy  = (Zeta*Py+EtaQy)*r1xZpE
        Wz  = (Zeta*Pz+EtaQz)*r1xZpE
        WPx = Wx-Px
        WPy = Wy-Py
        WPz = Wz-Pz
        WQx = Wx-Qx
        WQy = Wy-Qy
        WQz = Wz-Qz
        PQx = Px-Qx
        PQy = Py-Qy
        PQz = Pz-Qz
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE

        IF (T1<IB%Switch) THEN
          M=INT(T1*IB%Grid)
          T2=T1*T1
          T3=T2*T1
          TwoT=2.0D0*T1
          ET=IB%ET%D(0,M)+T1*IB%ET%D(1,M)+T2*IB%ET%D(2,M)+T3*IB%ET%D(3,M)
          G2=IB%GT%D(0,M)+T1*IB%GT%D(1,M)+T2*IB%GT%D(2,M)+T3*IB%GT%D(3,M)
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
        U(4)=R2
        U(2)=R3
        U(5)=QBx*U(1)+WQx*U(4)
        U(9)=QBy*U(1)+WQy*U(4)
        U(13)=QBz*U(1)+WQz*U(4)
        U(8)=QBx*U(4)+WQx*U(2)
        U(12)=QBy*U(4)+WQy*U(2)
        U(15)=QBz*U(4)+WQz*U(2)
        Ts0 = r1x2ZE*U(4)
        U(6)=PAx*U(5)+WPx*U(8)+Ts0
        U(11)=PAy*U(9)+WPy*U(12)+Ts0
        U(16)=PAz*U(13)+WPz*U(15)+Ts0
        U(2)=PAx*U(1)+WPx*U(4)
        U(3)=PAy*U(1)+WPy*U(4)
        U(4)=PAz*U(1)+WPz*U(4)
        U(7)=PAy*U(5)+WPy*U(8)
        U(8)=PAz*U(5)+WPz*U(8)
        U(10)=PAx*U(9)+WPx*U(12)
        U(12)=PAz*U(9)+WPz*U(12)
        U(14)=PAx*U(13)+WPx*U(15)
        U(15)=PAy*U(13)+WPy*U(15)

        C(I,1)=C(I,1)+U(2)*cq1
        C(I,2)=C(I,2)+U(3)*cq1
        C(I,3)=C(I,3)+U(4)*cq1
        C(I,5)=C(I,5)+U(6)
        C(I,6)=C(I,6)+U(7)
        C(I,7)=C(I,7)+U(8)
        C(I,9)=C(I,9)+U(10)
        C(I,10)=C(I,10)+U(11)
        C(I,11)=C(I,11)+U(12)
        C(I,13)=C(I,13)+U(14)
        C(I,14)=C(I,14)+U(15)
        C(I,15)=C(I,15)+U(16)

      END DO ! K
    END DO ! J
  END DO ! I

  ELSEIF(IntCode.EQ.02010301) THEN

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
        cp1  = PrmBufB(6,K)
        cp2  = PrmBufB(7,K)
        cp3  = PrmBufB(8,K)
        ZpE  = Zeta+Eta
        r1xZpE = 1.0D0/ZpE
        r1x2ZE = 0.5D0*r1xZpE
        r1xZ   = 1.0D0/Zeta
        r1x2Z  = 0.5D0*r1xZ
        rpxZ   = Eta*r1xZpE
        rpxE   = Zeta*r1xZpE
        Rkk = Up*Uq*DSQRT(r1xZpE)
        PAx = Px-Cx
        PAy = Py-Cy
        PAz = Pz-Cz
        Wx  = (Zeta*Px+EtaQx)*r1xZpE
        Wy  = (Zeta*Py+EtaQy)*r1xZpE
        Wz  = (Zeta*Pz+EtaQz)*r1xZpE
        WPx = Wx-Px
        WPy = Wy-Py
        WPz = Wz-Pz
        WQx = Wx-Qx
        WQy = Wy-Qy
        WQz = Wz-Qz
        PQx = Px-Qx
        PQy = Py-Qy
        PQz = Pz-Qz
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE

        IF (T1<IB%Switch) THEN
          M=INT(T1*IB%Grid)
          T2=T1*T1
          T3=T2*T1
          TwoT=2.0D0*T1
          ET=IB%ET%D(0,M)+T1*IB%ET%D(1,M)+T2*IB%ET%D(2,M)+T3*IB%ET%D(3,M)
          G2=IB%GT%D(0,M)+T1*IB%GT%D(1,M)+T2*IB%GT%D(2,M)+T3*IB%GT%D(3,M)
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
        U(4)=R2
        U(2)=R3
        U(5)=QBx*U(1)+WQx*U(4)
        U(9)=QBy*U(1)+WQy*U(4)
        U(13)=QBz*U(1)+WQz*U(4)
        U(8)=QBx*U(4)+WQx*U(2)
        U(12)=QBy*U(4)+WQy*U(2)
        U(15)=QBz*U(4)+WQz*U(2)
        Ts0 = r1x2ZE*U(4)
        U(6)=PAx*U(5)+WPx*U(8)+ Ts0
        U(11)=PAy*U(9)+WPy*U(12)+ Ts0
        U(16)=PAz*U(13)+WPz*U(15)+ Ts0
        U(2)=PAx*U(1)+WPx*U(4)
        U(3)=PAy*U(1)+WPy*U(4)
        U(4)=PAz*U(1)+WPz*U(4)
        U(7)=PAy*U(5)+WPy*U(8)
        U(8)=PAz*U(5)+WPz*U(8)
        U(10)=PAx*U(9)+WPx*U(12)
        U(12)=PAz*U(9)+WPz*U(12)
        U(14)=PAx*U(13)+WPx*U(15)
        U(15)=PAy*U(13)+WPy*U(15)

        C(I,1)=C(I,1)+U(5)*cp1
        C(I,2)=C(I,2)+U(6)
        C(I,3)=C(I,3)+U(7)
        C(I,4)=C(I,4)+U(8)
        C(I,5)=C(I,5)+U(9)*cp1
        C(I,6)=C(I,6)+U(10)
        C(I,7)=C(I,7)+U(11)
        C(I,8)=C(I,8)+U(12)
        C(I,9)=C(I,9)+U(13)*cp1
        C(I,10)=C(I,10)+U(14)
        C(I,11)=C(I,11)+U(15)
        C(I,12)=C(I,12)+U(16)

      END DO ! K
    END DO ! J
  END DO ! I

  ELSEIF(IntCode.EQ.03010301) THEN

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
        cp1  = PrmBufB(6,K)
        cp2  = PrmBufB(7,K)
        cp3  = PrmBufB(8,K)
        ZpE  = Zeta+Eta
        r1xZpE = 1.0D0/ZpE
        r1x2ZE = 0.5D0*r1xZpE
        r1xZ   = 1.0D0/Zeta
        r1x2Z  = 0.5D0*r1xZ
        rpxZ   = Eta*r1xZpE
        rpxE   = Zeta*r1xZpE
        Rkk = Up*Uq*DSQRT(r1xZpE)
        PAx = Px-Cx
        PAy = Py-Cy
        PAz = Pz-Cz
        Wx  = (Zeta*Px+EtaQx)*r1xZpE
        Wy  = (Zeta*Py+EtaQy)*r1xZpE
        Wz  = (Zeta*Pz+EtaQz)*r1xZpE
        WPx = Wx-Px
        WPy = Wy-Py
        WPz = Wz-Pz
        WQx = Wx-Qx
        WQy = Wy-Qy
        WQz = Wz-Qz
        PQx = Px-Qx
        PQy = Py-Qy
        PQz = Pz-Qz
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE

        IF (T1<IB%Switch) THEN
          M=INT(T1*IB%Grid)
          T2=T1*T1
          T3=T2*T1
          TwoT=2.0D0*T1
          ET=IB%ET%D(0,M)+T1*IB%ET%D(1,M)+T2*IB%ET%D(2,M)+T3*IB%ET%D(3,M)
          G2=IB%GT%D(0,M)+T1*IB%GT%D(1,M)+T2*IB%GT%D(2,M)+T3*IB%GT%D(3,M)
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
        U(4)=R2
        U(2)=R3
        U(5)=QBx*U(1)+WQx*U(4)
        U(9)=QBy*U(1)+WQy*U(4)
        U(13)=QBz*U(1)+WQz*U(4)
        U(8)=QBx*U(4)+WQx*U(2)
        U(12)=QBy*U(4)+WQy*U(2)
        U(15)=QBz*U(4)+WQz*U(2)
        Ts0 = r1x2ZE*U(4)
        U(6)=PAx*U(5)+WPx*U(8)+ Ts0
        U(11)=PAy*U(9)+WPy*U(12)+ Ts0
        U(16)=PAz*U(13)+WPz*U(15)+ Ts0
        U(2)=PAx*U(1)+WPx*U(4)
        U(3)=PAy*U(1)+WPy*U(4)
        U(4)=PAz*U(1)+WPz*U(4)
        U(7)=PAy*U(5)+WPy*U(8)
        U(8)=PAz*U(5)+WPz*U(8)
        U(10)=PAx*U(9)+WPx*U(12)
        U(12)=PAz*U(9)+WPz*U(12)
        U(14)=PAx*U(13)+WPx*U(15)
        U(15)=PAy*U(13)+WPy*U(15)

        C(I,1)=C(I,1)+U(6)
        C(I,2)=C(I,2)+U(7)
        C(I,3)=C(I,3)+U(8)
        C(I,5)=C(I,5)+U(10)
        C(I,6)=C(I,6)+U(11)
        C(I,7)=C(I,7)+U(12)
        C(I,9)=C(I,9)+U(14)
        C(I,10)=C(I,10)+U(15)
        C(I,11)=C(I,11)+U(16)

      END DO ! K
    END DO ! J
  END DO ! I

  ELSE
    CALL Halt('Illegal IntCode in Int2121')
  ENDIF

END SUBROUTINE Int2121
