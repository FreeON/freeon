SUBROUTINE Int2122(N,IntCode,CBra,CKet,DisBufB,PrmBufB,DB,IB,SB,C,U)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE GammaF3
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
  REAL(DOUBLE)  :: U(40),C(N,68)
  REAL(DOUBLE)  :: Cx,Cy,Cz
  REAL(DOUBLE)  :: Zeta,Up,Px,Py,Pz,cp1,cp2,cp3
  REAL(DOUBLE)  :: Eta, Uq,Qx,Qy,Qz,cq1,cq2,cq3
  REAL(DOUBLE)  :: EtaQx,EtaQy,EtaQz
  REAL(DOUBLE)  :: r1xZpE,ZpE,Rkk,r1xZ,r1x2Z,rpxZ
  REAL(DOUBLE)  :: PAx,PAy,PAz,PQx,PQy,PQz
  REAL(DOUBLE)  :: QBx,QBy,QBz,WQx,WQy,WQz
  REAL(DOUBLE)  :: Wx,Wy,Wz,WPx,WPy,WPz
  REAL(DOUBLE)  :: T1,T2,T3,T4,TwoT
  REAL(DOUBLE)  :: R1,R2,R3,R4,G0,G1,G2,G3,ET
  REAL(DOUBLE)  :: Ts0,Ts1,r1xE,r1x2E,rpxE,r1x2ZE
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  INTEGER       :: I,J,K,ME,MG,l
  INTEGER       :: I0,I1,I2

  IF(N.EQ.0) RETURN

  C=0.0D0

  Cx=DisBufB( 8)
  Cy=DisBufB( 9)
  Cz=DisBufB(10)

  IF(IntCode.EQ.02010202) THEN

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

        IF (T1<Gamma_Switch)THEN
          MG=AINT(T1*Gamma_Grid)
          ME=AINT(T1*Exp_Grid)
          T2=T1*T1
          T3=T2*T1
          T4=T2*T2
          TwoT=2.0D0*T1
          ET=Exp_0(ME)+T1*Exp_1(ME)+T2*Exp_2(ME)+T3*Exp_3(ME)+T4*Exp_4(ME)
          G3=F3_0(MG) +T1*F3_1(MG) +T2*F3_2(MG) +T3*F3_3(MG) +T4*F3_4(MG)
          G2=.200000000000000D+00*(TwoT*G3+ET)
          G1=.333333333333333D+00*(TwoT*G2+ET)
          G0=TwoT*G1+ET
          R1=Rkk*G0
          R2=Rkk*G1
          R3=Rkk*G2
          R4=Rkk*G3
        ELSE
          T2=1.0D0/T1
          T3=DSQRT(T2)
          R1=Rkk*IB%GammaA%D(1)*T3
          T3=T3*T2
          R2=Rkk*IB%GammaA%D(2)*T3
          T3=T3*T2
          R3=Rkk*IB%GammaA%D(3)*T3
          T3=T3*T2
          R4=Rkk*IB%GammaA%D(4)*T3
        ENDIF

        U(1)=R1
        U(4)=R2     
        U(2)=R3     
        U(3)=R4     
        U(5)=QBx*U(1)+WQx*U(4)
        U(9)=QBy*U(1)+WQy*U(4)
        U(13)=QBz*U(1)+WQz*U(4)
        U(8)=QBx*U(4)+WQx*U(2)
        U(12)=QBy*U(4)+WQy*U(2)
        U(15)=QBz*U(4)+WQz*U(2)
        U(6)=QBx*U(2)+WQx*U(3)
        U(7)=QBy*U(2)+WQy*U(3)
        U(10)=QBz*U(2)+WQz*U(3)
        Ts0=r1x2E*(U(1)-rpxE*U(4))
        U(17)=QBx*U(5)+WQx*U(8)+ Ts0
        U(21)=QBx*U(9)+WQx*U(12)
        U(25)=QBy*U(9)+WQy*U(12)+ Ts0
        U(29)=QBx*U(13)+WQx*U(15)
        U(33)=QBy*U(13)+WQy*U(15)
        U(37)=QBz*U(13)+WQz*U(15)+ Ts0
        Ts1=r1x2E*(U(4)-rpxE*U(2))
        U(20)=QBx*U(8)+WQx*U(6)+ Ts1
        U(24)=QBx*U(12)+WQx*U(7)
        U(28)=QBy*U(12)+WQy*U(7)+ Ts1
        U(31)=QBx*U(15)+WQx*U(10)
        U(34)=QBy*U(15)+WQy*U(10)
        U(39)=QBz*U(15)+WQz*U(10)+ Ts1
        U(6)=PAx*U(5)+WPx*U(8)+r1x2ZE*U(4)
        U(11)=PAy*U(9)+WPy*U(12)+r1x2ZE*U(4)
        U(16)=PAz*U(13)+WPz*U(15)+r1x2ZE*U(4)
        U(18)=PAx*U(17)+WPx*U(20)+2.0D0*r1x2ZE*U(8)
        U(22)=PAx*U(21)+WPx*U(24)+r1x2ZE*U(12)
        U(23)=PAy*U(21)+WPy*U(24)+r1x2ZE*U(8)
        U(27)=PAy*U(25)+WPy*U(28)+2.0D0*r1x2ZE*U(12)
        U(30)=PAx*U(29)+WPx*U(31)+r1x2ZE*U(15)
        U(32)=PAz*U(29)+WPz*U(31)+r1x2ZE*U(8)
        U(35)=PAy*U(33)+WPy*U(34)+r1x2ZE*U(15)
        U(36)=PAz*U(33)+WPz*U(34)+r1x2ZE*U(12)
        U(40)=PAz*U(37)+WPz*U(39)+2.0D0*r1x2ZE*U(15)
        U(2)=PAx*U(1)+WPx*U(4)
        U(3)=PAy*U(1)+WPy*U(4)
        U(4)=PAz*U(1)+WPz*U(4)
        U(7)=PAy*U(5)+WPy*U(8)
        U(8)=PAz*U(5)+WPz*U(8)
        U(10)=PAx*U(9)+WPx*U(12)
        U(12)=PAz*U(9)+WPz*U(12)
        U(14)=PAx*U(13)+WPx*U(15)
        U(15)=PAy*U(13)+WPy*U(15)
        U(19)=PAy*U(17)+WPy*U(20)
        U(20)=PAz*U(17)+WPz*U(20)
        U(24)=PAz*U(21)+WPz*U(24)
        U(26)=PAx*U(25)+WPx*U(28)
        U(28)=PAz*U(25)+WPz*U(28)
        U(31)=PAy*U(29)+WPy*U(31)
        U(34)=PAx*U(33)+WPx*U(34)
        U(38)=PAx*U(37)+WPx*U(39)
        U(39)=PAy*U(37)+WPy*U(39)

        C(I,1)=C(I,1)+U(1)*cp1*cq1
        C(I,2)=C(I,2)+U(2)*cq1
        C(I,3)=C(I,3)+U(3)*cq1
        C(I,4)=C(I,4)+U(4)*cq1
        C(I,5)=C(I,5)+U(5)*cp1*cq3
        C(I,6)=C(I,6)+U(6)*cq3
        C(I,7)=C(I,7)+U(7)*cq3
        C(I,8)=C(I,8)+U(8)*cq3
        C(I,9)=C(I,9)+U(9)*cp1*cq3
        C(I,10)=C(I,10)+U(10)*cq3
        C(I,11)=C(I,11)+U(11)*cq3
        C(I,12)=C(I,12)+U(12)*cq3
        C(I,13)=C(I,13)+U(13)*cp1*cq3
        C(I,14)=C(I,14)+U(14)*cq3
        C(I,15)=C(I,15)+U(15)*cq3
        C(I,16)=C(I,16)+U(16)*cq3
        C(I,17)=C(I,17)+U(5)*cp1*cq2
        C(I,18)=C(I,18)+U(6)*cq2
        C(I,19)=C(I,19)+U(7)*cq2
        C(I,20)=C(I,20)+U(8)*cq2
        C(I,21)=C(I,21)+U(17)*cp1
        C(I,22)=C(I,22)+U(18)
        C(I,23)=C(I,23)+U(19)
        C(I,24)=C(I,24)+U(20)
        C(I,25)=C(I,25)+U(1)*cp1*cq3
        C(I,26)=C(I,26)+U(2)*cq3
        C(I,27)=C(I,27)+U(3)*cq3
        C(I,28)=C(I,28)+U(4)*cq3
        C(I,29)=C(I,29)+U(5)*cp1
        C(I,30)=C(I,30)+U(6)
        C(I,31)=C(I,31)+U(7)
        C(I,32)=C(I,32)+U(8)
        C(I,33)=C(I,33)+U(9)*cp1*cq2
        C(I,34)=C(I,34)+U(10)*cq2
        C(I,35)=C(I,35)+U(11)*cq2
        C(I,36)=C(I,36)+U(12)*cq2
        C(I,37)=C(I,37)+U(21)*cp1
        C(I,38)=C(I,38)+U(22)
        C(I,39)=C(I,39)+U(23)
        C(I,40)=C(I,40)+U(24)
        C(I,41)=C(I,41)+U(25)*cp1
        C(I,42)=C(I,42)+U(26)
        C(I,43)=C(I,43)+U(27)
        C(I,44)=C(I,44)+U(28)
        C(I,45)=C(I,45)+U(9)*cp1
        C(I,46)=C(I,46)+U(10)
        C(I,47)=C(I,47)+U(11)
        C(I,48)=C(I,48)+U(12)
        C(I,49)=C(I,49)+U(13)*cp1*cq2
        C(I,50)=C(I,50)+U(14)*cq2
        C(I,51)=C(I,51)+U(15)*cq2
        C(I,52)=C(I,52)+U(16)*cq2
        C(I,53)=C(I,53)+U(29)*cp1
        C(I,54)=C(I,54)+U(30)
        C(I,55)=C(I,55)+U(31)
        C(I,56)=C(I,56)+U(32)
        C(I,57)=C(I,57)+U(33)*cp1
        C(I,58)=C(I,58)+U(34)
        C(I,59)=C(I,59)+U(35)
        C(I,60)=C(I,60)+U(36)
        C(I,61)=C(I,61)+U(37)*cp1
        C(I,62)=C(I,62)+U(38)
        C(I,63)=C(I,63)+U(39)
        C(I,64)=C(I,64)+U(40)
        C(I,65)=C(I,65)+U(13)*cp1
        C(I,66)=C(I,66)+U(14)
        C(I,67)=C(I,67)+U(15)
        C(I,68)=C(I,68)+U(16)

      END DO ! K
    END DO ! J
  END DO ! I

  ELSEIF(IntCode.EQ.03010302) THEN
  ELSEIF(IntCode.EQ.02010303) THEN
  ELSEIF(IntCode.EQ.03010303) THEN
  ELSEIF(IntCode.EQ.03010202) THEN
  ELSEIF(IntCode.EQ.02010302) THEN
  ELSE
    CALL Halt('Illegal IntCode in Int2122')
  ENDIF
END SUBROUTINE Int2122
