SUBROUTINE Int6611(N,IntCode,CBra,CKet,DisBufB,PrmBufB,DB,IB,SB,C,U)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE GammaF4
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
  REAL(DOUBLE)  :: U(48),C(N,31)
  REAL(DOUBLE)  :: Cx,Cy,Cz
  REAL(DOUBLE)  :: Zeta,Up,Px,Py,Pz,cp1,cp2,cp3
  REAL(DOUBLE)  :: Eta, Uq,Qx,Qy,Qz,cq1,cq2,cq3
  REAL(DOUBLE)  :: EtaQx,EtaQy,EtaQz
  REAL(DOUBLE)  :: r1xZpE,ZpE,Rkk,r1xZ,r1x2Z,rpxZ
  REAL(DOUBLE)  :: PAx,PAy,PAz,PQx,PQy,PQz
  REAL(DOUBLE)  :: QBx,QBy,QBz,WQx,WQy,WQz
  REAL(DOUBLE)  :: Wx,Wy,Wz,WPx,WPy,WPz
  REAL(DOUBLE)  :: T1,T2,T3,T4,TwoT
  REAL(DOUBLE)  :: R1,R2,R3,R4,R5,G0,G1,G2,G3,G4,ET
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

  IF(IntCode.EQ.06060101) THEN

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
          G4=F4_0(MG) +T1*F4_1(MG) +T2*F4_2(MG) +T3*F4_3(MG) +T4*F4_4(MG)
          G3=.142857142857143D+00*(TwoT*G4+ET)
          G2=.200000000000000D+00*(TwoT*G3+ET)
          G1=.333333333333333D+00*(TwoT*G2+ET)
          G0=TwoT*G1+ET
          R1=Rkk*G0
          R2=Rkk*G1
          R3=Rkk*G2
          R4=Rkk*G3
          R5=Rkk*G4
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
          T3=T3*T2
          R5=Rkk*IB%GammaA%D(5)*T3
        ENDIF

      U(1)=R1
      U(11)=R2
      U(12)=R3
      U(13)=R4
      U(5)=R5
      U(2)=PAx*U(1)+WPx*U(11)
      U(3)=PAy*U(1)+WPy*U(11)
      U(4)=PAz*U(1)+WPz*U(11)
      U(21)=PAx*U(11)+WPx*U(12)
      U(22)=PAy*U(11)+WPy*U(12)
      U(23)=PAz*U(11)+WPz*U(12)
      U(24)=PAx*U(12)+WPx*U(13)
      U(25)=PAy*U(12)+WPy*U(13)
      U(26)=PAz*U(12)+WPz*U(13)
      U(14)=PAx*U(13)+WPx*U(5)
      U(15)=PAy*U(13)+WPy*U(5)
      U(16)=PAz*U(13)+WPz*U(5)
      U(5)=PAx*U(2)+WPx*U(21)+r1x2Z*(U(1)-rpxZ*U(11))
      U(6)=PAx*U(3)+WPx*U(22)
      U(7)=PAy*U(3)+WPy*U(22)+r1x2Z*(U(1)-rpxZ*U(11))
      U(8)=PAx*U(4)+WPx*U(23)
      U(9)=PAy*U(4)+WPy*U(23)
      U(10)=PAz*U(4)+WPz*U(23)+r1x2Z*(U(1)-rpxZ*U(11))
      U(27)=PAx*U(21)+WPx*U(24)+r1x2Z*(U(11)-rpxZ*U(12))
      U(28)=PAx*U(22)+WPx*U(25)
      U(29)=PAy*U(22)+WPy*U(25)+r1x2Z*(U(11)-rpxZ*U(12))
      U(30)=PAx*U(23)+WPx*U(26)
      U(31)=PAy*U(23)+WPy*U(26)
      U(33)=PAz*U(23)+WPz*U(26)+r1x2Z*(U(11)-rpxZ*U(12))
      U(32)=PAx*U(24)+WPx*U(14)+r1x2Z*(U(12)-rpxZ*U(13))
      U(34)=PAx*U(25)+WPx*U(15)
      U(35)=PAy*U(25)+WPy*U(15)+r1x2Z*(U(12)-rpxZ*U(13))
      U(36)=PAx*U(26)+WPx*U(16)
      U(37)=PAy*U(26)+WPy*U(16)
      U(38)=PAz*U(26)+WPz*U(16)+r1x2Z*(U(12)-rpxZ*U(13))
      U(11)=PAx*U(5)+WPx*U(27)+2.0D0*r1x2Z*(U(2)-rpxZ*U(21))
      U(12)=PAx*U(6)+WPx*U(28)+r1x2Z*(U(3)-rpxZ*U(22))
      U(13)=PAx*U(7)+WPx*U(29)
      U(14)=PAy*U(7)+WPy*U(29)+2.0D0*r1x2Z*(U(3)-rpxZ*U(22))
      U(15)=PAx*U(8)+WPx*U(30)+r1x2Z*(U(4)-rpxZ*U(23))
      U(16)=PAx*U(9)+WPx*U(31)
      U(17)=PAy*U(9)+WPy*U(31)+r1x2Z*(U(4)-rpxZ*U(23))
      U(18)=PAx*U(10)+WPx*U(33)
      U(19)=PAy*U(10)+WPy*U(33)
      U(20)=PAz*U(10)+WPz*U(33)+2.0D0*r1x2Z*(U(4)-rpxZ*U(23))
      U(39)=PAx*U(27)+WPx*U(32)+2.0D0*r1x2Z*(U(21)-rpxZ*U(24))
      U(40)=PAx*U(28)+WPx*U(34)+r1x2Z*(U(22)-rpxZ*U(25))
      U(41)=PAx*U(29)+WPx*U(35)
      U(42)=PAy*U(29)+WPy*U(35)+2.0D0*r1x2Z*(U(22)-rpxZ*U(25))
      U(43)=PAx*U(30)+WPx*U(36)+r1x2Z*(U(23)-rpxZ*U(26))
      U(44)=PAx*U(31)+WPx*U(37)
      U(45)=PAy*U(31)+WPy*U(37)+r1x2Z*(U(23)-rpxZ*U(26))
      U(46)=PAx*U(33)+WPx*U(38)
      U(47)=PAy*U(33)+WPy*U(38)
      U(48)=PAz*U(33)+WPz*U(38)+2.0D0*r1x2Z*(U(23)-rpxZ*U(26))
      U(21)=PAx*U(11)+WPx*U(39)+3.0D0*r1x2Z*(U(5)-rpxZ*U(27))
      U(22)=PAx*U(12)+WPx*U(40)+2.0D0*r1x2Z*(U(6)-rpxZ*U(28))
      U(23)=PAx*U(13)+WPx*U(41)+r1x2Z*(U(7)-rpxZ*U(29))
      U(25)=PAy*U(14)+WPy*U(42)+3.0D0*r1x2Z*(U(7)-rpxZ*U(29))
      U(26)=PAx*U(15)+WPx*U(43)+2.0D0*r1x2Z*(U(8)-rpxZ*U(30))
      U(27)=PAx*U(16)+WPx*U(44)+r1x2Z*(U(9)-rpxZ*U(31))
      U(29)=PAy*U(17)+WPy*U(45)+2.0D0*r1x2Z*(U(9)-rpxZ*U(31))
      U(30)=PAx*U(18)+WPx*U(46)+r1x2Z*(U(10)-rpxZ*U(33))
      U(32)=PAy*U(19)+WPy*U(47)+r1x2Z*(U(10)-rpxZ*U(33))
      U(35)=PAz*U(20)+WPz*U(48)+3.0D0*r1x2Z*(U(10)-rpxZ*U(33))
      U(24)=PAx*U(14)+WPx*U(42)
      U(28)=PAx*U(17)+WPx*U(45)
      U(31)=PAx*U(19)+WPx*U(47)
      U(33)=PAx*U(20)+WPx*U(48)
      U(34)=PAy*U(20)+WPy*U(48)

      C(I,1)=C(I,1)+U(5)
      C(I,2)=C(I,2)+U(6)
      C(I,3)=C(I,3)+U(7)
      C(I,4)=C(I,4)+U(8)
      C(I,5)=C(I,5)+U(9)
      C(I,6)=C(I,6)+U(10)
      C(I,7)=C(I,7)+U(11)
      C(I,8)=C(I,8)+U(12)
      C(I,9)=C(I,9)+U(13)
      C(I,10)=C(I,10)+U(14)
      C(I,11)=C(I,11)+U(15)
      C(I,12)=C(I,12)+U(16)
      C(I,13)=C(I,13)+U(17)
      C(I,14)=C(I,14)+U(18)
      C(I,15)=C(I,15)+U(19)
      C(I,16)=C(I,16)+U(20)
      C(I,17)=C(I,17)+U(21)
      C(I,18)=C(I,18)+U(22)
      C(I,19)=C(I,19)+U(23)
      C(I,20)=C(I,20)+U(24)
      C(I,21)=C(I,21)+U(25)
      C(I,22)=C(I,22)+U(26)
      C(I,23)=C(I,23)+U(27)
      C(I,24)=C(I,24)+U(28)
      C(I,25)=C(I,25)+U(29)
      C(I,26)=C(I,26)+U(30)
      C(I,27)=C(I,27)+U(31)
      C(I,28)=C(I,28)+U(32)
      C(I,29)=C(I,29)+U(33)
      C(I,30)=C(I,30)+U(34)
      C(I,31)=C(I,31)+U(35)

      END DO ! K
    END DO ! J
  END DO ! I

  ELSE
    WRITE(*,*) "IntCode = ",IntCode
    CALL Halt('Illegal IntCode in Int6611')
  ENDIF
END SUBROUTINE Int6611
