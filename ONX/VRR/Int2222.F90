SUBROUTINE Int2222(N,IntCode,CBra,CKet,DisBufB,PrmBufB,DB,IB,SB,C,U)
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
  REAL(DOUBLE)  :: U(163),C(N,289)
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

  IF(IntCode.EQ.02020202) THEN

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
      U(6)=R2
      U(5)=R3
      U(2)=R4
      U(3)=R5
      U(18)=QBx*U(1)+WQx*U(6)
      U(35)=QBy*U(1)+WQy*U(6)
      U(52)=QBz*U(1)+WQz*U(6)
      U(8)=QBx*U(6)+WQx*U(5)
      U(9)=QBy*U(6)+WQy*U(5)
      U(11)=QBz*U(6)+WQz*U(5)
      U(7)=QBx*U(5)+WQx*U(2)
      U(10)=QBy*U(5)+WQy*U(2)
      U(12)=QBz*U(5)+WQz*U(2)
      U(4)=QBx*U(2)+WQx*U(3)
      U(13)=QBy*U(2)+WQy*U(3)
      U(14)=QBz*U(2)+WQz*U(3)
      U(69)=QBx*U(18)+WQx*U(8)+r1x2E*(U(1)-rpxE*U(6))
      U(86)=QBx*U(35)+WQx*U(9)
      U(103)=QBy*U(35)+WQy*U(9)+r1x2E*(U(1)-rpxE*U(6))
      U(120)=QBx*U(52)+WQx*U(11)
      U(137)=QBy*U(52)+WQy*U(11)
      U(154)=QBz*U(52)+WQz*U(11)+r1x2E*(U(1)-rpxE*U(6))
      U(15)=QBx*U(8)+WQx*U(7)+r1x2E*(U(6)-rpxE*U(5))
      U(16)=QBx*U(9)+WQx*U(10)
      U(17)=QBy*U(9)+WQy*U(10)+r1x2E*(U(6)-rpxE*U(5))
      U(26)=QBx*U(11)+WQx*U(12)
      U(28)=QBy*U(11)+WQy*U(12)
      U(29)=QBz*U(11)+WQz*U(12)+r1x2E*(U(6)-rpxE*U(5))
      U(22)=QBx*U(7)+WQx*U(4)+r1x2E*(U(5)-rpxE*U(2))
      U(23)=QBx*U(10)+WQx*U(13)
      U(24)=QBy*U(10)+WQy*U(13)+r1x2E*(U(5)-rpxE*U(2))
      U(25)=QBx*U(12)+WQx*U(14)
      U(27)=QBy*U(12)+WQy*U(14)
      U(30)=QBz*U(12)+WQz*U(14)+r1x2E*(U(5)-rpxE*U(2))
      U(2)=PAx*U(1)+WPx*U(6)
      U(3)=PAy*U(1)+WPy*U(6)
      U(4)=PAz*U(1)+WPz*U(6)
      U(31)=PAx*U(6)+WPx*U(5)
      U(32)=PAy*U(6)+WPy*U(5)
      U(33)=PAz*U(6)+WPz*U(5)
      U(19)=PAx*U(18)+WPx*U(8)+r1x2ZE*U(6)
      U(20)=PAy*U(18)+WPy*U(8)
      U(21)=PAz*U(18)+WPz*U(8)
      U(36)=PAx*U(35)+WPx*U(9)
      U(37)=PAy*U(35)+WPy*U(9)+r1x2ZE*U(6)
      U(38)=PAz*U(35)+WPz*U(9)
      U(53)=PAx*U(52)+WPx*U(11)
      U(54)=PAy*U(52)+WPy*U(11)
      U(55)=PAz*U(52)+WPz*U(11)+r1x2ZE*U(6)
      U(34)=PAx*U(8)+WPx*U(7)+r1x2ZE*U(5)
      U(39)=PAy*U(8)+WPy*U(7)
      U(40)=PAz*U(8)+WPz*U(7)
      U(42)=PAx*U(9)+WPx*U(10)
      U(45)=PAy*U(9)+WPy*U(10)+r1x2ZE*U(5)
      U(46)=PAz*U(9)+WPz*U(10)
      U(47)=PAx*U(11)+WPx*U(12)
      U(48)=PAy*U(11)+WPy*U(12)
      U(49)=PAz*U(11)+WPz*U(12)+r1x2ZE*U(5)
      U(70)=PAx*U(69)+WPx*U(15)+2.0D0*r1x2ZE*U(8)
      U(71)=PAy*U(69)+WPy*U(15)
      U(72)=PAz*U(69)+WPz*U(15)
      U(87)=PAx*U(86)+WPx*U(16)+r1x2ZE*U(9)
      U(88)=PAy*U(86)+WPy*U(16)+r1x2ZE*U(8)
      U(89)=PAz*U(86)+WPz*U(16)
      U(104)=PAx*U(103)+WPx*U(17)
      U(105)=PAy*U(103)+WPy*U(17)+2.0D0*r1x2ZE*U(9)
      U(106)=PAz*U(103)+WPz*U(17)
      U(121)=PAx*U(120)+WPx*U(26)+r1x2ZE*U(11)
      U(122)=PAy*U(120)+WPy*U(26)
      U(123)=PAz*U(120)+WPz*U(26)+r1x2ZE*U(8)
      U(138)=PAx*U(137)+WPx*U(28)
      U(139)=PAy*U(137)+WPy*U(28)+r1x2ZE*U(11)
      U(140)=PAz*U(137)+WPz*U(28)+r1x2ZE*U(9)
      U(155)=PAx*U(154)+WPx*U(29)
      U(156)=PAy*U(154)+WPy*U(29)
      U(157)=PAz*U(154)+WPz*U(29)+2.0D0*r1x2ZE*U(11)
      U(44)=PAx*U(15)+WPx*U(22)+2.0D0*r1x2ZE*U(7)
      U(50)=PAy*U(15)+WPy*U(22)
      U(51)=PAz*U(15)+WPz*U(22)
      U(56)=PAx*U(16)+WPx*U(23)+r1x2ZE*U(10)
      U(57)=PAy*U(16)+WPy*U(23)+r1x2ZE*U(7)
      U(59)=PAz*U(16)+WPz*U(23)
      U(60)=PAx*U(17)+WPx*U(24)
      U(62)=PAy*U(17)+WPy*U(24)+2.0D0*r1x2ZE*U(10)
      U(63)=PAz*U(17)+WPz*U(24)
      U(58)=PAx*U(26)+WPx*U(25)+r1x2ZE*U(12)
      U(64)=PAy*U(26)+WPy*U(25)
      U(65)=PAz*U(26)+WPz*U(25)+r1x2ZE*U(7)
      U(66)=PAx*U(28)+WPx*U(27)
      U(67)=PAy*U(28)+WPy*U(27)+r1x2ZE*U(12)
      U(68)=PAz*U(28)+WPz*U(27)+r1x2ZE*U(10)
      U(77)=PAx*U(29)+WPx*U(30)
      U(79)=PAy*U(29)+WPy*U(30)
      U(80)=PAz*U(29)+WPz*U(30)+2.0D0*r1x2ZE*U(12)
      U(22)=PAx*U(19)+WPx*U(34)+r1x2Z*(U(18)-rpxZ*U(8))+r1x2ZE*U(31)
      U(23)=PAx*U(20)+WPx*U(39)+r1x2ZE*U(32)
      U(25)=PAx*U(21)+WPx*U(40)+r1x2ZE*U(33)
      U(41)=PAy*U(37)+WPy*U(45)+r1x2Z*(U(35)-rpxZ*U(9))+r1x2ZE*U(32)
      U(43)=PAy*U(38)+WPy*U(46)+r1x2ZE*U(33)
      U(61)=PAz*U(55)+WPz*U(49)+r1x2Z*(U(52)-rpxZ*U(11))+r1x2ZE*U(33)
      U(73)=PAx*U(70)+WPx*U(44)+r1x2Z*(U(69)-rpxZ*U(15))+2.0D0*r1x2ZE*U(34)
      U(74)=PAx*U(71)+WPx*U(50)+2.0D0*r1x2ZE*U(39)
      U(76)=PAx*U(72)+WPx*U(51)+2.0D0*r1x2ZE*U(40)
      U(90)=PAx*U(87)+WPx*U(56)+r1x2Z*(U(86)-rpxZ*U(16))+r1x2ZE*U(42)
      U(91)=PAx*U(88)+WPx*U(57)+r1x2ZE*U(45)
      U(92)=PAy*U(88)+WPy*U(57)+r1x2Z*(U(86)-rpxZ*U(16))+r1x2ZE*U(39)
      U(93)=PAx*U(89)+WPx*U(59)+r1x2ZE*U(46)
      U(94)=PAy*U(89)+WPy*U(59)+r1x2ZE*U(40)
      U(109)=PAy*U(105)+WPy*U(62)+r1x2Z*(U(103)-rpxZ*U(17))+2.0D0*r1x2ZE*U(45)
      U(111)=PAy*U(106)+WPy*U(63)+2.0D0*r1x2ZE*U(46)
      U(124)=PAx*U(121)+WPx*U(58)+r1x2Z*(U(120)-rpxZ*U(26))+r1x2ZE*U(47)
      U(125)=PAx*U(122)+WPx*U(64)+r1x2ZE*U(48)
      U(127)=PAx*U(123)+WPx*U(65)+r1x2ZE*U(49)
      U(129)=PAz*U(123)+WPz*U(65)+r1x2Z*(U(120)-rpxZ*U(26))+r1x2ZE*U(40)
      U(143)=PAy*U(139)+WPy*U(67)+r1x2Z*(U(137)-rpxZ*U(28))+r1x2ZE*U(48)
      U(145)=PAy*U(140)+WPy*U(68)+r1x2ZE*U(49)
      U(146)=PAz*U(140)+WPz*U(68)+r1x2Z*(U(137)-rpxZ*U(28))+r1x2ZE*U(46)
      U(163)=PAz*U(157)+WPz*U(80)+r1x2Z*(U(154)-rpxZ*U(29))+2.0D0*r1x2ZE*U(49)
      U(5)=PAx*U(2)+WPx*U(31)+r1x2Z*(U(1)-rpxZ*U(6))
      U(7)=PAy*U(3)+WPy*U(32)+r1x2Z*(U(1)-rpxZ*U(6))
      U(10)=PAz*U(4)+WPz*U(33)+r1x2Z*(U(1)-rpxZ*U(6))
      U(24)=PAy*U(20)+WPy*U(39)+r1x2Z*(U(18)-rpxZ*U(8))
      U(27)=PAz*U(21)+WPz*U(40)+r1x2Z*(U(18)-rpxZ*U(8))
      U(39)=PAx*U(36)+WPx*U(42)+r1x2Z*(U(35)-rpxZ*U(9))
      U(44)=PAz*U(38)+WPz*U(46)+r1x2Z*(U(35)-rpxZ*U(9))
      U(56)=PAx*U(53)+WPx*U(47)+r1x2Z*(U(52)-rpxZ*U(11))
      U(58)=PAy*U(54)+WPy*U(48)+r1x2Z*(U(52)-rpxZ*U(11))
      U(75)=PAy*U(71)+WPy*U(50)+r1x2Z*(U(69)-rpxZ*U(15))
      U(78)=PAz*U(72)+WPz*U(51)+r1x2Z*(U(69)-rpxZ*U(15))
      U(95)=PAz*U(89)+WPz*U(59)+r1x2Z*(U(86)-rpxZ*U(16))
      U(107)=PAx*U(104)+WPx*U(60)+r1x2Z*(U(103)-rpxZ*U(17))
      U(112)=PAz*U(106)+WPz*U(63)+r1x2Z*(U(103)-rpxZ*U(17))
      U(126)=PAy*U(122)+WPy*U(64)+r1x2Z*(U(120)-rpxZ*U(26))
      U(141)=PAx*U(138)+WPx*U(66)+r1x2Z*(U(137)-rpxZ*U(28))
      U(158)=PAx*U(155)+WPx*U(77)+r1x2Z*(U(154)-rpxZ*U(29))
      U(160)=PAy*U(156)+WPy*U(79)+r1x2Z*(U(154)-rpxZ*U(29))
      U(6)=PAx*U(3)+WPx*U(32)
      U(8)=PAx*U(4)+WPx*U(33)
      U(9)=PAy*U(4)+WPy*U(33)
      U(26)=PAy*U(21)+WPy*U(40)
      U(40)=PAx*U(37)+WPx*U(45)
      U(42)=PAx*U(38)+WPx*U(46)
      U(57)=PAx*U(54)+WPx*U(48)
      U(59)=PAx*U(55)+WPx*U(49)
      U(60)=PAy*U(55)+WPy*U(49)
      U(77)=PAy*U(72)+WPy*U(51)
      U(108)=PAx*U(105)+WPx*U(62)
      U(110)=PAx*U(106)+WPx*U(63)
      U(128)=PAy*U(123)+WPy*U(65)
      U(142)=PAx*U(139)+WPx*U(67)
      U(144)=PAx*U(140)+WPx*U(68)
      U(159)=PAx*U(156)+WPx*U(79)
      U(161)=PAx*U(157)+WPx*U(80)
      U(162)=PAy*U(157)+WPy*U(80)

      C(I,1)=C(I,1)+U(1)*cp1*cq1
      C(I,2)=C(I,2)+U(2)*cp3*cq1
      C(I,3)=C(I,3)+U(3)*cp3*cq1
      C(I,4)=C(I,4)+U(4)*cp3*cq1
      C(I,5)=C(I,5)+U(2)*cp2*cq1
      C(I,6)=C(I,6)+U(5)*cq1
      C(I,7)=C(I,7)+U(1)*cp3*cq1
      C(I,8)=C(I,8)+U(2)*cq1
      C(I,9)=C(I,9)+U(3)*cp2*cq1
      C(I,10)=C(I,10)+U(6)*cq1
      C(I,11)=C(I,11)+U(7)*cq1
      C(I,12)=C(I,12)+U(3)*cq1
      C(I,13)=C(I,13)+U(4)*cp2*cq1
      C(I,14)=C(I,14)+U(8)*cq1
      C(I,15)=C(I,15)+U(9)*cq1
      C(I,16)=C(I,16)+U(10)*cq1
      C(I,17)=C(I,17)+U(4)*cq1
      C(I,18)=C(I,18)+U(18)*cp1*cq3
      C(I,19)=C(I,19)+U(19)*cp3*cq3
      C(I,20)=C(I,20)+U(20)*cp3*cq3
      C(I,21)=C(I,21)+U(21)*cp3*cq3
      C(I,22)=C(I,22)+U(19)*cp2*cq3
      C(I,23)=C(I,23)+U(22)*cq3
      C(I,24)=C(I,24)+U(18)*cp3*cq3
      C(I,25)=C(I,25)+U(19)*cq3
      C(I,26)=C(I,26)+U(20)*cp2*cq3
      C(I,27)=C(I,27)+U(23)*cq3
      C(I,28)=C(I,28)+U(24)*cq3
      C(I,29)=C(I,29)+U(20)*cq3
      C(I,30)=C(I,30)+U(21)*cp2*cq3
      C(I,31)=C(I,31)+U(25)*cq3
      C(I,32)=C(I,32)+U(26)*cq3
      C(I,33)=C(I,33)+U(27)*cq3
      C(I,34)=C(I,34)+U(21)*cq3
      C(I,35)=C(I,35)+U(35)*cp1*cq3
      C(I,36)=C(I,36)+U(36)*cp3*cq3
      C(I,37)=C(I,37)+U(37)*cp3*cq3
      C(I,38)=C(I,38)+U(38)*cp3*cq3
      C(I,39)=C(I,39)+U(36)*cp2*cq3
      C(I,40)=C(I,40)+U(39)*cq3
      C(I,41)=C(I,41)+U(35)*cp3*cq3
      C(I,42)=C(I,42)+U(36)*cq3
      C(I,43)=C(I,43)+U(37)*cp2*cq3
      C(I,44)=C(I,44)+U(40)*cq3
      C(I,45)=C(I,45)+U(41)*cq3
      C(I,46)=C(I,46)+U(37)*cq3
      C(I,47)=C(I,47)+U(38)*cp2*cq3
      C(I,48)=C(I,48)+U(42)*cq3
      C(I,49)=C(I,49)+U(43)*cq3
      C(I,50)=C(I,50)+U(44)*cq3
      C(I,51)=C(I,51)+U(38)*cq3
      C(I,52)=C(I,52)+U(52)*cp1*cq3
      C(I,53)=C(I,53)+U(53)*cp3*cq3
      C(I,54)=C(I,54)+U(54)*cp3*cq3
      C(I,55)=C(I,55)+U(55)*cp3*cq3
      C(I,56)=C(I,56)+U(53)*cp2*cq3
      C(I,57)=C(I,57)+U(56)*cq3
      C(I,58)=C(I,58)+U(52)*cp3*cq3
      C(I,59)=C(I,59)+U(53)*cq3
      C(I,60)=C(I,60)+U(54)*cp2*cq3
      C(I,61)=C(I,61)+U(57)*cq3
      C(I,62)=C(I,62)+U(58)*cq3
      C(I,63)=C(I,63)+U(54)*cq3
      C(I,64)=C(I,64)+U(55)*cp2*cq3
      C(I,65)=C(I,65)+U(59)*cq3
      C(I,66)=C(I,66)+U(60)*cq3
      C(I,67)=C(I,67)+U(61)*cq3
      C(I,68)=C(I,68)+U(55)*cq3
      C(I,69)=C(I,69)+U(18)*cp1*cq2
      C(I,70)=C(I,70)+U(19)*cp3*cq2
      C(I,71)=C(I,71)+U(20)*cp3*cq2
      C(I,72)=C(I,72)+U(21)*cp3*cq2
      C(I,73)=C(I,73)+U(19)*cp2*cq2
      C(I,74)=C(I,74)+U(22)*cq2
      C(I,75)=C(I,75)+U(18)*cp3*cq2
      C(I,76)=C(I,76)+U(19)*cq2
      C(I,77)=C(I,77)+U(20)*cp2*cq2
      C(I,78)=C(I,78)+U(23)*cq2
      C(I,79)=C(I,79)+U(24)*cq2
      C(I,80)=C(I,80)+U(20)*cq2
      C(I,81)=C(I,81)+U(21)*cp2*cq2
      C(I,82)=C(I,82)+U(25)*cq2
      C(I,83)=C(I,83)+U(26)*cq2
      C(I,84)=C(I,84)+U(27)*cq2
      C(I,85)=C(I,85)+U(21)*cq2
      C(I,86)=C(I,86)+U(69)*cp1
      C(I,87)=C(I,87)+U(70)*cp3
      C(I,88)=C(I,88)+U(71)*cp3
      C(I,89)=C(I,89)+U(72)*cp3
      C(I,90)=C(I,90)+U(70)*cp2
      C(I,91)=C(I,91)+U(73)
      C(I,92)=C(I,92)+U(69)*cp3
      C(I,93)=C(I,93)+U(70)
      C(I,94)=C(I,94)+U(71)*cp2
      C(I,95)=C(I,95)+U(74)
      C(I,96)=C(I,96)+U(75)
      C(I,97)=C(I,97)+U(71)
      C(I,98)=C(I,98)+U(72)*cp2
      C(I,99)=C(I,99)+U(76)
      C(I,100)=C(I,100)+U(77)
      C(I,101)=C(I,101)+U(78)
      C(I,102)=C(I,102)+U(72)
      C(I,103)=C(I,103)+U(1)*cp1*cq3
      C(I,104)=C(I,104)+U(2)*cp3*cq3
      C(I,105)=C(I,105)+U(3)*cp3*cq3
      C(I,106)=C(I,106)+U(4)*cp3*cq3
      C(I,107)=C(I,107)+U(2)*cp2*cq3
      C(I,108)=C(I,108)+U(5)*cq3
      C(I,109)=C(I,109)+U(1)*cp3*cq3
      C(I,110)=C(I,110)+U(2)*cq3
      C(I,111)=C(I,111)+U(3)*cp2*cq3
      C(I,112)=C(I,112)+U(6)*cq3
      C(I,113)=C(I,113)+U(7)*cq3
      C(I,114)=C(I,114)+U(3)*cq3
      C(I,115)=C(I,115)+U(4)*cp2*cq3
      C(I,116)=C(I,116)+U(8)*cq3
      C(I,117)=C(I,117)+U(9)*cq3
      C(I,118)=C(I,118)+U(10)*cq3
      C(I,119)=C(I,119)+U(4)*cq3
      C(I,120)=C(I,120)+U(18)*cp1
      C(I,121)=C(I,121)+U(19)*cp3
      C(I,122)=C(I,122)+U(20)*cp3
      C(I,123)=C(I,123)+U(21)*cp3
      C(I,124)=C(I,124)+U(19)*cp2
      C(I,125)=C(I,125)+U(22)
      C(I,126)=C(I,126)+U(18)*cp3
      C(I,127)=C(I,127)+U(19)
      C(I,128)=C(I,128)+U(20)*cp2
      C(I,129)=C(I,129)+U(23)
      C(I,130)=C(I,130)+U(24)
      C(I,131)=C(I,131)+U(20)
      C(I,132)=C(I,132)+U(21)*cp2
      C(I,133)=C(I,133)+U(25)
      C(I,134)=C(I,134)+U(26)
      C(I,135)=C(I,135)+U(27)
      C(I,136)=C(I,136)+U(21)
      C(I,137)=C(I,137)+U(35)*cp1*cq2
      C(I,138)=C(I,138)+U(36)*cp3*cq2
      C(I,139)=C(I,139)+U(37)*cp3*cq2
      C(I,140)=C(I,140)+U(38)*cp3*cq2
      C(I,141)=C(I,141)+U(36)*cp2*cq2
      C(I,142)=C(I,142)+U(39)*cq2
      C(I,143)=C(I,143)+U(35)*cp3*cq2
      C(I,144)=C(I,144)+U(36)*cq2
      C(I,145)=C(I,145)+U(37)*cp2*cq2
      C(I,146)=C(I,146)+U(40)*cq2
      C(I,147)=C(I,147)+U(41)*cq2
      C(I,148)=C(I,148)+U(37)*cq2
      C(I,149)=C(I,149)+U(38)*cp2*cq2
      C(I,150)=C(I,150)+U(42)*cq2
      C(I,151)=C(I,151)+U(43)*cq2
      C(I,152)=C(I,152)+U(44)*cq2
      C(I,153)=C(I,153)+U(38)*cq2
      C(I,154)=C(I,154)+U(86)*cp1
      C(I,155)=C(I,155)+U(87)*cp3
      C(I,156)=C(I,156)+U(88)*cp3
      C(I,157)=C(I,157)+U(89)*cp3
      C(I,158)=C(I,158)+U(87)*cp2
      C(I,159)=C(I,159)+U(90)
      C(I,160)=C(I,160)+U(86)*cp3
      C(I,161)=C(I,161)+U(87)
      C(I,162)=C(I,162)+U(88)*cp2
      C(I,163)=C(I,163)+U(91)
      C(I,164)=C(I,164)+U(92)
      C(I,165)=C(I,165)+U(88)
      C(I,166)=C(I,166)+U(89)*cp2
      C(I,167)=C(I,167)+U(93)
      C(I,168)=C(I,168)+U(94)
      C(I,169)=C(I,169)+U(95)
      C(I,170)=C(I,170)+U(89)
      C(I,171)=C(I,171)+U(103)*cp1
      C(I,172)=C(I,172)+U(104)*cp3
      C(I,173)=C(I,173)+U(105)*cp3
      C(I,174)=C(I,174)+U(106)*cp3
      C(I,175)=C(I,175)+U(104)*cp2
      C(I,176)=C(I,176)+U(107)
      C(I,177)=C(I,177)+U(103)*cp3
      C(I,178)=C(I,178)+U(104)
      C(I,179)=C(I,179)+U(105)*cp2
      C(I,180)=C(I,180)+U(108)
      C(I,181)=C(I,181)+U(109)
      C(I,182)=C(I,182)+U(105)
      C(I,183)=C(I,183)+U(106)*cp2
      C(I,184)=C(I,184)+U(110)
      C(I,185)=C(I,185)+U(111)
      C(I,186)=C(I,186)+U(112)
      C(I,187)=C(I,187)+U(106)
      C(I,188)=C(I,188)+U(35)*cp1
      C(I,189)=C(I,189)+U(36)*cp3
      C(I,190)=C(I,190)+U(37)*cp3
      C(I,191)=C(I,191)+U(38)*cp3
      C(I,192)=C(I,192)+U(36)*cp2
      C(I,193)=C(I,193)+U(39)
      C(I,194)=C(I,194)+U(35)*cp3
      C(I,195)=C(I,195)+U(36)
      C(I,196)=C(I,196)+U(37)*cp2
      C(I,197)=C(I,197)+U(40)
      C(I,198)=C(I,198)+U(41)
      C(I,199)=C(I,199)+U(37)
      C(I,200)=C(I,200)+U(38)*cp2
      C(I,201)=C(I,201)+U(42)
      C(I,202)=C(I,202)+U(43)
      C(I,203)=C(I,203)+U(44)
      C(I,204)=C(I,204)+U(38)
      C(I,205)=C(I,205)+U(52)*cp1*cq2
      C(I,206)=C(I,206)+U(53)*cp3*cq2
      C(I,207)=C(I,207)+U(54)*cp3*cq2
      C(I,208)=C(I,208)+U(55)*cp3*cq2
      C(I,209)=C(I,209)+U(53)*cp2*cq2
      C(I,210)=C(I,210)+U(56)*cq2
      C(I,211)=C(I,211)+U(52)*cp3*cq2
      C(I,212)=C(I,212)+U(53)*cq2
      C(I,213)=C(I,213)+U(54)*cp2*cq2
      C(I,214)=C(I,214)+U(57)*cq2
      C(I,215)=C(I,215)+U(58)*cq2
      C(I,216)=C(I,216)+U(54)*cq2
      C(I,217)=C(I,217)+U(55)*cp2*cq2
      C(I,218)=C(I,218)+U(59)*cq2
      C(I,219)=C(I,219)+U(60)*cq2
      C(I,220)=C(I,220)+U(61)*cq2
      C(I,221)=C(I,221)+U(55)*cq2
      C(I,222)=C(I,222)+U(120)*cp1
      C(I,223)=C(I,223)+U(121)*cp3
      C(I,224)=C(I,224)+U(122)*cp3
      C(I,225)=C(I,225)+U(123)*cp3
      C(I,226)=C(I,226)+U(121)*cp2
      C(I,227)=C(I,227)+U(124)
      C(I,228)=C(I,228)+U(120)*cp3
      C(I,229)=C(I,229)+U(121)
      C(I,230)=C(I,230)+U(122)*cp2
      C(I,231)=C(I,231)+U(125)
      C(I,232)=C(I,232)+U(126)
      C(I,233)=C(I,233)+U(122)
      C(I,234)=C(I,234)+U(123)*cp2
      C(I,235)=C(I,235)+U(127)
      C(I,236)=C(I,236)+U(128)
      C(I,237)=C(I,237)+U(129)
      C(I,238)=C(I,238)+U(123)
      C(I,239)=C(I,239)+U(137)*cp1
      C(I,240)=C(I,240)+U(138)*cp3
      C(I,241)=C(I,241)+U(139)*cp3
      C(I,242)=C(I,242)+U(140)*cp3
      C(I,243)=C(I,243)+U(138)*cp2
      C(I,244)=C(I,244)+U(141)
      C(I,245)=C(I,245)+U(137)*cp3
      C(I,246)=C(I,246)+U(138)
      C(I,247)=C(I,247)+U(139)*cp2
      C(I,248)=C(I,248)+U(142)
      C(I,249)=C(I,249)+U(143)
      C(I,250)=C(I,250)+U(139)
      C(I,251)=C(I,251)+U(140)*cp2
      C(I,252)=C(I,252)+U(144)
      C(I,253)=C(I,253)+U(145)
      C(I,254)=C(I,254)+U(146)
      C(I,255)=C(I,255)+U(140)
      C(I,256)=C(I,256)+U(154)*cp1
      C(I,257)=C(I,257)+U(155)*cp3
      C(I,258)=C(I,258)+U(156)*cp3
      C(I,259)=C(I,259)+U(157)*cp3
      C(I,260)=C(I,260)+U(155)*cp2
      C(I,261)=C(I,261)+U(158)
      C(I,262)=C(I,262)+U(154)*cp3
      C(I,263)=C(I,263)+U(155)
      C(I,264)=C(I,264)+U(156)*cp2
      C(I,265)=C(I,265)+U(159)
      C(I,266)=C(I,266)+U(160)
      C(I,267)=C(I,267)+U(156)
      C(I,268)=C(I,268)+U(157)*cp2
      C(I,269)=C(I,269)+U(161)
      C(I,270)=C(I,270)+U(162)
      C(I,271)=C(I,271)+U(163)
      C(I,272)=C(I,272)+U(157)
      C(I,273)=C(I,273)+U(52)*cp1
      C(I,274)=C(I,274)+U(53)*cp3
      C(I,275)=C(I,275)+U(54)*cp3
      C(I,276)=C(I,276)+U(55)*cp3
      C(I,277)=C(I,277)+U(53)*cp2
      C(I,278)=C(I,278)+U(56)
      C(I,279)=C(I,279)+U(52)*cp3
      C(I,280)=C(I,280)+U(53)
      C(I,281)=C(I,281)+U(54)*cp2
      C(I,282)=C(I,282)+U(57)
      C(I,283)=C(I,283)+U(58)
      C(I,284)=C(I,284)+U(54)
      C(I,285)=C(I,285)+U(55)*cp2
      C(I,286)=C(I,286)+U(59)
      C(I,287)=C(I,287)+U(60)
      C(I,288)=C(I,288)+U(61)
      C(I,289)=C(I,289)+U(55)

      END DO ! K
    END DO ! J
  END DO ! I

  ELSEIF(IntCode.EQ.03020202) THEN
  ELSEIF(IntCode.EQ.02020302) THEN
  ELSEIF(IntCode.EQ.03030202) THEN
  ELSEIF(IntCode.EQ.02020303) THEN
  ELSEIF(IntCode.EQ.03020302) THEN
  ELSEIF(IntCode.EQ.03030302) THEN
  ELSEIF(IntCode.EQ.03020303) THEN
  ELSEIF(IntCode.EQ.03030303) THEN
  ELSE
    CALL Halt('Illegal IntCode in Int2222')
  ENDIF
END SUBROUTINE Int2222
