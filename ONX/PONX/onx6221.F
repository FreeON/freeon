      SUBROUTINE ONX6221(N,iBgn,iEnd,IndA,BfnBufB,BfnBufK,KBra,KKet,
     &                   SwhB,SwhK,RAC,RBD,
     &                   PrmBufB,PrmBufK,U,C)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 U(97),C(N,94)
      REAL*8 F4_0(0:486),F4_1(0:486),F4_2(0:486),F4_3(0:486)
      REAL*8 EX_0(0:486),EX_1(0:486),EX_2(0:486),EX_3(0:486)
      REAL*8 PrmBufB(MXB,KBra+MInfo,*)
      REAL*8 PrmBufK(MXB,KKet+MInfo,*)
      REAL*8 RAC(3,*),RBD(3,*)
      INTEGER IndA(0:*)
      INTEGER BfnBufB(*),BfnBufK(*)
      INCLUDE "bkonx.inc"
      INCLUDE "Gamma_A.inc"
      INCLUDE "Gamma_4.inc"
      INCLUDE "Exponential.inc"

      CALL RLoader(C,N*94,0.0D0)

      I0=0
      DO 4000 I=iBgn,iEnd
         Is=I-iBgn+1
         Ixa = BfnBufB(I)
         DO 3000 J=1,IndA(I)
            Ixb = BfnBufK(J)
            IF(Ixa.GE.Ixb) THEN
               I0=I0+1
               rDx=Dx+SwhK*RBD(1,J)
               rDy=Dy+SwhK*RBD(2,J)
               rDz=Dz+SwhK*RBD(3,J)
               rCx=Cx+SwhB*RAC(1,Is)
               rCy=Cy+SwhB*RAC(2,Is)
               rCz=Cz+SwhB*RAC(3,Is)
               DO 2000 IKK=1,KKet
                  Eta  = PrmBufK(1,IKK,J)
                  Qx   = PrmBufK(2,IKK,J)
                  Qy   = PrmBufK(3,IKK,J)
                  Qz   = PrmBufK(4,IKK,J)
                  Uq   = PrmBufK(5,IKK,J)
                  cq1  = PrmBufK(6,IKK,J)
                  schK = PrmBufK(9,IKK,J)
                  TTsh = Tsh/(Dcd*schK)
                  r1xe=1.0D0/Eta
                  r1x2e=0.5D0*r1xe
                  EtaQx=Eta*Qx
                  EtaQy=Eta*Qy
                  EtaQz=Eta*Qz
                  QBx = Qx-rDx
                  QBy = Qy-rDy
                  QBz = Qz-rDz
                  DO 1000 IKB=1,KBra
                     Zeta = PrmBufB(1,IKB,I)
                     Px   = PrmBufB(2,IKB,I)
                     Py   = PrmBufB(3,IKB,I)
                     Pz   = PrmBufB(4,IKB,I)
                     Up   = PrmBufB(5,IKB,I)
                     cp1  = PrmBufB(6,IKB,I)
                     schB = PrmBufB(9,IKB,I)
                     IF(SchB.LE.TTsh) GOTO 1010
                     ZpE=Zeta+Eta
                     r1xZpE=1.0D0/ZpE
                     r1x2ZE=0.5D0*r1xZpE
                     r1xZ=1.0D0/Zeta
                     r1x2Z=0.5D0*r1xZ
                     rpxZ=Eta*r1xZpE
                     rpxE=Zeta*r1xZpE
                     Rkk=Up*Uq*DSQRT(r1xZpE)
                     PAx = Px-rCx
                     PAy = Py-rCy
                     PAz = Pz-rCz
                     Wx  = (Zeta*Px+EtaQx)*r1xZpE
                     Wy  = (Zeta*Py+EtaQy)*r1xZpE
                     Wz  = (Zeta*Pz+EtaQz)*r1xZpE
                     WPx = Wx - Px
                     WPy = Wy - Py
                     WPz = Wz - Pz
                     WQx = Wx - Qx
                     WQy = Wy - Qy
                     WQz = Wz - Qz
                     PQx = Px - Qx
                     PQy = Py - Qy
                     PQz = Pz - Qz
                     T=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
                     IF(T.LT.Switch) THEN
                        m=INT(T*Grid)
                        T2=T*T
                        T3=T2*T
                        G4=F4_0(m)+T*F4_1(m)+T2*F4_2(m)+T3*F4_3(m)
                        ET=EX_0(m)+T*EX_1(m)+T2*EX_2(m)+T3*EX_3(m)
                        TwoT=2.0D0*T
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
                        T1=1.0D0/T
                        TS=DSQRT(T1)
                        R1=Rkk*F0Asymp*Ts
                        Ts=Ts*T1
                        R2=Rkk*F1Asymp*Ts
                        Ts=Ts*T1
                        R3=Rkk*F2Asymp*Ts
                        Ts=Ts*T1
                        R4=Rkk*F3Asymp*Ts
                        Ts=Ts*T1
                        R5=Rkk*F4Asymp*Ts
                     ENDIF

      U(1)=R1
      U(11)=R2
      U(12)=R3
      U(5)=R4
      U(2)=R5
      U(25)=QBx*U(1)+WQx*U(11)
      U(49)=QBy*U(1)+WQy*U(11)
      U(73)=QBz*U(1)+WQz*U(11)
      U(13)=QBx*U(11)+WQx*U(12)
      U(14)=QBy*U(11)+WQy*U(12)
      U(15)=QBz*U(11)+WQz*U(12)
      U(16)=QBx*U(12)+WQx*U(5)
      U(17)=QBy*U(12)+WQy*U(5)
      U(18)=QBz*U(12)+WQz*U(5)
      U(6)=QBx*U(5)+WQx*U(2)
      U(7)=QBy*U(5)+WQy*U(2)
      U(8)=QBz*U(5)+WQz*U(2)
      U(2)=PAx*U(1)+WPx*U(11)
      U(3)=PAy*U(1)+WPy*U(11)
      U(4)=PAz*U(1)+WPz*U(11)
      U(19)=PAx*U(11)+WPx*U(12)
      U(20)=PAy*U(11)+WPy*U(12)
      U(21)=PAz*U(11)+WPz*U(12)
      U(22)=PAx*U(12)+WPx*U(5)
      U(23)=PAy*U(12)+WPy*U(5)
      U(24)=PAz*U(12)+WPz*U(5)
      U(26)=PAx*U(25)+WPx*U(13)+
     &  r1x2ZE*U(11)
      U(27)=PAy*U(25)+WPy*U(13)
      U(28)=PAz*U(25)+WPz*U(13)
      U(50)=PAx*U(49)+WPx*U(14)
      U(51)=PAy*U(49)+WPy*U(14)+
     &  r1x2ZE*U(11)
      U(52)=PAz*U(49)+WPz*U(14)
      U(74)=PAx*U(73)+WPx*U(15)
      U(75)=PAy*U(73)+WPy*U(15)
      U(76)=PAz*U(73)+WPz*U(15)+
     &  r1x2ZE*U(11)
      U(35)=PAx*U(13)+WPx*U(16)+
     &  r1x2ZE*U(12)
      U(38)=PAy*U(13)+WPy*U(16)
      U(43)=PAz*U(13)+WPz*U(16)
      U(45)=PAx*U(14)+WPx*U(17)
      U(46)=PAy*U(14)+WPy*U(17)+
     &  r1x2ZE*U(12)
      U(47)=PAz*U(14)+WPz*U(17)
      U(48)=PAx*U(15)+WPx*U(18)
      U(61)=PAy*U(15)+WPy*U(18)
      U(64)=PAz*U(15)+WPz*U(18)+
     &  r1x2ZE*U(12)
      U(36)=PAx*U(16)+WPx*U(6)+
     &  r1x2ZE*U(5)
      U(37)=PAy*U(16)+WPy*U(6)
      U(39)=PAz*U(16)+WPz*U(6)
      U(40)=PAx*U(17)+WPx*U(7)
      U(41)=PAy*U(17)+WPy*U(7)+
     &  r1x2ZE*U(5)
      U(42)=PAz*U(17)+WPz*U(7)
      U(44)=PAx*U(18)+WPx*U(8)
      U(59)=PAy*U(18)+WPy*U(8)
      U(60)=PAz*U(18)+WPz*U(8)+
     &  r1x2ZE*U(5)
      U(5)=PAx*U(2)+WPx*U(19)+
     &  r1x2Z*(U(1)-rpxZ*U(11))
      U(6)=PAx*U(3)+WPx*U(20)
      U(7)=PAy*U(3)+WPy*U(20)+
     &  r1x2Z*(U(1)-rpxZ*U(11))
      U(8)=PAx*U(4)+WPx*U(21)
      U(9)=PAy*U(4)+WPy*U(21)
      U(10)=PAz*U(4)+WPz*U(21)+
     &  r1x2Z*(U(1)-rpxZ*U(11))
      U(63)=PAx*U(19)+WPx*U(22)+
     &  r1x2Z*(U(11)-rpxZ*U(12))
      U(66)=PAx*U(20)+WPx*U(23)
      U(69)=PAy*U(20)+WPy*U(23)+
     &  r1x2Z*(U(11)-rpxZ*U(12))
      U(68)=PAx*U(21)+WPx*U(24)
      U(70)=PAy*U(21)+WPy*U(24)
      U(71)=PAz*U(21)+WPz*U(24)+
     &  r1x2Z*(U(11)-rpxZ*U(12))
      U(29)=PAx*U(26)+WPx*U(35)+
     &  r1x2Z*(U(25)-rpxZ*U(13))+
     &  r1x2ZE*U(19)
      U(30)=PAx*U(27)+WPx*U(38)+
     &  r1x2ZE*U(20)
      U(31)=PAy*U(27)+WPy*U(38)+
     &  r1x2Z*(U(25)-rpxZ*U(13))
      U(32)=PAx*U(28)+WPx*U(43)+
     &  r1x2ZE*U(21)
      U(33)=PAy*U(28)+WPy*U(43)
      U(34)=PAz*U(28)+WPz*U(43)+
     &  r1x2Z*(U(25)-rpxZ*U(13))
      U(53)=PAx*U(50)+WPx*U(45)+
     &  r1x2Z*(U(49)-rpxZ*U(14))
      U(54)=PAx*U(51)+WPx*U(46)
      U(55)=PAy*U(51)+WPy*U(46)+
     &  r1x2Z*(U(49)-rpxZ*U(14))+
     &  r1x2ZE*U(20)
      U(56)=PAx*U(52)+WPx*U(47)
      U(57)=PAy*U(52)+WPy*U(47)+
     &  r1x2ZE*U(21)
      U(58)=PAz*U(52)+WPz*U(47)+
     &  r1x2Z*(U(49)-rpxZ*U(14))
      U(77)=PAx*U(74)+WPx*U(48)+
     &  r1x2Z*(U(73)-rpxZ*U(15))
      U(78)=PAx*U(75)+WPx*U(61)
      U(79)=PAy*U(75)+WPy*U(61)+
     &  r1x2Z*(U(73)-rpxZ*U(15))
      U(80)=PAx*U(76)+WPx*U(64)
      U(81)=PAy*U(76)+WPy*U(64)
      U(82)=PAz*U(76)+WPz*U(64)+
     &  r1x2Z*(U(73)-rpxZ*U(15))+
     &  r1x2ZE*U(21)
      U(62)=PAx*U(35)+WPx*U(36)+
     &  r1x2Z*(U(13)-rpxZ*U(16))+
     &  r1x2ZE*U(22)
      U(65)=PAx*U(38)+WPx*U(37)+
     &  r1x2ZE*U(23)
      U(72)=PAy*U(38)+WPy*U(37)+
     &  r1x2Z*(U(13)-rpxZ*U(16))
      U(67)=PAx*U(43)+WPx*U(39)+
     &  r1x2ZE*U(24)
      U(83)=PAy*U(43)+WPy*U(39)
      U(85)=PAz*U(43)+WPz*U(39)+
     &  r1x2Z*(U(13)-rpxZ*U(16))
      U(84)=PAx*U(45)+WPx*U(40)+
     &  r1x2Z*(U(14)-rpxZ*U(17))
      U(86)=PAx*U(46)+WPx*U(41)
      U(88)=PAy*U(46)+WPy*U(41)+
     &  r1x2Z*(U(14)-rpxZ*U(17))+
     &  r1x2ZE*U(23)
      U(87)=PAx*U(47)+WPx*U(42)
      U(90)=PAy*U(47)+WPy*U(42)+
     &  r1x2ZE*U(24)
      U(91)=PAz*U(47)+WPz*U(42)+
     &  r1x2Z*(U(14)-rpxZ*U(17))
      U(89)=PAx*U(48)+WPx*U(44)+
     &  r1x2Z*(U(15)-rpxZ*U(18))
      U(93)=PAx*U(61)+WPx*U(59)
      U(94)=PAy*U(61)+WPy*U(59)+
     &  r1x2Z*(U(15)-rpxZ*U(18))
      U(95)=PAx*U(64)+WPx*U(60)
      U(96)=PAy*U(64)+WPy*U(60)
      U(97)=PAz*U(64)+WPz*U(60)+
     &  r1x2Z*(U(15)-rpxZ*U(18))+
     &  r1x2ZE*U(24)
      U(35)=PAx*U(29)+WPx*U(62)+
     &  2.0D0*r1x2Z*(U(26)-rpxZ*U(35))+
     &  r1x2ZE*U(63)
      U(36)=PAx*U(30)+WPx*U(65)+
     &  r1x2Z*(U(27)-rpxZ*U(38))+
     &  r1x2ZE*U(66)
      U(37)=PAx*U(31)+WPx*U(72)+
     &  r1x2ZE*U(69)
      U(39)=PAx*U(32)+WPx*U(67)+
     &  r1x2Z*(U(28)-rpxZ*U(43))+
     &  r1x2ZE*U(68)
      U(40)=PAx*U(33)+WPx*U(83)+
     &  r1x2ZE*U(70)
      U(42)=PAx*U(34)+WPx*U(85)+
     &  r1x2ZE*U(71)
      U(62)=PAy*U(55)+WPy*U(88)+
     &  2.0D0*r1x2Z*(U(51)-rpxZ*U(46))+
     &  r1x2ZE*U(69)
      U(65)=PAy*U(57)+WPy*U(90)+
     &  r1x2Z*(U(52)-rpxZ*U(47))+
     &  r1x2ZE*U(70)
      U(67)=PAy*U(58)+WPy*U(91)+
     &  r1x2ZE*U(71)
      U(92)=PAz*U(82)+WPz*U(97)+
     &  2.0D0*r1x2Z*(U(76)-rpxZ*U(64))+
     &  r1x2ZE*U(71)
      U(11)=PAx*U(5)+WPx*U(63)+
     &  2.0D0*r1x2Z*(U(2)-rpxZ*U(19))
      U(12)=PAx*U(6)+WPx*U(66)+
     &  r1x2Z*(U(3)-rpxZ*U(20))
      U(14)=PAy*U(7)+WPy*U(69)+
     &  2.0D0*r1x2Z*(U(3)-rpxZ*U(20))
      U(15)=PAx*U(8)+WPx*U(68)+
     &  r1x2Z*(U(4)-rpxZ*U(21))
      U(17)=PAy*U(9)+WPy*U(70)+
     &  r1x2Z*(U(4)-rpxZ*U(21))
      U(20)=PAz*U(10)+WPz*U(71)+
     &  2.0D0*r1x2Z*(U(4)-rpxZ*U(21))
      U(38)=PAy*U(31)+WPy*U(72)+
     &  2.0D0*r1x2Z*(U(27)-rpxZ*U(38))
      U(41)=PAy*U(33)+WPy*U(83)+
     &  r1x2Z*(U(28)-rpxZ*U(43))
      U(44)=PAz*U(34)+WPz*U(85)+
     &  2.0D0*r1x2Z*(U(28)-rpxZ*U(43))
      U(59)=PAx*U(53)+WPx*U(84)+
     &  2.0D0*r1x2Z*(U(50)-rpxZ*U(45))
      U(60)=PAx*U(54)+WPx*U(86)+
     &  r1x2Z*(U(51)-rpxZ*U(46))
      U(63)=PAx*U(56)+WPx*U(87)+
     &  r1x2Z*(U(52)-rpxZ*U(47))
      U(68)=PAz*U(58)+WPz*U(91)+
     &  2.0D0*r1x2Z*(U(52)-rpxZ*U(47))
      U(83)=PAx*U(77)+WPx*U(89)+
     &  2.0D0*r1x2Z*(U(74)-rpxZ*U(48))
      U(84)=PAx*U(78)+WPx*U(93)+
     &  r1x2Z*(U(75)-rpxZ*U(61))
      U(86)=PAy*U(79)+WPy*U(94)+
     &  2.0D0*r1x2Z*(U(75)-rpxZ*U(61))
      U(87)=PAx*U(80)+WPx*U(95)+
     &  r1x2Z*(U(76)-rpxZ*U(64))
      U(89)=PAy*U(81)+WPy*U(96)+
     &  r1x2Z*(U(76)-rpxZ*U(64))
      U(13)=PAx*U(7)+WPx*U(69)
      U(16)=PAx*U(9)+WPx*U(70)
      U(18)=PAx*U(10)+WPx*U(71)
      U(19)=PAy*U(10)+WPy*U(71)
      U(43)=PAy*U(34)+WPy*U(85)
      U(61)=PAx*U(55)+WPx*U(88)
      U(64)=PAx*U(57)+WPx*U(90)
      U(66)=PAx*U(58)+WPx*U(91)
      U(85)=PAx*U(79)+WPx*U(94)
      U(88)=PAx*U(81)+WPx*U(96)
      U(90)=PAx*U(82)+WPx*U(97)
      U(91)=PAy*U(82)+WPy*U(97)

      C(I0,1)=C(I0,1)+U(5)*cp1*cq1
      C(I0,2)=C(I0,2)+U(5)*cq1
      C(I0,3)=C(I0,3)+U(6)*cq1
      C(I0,4)=C(I0,4)+U(7)*cq1
      C(I0,5)=C(I0,5)+U(6)*cp1*cq1
      C(I0,6)=C(I0,6)+U(8)*cq1
      C(I0,7)=C(I0,7)+U(9)*cq1
      C(I0,8)=C(I0,8)+U(10)*cq1
      C(I0,9)=C(I0,9)+U(7)*cp1*cq1
      C(I0,10)=C(I0,10)+U(11)*cq1
      C(I0,11)=C(I0,11)+U(12)*cq1
      C(I0,12)=C(I0,12)+U(13)*cq1
      C(I0,13)=C(I0,13)+U(8)*cp1*cq1
      C(I0,14)=C(I0,14)+U(14)*cq1
      C(I0,15)=C(I0,15)+U(15)*cq1
      C(I0,16)=C(I0,16)+U(16)*cq1
      C(I0,17)=C(I0,17)+U(9)*cp1*cq1
      C(I0,18)=C(I0,18)+U(17)*cq1
      C(I0,19)=C(I0,19)+U(18)*cq1
      C(I0,20)=C(I0,20)+U(19)*cq1
      C(I0,21)=C(I0,21)+U(10)*cp1*cq1
      C(I0,22)=C(I0,22)+U(20)*cq1
      C(I0,25)=C(I0,25)+U(29)*cp1
      C(I0,26)=C(I0,26)+U(29)
      C(I0,27)=C(I0,27)+U(30)
      C(I0,28)=C(I0,28)+U(31)
      C(I0,29)=C(I0,29)+U(30)*cp1
      C(I0,30)=C(I0,30)+U(32)
      C(I0,31)=C(I0,31)+U(33)
      C(I0,32)=C(I0,32)+U(34)
      C(I0,33)=C(I0,33)+U(31)*cp1
      C(I0,34)=C(I0,34)+U(35)
      C(I0,35)=C(I0,35)+U(36)
      C(I0,36)=C(I0,36)+U(37)
      C(I0,37)=C(I0,37)+U(32)*cp1
      C(I0,38)=C(I0,38)+U(38)
      C(I0,39)=C(I0,39)+U(39)
      C(I0,40)=C(I0,40)+U(40)
      C(I0,41)=C(I0,41)+U(33)*cp1
      C(I0,42)=C(I0,42)+U(41)
      C(I0,43)=C(I0,43)+U(42)
      C(I0,44)=C(I0,44)+U(43)
      C(I0,45)=C(I0,45)+U(34)*cp1
      C(I0,46)=C(I0,46)+U(44)
      C(I0,49)=C(I0,49)+U(53)*cp1
      C(I0,50)=C(I0,50)+U(53)
      C(I0,51)=C(I0,51)+U(54)
      C(I0,52)=C(I0,52)+U(55)
      C(I0,53)=C(I0,53)+U(54)*cp1
      C(I0,54)=C(I0,54)+U(56)
      C(I0,55)=C(I0,55)+U(57)
      C(I0,56)=C(I0,56)+U(58)
      C(I0,57)=C(I0,57)+U(55)*cp1
      C(I0,58)=C(I0,58)+U(59)
      C(I0,59)=C(I0,59)+U(60)
      C(I0,60)=C(I0,60)+U(61)
      C(I0,61)=C(I0,61)+U(56)*cp1
      C(I0,62)=C(I0,62)+U(62)
      C(I0,63)=C(I0,63)+U(63)
      C(I0,64)=C(I0,64)+U(64)
      C(I0,65)=C(I0,65)+U(57)*cp1
      C(I0,66)=C(I0,66)+U(65)
      C(I0,67)=C(I0,67)+U(66)
      C(I0,68)=C(I0,68)+U(67)
      C(I0,69)=C(I0,69)+U(58)*cp1
      C(I0,70)=C(I0,70)+U(68)
      C(I0,73)=C(I0,73)+U(77)*cp1
      C(I0,74)=C(I0,74)+U(77)
      C(I0,75)=C(I0,75)+U(78)
      C(I0,76)=C(I0,76)+U(79)
      C(I0,77)=C(I0,77)+U(78)*cp1
      C(I0,78)=C(I0,78)+U(80)
      C(I0,79)=C(I0,79)+U(81)
      C(I0,80)=C(I0,80)+U(82)
      C(I0,81)=C(I0,81)+U(79)*cp1
      C(I0,82)=C(I0,82)+U(83)
      C(I0,83)=C(I0,83)+U(84)
      C(I0,84)=C(I0,84)+U(85)
      C(I0,85)=C(I0,85)+U(80)*cp1
      C(I0,86)=C(I0,86)+U(86)
      C(I0,87)=C(I0,87)+U(87)
      C(I0,88)=C(I0,88)+U(88)
      C(I0,89)=C(I0,89)+U(81)*cp1
      C(I0,90)=C(I0,90)+U(89)
      C(I0,91)=C(I0,91)+U(90)
      C(I0,92)=C(I0,92)+U(91)
      C(I0,93)=C(I0,93)+U(82)*cp1
      C(I0,94)=C(I0,94)+U(92)


 1000             CONTINUE
 1010             CONTINUE
 2000          CONTINUE
            ENDIF
 3000    CONTINUE
 4000 CONTINUE 
      RETURN
      END
