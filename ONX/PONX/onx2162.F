      SUBROUTINE ONX2162(N,iBgn,iEnd,IndA,BfnBufB,BfnBufK,KBra,KKet,
     &                   SwhB,SwhK,RAC,RBD,
     &                   PrmBufB,PrmBufK,U,C)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 U(80),C(N,88)
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

      CALL RLoader(C,N*88,0.0D0)

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
      U(4)=R2
      U(2)=R3
      U(3)=R4
      U(6)=R5
      U(5)=QBx*U(1)+WQx*U(4)
      U(9)=QBy*U(1)+WQy*U(4)
      U(13)=QBz*U(1)+WQz*U(4)
      U(8)=QBx*U(4)+WQx*U(2)
      U(12)=QBy*U(4)+WQy*U(2)
      U(15)=QBz*U(4)+WQz*U(2)
      U(7)=QBx*U(2)+WQx*U(3)
      U(10)=QBy*U(2)+WQy*U(3)
      U(11)=QBz*U(2)+WQz*U(3)
      U(14)=QBx*U(3)+WQx*U(6)
      U(16)=QBy*U(3)+WQy*U(6)
      U(18)=QBz*U(3)+WQz*U(6)
      U(17)=QBx*U(5)+WQx*U(8)+
     &  r1x2E*(U(1)-rpxE*U(4))
      U(21)=QBx*U(9)+WQx*U(12)
      U(25)=QBy*U(9)+WQy*U(12)+
     &  r1x2E*(U(1)-rpxE*U(4))
      U(29)=QBx*U(13)+WQx*U(15)
      U(33)=QBy*U(13)+WQy*U(15)
      U(37)=QBz*U(13)+WQz*U(15)+
     &  r1x2E*(U(1)-rpxE*U(4))
      U(20)=QBx*U(8)+WQx*U(7)+
     &  r1x2E*(U(4)-rpxE*U(2))
      U(24)=QBx*U(12)+WQx*U(10)
      U(28)=QBy*U(12)+WQy*U(10)+
     &  r1x2E*(U(4)-rpxE*U(2))
      U(31)=QBx*U(15)+WQx*U(11)
      U(34)=QBy*U(15)+WQy*U(11)
      U(39)=QBz*U(15)+WQz*U(11)+
     &  r1x2E*(U(4)-rpxE*U(2))
      U(19)=QBx*U(7)+WQx*U(14)+
     &  r1x2E*(U(2)-rpxE*U(3))
      U(22)=QBx*U(10)+WQx*U(16)
      U(23)=QBy*U(10)+WQy*U(16)+
     &  r1x2E*(U(2)-rpxE*U(3))
      U(26)=QBx*U(11)+WQx*U(18)
      U(27)=QBy*U(11)+WQy*U(18)
      U(30)=QBz*U(11)+WQz*U(18)+
     &  r1x2E*(U(2)-rpxE*U(3))
      U(41)=QBx*U(17)+WQx*U(20)+
     &  2.0D0*r1x2E*(U(5)-rpxE*U(8))
      U(45)=QBx*U(21)+WQx*U(24)+
     &  r1x2E*(U(9)-rpxE*U(12))
      U(49)=QBx*U(25)+WQx*U(28)
      U(53)=QBy*U(25)+WQy*U(28)+
     &  2.0D0*r1x2E*(U(9)-rpxE*U(12))
      U(57)=QBx*U(29)+WQx*U(31)+
     &  r1x2E*(U(13)-rpxE*U(15))
      U(61)=QBx*U(33)+WQx*U(34)
      U(65)=QBy*U(33)+WQy*U(34)+
     &  r1x2E*(U(13)-rpxE*U(15))
      U(69)=QBx*U(37)+WQx*U(39)
      U(73)=QBy*U(37)+WQy*U(39)
      U(77)=QBz*U(37)+WQz*U(39)+
     &  2.0D0*r1x2E*(U(13)-rpxE*U(15))
      U(44)=QBx*U(20)+WQx*U(19)+
     &  2.0D0*r1x2E*(U(8)-rpxE*U(7))
      U(48)=QBx*U(24)+WQx*U(22)+
     &  r1x2E*(U(12)-rpxE*U(10))
      U(52)=QBx*U(28)+WQx*U(23)
      U(56)=QBy*U(28)+WQy*U(23)+
     &  2.0D0*r1x2E*(U(12)-rpxE*U(10))
      U(59)=QBx*U(31)+WQx*U(26)+
     &  r1x2E*(U(15)-rpxE*U(11))
      U(38)=QBx*U(34)+WQx*U(27)
      U(66)=QBy*U(34)+WQy*U(27)+
     &  r1x2E*(U(15)-rpxE*U(11))
      U(71)=QBx*U(39)+WQx*U(30)
      U(74)=QBy*U(39)+WQy*U(30)
      U(79)=QBz*U(39)+WQz*U(30)+
     &  2.0D0*r1x2E*(U(15)-rpxE*U(11))
      U(6)=PAx*U(5)+WPx*U(8)+
     &  r1x2ZE*U(4)
      U(11)=PAy*U(9)+WPy*U(12)+
     &  r1x2ZE*U(4)
      U(16)=PAz*U(13)+WPz*U(15)+
     &  r1x2ZE*U(4)
      U(18)=PAx*U(17)+WPx*U(20)+
     &  2.0D0*r1x2ZE*U(8)
      U(22)=PAx*U(21)+WPx*U(24)+
     &  r1x2ZE*U(12)
      U(23)=PAy*U(21)+WPy*U(24)+
     &  r1x2ZE*U(8)
      U(27)=PAy*U(25)+WPy*U(28)+
     &  2.0D0*r1x2ZE*U(12)
      U(30)=PAx*U(29)+WPx*U(31)+
     &  r1x2ZE*U(15)
      U(32)=PAz*U(29)+WPz*U(31)+
     &  r1x2ZE*U(8)
      U(35)=PAy*U(33)+WPy*U(34)+
     &  r1x2ZE*U(15)
      U(36)=PAz*U(33)+WPz*U(34)+
     &  r1x2ZE*U(12)
      U(40)=PAz*U(37)+WPz*U(39)+
     &  2.0D0*r1x2ZE*U(15)
      U(42)=PAx*U(41)+WPx*U(44)+
     &  3.0D0*r1x2ZE*U(20)
      U(46)=PAx*U(45)+WPx*U(48)+
     &  2.0D0*r1x2ZE*U(24)
      U(47)=PAy*U(45)+WPy*U(48)+
     &  r1x2ZE*U(20)
      U(50)=PAx*U(49)+WPx*U(52)+
     &  r1x2ZE*U(28)
      U(51)=PAy*U(49)+WPy*U(52)+
     &  2.0D0*r1x2ZE*U(24)
      U(55)=PAy*U(53)+WPy*U(56)+
     &  3.0D0*r1x2ZE*U(28)
      U(58)=PAx*U(57)+WPx*U(59)+
     &  2.0D0*r1x2ZE*U(31)
      U(60)=PAz*U(57)+WPz*U(59)+
     &  r1x2ZE*U(20)
      U(62)=PAx*U(61)+WPx*U(38)+
     &  r1x2ZE*U(34)
      U(63)=PAy*U(61)+WPy*U(38)+
     &  r1x2ZE*U(31)
      U(64)=PAz*U(61)+WPz*U(38)+
     &  r1x2ZE*U(24)
      U(67)=PAy*U(65)+WPy*U(66)+
     &  2.0D0*r1x2ZE*U(34)
      U(68)=PAz*U(65)+WPz*U(66)+
     &  r1x2ZE*U(28)
      U(70)=PAx*U(69)+WPx*U(71)+
     &  r1x2ZE*U(39)
      U(72)=PAz*U(69)+WPz*U(71)+
     &  2.0D0*r1x2ZE*U(31)
      U(75)=PAy*U(73)+WPy*U(74)+
     &  r1x2ZE*U(39)
      U(76)=PAz*U(73)+WPz*U(74)+
     &  2.0D0*r1x2ZE*U(34)
      U(80)=PAz*U(77)+WPz*U(79)+
     &  3.0D0*r1x2ZE*U(39)
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
      U(43)=PAy*U(41)+WPy*U(44)
      U(44)=PAz*U(41)+WPz*U(44)
      U(48)=PAz*U(45)+WPz*U(48)
      U(52)=PAz*U(49)+WPz*U(52)
      U(54)=PAx*U(53)+WPx*U(56)
      U(56)=PAz*U(53)+WPz*U(56)
      U(59)=PAy*U(57)+WPy*U(59)
      U(66)=PAx*U(65)+WPx*U(66)
      U(71)=PAy*U(69)+WPy*U(71)
      U(74)=PAx*U(73)+WPx*U(74)
      U(78)=PAx*U(77)+WPx*U(79)
      U(79)=PAy*U(77)+WPy*U(79)

      C(I0,1)=C(I0,1)+U(17)*cp1*cq1
      C(I0,2)=C(I0,2)+U(18)*cq1
      C(I0,3)=C(I0,3)+U(19)*cq1
      C(I0,4)=C(I0,4)+U(20)*cq1
      C(I0,5)=C(I0,5)+U(17)*cp1
      C(I0,6)=C(I0,6)+U(18)
      C(I0,7)=C(I0,7)+U(19)
      C(I0,8)=C(I0,8)+U(20)
      C(I0,9)=C(I0,9)+U(21)*cp1
      C(I0,10)=C(I0,10)+U(22)
      C(I0,11)=C(I0,11)+U(23)
      C(I0,12)=C(I0,12)+U(24)
      C(I0,13)=C(I0,13)+U(25)*cp1
      C(I0,14)=C(I0,14)+U(26)
      C(I0,15)=C(I0,15)+U(27)
      C(I0,16)=C(I0,16)+U(28)
      C(I0,17)=C(I0,17)+U(21)*cp1*cq1
      C(I0,18)=C(I0,18)+U(22)*cq1
      C(I0,19)=C(I0,19)+U(23)*cq1
      C(I0,20)=C(I0,20)+U(24)*cq1
      C(I0,21)=C(I0,21)+U(29)*cp1
      C(I0,22)=C(I0,22)+U(30)
      C(I0,23)=C(I0,23)+U(31)
      C(I0,24)=C(I0,24)+U(32)
      C(I0,25)=C(I0,25)+U(33)*cp1
      C(I0,26)=C(I0,26)+U(34)
      C(I0,27)=C(I0,27)+U(35)
      C(I0,28)=C(I0,28)+U(36)
      C(I0,29)=C(I0,29)+U(37)*cp1
      C(I0,30)=C(I0,30)+U(38)
      C(I0,31)=C(I0,31)+U(39)
      C(I0,32)=C(I0,32)+U(40)
      C(I0,33)=C(I0,33)+U(25)*cp1*cq1
      C(I0,34)=C(I0,34)+U(26)*cq1
      C(I0,35)=C(I0,35)+U(27)*cq1
      C(I0,36)=C(I0,36)+U(28)*cq1
      C(I0,37)=C(I0,37)+U(41)*cp1
      C(I0,38)=C(I0,38)+U(42)
      C(I0,39)=C(I0,39)+U(43)
      C(I0,40)=C(I0,40)+U(44)
      C(I0,41)=C(I0,41)+U(45)*cp1
      C(I0,42)=C(I0,42)+U(46)
      C(I0,43)=C(I0,43)+U(47)
      C(I0,44)=C(I0,44)+U(48)
      C(I0,45)=C(I0,45)+U(49)*cp1
      C(I0,46)=C(I0,46)+U(50)
      C(I0,47)=C(I0,47)+U(51)
      C(I0,48)=C(I0,48)+U(52)
      C(I0,49)=C(I0,49)+U(29)*cp1*cq1
      C(I0,50)=C(I0,50)+U(30)*cq1
      C(I0,51)=C(I0,51)+U(31)*cq1
      C(I0,52)=C(I0,52)+U(32)*cq1
      C(I0,53)=C(I0,53)+U(53)*cp1
      C(I0,54)=C(I0,54)+U(54)
      C(I0,55)=C(I0,55)+U(55)
      C(I0,56)=C(I0,56)+U(56)
      C(I0,57)=C(I0,57)+U(57)*cp1
      C(I0,58)=C(I0,58)+U(58)
      C(I0,59)=C(I0,59)+U(59)
      C(I0,60)=C(I0,60)+U(60)
      C(I0,61)=C(I0,61)+U(61)*cp1
      C(I0,62)=C(I0,62)+U(62)
      C(I0,63)=C(I0,63)+U(63)
      C(I0,64)=C(I0,64)+U(64)
      C(I0,65)=C(I0,65)+U(33)*cp1*cq1
      C(I0,66)=C(I0,66)+U(34)*cq1
      C(I0,67)=C(I0,67)+U(35)*cq1
      C(I0,68)=C(I0,68)+U(36)*cq1
      C(I0,69)=C(I0,69)+U(65)*cp1
      C(I0,70)=C(I0,70)+U(66)
      C(I0,71)=C(I0,71)+U(67)
      C(I0,72)=C(I0,72)+U(68)
      C(I0,73)=C(I0,73)+U(69)*cp1
      C(I0,74)=C(I0,74)+U(70)
      C(I0,75)=C(I0,75)+U(71)
      C(I0,76)=C(I0,76)+U(72)
      C(I0,77)=C(I0,77)+U(73)*cp1
      C(I0,78)=C(I0,78)+U(74)
      C(I0,79)=C(I0,79)+U(75)
      C(I0,80)=C(I0,80)+U(76)
      C(I0,81)=C(I0,81)+U(37)*cp1*cq1
      C(I0,82)=C(I0,82)+U(38)*cq1
      C(I0,83)=C(I0,83)+U(39)*cq1
      C(I0,84)=C(I0,84)+U(40)*cq1
      C(I0,85)=C(I0,85)+U(77)*cp1
      C(I0,86)=C(I0,86)+U(78)
      C(I0,87)=C(I0,87)+U(79)
      C(I0,88)=C(I0,88)+U(80)

 1000             CONTINUE
 1010             CONTINUE
 2000          CONTINUE
            ENDIF
 3000    CONTINUE
 4000 CONTINUE 
      RETURN
      END
