      SUBROUTINE ONX2222(N,iBgn,iEnd,IndA,BfnBufB,BfnBufK,KBra,KKet,
     &                   SwhB,SwhK,RAC,RBD,
     &                   PrmBufB,PrmBufK,U,C)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 U(163),C(N,289)
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

      CALL RLoader(C,N*289,0.0D0)
 
      Tsh0=Tsh/Dcd
      I0=0
      DO 4000 I=iBgn,iEnd
         Is=I-iBgn+1
         Ixa = BfnBufB(I)
         Tsh1=Tsh0/PrmBufB(9,1,I)
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
                  schK = PrmBufK(9,IKK,J)
                  IF(schK.LE.Tsh1) GOTO 2010
                  Tsh2 = Tsh0/schK
                  Eta  = PrmBufK(1,IKK,J)
                  Qx   = PrmBufK(2,IKK,J)
                  Qy   = PrmBufK(3,IKK,J)
                  Qz   = PrmBufK(4,IKK,J)
                  Uq   = PrmBufK(5,IKK,J)
                  cq1  = PrmBufK(6,IKK,J)
                  cq2  = PrmBufK(7,IKK,J)
                  cq3  = PrmBufK(8,IKK,J)
                  r1xe=1.0D0/Eta
                  r1x2e=0.5D0*r1xe
                  EtaQx=Eta*Qx
                  EtaQy=Eta*Qy
                  EtaQz=Eta*Qz
                  QBx = Qx-rDx
                  QBy = Qy-rDy
                  QBz = Qz-rDz

                  DO 1000 IKB=1,KBra
                     IF(PrmBufB(9,IKB,I).LE.Tsh2) GOTO 1010
                     Zeta = PrmBufB(1,IKB,I)
                     Px   = PrmBufB(2,IKB,I)
                     Py   = PrmBufB(3,IKB,I)
                     Pz   = PrmBufB(4,IKB,I)
                     Up   = PrmBufB(5,IKB,I)
                     cp1  = PrmBufB(6,IKB,I)
                     cp2  = PrmBufB(7,IKB,I)
                     cp3  = PrmBufB(8,IKB,I)
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
      U(69)=QBx*U(18)+WQx*U(8)+                                        
     &  r1x2E*(U(1)-rpxE*U(6))                                           
      U(86)=QBx*U(35)+WQx*U(9)                                         
      U(103)=QBy*U(35)+WQy*U(9)+                                       
     &  r1x2E*(U(1)-rpxE*U(6))                                           
      U(120)=QBx*U(52)+WQx*U(11)                                       
      U(137)=QBy*U(52)+WQy*U(11)                                       
      U(154)=QBz*U(52)+WQz*U(11)+                                      
     &  r1x2E*(U(1)-rpxE*U(6))                                           
      U(15)=QBx*U(8)+WQx*U(7)+                                         
     &  r1x2E*(U(6)-rpxE*U(5))                                           
      U(16)=QBx*U(9)+WQx*U(10)                                         
      U(17)=QBy*U(9)+WQy*U(10)+                                        
     &  r1x2E*(U(6)-rpxE*U(5))                                           
      U(26)=QBx*U(11)+WQx*U(12)                                        
      U(28)=QBy*U(11)+WQy*U(12)                                        
      U(29)=QBz*U(11)+WQz*U(12)+                                       
     &  r1x2E*(U(6)-rpxE*U(5))                                           
      U(22)=QBx*U(7)+WQx*U(4)+                                         
     &  r1x2E*(U(5)-rpxE*U(2))                                           
      U(23)=QBx*U(10)+WQx*U(13)                                        
      U(24)=QBy*U(10)+WQy*U(13)+                                       
     &  r1x2E*(U(5)-rpxE*U(2))                                           
      U(25)=QBx*U(12)+WQx*U(14)                                        
      U(27)=QBy*U(12)+WQy*U(14)                                        
      U(30)=QBz*U(12)+WQz*U(14)+                                       
     &  r1x2E*(U(5)-rpxE*U(2))                                           
      U(2)=PAx*U(1)+WPx*U(6)                                           
      U(3)=PAy*U(1)+WPy*U(6)                                           
      U(4)=PAz*U(1)+WPz*U(6)                                           
      U(31)=PAx*U(6)+WPx*U(5)                                          
      U(32)=PAy*U(6)+WPy*U(5)                                          
      U(33)=PAz*U(6)+WPz*U(5)                                          
      U(19)=PAx*U(18)+WPx*U(8)+                                        
     &  r1x2ZE*U(6)                                                      
      U(20)=PAy*U(18)+WPy*U(8)                                         
      U(21)=PAz*U(18)+WPz*U(8)                                         
      U(36)=PAx*U(35)+WPx*U(9)                                         
      U(37)=PAy*U(35)+WPy*U(9)+                                        
     &  r1x2ZE*U(6)                                                      
      U(38)=PAz*U(35)+WPz*U(9)                                         
      U(53)=PAx*U(52)+WPx*U(11)                                        
      U(54)=PAy*U(52)+WPy*U(11)                                        
      U(55)=PAz*U(52)+WPz*U(11)+                                       
     &  r1x2ZE*U(6)                                                      
      U(34)=PAx*U(8)+WPx*U(7)+                                         
     &  r1x2ZE*U(5)                                                      
      U(39)=PAy*U(8)+WPy*U(7)                                          
      U(40)=PAz*U(8)+WPz*U(7)                                          
      U(42)=PAx*U(9)+WPx*U(10)                                         
      U(45)=PAy*U(9)+WPy*U(10)+                                        
     &  r1x2ZE*U(5)                                                      
      U(46)=PAz*U(9)+WPz*U(10)                                         
      U(47)=PAx*U(11)+WPx*U(12)                                        
      U(48)=PAy*U(11)+WPy*U(12)                                        
      U(49)=PAz*U(11)+WPz*U(12)+                                       
     &  r1x2ZE*U(5)                                                      
      U(70)=PAx*U(69)+WPx*U(15)+                                       
     &  2.0D0*r1x2ZE*U(8)                                                
      U(71)=PAy*U(69)+WPy*U(15)                                        
      U(72)=PAz*U(69)+WPz*U(15)                                        
      U(87)=PAx*U(86)+WPx*U(16)+                                       
     &  r1x2ZE*U(9)                                                      
      U(88)=PAy*U(86)+WPy*U(16)+                                       
     &  r1x2ZE*U(8)                                                      
      U(89)=PAz*U(86)+WPz*U(16)                                        
      U(104)=PAx*U(103)+WPx*U(17)                                      
      U(105)=PAy*U(103)+WPy*U(17)+                                     
     &  2.0D0*r1x2ZE*U(9)                                                
      U(106)=PAz*U(103)+WPz*U(17)                                      
      U(121)=PAx*U(120)+WPx*U(26)+                                     
     &  r1x2ZE*U(11)                                                     
      U(122)=PAy*U(120)+WPy*U(26)                                      
      U(123)=PAz*U(120)+WPz*U(26)+                                     
     &  r1x2ZE*U(8)                                                      
      U(138)=PAx*U(137)+WPx*U(28)                                      
      U(139)=PAy*U(137)+WPy*U(28)+                                     
     &  r1x2ZE*U(11)                                                     
      U(140)=PAz*U(137)+WPz*U(28)+                                     
     &  r1x2ZE*U(9)                                                      
      U(155)=PAx*U(154)+WPx*U(29)                                      
      U(156)=PAy*U(154)+WPy*U(29)                                      
      U(157)=PAz*U(154)+WPz*U(29)+                                     
     &  2.0D0*r1x2ZE*U(11)                                               
      U(44)=PAx*U(15)+WPx*U(22)+                                       
     &  2.0D0*r1x2ZE*U(7)                                                
      U(50)=PAy*U(15)+WPy*U(22)                                        
      U(51)=PAz*U(15)+WPz*U(22)                                        
      U(56)=PAx*U(16)+WPx*U(23)+                                       
     &  r1x2ZE*U(10)                                                     
      U(57)=PAy*U(16)+WPy*U(23)+                                       
     &  r1x2ZE*U(7)                                                      
      U(59)=PAz*U(16)+WPz*U(23)                                        
      U(60)=PAx*U(17)+WPx*U(24)                                        
      U(62)=PAy*U(17)+WPy*U(24)+                                       
     &  2.0D0*r1x2ZE*U(10)                                               
      U(63)=PAz*U(17)+WPz*U(24)                                        
      U(58)=PAx*U(26)+WPx*U(25)+                                       
     &  r1x2ZE*U(12)                                                     
      U(64)=PAy*U(26)+WPy*U(25)                                        
      U(65)=PAz*U(26)+WPz*U(25)+                                       
     &  r1x2ZE*U(7)                                                      
      U(66)=PAx*U(28)+WPx*U(27)                                        
      U(67)=PAy*U(28)+WPy*U(27)+                                       
     &  r1x2ZE*U(12)                                                     
      U(68)=PAz*U(28)+WPz*U(27)+                                       
     &  r1x2ZE*U(10)                                                     
      U(77)=PAx*U(29)+WPx*U(30)                                        
      U(79)=PAy*U(29)+WPy*U(30)                                        
      U(80)=PAz*U(29)+WPz*U(30)+                                       
     &  2.0D0*r1x2ZE*U(12)                                               
      U(22)=PAx*U(19)+WPx*U(34)+                                       
     &  r1x2Z*(U(18)-rpxZ*U(8))+                                         
     &  r1x2ZE*U(31)                                                     
      U(23)=PAx*U(20)+WPx*U(39)+                                       
     &  r1x2ZE*U(32)                                                     
      U(25)=PAx*U(21)+WPx*U(40)+                                       
     &  r1x2ZE*U(33)                                                     
      U(41)=PAy*U(37)+WPy*U(45)+                                       
     &  r1x2Z*(U(35)-rpxZ*U(9))+                                         
     &  r1x2ZE*U(32)                                                     
      U(43)=PAy*U(38)+WPy*U(46)+                                       
     &  r1x2ZE*U(33)                                                     
      U(61)=PAz*U(55)+WPz*U(49)+                                       
     &  r1x2Z*(U(52)-rpxZ*U(11))+                                        
     &  r1x2ZE*U(33)                                                     
      U(73)=PAx*U(70)+WPx*U(44)+                                       
     &  r1x2Z*(U(69)-rpxZ*U(15))+                                        
     &  2.0D0*r1x2ZE*U(34)                                               
      U(74)=PAx*U(71)+WPx*U(50)+                                       
     &  2.0D0*r1x2ZE*U(39)                                               
      U(76)=PAx*U(72)+WPx*U(51)+                                       
     &  2.0D0*r1x2ZE*U(40)                                               
      U(90)=PAx*U(87)+WPx*U(56)+                                       
     &  r1x2Z*(U(86)-rpxZ*U(16))+                                        
     &  r1x2ZE*U(42)                                                     
      U(91)=PAx*U(88)+WPx*U(57)+                                       
     &  r1x2ZE*U(45)                                                     
      U(92)=PAy*U(88)+WPy*U(57)+                                       
     &  r1x2Z*(U(86)-rpxZ*U(16))+                                        
     &  r1x2ZE*U(39)                                                     
      U(93)=PAx*U(89)+WPx*U(59)+                                       
     &  r1x2ZE*U(46)                                                     
      U(94)=PAy*U(89)+WPy*U(59)+                                       
     &  r1x2ZE*U(40)                                                     
      U(109)=PAy*U(105)+WPy*U(62)+                                     
     &  r1x2Z*(U(103)-rpxZ*U(17))+                                       
     &  2.0D0*r1x2ZE*U(45)                                               
      U(111)=PAy*U(106)+WPy*U(63)+                                     
     &  2.0D0*r1x2ZE*U(46)                                               
      U(124)=PAx*U(121)+WPx*U(58)+                                     
     &  r1x2Z*(U(120)-rpxZ*U(26))+                                       
     &  r1x2ZE*U(47)                                                     
      U(125)=PAx*U(122)+WPx*U(64)+                                     
     &  r1x2ZE*U(48)                                                     
      U(127)=PAx*U(123)+WPx*U(65)+                                     
     &  r1x2ZE*U(49)                                                     
      U(129)=PAz*U(123)+WPz*U(65)+                                     
     &  r1x2Z*(U(120)-rpxZ*U(26))+                                       
     &  r1x2ZE*U(40)                                                     
      U(143)=PAy*U(139)+WPy*U(67)+                                     
     &  r1x2Z*(U(137)-rpxZ*U(28))+                                       
     &  r1x2ZE*U(48)                                                     
      U(145)=PAy*U(140)+WPy*U(68)+                                     
     &  r1x2ZE*U(49)                                                     
      U(146)=PAz*U(140)+WPz*U(68)+                                     
     &  r1x2Z*(U(137)-rpxZ*U(28))+                                       
     &  r1x2ZE*U(46)                                                     
      U(163)=PAz*U(157)+WPz*U(80)+                                     
     &  r1x2Z*(U(154)-rpxZ*U(29))+                                       
     &  2.0D0*r1x2ZE*U(49)                                               
      U(5)=PAx*U(2)+WPx*U(31)+                                         
     &  r1x2Z*(U(1)-rpxZ*U(6))                                           
      U(7)=PAy*U(3)+WPy*U(32)+                                         
     &  r1x2Z*(U(1)-rpxZ*U(6))                                           
      U(10)=PAz*U(4)+WPz*U(33)+                                        
     &  r1x2Z*(U(1)-rpxZ*U(6))                                           
      U(24)=PAy*U(20)+WPy*U(39)+                                       
     &  r1x2Z*(U(18)-rpxZ*U(8))                                          
      U(27)=PAz*U(21)+WPz*U(40)+                                       
     &  r1x2Z*(U(18)-rpxZ*U(8))                                          
      U(39)=PAx*U(36)+WPx*U(42)+                                       
     &  r1x2Z*(U(35)-rpxZ*U(9))                                          
      U(44)=PAz*U(38)+WPz*U(46)+                                       
     &  r1x2Z*(U(35)-rpxZ*U(9))                                          
      U(56)=PAx*U(53)+WPx*U(47)+                                       
     &  r1x2Z*(U(52)-rpxZ*U(11))                                         
      U(58)=PAy*U(54)+WPy*U(48)+                                       
     &  r1x2Z*(U(52)-rpxZ*U(11))                                         
      U(75)=PAy*U(71)+WPy*U(50)+                                       
     &  r1x2Z*(U(69)-rpxZ*U(15))                                         
      U(78)=PAz*U(72)+WPz*U(51)+                                       
     &  r1x2Z*(U(69)-rpxZ*U(15))                                         
      U(95)=PAz*U(89)+WPz*U(59)+                                       
     &  r1x2Z*(U(86)-rpxZ*U(16))                                         
      U(107)=PAx*U(104)+WPx*U(60)+                                     
     &  r1x2Z*(U(103)-rpxZ*U(17))                                        
      U(112)=PAz*U(106)+WPz*U(63)+                                     
     &  r1x2Z*(U(103)-rpxZ*U(17))                                        
      U(126)=PAy*U(122)+WPy*U(64)+                                     
     &  r1x2Z*(U(120)-rpxZ*U(26))                                        
      U(141)=PAx*U(138)+WPx*U(66)+                                     
     &  r1x2Z*(U(137)-rpxZ*U(28))                                        
      U(158)=PAx*U(155)+WPx*U(77)+                                     
     &  r1x2Z*(U(154)-rpxZ*U(29))                                        
      U(160)=PAy*U(156)+WPy*U(79)+                                     
     &  r1x2Z*(U(154)-rpxZ*U(29))                                        
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

      C(I0,1)=C(I0,1)+U(1)*cp1*cq1                                     
      C(I0,2)=C(I0,2)+U(2)*cp3*cq1                                     
      C(I0,3)=C(I0,3)+U(3)*cp3*cq1                                     
      C(I0,4)=C(I0,4)+U(4)*cp3*cq1                                     
      C(I0,5)=C(I0,5)+U(2)*cp2*cq1                                     
      C(I0,6)=C(I0,6)+U(5)*cq1                                         
      C(I0,7)=C(I0,7)+U(1)*cp3*cq1                                     
      C(I0,8)=C(I0,8)+U(2)*cq1                                         
      C(I0,9)=C(I0,9)+U(3)*cp2*cq1                                     
      C(I0,10)=C(I0,10)+U(6)*cq1                                       
      C(I0,11)=C(I0,11)+U(7)*cq1                                       
      C(I0,12)=C(I0,12)+U(3)*cq1                                       
      C(I0,13)=C(I0,13)+U(4)*cp2*cq1                                   
      C(I0,14)=C(I0,14)+U(8)*cq1                                       
      C(I0,15)=C(I0,15)+U(9)*cq1                                       
      C(I0,16)=C(I0,16)+U(10)*cq1                                      
      C(I0,17)=C(I0,17)+U(4)*cq1                                       
      C(I0,18)=C(I0,18)+U(18)*cp1*cq3                                  
      C(I0,19)=C(I0,19)+U(19)*cp3*cq3                                  
      C(I0,20)=C(I0,20)+U(20)*cp3*cq3                                  
      C(I0,21)=C(I0,21)+U(21)*cp3*cq3                                  
      C(I0,22)=C(I0,22)+U(19)*cp2*cq3                                  
      C(I0,23)=C(I0,23)+U(22)*cq3                                      
      C(I0,24)=C(I0,24)+U(18)*cp3*cq3                                  
      C(I0,25)=C(I0,25)+U(19)*cq3                                      
      C(I0,26)=C(I0,26)+U(20)*cp2*cq3                                  
      C(I0,27)=C(I0,27)+U(23)*cq3                                      
      C(I0,28)=C(I0,28)+U(24)*cq3                                      
      C(I0,29)=C(I0,29)+U(20)*cq3                                      
      C(I0,30)=C(I0,30)+U(21)*cp2*cq3                                  
      C(I0,31)=C(I0,31)+U(25)*cq3                                      
      C(I0,32)=C(I0,32)+U(26)*cq3                                      
      C(I0,33)=C(I0,33)+U(27)*cq3                                      
      C(I0,34)=C(I0,34)+U(21)*cq3                                      
      C(I0,35)=C(I0,35)+U(35)*cp1*cq3                                  
      C(I0,36)=C(I0,36)+U(36)*cp3*cq3                                  
      C(I0,37)=C(I0,37)+U(37)*cp3*cq3                                  
      C(I0,38)=C(I0,38)+U(38)*cp3*cq3                                  
      C(I0,39)=C(I0,39)+U(36)*cp2*cq3                                  
      C(I0,40)=C(I0,40)+U(39)*cq3                                      
      C(I0,41)=C(I0,41)+U(35)*cp3*cq3                                  
      C(I0,42)=C(I0,42)+U(36)*cq3                                      
      C(I0,43)=C(I0,43)+U(37)*cp2*cq3                                  
      C(I0,44)=C(I0,44)+U(40)*cq3                                      
      C(I0,45)=C(I0,45)+U(41)*cq3                                      
      C(I0,46)=C(I0,46)+U(37)*cq3                                      
      C(I0,47)=C(I0,47)+U(38)*cp2*cq3                                  
      C(I0,48)=C(I0,48)+U(42)*cq3                                      
      C(I0,49)=C(I0,49)+U(43)*cq3                                      
      C(I0,50)=C(I0,50)+U(44)*cq3                                      
      C(I0,51)=C(I0,51)+U(38)*cq3                                      
      C(I0,52)=C(I0,52)+U(52)*cp1*cq3                                  
      C(I0,53)=C(I0,53)+U(53)*cp3*cq3                                  
      C(I0,54)=C(I0,54)+U(54)*cp3*cq3                                  
      C(I0,55)=C(I0,55)+U(55)*cp3*cq3                                  
      C(I0,56)=C(I0,56)+U(53)*cp2*cq3                                  
      C(I0,57)=C(I0,57)+U(56)*cq3                                      
      C(I0,58)=C(I0,58)+U(52)*cp3*cq3                                  
      C(I0,59)=C(I0,59)+U(53)*cq3                                      
      C(I0,60)=C(I0,60)+U(54)*cp2*cq3                                  
      C(I0,61)=C(I0,61)+U(57)*cq3                                      
      C(I0,62)=C(I0,62)+U(58)*cq3                                      
      C(I0,63)=C(I0,63)+U(54)*cq3                                      
      C(I0,64)=C(I0,64)+U(55)*cp2*cq3                                  
      C(I0,65)=C(I0,65)+U(59)*cq3                                      
      C(I0,66)=C(I0,66)+U(60)*cq3                                      
      C(I0,67)=C(I0,67)+U(61)*cq3                                      
      C(I0,68)=C(I0,68)+U(55)*cq3                                      
      C(I0,69)=C(I0,69)+U(18)*cp1*cq2                                  
      C(I0,70)=C(I0,70)+U(19)*cp3*cq2                                  
      C(I0,71)=C(I0,71)+U(20)*cp3*cq2                                  
      C(I0,72)=C(I0,72)+U(21)*cp3*cq2                                  
      C(I0,73)=C(I0,73)+U(19)*cp2*cq2                                  
      C(I0,74)=C(I0,74)+U(22)*cq2                                      
      C(I0,75)=C(I0,75)+U(18)*cp3*cq2                                  
      C(I0,76)=C(I0,76)+U(19)*cq2                                      
      C(I0,77)=C(I0,77)+U(20)*cp2*cq2                                  
      C(I0,78)=C(I0,78)+U(23)*cq2                                      
      C(I0,79)=C(I0,79)+U(24)*cq2                                      
      C(I0,80)=C(I0,80)+U(20)*cq2                                      
      C(I0,81)=C(I0,81)+U(21)*cp2*cq2                                  
      C(I0,82)=C(I0,82)+U(25)*cq2                                      
      C(I0,83)=C(I0,83)+U(26)*cq2                                      
      C(I0,84)=C(I0,84)+U(27)*cq2                                      
      C(I0,85)=C(I0,85)+U(21)*cq2                                      
      C(I0,86)=C(I0,86)+U(69)*cp1                                      
      C(I0,87)=C(I0,87)+U(70)*cp3                                      
      C(I0,88)=C(I0,88)+U(71)*cp3                                      
      C(I0,89)=C(I0,89)+U(72)*cp3                                      
      C(I0,90)=C(I0,90)+U(70)*cp2                                      
      C(I0,91)=C(I0,91)+U(73)                                          
      C(I0,92)=C(I0,92)+U(69)*cp3                                      
      C(I0,93)=C(I0,93)+U(70)                                          
      C(I0,94)=C(I0,94)+U(71)*cp2                                      
      C(I0,95)=C(I0,95)+U(74)                                          
      C(I0,96)=C(I0,96)+U(75)                                          
      C(I0,97)=C(I0,97)+U(71)                                          
      C(I0,98)=C(I0,98)+U(72)*cp2                                      
      C(I0,99)=C(I0,99)+U(76)                                          
      C(I0,100)=C(I0,100)+U(77)                                        
      C(I0,101)=C(I0,101)+U(78)                                        
      C(I0,102)=C(I0,102)+U(72)                                        
      C(I0,103)=C(I0,103)+U(1)*cp1*cq3                                 
      C(I0,104)=C(I0,104)+U(2)*cp3*cq3                                 
      C(I0,105)=C(I0,105)+U(3)*cp3*cq3                                 
      C(I0,106)=C(I0,106)+U(4)*cp3*cq3                                 
      C(I0,107)=C(I0,107)+U(2)*cp2*cq3                                 
      C(I0,108)=C(I0,108)+U(5)*cq3                                     
      C(I0,109)=C(I0,109)+U(1)*cp3*cq3                                 
      C(I0,110)=C(I0,110)+U(2)*cq3                                     
      C(I0,111)=C(I0,111)+U(3)*cp2*cq3                                 
      C(I0,112)=C(I0,112)+U(6)*cq3                                     
      C(I0,113)=C(I0,113)+U(7)*cq3                                     
      C(I0,114)=C(I0,114)+U(3)*cq3                                     
      C(I0,115)=C(I0,115)+U(4)*cp2*cq3                                 
      C(I0,116)=C(I0,116)+U(8)*cq3                                     
      C(I0,117)=C(I0,117)+U(9)*cq3                                     
      C(I0,118)=C(I0,118)+U(10)*cq3                                    
      C(I0,119)=C(I0,119)+U(4)*cq3                                     
      C(I0,120)=C(I0,120)+U(18)*cp1                                    
      C(I0,121)=C(I0,121)+U(19)*cp3                                    
      C(I0,122)=C(I0,122)+U(20)*cp3                                    
      C(I0,123)=C(I0,123)+U(21)*cp3                                    
      C(I0,124)=C(I0,124)+U(19)*cp2                                    
      C(I0,125)=C(I0,125)+U(22)                                        
      C(I0,126)=C(I0,126)+U(18)*cp3                                    
      C(I0,127)=C(I0,127)+U(19)                                        
      C(I0,128)=C(I0,128)+U(20)*cp2                                    
      C(I0,129)=C(I0,129)+U(23)                                        
      C(I0,130)=C(I0,130)+U(24)                                        
      C(I0,131)=C(I0,131)+U(20)                                        
      C(I0,132)=C(I0,132)+U(21)*cp2                                    
      C(I0,133)=C(I0,133)+U(25)                                        
      C(I0,134)=C(I0,134)+U(26)                                        
      C(I0,135)=C(I0,135)+U(27)                                        
      C(I0,136)=C(I0,136)+U(21)                                        
      C(I0,137)=C(I0,137)+U(35)*cp1*cq2                                
      C(I0,138)=C(I0,138)+U(36)*cp3*cq2                                
      C(I0,139)=C(I0,139)+U(37)*cp3*cq2                                
      C(I0,140)=C(I0,140)+U(38)*cp3*cq2                                
      C(I0,141)=C(I0,141)+U(36)*cp2*cq2                                
      C(I0,142)=C(I0,142)+U(39)*cq2                                    
      C(I0,143)=C(I0,143)+U(35)*cp3*cq2                                
      C(I0,144)=C(I0,144)+U(36)*cq2                                    
      C(I0,145)=C(I0,145)+U(37)*cp2*cq2                                
      C(I0,146)=C(I0,146)+U(40)*cq2                                    
      C(I0,147)=C(I0,147)+U(41)*cq2                                    
      C(I0,148)=C(I0,148)+U(37)*cq2                                    
      C(I0,149)=C(I0,149)+U(38)*cp2*cq2                                
      C(I0,150)=C(I0,150)+U(42)*cq2                                    
      C(I0,151)=C(I0,151)+U(43)*cq2                                    
      C(I0,152)=C(I0,152)+U(44)*cq2                                    
      C(I0,153)=C(I0,153)+U(38)*cq2                                    
      C(I0,154)=C(I0,154)+U(86)*cp1                                    
      C(I0,155)=C(I0,155)+U(87)*cp3                                    
      C(I0,156)=C(I0,156)+U(88)*cp3                                    
      C(I0,157)=C(I0,157)+U(89)*cp3                                    
      C(I0,158)=C(I0,158)+U(87)*cp2                                    
      C(I0,159)=C(I0,159)+U(90)                                        
      C(I0,160)=C(I0,160)+U(86)*cp3                                    
      C(I0,161)=C(I0,161)+U(87)                                        
      C(I0,162)=C(I0,162)+U(88)*cp2                                    
      C(I0,163)=C(I0,163)+U(91)                                        
      C(I0,164)=C(I0,164)+U(92)                                        
      C(I0,165)=C(I0,165)+U(88)                                        
      C(I0,166)=C(I0,166)+U(89)*cp2                                    
      C(I0,167)=C(I0,167)+U(93)                                        
      C(I0,168)=C(I0,168)+U(94)                                        
      C(I0,169)=C(I0,169)+U(95)                                        
      C(I0,170)=C(I0,170)+U(89)                                        
      C(I0,171)=C(I0,171)+U(103)*cp1                                   
      C(I0,172)=C(I0,172)+U(104)*cp3                                   
      C(I0,173)=C(I0,173)+U(105)*cp3                                   
      C(I0,174)=C(I0,174)+U(106)*cp3                                   
      C(I0,175)=C(I0,175)+U(104)*cp2                                   
      C(I0,176)=C(I0,176)+U(107)                                       
      C(I0,177)=C(I0,177)+U(103)*cp3                                   
      C(I0,178)=C(I0,178)+U(104)                                       
      C(I0,179)=C(I0,179)+U(105)*cp2                                   
      C(I0,180)=C(I0,180)+U(108)                                       
      C(I0,181)=C(I0,181)+U(109)                                       
      C(I0,182)=C(I0,182)+U(105)                                       
      C(I0,183)=C(I0,183)+U(106)*cp2                                   
      C(I0,184)=C(I0,184)+U(110)                                       
      C(I0,185)=C(I0,185)+U(111)                                       
      C(I0,186)=C(I0,186)+U(112)                                       
      C(I0,187)=C(I0,187)+U(106)                                       
      C(I0,188)=C(I0,188)+U(35)*cp1                                    
      C(I0,189)=C(I0,189)+U(36)*cp3                                    
      C(I0,190)=C(I0,190)+U(37)*cp3                                    
      C(I0,191)=C(I0,191)+U(38)*cp3                                    
      C(I0,192)=C(I0,192)+U(36)*cp2                                    
      C(I0,193)=C(I0,193)+U(39)                                        
      C(I0,194)=C(I0,194)+U(35)*cp3                                    
      C(I0,195)=C(I0,195)+U(36)                                        
      C(I0,196)=C(I0,196)+U(37)*cp2                                    
      C(I0,197)=C(I0,197)+U(40)                                        
      C(I0,198)=C(I0,198)+U(41)                                        
      C(I0,199)=C(I0,199)+U(37)                                        
      C(I0,200)=C(I0,200)+U(38)*cp2                                    
      C(I0,201)=C(I0,201)+U(42)                                        
      C(I0,202)=C(I0,202)+U(43)                                        
      C(I0,203)=C(I0,203)+U(44)                                        
      C(I0,204)=C(I0,204)+U(38)                                        
      C(I0,205)=C(I0,205)+U(52)*cp1*cq2                                
      C(I0,206)=C(I0,206)+U(53)*cp3*cq2                                
      C(I0,207)=C(I0,207)+U(54)*cp3*cq2                                
      C(I0,208)=C(I0,208)+U(55)*cp3*cq2                                
      C(I0,209)=C(I0,209)+U(53)*cp2*cq2                                
      C(I0,210)=C(I0,210)+U(56)*cq2                                    
      C(I0,211)=C(I0,211)+U(52)*cp3*cq2                                
      C(I0,212)=C(I0,212)+U(53)*cq2                                    
      C(I0,213)=C(I0,213)+U(54)*cp2*cq2                                
      C(I0,214)=C(I0,214)+U(57)*cq2                                    
      C(I0,215)=C(I0,215)+U(58)*cq2                                    
      C(I0,216)=C(I0,216)+U(54)*cq2                                    
      C(I0,217)=C(I0,217)+U(55)*cp2*cq2                                
      C(I0,218)=C(I0,218)+U(59)*cq2                                    
      C(I0,219)=C(I0,219)+U(60)*cq2                                    
      C(I0,220)=C(I0,220)+U(61)*cq2                                    
      C(I0,221)=C(I0,221)+U(55)*cq2                                    
      C(I0,222)=C(I0,222)+U(120)*cp1                                   
      C(I0,223)=C(I0,223)+U(121)*cp3                                   
      C(I0,224)=C(I0,224)+U(122)*cp3                                   
      C(I0,225)=C(I0,225)+U(123)*cp3                                   
      C(I0,226)=C(I0,226)+U(121)*cp2                                   
      C(I0,227)=C(I0,227)+U(124)                                       
      C(I0,228)=C(I0,228)+U(120)*cp3                                   
      C(I0,229)=C(I0,229)+U(121)                                       
      C(I0,230)=C(I0,230)+U(122)*cp2                                   
      C(I0,231)=C(I0,231)+U(125)                                       
      C(I0,232)=C(I0,232)+U(126)                                       
      C(I0,233)=C(I0,233)+U(122)                                       
      C(I0,234)=C(I0,234)+U(123)*cp2                                   
      C(I0,235)=C(I0,235)+U(127)                                       
      C(I0,236)=C(I0,236)+U(128)                                       
      C(I0,237)=C(I0,237)+U(129)                                       
      C(I0,238)=C(I0,238)+U(123)                                       
      C(I0,239)=C(I0,239)+U(137)*cp1                                   
      C(I0,240)=C(I0,240)+U(138)*cp3                                   
      C(I0,241)=C(I0,241)+U(139)*cp3                                   
      C(I0,242)=C(I0,242)+U(140)*cp3                                   
      C(I0,243)=C(I0,243)+U(138)*cp2                                   
      C(I0,244)=C(I0,244)+U(141)                                       
      C(I0,245)=C(I0,245)+U(137)*cp3                                   
      C(I0,246)=C(I0,246)+U(138)                                       
      C(I0,247)=C(I0,247)+U(139)*cp2                                   
      C(I0,248)=C(I0,248)+U(142)                                       
      C(I0,249)=C(I0,249)+U(143)                                       
      C(I0,250)=C(I0,250)+U(139)                                       
      C(I0,251)=C(I0,251)+U(140)*cp2                                   
      C(I0,252)=C(I0,252)+U(144)                                       
      C(I0,253)=C(I0,253)+U(145)                                       
      C(I0,254)=C(I0,254)+U(146)                                       
      C(I0,255)=C(I0,255)+U(140)                                       
      C(I0,256)=C(I0,256)+U(154)*cp1                                   
      C(I0,257)=C(I0,257)+U(155)*cp3                                   
      C(I0,258)=C(I0,258)+U(156)*cp3                                   
      C(I0,259)=C(I0,259)+U(157)*cp3                                   
      C(I0,260)=C(I0,260)+U(155)*cp2                                   
      C(I0,261)=C(I0,261)+U(158)                                       
      C(I0,262)=C(I0,262)+U(154)*cp3                                   
      C(I0,263)=C(I0,263)+U(155)                                       
      C(I0,264)=C(I0,264)+U(156)*cp2                                   
      C(I0,265)=C(I0,265)+U(159)                                       
      C(I0,266)=C(I0,266)+U(160)                                       
      C(I0,267)=C(I0,267)+U(156)                                       
      C(I0,268)=C(I0,268)+U(157)*cp2                                   
      C(I0,269)=C(I0,269)+U(161)                                       
      C(I0,270)=C(I0,270)+U(162)                                       
      C(I0,271)=C(I0,271)+U(163)                                       
      C(I0,272)=C(I0,272)+U(157)                                       
      C(I0,273)=C(I0,273)+U(52)*cp1                                    
      C(I0,274)=C(I0,274)+U(53)*cp3                                    
      C(I0,275)=C(I0,275)+U(54)*cp3                                    
      C(I0,276)=C(I0,276)+U(55)*cp3                                    
      C(I0,277)=C(I0,277)+U(53)*cp2                                    
      C(I0,278)=C(I0,278)+U(56)                                        
      C(I0,279)=C(I0,279)+U(52)*cp3                                    
      C(I0,280)=C(I0,280)+U(53)                                        
      C(I0,281)=C(I0,281)+U(54)*cp2                                    
      C(I0,282)=C(I0,282)+U(57)                                        
      C(I0,283)=C(I0,283)+U(58)                                        
      C(I0,284)=C(I0,284)+U(54)                                        
      C(I0,285)=C(I0,285)+U(55)*cp2                                    
      C(I0,286)=C(I0,286)+U(59)                                        
      C(I0,287)=C(I0,287)+U(60)                                        
      C(I0,288)=C(I0,288)+U(61)                                        
      C(I0,289)=C(I0,289)+U(55)                                        

      
 1000             CONTINUE
 1010             CONTINUE
 2000          CONTINUE
 2010          CONTINUE

            ENDIF
 3000    CONTINUE
 4000 CONTINUE 
      RETURN
      END
