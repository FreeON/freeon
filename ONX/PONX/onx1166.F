      SUBROUTINE ONX1166(N,iBgn,iEnd,IndA,BfnBufB,BfnBufK,KBra,KKet,
     &                   SwhB,SwhK,RAC,RBD,
     &                   PrmBufB,PrmBufK,U,C)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 U(48),C(N,31)
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

      CALL RLoader(C,N*31,0.0D0)

      I0=0
      DO 4000 I=iBgn,iEnd
         Ixa = BfnBufB(I)
         DO 3000 J=1,IndA(I)
            Ixb = BfnBufK(J)
            IF(Ixa.GE.Ixb) THEN
               I0=I0+1
               rDx=Dx+SwhK*RBD(1,J)
               rDy=Dy+SwhK*RBD(2,J)
               rDz=Dz+SwhK*RBD(3,J)
               DO 2000 IKK=1,KKet
                  Eta  = PrmBufK(1,IKK,J)
                  Qx   = PrmBufK(2,IKK,J)
                  Qy   = PrmBufK(3,IKK,J)
                  Qz   = PrmBufK(4,IKK,J)
                  Uq   = PrmBufK(5,IKK,J)
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
      U(13)=R4
      U(5)=R5
      U(2)=QBx*U(1)+WQx*U(11)
      U(3)=QBy*U(1)+WQy*U(11)
      U(4)=QBz*U(1)+WQz*U(11)
      U(21)=QBx*U(11)+WQx*U(12)
      U(22)=QBy*U(11)+WQy*U(12)
      U(23)=QBz*U(11)+WQz*U(12)
      U(24)=QBx*U(12)+WQx*U(13)
      U(25)=QBy*U(12)+WQy*U(13)
      U(26)=QBz*U(12)+WQz*U(13)
      U(14)=QBx*U(13)+WQx*U(5)
      U(15)=QBy*U(13)+WQy*U(5)
      U(16)=QBz*U(13)+WQz*U(5)
      Ts0=r1x2E*(U(1)-rpxE*U(11))
      U(5)=QBx*U(2)+WQx*U(21)+ Ts0
      U(6)=QBx*U(3)+WQx*U(22)
      U(7)=QBy*U(3)+WQy*U(22)+ Ts0
      U(8)=QBx*U(4)+WQx*U(23)
      U(9)=QBy*U(4)+WQy*U(23)
      U(10)=QBz*U(4)+WQz*U(23)+ Ts0
      Ts1=r1x2E*(U(11)-rpxE*U(12))
      U(27)=QBx*U(21)+WQx*U(24)+ Ts1
      U(28)=QBx*U(22)+WQx*U(25)
      U(29)=QBy*U(22)+WQy*U(25)+ Ts1
      U(30)=QBx*U(23)+WQx*U(26)
      U(31)=QBy*U(23)+WQy*U(26)
      U(33)=QBz*U(23)+WQz*U(26)+ Ts1
      Ts2=r1x2E*(U(12)-rpxE*U(13))
      U(32)=QBx*U(24)+WQx*U(14)+ Ts2
      U(34)=QBx*U(25)+WQx*U(15)
      U(35)=QBy*U(25)+WQy*U(15)+ Ts2
      U(36)=QBx*U(26)+WQx*U(16)
      U(37)=QBy*U(26)+WQy*U(16)
      U(38)=QBz*U(26)+WQz*U(16)+ Ts2
      U(11)=QBx*U(5)+WQx*U(27)+
     &  2.0D0*r1x2E*(U(2)-rpxE*U(21))
      U(12)=QBx*U(6)+WQx*U(28)+
     &  r1x2E*(U(3)-rpxE*U(22))
      U(13)=QBx*U(7)+WQx*U(29)
      U(14)=QBy*U(7)+WQy*U(29)+
     &  2.0D0*r1x2E*(U(3)-rpxE*U(22))
      U(15)=QBx*U(8)+WQx*U(30)+
     &  r1x2E*(U(4)-rpxE*U(23))
      U(16)=QBx*U(9)+WQx*U(31)
      U(17)=QBy*U(9)+WQy*U(31)+
     &  r1x2E*(U(4)-rpxE*U(23))
      U(18)=QBx*U(10)+WQx*U(33)
      U(19)=QBy*U(10)+WQy*U(33)
      U(20)=QBz*U(10)+WQz*U(33)+
     &  2.0D0*r1x2E*(U(4)-rpxE*U(23))
      U(39)=QBx*U(27)+WQx*U(32)+
     &  2.0D0*r1x2E*(U(21)-rpxE*U(24))
      U(40)=QBx*U(28)+WQx*U(34)+
     &  r1x2E*(U(22)-rpxE*U(25))
      U(41)=QBx*U(29)+WQx*U(35)
      U(42)=QBy*U(29)+WQy*U(35)+
     &  2.0D0*r1x2E*(U(22)-rpxE*U(25))
      U(43)=QBx*U(30)+WQx*U(36)+
     &  r1x2E*(U(23)-rpxE*U(26))
      U(44)=QBx*U(31)+WQx*U(37)
      U(45)=QBy*U(31)+WQy*U(37)+
     &  r1x2E*(U(23)-rpxE*U(26))
      U(46)=QBx*U(33)+WQx*U(38)
      U(47)=QBy*U(33)+WQy*U(38)
      U(48)=QBz*U(33)+WQz*U(38)+
     &  2.0D0*r1x2E*(U(23)-rpxE*U(26))
      U(21)=QBx*U(11)+WQx*U(39)+
     &  3.0D0*r1x2E*(U(5)-rpxE*U(27))
      U(22)=QBx*U(12)+WQx*U(40)+
     &  2.0D0*r1x2E*(U(6)-rpxE*U(28))
      U(23)=QBx*U(13)+WQx*U(41)+
     &  r1x2E*(U(7)-rpxE*U(29))
      U(25)=QBy*U(14)+WQy*U(42)+
     &  3.0D0*r1x2E*(U(7)-rpxE*U(29))
      U(26)=QBx*U(15)+WQx*U(43)+
     &  2.0D0*r1x2E*(U(8)-rpxE*U(30))
      U(27)=QBx*U(16)+WQx*U(44)+
     &  r1x2E*(U(9)-rpxE*U(31))
      U(29)=QBy*U(17)+WQy*U(45)+
     &  2.0D0*r1x2E*(U(9)-rpxE*U(31))
      U(30)=QBx*U(18)+WQx*U(46)+
     &  r1x2E*(U(10)-rpxE*U(33))
      U(32)=QBy*U(19)+WQy*U(47)+
     &  r1x2E*(U(10)-rpxE*U(33))
      U(35)=QBz*U(20)+WQz*U(48)+
     &  3.0D0*r1x2E*(U(10)-rpxE*U(33))
      U(24)=QBx*U(14)+WQx*U(42)
      U(28)=QBx*U(17)+WQx*U(45)
      U(31)=QBx*U(19)+WQx*U(47)
      U(33)=QBx*U(20)+WQx*U(48)
      U(34)=QBy*U(20)+WQy*U(48)

      C(I0,1)=C(I0,1)+U(5)                                             
      C(I0,2)=C(I0,2)+U(6)                                             
      C(I0,3)=C(I0,3)+U(7)                                             
      C(I0,4)=C(I0,4)+U(8)                                             
      C(I0,5)=C(I0,5)+U(9)                                             
      C(I0,6)=C(I0,6)+U(10)                                            
      C(I0,7)=C(I0,7)+U(11)                                            
      C(I0,8)=C(I0,8)+U(12)                                            
      C(I0,9)=C(I0,9)+U(13)                                            
      C(I0,10)=C(I0,10)+U(14)                                          
      C(I0,11)=C(I0,11)+U(15)                                          
      C(I0,12)=C(I0,12)+U(16)                                          
      C(I0,13)=C(I0,13)+U(17)                                          
      C(I0,14)=C(I0,14)+U(18)                                          
      C(I0,15)=C(I0,15)+U(19)                                          
      C(I0,16)=C(I0,16)+U(20)                                          
      C(I0,17)=C(I0,17)+U(21)                                          
      C(I0,18)=C(I0,18)+U(22)                                          
      C(I0,19)=C(I0,19)+U(23)                                          
      C(I0,20)=C(I0,20)+U(24)                                          
      C(I0,21)=C(I0,21)+U(25)                                          
      C(I0,22)=C(I0,22)+U(26)                                          
      C(I0,23)=C(I0,23)+U(27)                                          
      C(I0,24)=C(I0,24)+U(28)                                          
      C(I0,25)=C(I0,25)+U(29)                                          
      C(I0,26)=C(I0,26)+U(30)                                          
      C(I0,27)=C(I0,27)+U(31)                                          
      C(I0,28)=C(I0,28)+U(32)                                          
      C(I0,29)=C(I0,29)+U(33)                                          
      C(I0,30)=C(I0,30)+U(34)                                          
      C(I0,31)=C(I0,31)+U(35)                           

 1000             CONTINUE
 1010             CONTINUE
 2000          CONTINUE
            ENDIF
 3000    CONTINUE
 4000 CONTINUE 
      RETURN
      END
