      SUBROUTINE ONX6611(N,iBgn,iEnd,IndA,BfnBufB,BfnBufK,KBra,KKet,
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
         Is=I-iBgn+1
         Ixa = BfnBufB(I)
         DO 3000 J=1,IndA(I)
            Ixb = BfnBufK(J)
            IF(Ixa.GE.Ixb) THEN
               I0=I0+1
               rCx=Cx+SwhB*RAC(1,Is)
               rCy=Cy+SwhB*RAC(2,Is)
               rCz=Cz+SwhB*RAC(3,Is)
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
      U(5)=PAx*U(2)+WPx*U(21)+
     &  r1x2Z*(U(1)-rpxZ*U(11))
      U(6)=PAx*U(3)+WPx*U(22)
      U(7)=PAy*U(3)+WPy*U(22)+
     &  r1x2Z*(U(1)-rpxZ*U(11))
      U(8)=PAx*U(4)+WPx*U(23)
      U(9)=PAy*U(4)+WPy*U(23)
      U(10)=PAz*U(4)+WPz*U(23)+
     &  r1x2Z*(U(1)-rpxZ*U(11))
      U(27)=PAx*U(21)+WPx*U(24)+
     &  r1x2Z*(U(11)-rpxZ*U(12))
      U(28)=PAx*U(22)+WPx*U(25)
      U(29)=PAy*U(22)+WPy*U(25)+
     &  r1x2Z*(U(11)-rpxZ*U(12))
      U(30)=PAx*U(23)+WPx*U(26)
      U(31)=PAy*U(23)+WPy*U(26)
      U(33)=PAz*U(23)+WPz*U(26)+
     &  r1x2Z*(U(11)-rpxZ*U(12))
      U(32)=PAx*U(24)+WPx*U(14)+
     &  r1x2Z*(U(12)-rpxZ*U(13))
      U(34)=PAx*U(25)+WPx*U(15)
      U(35)=PAy*U(25)+WPy*U(15)+
     &  r1x2Z*(U(12)-rpxZ*U(13))
      U(36)=PAx*U(26)+WPx*U(16)
      U(37)=PAy*U(26)+WPy*U(16)
      U(38)=PAz*U(26)+WPz*U(16)+
     &  r1x2Z*(U(12)-rpxZ*U(13))
      U(11)=PAx*U(5)+WPx*U(27)+
     &  2.0D0*r1x2Z*(U(2)-rpxZ*U(21))
      U(12)=PAx*U(6)+WPx*U(28)+
     &  r1x2Z*(U(3)-rpxZ*U(22))
      U(13)=PAx*U(7)+WPx*U(29)
      U(14)=PAy*U(7)+WPy*U(29)+
     &  2.0D0*r1x2Z*(U(3)-rpxZ*U(22))
      U(15)=PAx*U(8)+WPx*U(30)+
     &  r1x2Z*(U(4)-rpxZ*U(23))
      U(16)=PAx*U(9)+WPx*U(31)
      U(17)=PAy*U(9)+WPy*U(31)+
     &  r1x2Z*(U(4)-rpxZ*U(23))
      U(18)=PAx*U(10)+WPx*U(33)
      U(19)=PAy*U(10)+WPy*U(33)
      U(20)=PAz*U(10)+WPz*U(33)+
     &  2.0D0*r1x2Z*(U(4)-rpxZ*U(23))
      U(39)=PAx*U(27)+WPx*U(32)+
     &  2.0D0*r1x2Z*(U(21)-rpxZ*U(24))
      U(40)=PAx*U(28)+WPx*U(34)+
     &  r1x2Z*(U(22)-rpxZ*U(25))
      U(41)=PAx*U(29)+WPx*U(35)
      U(42)=PAy*U(29)+WPy*U(35)+
     &  2.0D0*r1x2Z*(U(22)-rpxZ*U(25))
      U(43)=PAx*U(30)+WPx*U(36)+
     &  r1x2Z*(U(23)-rpxZ*U(26))
      U(44)=PAx*U(31)+WPx*U(37)
      U(45)=PAy*U(31)+WPy*U(37)+
     &  r1x2Z*(U(23)-rpxZ*U(26))
      U(46)=PAx*U(33)+WPx*U(38)
      U(47)=PAy*U(33)+WPy*U(38)
      U(48)=PAz*U(33)+WPz*U(38)+
     &  2.0D0*r1x2Z*(U(23)-rpxZ*U(26))
      U(21)=PAx*U(11)+WPx*U(39)+
     &  3.0D0*r1x2Z*(U(5)-rpxZ*U(27))
      U(22)=PAx*U(12)+WPx*U(40)+
     &  2.0D0*r1x2Z*(U(6)-rpxZ*U(28))
      U(23)=PAx*U(13)+WPx*U(41)+
     &  r1x2Z*(U(7)-rpxZ*U(29))
      U(25)=PAy*U(14)+WPy*U(42)+
     &  3.0D0*r1x2Z*(U(7)-rpxZ*U(29))
      U(26)=PAx*U(15)+WPx*U(43)+
     &  2.0D0*r1x2Z*(U(8)-rpxZ*U(30))
      U(27)=PAx*U(16)+WPx*U(44)+
     &  r1x2Z*(U(9)-rpxZ*U(31))
      U(29)=PAy*U(17)+WPy*U(45)+
     &  2.0D0*r1x2Z*(U(9)-rpxZ*U(31))
      U(30)=PAx*U(18)+WPx*U(46)+
     &  r1x2Z*(U(10)-rpxZ*U(33))
      U(32)=PAy*U(19)+WPy*U(47)+
     &  r1x2Z*(U(10)-rpxZ*U(33))
      U(35)=PAz*U(20)+WPz*U(48)+
     &  3.0D0*r1x2Z*(U(10)-rpxZ*U(33))
      U(24)=PAx*U(14)+WPx*U(42)
      U(28)=PAx*U(17)+WPx*U(45)
      U(31)=PAx*U(19)+WPx*U(47)
      U(33)=PAx*U(20)+WPx*U(48)
      U(34)=PAy*U(20)+WPy*U(48)

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
