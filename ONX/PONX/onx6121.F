      SUBROUTINE ONX6121(N,iBgn,iEnd,IndA,BfnBufB,BfnBufK,KBra,KKet,
     &                   SwhB,SwhK,RAC,RBD,
     &                   PrmBufB,PrmBufK,U,C)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 U(61),C(N,57)
      REAL*8 F3_0(0:486),F3_1(0:486),F3_2(0:486),F3_3(0:486)
      REAL*8 EX_0(0:486),EX_1(0:486),EX_2(0:486),EX_3(0:486)
      REAL*8 PrmBufB(MXB,KBra+MInfo,*)
      REAL*8 PrmBufK(MXB,KKet+MInfo,*)
      REAL*8 RAC(3,*),RBD(3,*)
      INTEGER IndA(0:*)
      INTEGER BfnBufB(*),BfnBufK(*)
      INCLUDE "bkonx.inc"
      INCLUDE "Gamma_A.inc"
      INCLUDE "Gamma_3.inc"
      INCLUDE "Exponential.inc"

      CALL RLoader(C,N*57,0.0D0)

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
                        G3=F3_0(m)+T*F3_1(m)+T2*F3_2(m)+T3*F3_3(m)
                        ET=EX_0(m)+T*EX_1(m)+T2*EX_2(m)+T3*EX_3(m)
                        TwoT=2.0D0*T
                        G2=.200000000000000D+00*(TwoT*G3+ET)
                        G1=.333333333333333D+00*(TwoT*G2+ET)
                        G0=TwoT*G1+ET
                        R1=Rkk*G0
                        R2=Rkk*G1
                        R3=Rkk*G2
                        R4=Rkk*G3
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
                     ENDIF

      U(1)=R1
      U(6)=R2
      U(5)=R3
      U(2)=R4
      U(18)=QBx*U(1)+WQx*U(6)
      U(35)=QBy*U(1)+WQy*U(6)
      U(52)=QBz*U(1)+WQz*U(6)
      U(8)=QBx*U(6)+WQx*U(5)
      U(9)=QBy*U(6)+WQy*U(5)
      U(11)=QBz*U(6)+WQz*U(5)
      U(7)=QBx*U(5)+WQx*U(2)
      U(10)=QBy*U(5)+WQy*U(2)
      U(12)=QBz*U(5)+WQz*U(2)
      U(2)=PAx*U(1)+WPx*U(6)
      U(3)=PAy*U(1)+WPy*U(6)
      U(4)=PAz*U(1)+WPz*U(6)
      U(13)=PAx*U(6)+WPx*U(5)
      U(14)=PAy*U(6)+WPy*U(5)
      U(15)=PAz*U(6)+WPz*U(5)
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
      U(16)=PAx*U(8)+WPx*U(7)+
     &  r1x2ZE*U(5)
      U(17)=PAy*U(8)+WPy*U(7)
      U(26)=PAz*U(8)+WPz*U(7)
      U(28)=PAx*U(9)+WPx*U(10)
      U(29)=PAy*U(9)+WPy*U(10)+
     &  r1x2ZE*U(5)
      U(30)=PAz*U(9)+WPz*U(10)
      U(31)=PAx*U(11)+WPx*U(12)
      U(32)=PAy*U(11)+WPy*U(12)
      U(33)=PAz*U(11)+WPz*U(12)+
     &  r1x2ZE*U(5)
      U(22)=PAx*U(19)+WPx*U(16)+
     &  r1x2Z*(U(18)-rpxZ*U(8))+
     &  r1x2ZE*U(13)
      U(23)=PAx*U(20)+WPx*U(17)+
     &  r1x2ZE*U(14)
      U(25)=PAx*U(21)+WPx*U(26)+
     &  r1x2ZE*U(15)
      U(41)=PAy*U(37)+WPy*U(29)+
     &  r1x2Z*(U(35)-rpxZ*U(9))+
     &  r1x2ZE*U(14)
      U(43)=PAy*U(38)+WPy*U(30)+
     &  r1x2ZE*U(15)
      U(61)=PAz*U(55)+WPz*U(33)+
     &  r1x2Z*(U(52)-rpxZ*U(11))+
     &  r1x2ZE*U(15)
      U(5)=PAx*U(2)+WPx*U(13)+
     &  r1x2Z*(U(1)-rpxZ*U(6))
      U(7)=PAy*U(3)+WPy*U(14)+
     &  r1x2Z*(U(1)-rpxZ*U(6))
      U(10)=PAz*U(4)+WPz*U(15)+
     &  r1x2Z*(U(1)-rpxZ*U(6))
      U(24)=PAy*U(20)+WPy*U(17)+
     &  r1x2Z*(U(18)-rpxZ*U(8))
      U(27)=PAz*U(21)+WPz*U(26)+
     &  r1x2Z*(U(18)-rpxZ*U(8))
      U(39)=PAx*U(36)+WPx*U(28)+
     &  r1x2Z*(U(35)-rpxZ*U(9))
      U(44)=PAz*U(38)+WPz*U(30)+
     &  r1x2Z*(U(35)-rpxZ*U(9))
      U(56)=PAx*U(53)+WPx*U(31)+
     &  r1x2Z*(U(52)-rpxZ*U(11))
      U(58)=PAy*U(54)+WPy*U(32)+
     &  r1x2Z*(U(52)-rpxZ*U(11))
      U(6)=PAx*U(3)+WPx*U(14)
      U(8)=PAx*U(4)+WPx*U(15)
      U(9)=PAy*U(4)+WPy*U(15)
      U(26)=PAy*U(21)+WPy*U(26)
      U(40)=PAx*U(37)+WPx*U(29)
      U(42)=PAx*U(38)+WPx*U(30)
      U(57)=PAx*U(54)+WPx*U(32)
      U(59)=PAx*U(55)+WPx*U(33)
      U(60)=PAy*U(55)+WPy*U(33)

      C(I0,1)=C(I0,1)+U(5)*cq1                                         
      C(I0,2)=C(I0,2)+U(6)*cq1                                         
      C(I0,3)=C(I0,3)+U(7)*cq1                                         
      C(I0,4)=C(I0,4)+U(8)*cq1                                         
      C(I0,5)=C(I0,5)+U(9)*cq1                                         
      C(I0,6)=C(I0,6)+U(10)*cq1                                        
      C(I0,18)=C(I0,18)+U(22)                                          
      C(I0,19)=C(I0,19)+U(23)                                          
      C(I0,20)=C(I0,20)+U(24)                                          
      C(I0,21)=C(I0,21)+U(25)                                          
      C(I0,22)=C(I0,22)+U(26)                                          
      C(I0,23)=C(I0,23)+U(27)                                          
      C(I0,35)=C(I0,35)+U(39)                                          
      C(I0,36)=C(I0,36)+U(40)                                          
      C(I0,37)=C(I0,37)+U(41)                                          
      C(I0,38)=C(I0,38)+U(42)                                          
      C(I0,39)=C(I0,39)+U(43)                                          
      C(I0,40)=C(I0,40)+U(44)                                          
      C(I0,52)=C(I0,52)+U(56)                                          
      C(I0,53)=C(I0,53)+U(57)                                          
      C(I0,54)=C(I0,54)+U(58)                                          
      C(I0,55)=C(I0,55)+U(59)                                          
      C(I0,56)=C(I0,56)+U(60)                                          
      C(I0,57)=C(I0,57)+U(61)                         

 1000             CONTINUE
 1010             CONTINUE
 2000          CONTINUE
            ENDIF
 3000    CONTINUE
 4000 CONTINUE 
      RETURN
      END
