      SUBROUTINE ONX2161(N,iBgn,iEnd,IndA,BfnBufB,BfnBufK,KBra,KKet,
     &                   SwhB,SwhK,RAC,RBD,
     &                   PrmBufB,PrmBufK,U,C)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 U(40),C(N,24)
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

      CALL RLoader(C,N*24,0.0D0)

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
      U(17)=QBx*U(5)+WQx*U(8)+                                         
     &  r1x2E*(U(1)-rpxE*U(4))                                           
      U(21)=QBx*U(9)+WQx*U(12)                                         
      U(25)=QBy*U(9)+WQy*U(12)+                                        
     &  r1x2E*(U(1)-rpxE*U(4))                                           
      U(29)=QBx*U(13)+WQx*U(15)                                        
      U(33)=QBy*U(13)+WQy*U(15)                                        
      U(37)=QBz*U(13)+WQz*U(15)+                                       
     &  r1x2E*(U(1)-rpxE*U(4))                                           
      U(20)=QBx*U(8)+WQx*U(6)+                                         
     &  r1x2E*(U(4)-rpxE*U(2))                                           
      U(24)=QBx*U(12)+WQx*U(7)                                         
      U(28)=QBy*U(12)+WQy*U(7)+                                        
     &  r1x2E*(U(4)-rpxE*U(2))                                           
      U(31)=QBx*U(15)+WQx*U(10)                                        
      U(34)=QBy*U(15)+WQy*U(10)                                        
      U(39)=QBz*U(15)+WQz*U(10)+                                       
     &  r1x2E*(U(4)-rpxE*U(2))                                           
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

      C(I0,1)=C(I0,1)+U(17)*cp1                                        
      C(I0,2)=C(I0,2)+U(18)                                            
      C(I0,3)=C(I0,3)+U(19)                                            
      C(I0,4)=C(I0,4)+U(20)                                            
      C(I0,5)=C(I0,5)+U(21)*cp1                                        
      C(I0,6)=C(I0,6)+U(22)                                            
      C(I0,7)=C(I0,7)+U(23)                                            
      C(I0,8)=C(I0,8)+U(24)                                            
      C(I0,9)=C(I0,9)+U(25)*cp1                                        
      C(I0,10)=C(I0,10)+U(26)                                          
      C(I0,11)=C(I0,11)+U(27)                                          
      C(I0,12)=C(I0,12)+U(28)                                          
      C(I0,13)=C(I0,13)+U(29)*cp1                                      
      C(I0,14)=C(I0,14)+U(30)                                          
      C(I0,15)=C(I0,15)+U(31)                                          
      C(I0,16)=C(I0,16)+U(32)                                          
      C(I0,17)=C(I0,17)+U(33)*cp1                                      
      C(I0,18)=C(I0,18)+U(34)                                          
      C(I0,19)=C(I0,19)+U(35)                                          
      C(I0,20)=C(I0,20)+U(36)                                          
      C(I0,21)=C(I0,21)+U(37)*cp1                                      
      C(I0,22)=C(I0,22)+U(38)                                          
      C(I0,23)=C(I0,23)+U(39)                                          
      C(I0,24)=C(I0,24)+U(40)               

 1000             CONTINUE
 1010             CONTINUE
 2000          CONTINUE
            ENDIF
 3000    CONTINUE
 4000 CONTINUE 
      RETURN
      END
