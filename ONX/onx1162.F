      SUBROUTINE ONX1162(N,iBgn,iEnd,IndA,BfnBufB,BfnBufK,KBra,KKet,
     &                   SwhB,SwhK,RAC,RBD,
     &                   PrmBufB,PrmBufK,U,C)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 U(24),C(N,22)
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

      CALL RLoader(C,N*22,0.0D0)

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
      U(11)=R2                                                         
      U(12)=R3                                                         
      U(5)=R4                                                          
      U(2)=QBx*U(1)+WQx*U(11)                                          
      U(3)=QBy*U(1)+WQy*U(11)                                          
      U(4)=QBz*U(1)+WQz*U(11)                                          
      U(13)=QBx*U(11)+WQx*U(12)                                        
      U(14)=QBy*U(11)+WQy*U(12)                                        
      U(16)=QBz*U(11)+WQz*U(12)                                        
      U(15)=QBx*U(12)+WQx*U(5)                                         
      U(17)=QBy*U(12)+WQy*U(5)                                         
      U(18)=QBz*U(12)+WQz*U(5)                                         
      Ts0=r1x2E*(U(1)-rpxE*U(11))
      U(5)=QBx*U(2)+WQx*U(13)+ Ts0
      U(6)=QBx*U(3)+WQx*U(14)                                          
      U(7)=QBy*U(3)+WQy*U(14)+ Ts0
      U(8)=QBx*U(4)+WQx*U(16)                                          
      U(9)=QBy*U(4)+WQy*U(16)                                          
      U(10)=QBz*U(4)+WQz*U(16)+ Ts0
      Ts1=r1x2E*(U(11)-rpxE*U(12))
      U(19)=QBx*U(13)+WQx*U(15)+ Ts1
      U(20)=QBx*U(14)+WQx*U(17)                                        
      U(21)=QBy*U(14)+WQy*U(17)+ Ts1
      U(22)=QBx*U(16)+WQx*U(18)                                        
      U(23)=QBy*U(16)+WQy*U(18)                                        
      U(24)=QBz*U(16)+WQz*U(18)+ Ts1
      U(11)=QBx*U(5)+WQx*U(19)+                                        
     &  2.0D0*r1x2E*(U(2)-rpxE*U(13))                                    
      U(12)=QBx*U(6)+WQx*U(20)+                                        
     &  r1x2E*(U(3)-rpxE*U(14))                                          
      U(14)=QBy*U(7)+WQy*U(21)+                                        
     &  2.0D0*r1x2E*(U(3)-rpxE*U(14))                                    
      U(15)=QBx*U(8)+WQx*U(22)+                                        
     &  r1x2E*(U(4)-rpxE*U(16))                                          
      U(17)=QBy*U(9)+WQy*U(23)+                                        
     &  r1x2E*(U(4)-rpxE*U(16))                                          
      U(20)=QBz*U(10)+WQz*U(24)+                                       
     &  2.0D0*r1x2E*(U(4)-rpxE*U(16))                                    
      U(13)=QBx*U(7)+WQx*U(21)                                         
      U(16)=QBx*U(9)+WQx*U(23)                                         
      U(18)=QBx*U(10)+WQx*U(24)                                        
      U(19)=QBy*U(10)+WQy*U(24)                    

      C(I0,1)=C(I0,1)+U(5)*cq1                                         
      C(I0,2)=C(I0,2)+U(5)                                             
      C(I0,3)=C(I0,3)+U(6)                                             
      C(I0,4)=C(I0,4)+U(7)                                             
      C(I0,5)=C(I0,5)+U(6)*cq1                                         
      C(I0,6)=C(I0,6)+U(8)                                             
      C(I0,7)=C(I0,7)+U(9)                                             
      C(I0,8)=C(I0,8)+U(10)                                            
      C(I0,9)=C(I0,9)+U(7)*cq1                                         
      C(I0,10)=C(I0,10)+U(11)                                          
      C(I0,11)=C(I0,11)+U(12)                                          
      C(I0,12)=C(I0,12)+U(13)                                          
      C(I0,13)=C(I0,13)+U(8)*cq1                                       
      C(I0,14)=C(I0,14)+U(14)                                          
      C(I0,15)=C(I0,15)+U(15)                                          
      C(I0,16)=C(I0,16)+U(16)                                          
      C(I0,17)=C(I0,17)+U(9)*cq1                                       
      C(I0,18)=C(I0,18)+U(17)                                          
      C(I0,19)=C(I0,19)+U(18)                                          
      C(I0,20)=C(I0,20)+U(19)                                          
      C(I0,21)=C(I0,21)+U(10)*cq1                                      
      C(I0,22)=C(I0,22)+U(20)                      

 1000             CONTINUE
 1010             CONTINUE
 2000          CONTINUE
            ENDIF
 3000    CONTINUE
 4000 CONTINUE 
      RETURN
      END
