      SUBROUTINE ONX6211(N,iBgn,iEnd,IndA,BfnBufB,BfnBufK,KBra,KKet,
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
      U(11)=R2                                                         
      U(12)=R3                                                         
      U(5)=R4                                                          
      U(2)=PAx*U(1)+WPx*U(11)                                          
      U(3)=PAy*U(1)+WPy*U(11)                                          
      U(4)=PAz*U(1)+WPz*U(11)                                          
      U(13)=PAx*U(11)+WPx*U(12)                                        
      U(14)=PAy*U(11)+WPy*U(12)                                        
      U(16)=PAz*U(11)+WPz*U(12)                                        
      U(15)=PAx*U(12)+WPx*U(5)                                         
      U(17)=PAy*U(12)+WPy*U(5)                                         
      U(18)=PAz*U(12)+WPz*U(5)                                         
      U(5)=PAx*U(2)+WPx*U(13)+                                         
     &  r1x2Z*(U(1)-rpxZ*U(11))                                          
      U(6)=PAx*U(3)+WPx*U(14)                                          
      U(7)=PAy*U(3)+WPy*U(14)+                                         
     &  r1x2Z*(U(1)-rpxZ*U(11))                                          
      U(8)=PAx*U(4)+WPx*U(16)                                          
      U(9)=PAy*U(4)+WPy*U(16)                                          
      U(10)=PAz*U(4)+WPz*U(16)+                                        
     &  r1x2Z*(U(1)-rpxZ*U(11))                                          
      U(19)=PAx*U(13)+WPx*U(15)+                                       
     &  r1x2Z*(U(11)-rpxZ*U(12))                                         
      U(20)=PAx*U(14)+WPx*U(17)                                        
      U(21)=PAy*U(14)+WPy*U(17)+                                       
     &  r1x2Z*(U(11)-rpxZ*U(12))                                         
      U(22)=PAx*U(16)+WPx*U(18)                                        
      U(23)=PAy*U(16)+WPy*U(18)                                        
      U(24)=PAz*U(16)+WPz*U(18)+                                       
     &  r1x2Z*(U(11)-rpxZ*U(12))                                         
      U(11)=PAx*U(5)+WPx*U(19)+                                        
     &  2.0D0*r1x2Z*(U(2)-rpxZ*U(13))                                    
      U(12)=PAx*U(6)+WPx*U(20)+                                        
     &  r1x2Z*(U(3)-rpxZ*U(14))                                          
      U(14)=PAy*U(7)+WPy*U(21)+                                        
     &  2.0D0*r1x2Z*(U(3)-rpxZ*U(14))                                    
      U(15)=PAx*U(8)+WPx*U(22)+                                        
     &  r1x2Z*(U(4)-rpxZ*U(16))                                          
      U(17)=PAy*U(9)+WPy*U(23)+                                        
     &  r1x2Z*(U(4)-rpxZ*U(16))                                          
      U(20)=PAz*U(10)+WPz*U(24)+                                       
     &  2.0D0*r1x2Z*(U(4)-rpxZ*U(16))                                    
      U(13)=PAx*U(7)+WPx*U(21)                                         
      U(16)=PAx*U(9)+WPx*U(23)                                         
      U(18)=PAx*U(10)+WPx*U(24)                                        
      U(19)=PAy*U(10)+WPy*U(24)                

      C(I0,1)=C(I0,1)+U(5)*cp1                                         
      C(I0,2)=C(I0,2)+U(5)                                             
      C(I0,3)=C(I0,3)+U(6)                                             
      C(I0,4)=C(I0,4)+U(7)                                             
      C(I0,5)=C(I0,5)+U(6)*cp1                                         
      C(I0,6)=C(I0,6)+U(8)                                             
      C(I0,7)=C(I0,7)+U(9)                                             
      C(I0,8)=C(I0,8)+U(10)                                            
      C(I0,9)=C(I0,9)+U(7)*cp1                                         
      C(I0,10)=C(I0,10)+U(11)                                          
      C(I0,11)=C(I0,11)+U(12)                                          
      C(I0,12)=C(I0,12)+U(13)                                          
      C(I0,13)=C(I0,13)+U(8)*cp1                                       
      C(I0,14)=C(I0,14)+U(14)                                          
      C(I0,15)=C(I0,15)+U(15)                                          
      C(I0,16)=C(I0,16)+U(16)                                          
      C(I0,17)=C(I0,17)+U(9)*cp1                                       
      C(I0,18)=C(I0,18)+U(17)                                          
      C(I0,19)=C(I0,19)+U(18)                                          
      C(I0,20)=C(I0,20)+U(19)                                          
      C(I0,21)=C(I0,21)+U(10)*cp1                                      
      C(I0,22)=C(I0,22)+U(20)           

 1000             CONTINUE
 1010             CONTINUE
 2000          CONTINUE
            ENDIF
 3000    CONTINUE
 4000 CONTINUE 
      RETURN
      END
