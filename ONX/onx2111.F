      SUBROUTINE ONX2111(N,iBgn,iEnd,IndA,BfnBufB,BfnBufK,KBra,KKet,
     &                   SwhB,SwhK,RAC,RBD,
     &                   PrmBufB,PrmBufK,U,C)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 U(4),C(N,4)
      REAL*8 F1_0(0:486),F1_1(0:486),F1_2(0:486),F1_3(0:486)
      REAL*8 EX_0(0:486),EX_1(0:486),EX_2(0:486),EX_3(0:486)
      REAL*8 PrmBufB(MXB,KBra+MInfo,*)
      REAL*8 PrmBufK(MXB,KKet+MInfo,*)
      REAL*8 RAC(3,*),RBD(3,*)
      INTEGER IndA(0:*)
      INTEGER BfnBufB(*),BfnBufK(*)
      INCLUDE "bkonx.inc"
      INCLUDE "Gamma_A.inc"
      INCLUDE "Gamma_1.inc"
      INCLUDE "Exponential.inc"

      CALL RLoader(C,N*4,0.0D0)

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
                  EtaQx=Eta*Qx
                  EtaQy=Eta*Qy
                  EtaQz=Eta*Qz

                  DO 1000 IKB=1,KBra
                     IF(PrmBufB(9,IKB,I).LE.TTsh) GOTO 1010
                     Zeta = PrmBufB(1,IKB,I)
                     Px   = PrmBufB(2,IKB,I)
                     Py   = PrmBufB(3,IKB,I)
                     Pz   = PrmBufB(4,IKB,I)
                     Up   = PrmBufB(5,IKB,I)
                     cp1  = PrmBufB(6,IKB,I)
                     ZpE=Zeta+Eta
                     r1xZpE=1.0D0/ZpE
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
                     PQx = Px - Qx
                     PQy = Py - Qy
                     PQz = Pz - Qz
                     T=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
                     IF(T.LT.Switch) THEN
                        m=INT(T*Grid)
                        T2=T*T
                        T3=T2*T
                        G1=F1_0(m)+T*F1_1(m)+T2*F1_2(m)+T3*F1_3(m)
                        ET=EX_0(m)+T*EX_1(m)+T2*EX_2(m)+T3*EX_3(m)
                        TwoT=2.0D0*T
                        G0=TwoT*G1+ET
                        R1=Rkk*G0
                        R2=Rkk*G1
                     ELSE
                        T1=1.0D0/T
                        TS=DSQRT(T1)
                        R1=Rkk*F0Asymp*Ts
                        Ts=Ts*T1
                        R2=Rkk*F1Asymp*Ts
                     ENDIF



      U(1)=R1                                                          
      U(4)=R2                                                          
      U(2)=PAx*U(1)+WPx*U(4)                                           
      U(3)=PAy*U(1)+WPy*U(4)                                           
      U(4)=PAz*U(1)+WPz*U(4)         

      C(I0,1)=C(I0,1)+U(1)*cp1                                         
      C(I0,2)=C(I0,2)+U(2)                                             
      C(I0,3)=C(I0,3)+U(3)                                             
      C(I0,4)=C(I0,4)+U(4)           

 1000             CONTINUE
 1010             CONTINUE
 2000          CONTINUE
            ENDIF
 3000    CONTINUE
 4000 CONTINUE 
      RETURN
      END
