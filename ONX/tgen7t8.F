      SUBROUTINE TGen7t8(Ltot,KBra,KKet,N,Ind,
     >                   SwhB,SwhK,
     >                   PrmBufB,PrmBufK,
     >                   CBra,CKet,WR,WZ,RAC,RBD,R)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 PrmBufB(MXB,KBra+MInfo,*)
      REAL*8 PrmBufK(MXB,KKet+MInfo,*)
      REAL*8 CBra(N,KBra,*),CKet(N,KKet,*)
      REAL*8 WR(12,*),WZ(5,*)
      REAL*8 RAC(3,*),RBD(3,*)
      REAL*8 R(*)
      REAL*8 F7_0(0:486),F7_1(0:486),F7_2(0:486),F7_3(0:486)
      REAL*8 F8_0(0:486),F8_1(0:486),F8_2(0:486),F8_3(0:486)
      REAL*8 EX_0(0:486),EX_1(0:486),EX_2(0:486),EX_3(0:486)
      REAL*8 Half
      INTEGER Ind(2,*)
      PARAMETER (Half=0.5D0)
      INCLUDE "bkonx.inc"
      INCLUDE "Gamma_A.inc"
      INCLUDE "Gamma_7.inc"
      INCLUDE "Gamma_8.inc"
      INCLUDE "Exponential.inc"

      IF(LTot.EQ.7) THEN
         NBK =N*KBra*KKet
         NBK2=NBK*2
         NBK3=NBK*3
         NBK4=NBK*4
         NBK5=NBK*5
         NBK6=NBK*6
         NBK7=NBK*7
         I1=0
         DO IKK=1,KKet
            DO IKB=1,KBra
               DO I0=1,N
                  I=Ind(1,I0)
                  J=Ind(2,I0)
                  I1=I1+1
                  Zeta = PrmBufB(1,IKB,I)
                  Px   = PrmBufB(2,IKB,I)
                  Py   = PrmBufB(3,IKB,I)
                  Pz   = PrmBufB(4,IKB,I)
                  CBra(I0,IKB,1)=PrmBufB(6,IKB,I)
                  CBra(I0,IKB,2)=PrmBufB(7,IKB,I)
                  CBra(I0,IKB,3)=PrmBufB(8,IKB,I)
                  Eta  = PrmBufK(1,IKK,J)
                  Qx   = PrmBufK(2,IKK,J)
                  Qy   = PrmBufK(3,IKK,J)
                  Qz   = PrmBufK(4,IKK,J)
                  CKet(I0,IKK,1)=PrmBufK(6,IKK,J)
                  CKet(I0,IKK,2)=PrmBufK(7,IKK,J)
                  CKet(I0,IKK,3)=PrmBufK(8,IKK,J)
                  r1xZpE=1.0D0/(Zeta+Eta)
                  Rkk=PrmBufB(5,IKB,I)*PrmBufK(5,IKK,J)*DSQRT(r1xZpE)
                  t1=(Zeta*Px+Eta*Qx)*r1xZpE
                  t2=(Zeta*Py+Eta*Qy)*r1xZpE
                  t3=(Zeta*Pz+Eta*Qz)*r1xZpE
                  WR( 1,I1) = Px-Cx-SwhB*RAC(1,I)
                  WR( 2,I1) = Qx-Dx-SwhK*RBD(1,J)
                  WR( 3,I1) = Py-Cy-SwhB*RAC(2,I)
                  WR( 4,I1) = Qy-Dy-SwhK*RBD(2,J)
                  WR( 5,I1) = Pz-Cz-SwhB*RAC(3,I)
                  WR( 6,I1) = Qz-Dz-SwhK*RBD(3,J)
                  WR( 7,I1) = t1-Px
                  WR( 8,I1) = t1-Qx
                  WR( 9,I1) = t2-Py
                  WR(10,I1) = t2-Qy
                  WR(11,I1) = t3-Pz
                  WR(12,I1) = t3-Qz
                  WZ( 1,I1) = Half/Eta
                  WZ( 2,I1) = Half/Zeta
                  WZ( 3,I1) = Zeta*r1xZpE
                  WZ( 4,I1) = Eta*r1xZpE
                  WZ( 5,I1) = Half*r1xZpE
                  PQx = Px-Qx
                  PQy = Py-Qy
                  PQz = Pz-Qz
                  T   = (PQx*PQx+PQy*PQy+PQz*PQz)*Eta*Zeta*r1xZpE
                  IF(T.LT.Switch) THEN
                     m=INT(T*Grid)
                     T2=T*T
                     T3=T2*T
                     G7=F7_0(m)+T*F7_1(m)+T2*F7_2(m)+T3*F7_3(m)
                     ET=EX_0(m)+T*EX_1(m)+T2*EX_2(m)+T3*EX_3(m)
                     TwoT=2.0D0*T
                     G6=.769230769230769D-01*(TwoT*G7+ET)
                     G5=.909090909090909D-01*(TwoT*G6+ET)
                     G4=.111111111111111D+00*(TwoT*G5+ET)
                     G3=.142857142857143D+00*(TwoT*G4+ET)
                     G2=.200000000000000D+00*(TwoT*G3+ET)
                     G1=.333333333333333D+00*(TwoT*G2+ET)
                     G0=TwoT*G1+ET
                     R(I1     )=Rkk*G0
                     R(I1+NBK )=Rkk*G1
                     R(I1+NBK2)=Rkk*G2
                     R(I1+NBK3)=Rkk*G3
                     R(I1+NBK4)=Rkk*G4
                     R(I1+NBK5)=Rkk*G5
                     R(I1+NBK6)=Rkk*G6
                     R(I1+NBK7)=Rkk*G7
                  ELSE
                     T1=re1/T
                     TS=DSQRT(T1)
                     R(I1     )=Rkk*F0Asymp*TS
                     TS=TS*T1
                     R(I1+NBK )=Rkk*F1Asymp*TS
                     TS=TS*T1
                     R(I1+NBK2)=Rkk*F2Asymp*TS
                     TS=TS*T1
                     R(I1+NBK3)=Rkk*F3Asymp*TS
                     TS=TS*T1
                     R(I1+NBK4)=Rkk*F4Asymp*TS
                     TS=TS*T1
                     R(I1+NBK5)=Rkk*F5Asymp*TS
                     TS=TS*T1
                     R(I1+NBK6)=Rkk*F6Asymp*TS
                     TS=TS*T1
                     R(I1+NBK7)=Rkk*F7Asymp*TS
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ELSEIF(LTot.EQ.8) THEN
         NBK =N*KBra*KKet
         NBK2=NBK*2
         NBK3=NBK*3
         NBK4=NBK*4
         NBK5=NBK*5
         NBK6=NBK*6
         NBK7=NBK*7
         NBK8=NBK*8
         I1=0
         DO IKK=1,KKet
            DO IKB=1,KBra
               DO I0=1,N
                  I=Ind(1,I0)
                  J=Ind(2,I0)
                  I1=I1+1
                  Zeta = PrmBufB(1,IKB,I)
                  Px   = PrmBufB(2,IKB,I)
                  Py   = PrmBufB(3,IKB,I)
                  Pz   = PrmBufB(4,IKB,I)
                  CBra(I0,IKB,1)=PrmBufB(6,IKB,I)
                  CBra(I0,IKB,2)=PrmBufB(7,IKB,I)
                  CBra(I0,IKB,3)=PrmBufB(8,IKB,I)
                  Eta  = PrmBufK(1,IKK,J)
                  Qx   = PrmBufK(2,IKK,J)
                  Qy   = PrmBufK(3,IKK,J)
                  Qz   = PrmBufK(4,IKK,J)
                  CKet(I0,IKK,1)=PrmBufK(6,IKK,J)
                  CKet(I0,IKK,2)=PrmBufK(7,IKK,J)
                  CKet(I0,IKK,3)=PrmBufK(8,IKK,J)
                  r1xZpE=1.0D0/(Zeta+Eta)
                  Rkk=PrmBufB(5,IKB,I)*PrmBufK(5,IKK,J)*DSQRT(r1xZpE)
                  t1=(Zeta*Px+Eta*Qx)*r1xZpE
                  t2=(Zeta*Py+Eta*Qy)*r1xZpE
                  t3=(Zeta*Pz+Eta*Qz)*r1xZpE
                  WR( 1,I1) = Px-Cx-SwhB*RAC(1,I)
                  WR( 2,I1) = Qx-Dx-SwhK*RBD(1,J)
                  WR( 3,I1) = Py-Cy-SwhB*RAC(2,I)
                  WR( 4,I1) = Qy-Dy-SwhK*RBD(2,J)
                  WR( 5,I1) = Pz-Cz-SwhB*RAC(3,I)
                  WR( 6,I1) = Qz-Dz-SwhK*RBD(3,J)
                  WR( 7,I1) = t1-Px
                  WR( 8,I1) = t1-Qx
                  WR( 9,I1) = t2-Py
                  WR(10,I1) = t2-Qy
                  WR(11,I1) = t3-Pz
                  WR(12,I1) = t3-Qz
                  WZ( 1,I1) = Half/Eta
                  WZ( 2,I1) = Half/Zeta
                  WZ( 3,I1) = Zeta*r1xZpE
                  WZ( 4,I1) = Eta*r1xZpE
                  WZ( 5,I1) = Half*r1xZpE
                  PQx = Px-Qx
                  PQy = Py-Qy
                  PQz = Pz-Qz
                  T   = (PQx*PQx+PQy*PQy+PQz*PQz)*Eta*Zeta*r1xZpE
                  IF(T.LT.Switch) THEN
                     m=INT(T*Grid)
                     T2=T*T
                     T3=T2*T
                     G8=F8_0(m)+T*F8_1(m)+T2*F8_2(m)+T3*F8_3(m)
                     ET=EX_0(m)+T*EX_1(m)+T2*EX_2(m)+T3*EX_3(m)
                     TwoT=2.0D0*T
                     G7=.666666666666667D-01*(TwoT*G8+ET)
                     G6=.769230769230769D-01*(TwoT*G7+ET)
                     G5=.909090909090909D-01*(TwoT*G6+ET)
                     G4=.111111111111111D+00*(TwoT*G5+ET)
                     G3=.142857142857143D+00*(TwoT*G4+ET)
                     G2=.200000000000000D+00*(TwoT*G3+ET)
                     G1=.333333333333333D+00*(TwoT*G2+ET)
                     G0=TwoT*G1+ET
                     R(I1     )=Rkk*G0
                     R(I1+NBK )=Rkk*G1
                     R(I1+NBK2)=Rkk*G2
                     R(I1+NBK3)=Rkk*G3
                     R(I1+NBK4)=Rkk*G4
                     R(I1+NBK5)=Rkk*G5
                     R(I1+NBK6)=Rkk*G6
                     R(I1+NBK7)=Rkk*G7
                     R(I1+NBK8)=Rkk*G8
                  ELSE
                     T1=re1/T
                     TS=DSQRT(T1)
                     R(I1     )=Rkk*F0Asymp*TS
                     TS=TS*T1
                     R(I1+NBK )=Rkk*F1Asymp*TS
                     TS=TS*T1
                     R(I1+NBK2)=Rkk*F2Asymp*TS
                     TS=TS*T1
                     R(I1+NBK3)=Rkk*F3Asymp*TS
                     TS=TS*T1
                     R(I1+NBK4)=Rkk*F4Asymp*TS
                     TS=TS*T1
                     R(I1+NBK5)=Rkk*F5Asymp*TS
                     TS=TS*T1
                     R(I1+NBK6)=Rkk*F6Asymp*TS
                     TS=TS*T1
                     R(I1+NBK7)=Rkk*F7Asymp*TS
                     TS=TS*T1
                     R(I1+NBK8)=Rkk*F8Asymp*TS
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ELSE
         CALL MondoHalt('Ltot > 8 in TGen8')
      ENDIF
      RETURN
      END      

