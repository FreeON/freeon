      SUBROUTINE AuxIGen0(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K,IF0
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:0,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,TwoPi5x2
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "GammaGrid_77.Inc"
      INCLUDE "GammaDimensions_77.Inc"
      REAL*8 F0_0(0:Gamma_Mesh)
      REAL*8 F0_1(0:Gamma_Mesh)
      REAL*8 F0_2(0:Gamma_Mesh)
      REAL*8 F0_3(0:Gamma_Mesh)
      REAL*8 F0_4(0:Gamma_Mesh)
      REAL*8 F0_5(0:Gamma_Mesh)
      REAL*8 F0_6(0:Gamma_Mesh)
      INCLUDE "F0_77.Inc"
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=TwoPi5x2/(RTE*SQRT(RPE))
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         IF(T.LT.1D0)THEN
            J=AINT(T*Gamma_Grid)
            G0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >        +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))
         ELSEIF(T.LT.Gamma_Switch)THEN
            J=AINT(T*Gamma_Grid)
            G0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >        +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))
         ELSE
            SqrtT=DSQRT(T)
            OneOvT=One/T
            G0=SqrtPi/(Two*SqrtT)
         ENDIF
         o1=Upq
         AuxR(0,I)=o1*G0
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen1(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K,IF0,IF1
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:1,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,TwoPi5x2
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "GammaGrid_77.Inc"
      INCLUDE "GammaDimensions_77.Inc"
      REAL*8 F0_0(0:Gamma_Mesh)
      REAL*8 F0_1(0:Gamma_Mesh)
      REAL*8 F0_2(0:Gamma_Mesh)
      REAL*8 F0_3(0:Gamma_Mesh)
      REAL*8 F0_4(0:Gamma_Mesh)
      REAL*8 F0_5(0:Gamma_Mesh)
      REAL*8 F0_6(0:Gamma_Mesh)
      REAL*8 F1_0(0:Gamma_Mesh)
      REAL*8 F1_1(0:Gamma_Mesh)
      REAL*8 F1_2(0:Gamma_Mesh)
      REAL*8 F1_3(0:Gamma_Mesh)
      REAL*8 F1_4(0:Gamma_Mesh)
      REAL*8 F1_5(0:Gamma_Mesh)
      REAL*8 F1_6(0:Gamma_Mesh)
      INCLUDE "F0_77.Inc"
      INCLUDE "F1_77.Inc"
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=TwoPi5x2/(RTE*SQRT(RPE))
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         IF(T.LT.1D0)THEN
            J=AINT(T*Gamma_Grid)
            G1=F1_0(J)+T*(F1_1(J)+T*(F1_2(J)+T*(F1_3(J)
     >        +T*(F1_4(J)+T*(F1_5(J)+T*F1_6(J))))))
            ET=DEXP(-T)
            TwoT=2.0D0*T
            G0=TwoT*G1+ET
         ELSEIF(T.LT.Gamma_Switch)THEN
            J=AINT(T*Gamma_Grid)
            G0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >        +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))
            ET=DEXP(-T)
            OneOvr2T=5D-1/T
            G1=OneOvr2T*(0.1000000000000000D+01*G0-ET)
         ELSE
            SqrtT=DSQRT(T)
            OneOvT=One/T
            G0=SqrtPi/(Two*SqrtT)
            G1=G0*0.5000000000000000D+00*OneOvT
         ENDIF
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*G0
         o1=o2*o1
         AuxR(1,I)=o1*G1
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen2(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K,IF0,IF2
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:2,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,TwoPi5x2
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "GammaGrid_77.Inc"
      INCLUDE "GammaDimensions_77.Inc"
      REAL*8 F0_0(0:Gamma_Mesh)
      REAL*8 F0_1(0:Gamma_Mesh)
      REAL*8 F0_2(0:Gamma_Mesh)
      REAL*8 F0_3(0:Gamma_Mesh)
      REAL*8 F0_4(0:Gamma_Mesh)
      REAL*8 F0_5(0:Gamma_Mesh)
      REAL*8 F0_6(0:Gamma_Mesh)
      REAL*8 F2_0(0:Gamma_Mesh)
      REAL*8 F2_1(0:Gamma_Mesh)
      REAL*8 F2_2(0:Gamma_Mesh)
      REAL*8 F2_3(0:Gamma_Mesh)
      REAL*8 F2_4(0:Gamma_Mesh)
      REAL*8 F2_5(0:Gamma_Mesh)
      REAL*8 F2_6(0:Gamma_Mesh)
      INCLUDE "F0_77.Inc"
      INCLUDE "F2_77.Inc"
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=TwoPi5x2/(RTE*SQRT(RPE))
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         IF(T.LT.1D0)THEN
            J=AINT(T*Gamma_Grid)
            G2=F2_0(J)+T*(F2_1(J)+T*(F2_2(J)+T*(F2_3(J)
     >        +T*(F2_4(J)+T*(F2_5(J)+T*F2_6(J))))))
            ET=DEXP(-T)
            TwoT=2.0D0*T
            G1=0.3333333333333333D+00*(TwoT*G2+ET)
            G0=TwoT*G1+ET
         ELSEIF(T.LT.Gamma_Switch)THEN
            J=AINT(T*Gamma_Grid)
            G0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >        +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))
            ET=DEXP(-T)
            OneOvr2T=5D-1/T
            G1=OneOvr2T*(0.1000000000000000D+01*G0-ET)
            G2=OneOvr2T*(0.3000000000000000D+01*G1-ET)
         ELSE
            SqrtT=DSQRT(T)
            OneOvT=One/T
            G0=SqrtPi/(Two*SqrtT)
            G1=G0*0.5000000000000000D+00*OneOvT
            G2=G1*0.1500000000000000D+01*OneOvT
         ENDIF
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*G0
         o1=o2*o1
         AuxR(1,I)=o1*G1
         o1=o2*o1
         AuxR(2,I)=o1*G2
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen3(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K,IF0,IF3
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:3,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,TwoPi5x2
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "GammaGrid_77.Inc"
      INCLUDE "GammaDimensions_77.Inc"
      REAL*8 F0_0(0:Gamma_Mesh)
      REAL*8 F0_1(0:Gamma_Mesh)
      REAL*8 F0_2(0:Gamma_Mesh)
      REAL*8 F0_3(0:Gamma_Mesh)
      REAL*8 F0_4(0:Gamma_Mesh)
      REAL*8 F0_5(0:Gamma_Mesh)
      REAL*8 F0_6(0:Gamma_Mesh)
      REAL*8 F3_0(0:Gamma_Mesh)
      REAL*8 F3_1(0:Gamma_Mesh)
      REAL*8 F3_2(0:Gamma_Mesh)
      REAL*8 F3_3(0:Gamma_Mesh)
      REAL*8 F3_4(0:Gamma_Mesh)
      REAL*8 F3_5(0:Gamma_Mesh)
      REAL*8 F3_6(0:Gamma_Mesh)
      INCLUDE "F0_77.Inc"
      INCLUDE "F3_77.Inc"
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=TwoPi5x2/(RTE*SQRT(RPE))
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         IF(T.LT.1D0)THEN
            J=AINT(T*Gamma_Grid)
            G3=F3_0(J)+T*(F3_1(J)+T*(F3_2(J)+T*(F3_3(J)
     >        +T*(F3_4(J)+T*(F3_5(J)+T*F3_6(J))))))
            ET=DEXP(-T)
            TwoT=2.0D0*T
            G2=0.2000000000000000D+00*(TwoT*G3+ET)
            G1=0.3333333333333333D+00*(TwoT*G2+ET)
            G0=TwoT*G1+ET
         ELSEIF(T.LT.Gamma_Switch)THEN
            J=AINT(T*Gamma_Grid)
            G0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >        +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))
            ET=DEXP(-T)
            OneOvr2T=5D-1/T
            G1=OneOvr2T*(0.1000000000000000D+01*G0-ET)
            G2=OneOvr2T*(0.3000000000000000D+01*G1-ET)
            G3=OneOvr2T*(0.5000000000000000D+01*G2-ET)
         ELSE
            SqrtT=DSQRT(T)
            OneOvT=One/T
            G0=SqrtPi/(Two*SqrtT)
            G1=G0*0.5000000000000000D+00*OneOvT
            G2=G1*0.1500000000000000D+01*OneOvT
            G3=G2*0.2500000000000000D+01*OneOvT
         ENDIF
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*G0
         o1=o2*o1
         AuxR(1,I)=o1*G1
         o1=o2*o1
         AuxR(2,I)=o1*G2
         o1=o2*o1
         AuxR(3,I)=o1*G3
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen4(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K,IF0,IF4
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:4,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,TwoPi5x2
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "GammaGrid_77.Inc"
      INCLUDE "GammaDimensions_77.Inc"
      REAL*8 F0_0(0:Gamma_Mesh)
      REAL*8 F0_1(0:Gamma_Mesh)
      REAL*8 F0_2(0:Gamma_Mesh)
      REAL*8 F0_3(0:Gamma_Mesh)
      REAL*8 F0_4(0:Gamma_Mesh)
      REAL*8 F0_5(0:Gamma_Mesh)
      REAL*8 F0_6(0:Gamma_Mesh)
      REAL*8 F4_0(0:Gamma_Mesh)
      REAL*8 F4_1(0:Gamma_Mesh)
      REAL*8 F4_2(0:Gamma_Mesh)
      REAL*8 F4_3(0:Gamma_Mesh)
      REAL*8 F4_4(0:Gamma_Mesh)
      REAL*8 F4_5(0:Gamma_Mesh)
      REAL*8 F4_6(0:Gamma_Mesh)
      INCLUDE "F0_77.Inc"
      INCLUDE "F4_77.Inc"
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=TwoPi5x2/(RTE*SQRT(RPE))
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         IF(T.LT.1D0)THEN
            J=AINT(T*Gamma_Grid)
            G4=F4_0(J)+T*(F4_1(J)+T*(F4_2(J)+T*(F4_3(J)
     >        +T*(F4_4(J)+T*(F4_5(J)+T*F4_6(J))))))
            ET=DEXP(-T)
            TwoT=2.0D0*T
            G3=0.1428571428571428D+00*(TwoT*G4+ET)
            G2=0.2000000000000000D+00*(TwoT*G3+ET)
            G1=0.3333333333333333D+00*(TwoT*G2+ET)
            G0=TwoT*G1+ET
         ELSEIF(T.LT.Gamma_Switch)THEN
            J=AINT(T*Gamma_Grid)
            G0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >        +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))
            ET=DEXP(-T)
            OneOvr2T=5D-1/T
            G1=OneOvr2T*(0.1000000000000000D+01*G0-ET)
            G2=OneOvr2T*(0.3000000000000000D+01*G1-ET)
            G3=OneOvr2T*(0.5000000000000000D+01*G2-ET)
            G4=OneOvr2T*(0.7000000000000000D+01*G3-ET)
         ELSE
            SqrtT=DSQRT(T)
            OneOvT=One/T
            G0=SqrtPi/(Two*SqrtT)
            G1=G0*0.5000000000000000D+00*OneOvT
            G2=G1*0.1500000000000000D+01*OneOvT
            G3=G2*0.2500000000000000D+01*OneOvT
            G4=G3*0.3500000000000000D+01*OneOvT
         ENDIF
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*G0
         o1=o2*o1
         AuxR(1,I)=o1*G1
         o1=o2*o1
         AuxR(2,I)=o1*G2
         o1=o2*o1
         AuxR(3,I)=o1*G3
         o1=o2*o1
         AuxR(4,I)=o1*G4
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen5(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K,IF0,IF5
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:5,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,TwoPi5x2
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "GammaGrid_77.Inc"
      INCLUDE "GammaDimensions_77.Inc"
      REAL*8 F0_0(0:Gamma_Mesh)
      REAL*8 F0_1(0:Gamma_Mesh)
      REAL*8 F0_2(0:Gamma_Mesh)
      REAL*8 F0_3(0:Gamma_Mesh)
      REAL*8 F0_4(0:Gamma_Mesh)
      REAL*8 F0_5(0:Gamma_Mesh)
      REAL*8 F0_6(0:Gamma_Mesh)
      REAL*8 F5_0(0:Gamma_Mesh)
      REAL*8 F5_1(0:Gamma_Mesh)
      REAL*8 F5_2(0:Gamma_Mesh)
      REAL*8 F5_3(0:Gamma_Mesh)
      REAL*8 F5_4(0:Gamma_Mesh)
      REAL*8 F5_5(0:Gamma_Mesh)
      REAL*8 F5_6(0:Gamma_Mesh)
      INCLUDE "F0_77.Inc"
      INCLUDE "F5_77.Inc"
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=TwoPi5x2/(RTE*SQRT(RPE))
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         IF(T.LT.1D0)THEN
            J=AINT(T*Gamma_Grid)
            G5=F5_0(J)+T*(F5_1(J)+T*(F5_2(J)+T*(F5_3(J)
     >        +T*(F5_4(J)+T*(F5_5(J)+T*F5_6(J))))))
            ET=DEXP(-T)
            TwoT=2.0D0*T
            G4=0.1111111111111111D+00*(TwoT*G5+ET)
            G3=0.1428571428571428D+00*(TwoT*G4+ET)
            G2=0.2000000000000000D+00*(TwoT*G3+ET)
            G1=0.3333333333333333D+00*(TwoT*G2+ET)
            G0=TwoT*G1+ET
         ELSEIF(T.LT.Gamma_Switch)THEN
            J=AINT(T*Gamma_Grid)
            G0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >        +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))
            ET=DEXP(-T)
            OneOvr2T=5D-1/T
            G1=OneOvr2T*(0.1000000000000000D+01*G0-ET)
            G2=OneOvr2T*(0.3000000000000000D+01*G1-ET)
            G3=OneOvr2T*(0.5000000000000000D+01*G2-ET)
            G4=OneOvr2T*(0.7000000000000000D+01*G3-ET)
            G5=OneOvr2T*(0.9000000000000000D+01*G4-ET)
         ELSE
            SqrtT=DSQRT(T)
            OneOvT=One/T
            G0=SqrtPi/(Two*SqrtT)
            G1=G0*0.5000000000000000D+00*OneOvT
            G2=G1*0.1500000000000000D+01*OneOvT
            G3=G2*0.2500000000000000D+01*OneOvT
            G4=G3*0.3500000000000000D+01*OneOvT
            G5=G4*0.4500000000000000D+01*OneOvT
         ENDIF
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*G0
         o1=o2*o1
         AuxR(1,I)=o1*G1
         o1=o2*o1
         AuxR(2,I)=o1*G2
         o1=o2*o1
         AuxR(3,I)=o1*G3
         o1=o2*o1
         AuxR(4,I)=o1*G4
         o1=o2*o1
         AuxR(5,I)=o1*G5
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen6(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K,IF0,IF6
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:6,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,TwoPi5x2
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "GammaGrid_77.Inc"
      INCLUDE "GammaDimensions_77.Inc"
      REAL*8 F0_0(0:Gamma_Mesh)
      REAL*8 F0_1(0:Gamma_Mesh)
      REAL*8 F0_2(0:Gamma_Mesh)
      REAL*8 F0_3(0:Gamma_Mesh)
      REAL*8 F0_4(0:Gamma_Mesh)
      REAL*8 F0_5(0:Gamma_Mesh)
      REAL*8 F0_6(0:Gamma_Mesh)
      REAL*8 F6_0(0:Gamma_Mesh)
      REAL*8 F6_1(0:Gamma_Mesh)
      REAL*8 F6_2(0:Gamma_Mesh)
      REAL*8 F6_3(0:Gamma_Mesh)
      REAL*8 F6_4(0:Gamma_Mesh)
      REAL*8 F6_5(0:Gamma_Mesh)
      REAL*8 F6_6(0:Gamma_Mesh)
      INCLUDE "F0_77.Inc"
      INCLUDE "F6_77.Inc"
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=TwoPi5x2/(RTE*SQRT(RPE))
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         IF(T.LT.1D0)THEN
            J=AINT(T*Gamma_Grid)
            G6=F6_0(J)+T*(F6_1(J)+T*(F6_2(J)+T*(F6_3(J)
     >        +T*(F6_4(J)+T*(F6_5(J)+T*F6_6(J))))))
            ET=DEXP(-T)
            TwoT=2.0D0*T
            G5=0.9090909090909091D-01*(TwoT*G6+ET)
            G4=0.1111111111111111D+00*(TwoT*G5+ET)
            G3=0.1428571428571428D+00*(TwoT*G4+ET)
            G2=0.2000000000000000D+00*(TwoT*G3+ET)
            G1=0.3333333333333333D+00*(TwoT*G2+ET)
            G0=TwoT*G1+ET
         ELSEIF(T.LT.Gamma_Switch)THEN
            J=AINT(T*Gamma_Grid)
            G0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >        +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))
            ET=DEXP(-T)
            OneOvr2T=5D-1/T
            G1=OneOvr2T*(0.1000000000000000D+01*G0-ET)
            G2=OneOvr2T*(0.3000000000000000D+01*G1-ET)
            G3=OneOvr2T*(0.5000000000000000D+01*G2-ET)
            G4=OneOvr2T*(0.7000000000000000D+01*G3-ET)
            G5=OneOvr2T*(0.9000000000000000D+01*G4-ET)
            G6=OneOvr2T*(0.1100000000000000D+02*G5-ET)
         ELSE
            SqrtT=DSQRT(T)
            OneOvT=One/T
            G0=SqrtPi/(Two*SqrtT)
            G1=G0*0.5000000000000000D+00*OneOvT
            G2=G1*0.1500000000000000D+01*OneOvT
            G3=G2*0.2500000000000000D+01*OneOvT
            G4=G3*0.3500000000000000D+01*OneOvT
            G5=G4*0.4500000000000000D+01*OneOvT
            G6=G5*0.5500000000000000D+01*OneOvT
         ENDIF
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*G0
         o1=o2*o1
         AuxR(1,I)=o1*G1
         o1=o2*o1
         AuxR(2,I)=o1*G2
         o1=o2*o1
         AuxR(3,I)=o1*G3
         o1=o2*o1
         AuxR(4,I)=o1*G4
         o1=o2*o1
         AuxR(5,I)=o1*G5
         o1=o2*o1
         AuxR(6,I)=o1*G6
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen7(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K,IF0,IF7
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:7,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,TwoPi5x2
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "GammaGrid_77.Inc"
      INCLUDE "GammaDimensions_77.Inc"
      REAL*8 F0_0(0:Gamma_Mesh)
      REAL*8 F0_1(0:Gamma_Mesh)
      REAL*8 F0_2(0:Gamma_Mesh)
      REAL*8 F0_3(0:Gamma_Mesh)
      REAL*8 F0_4(0:Gamma_Mesh)
      REAL*8 F0_5(0:Gamma_Mesh)
      REAL*8 F0_6(0:Gamma_Mesh)
      REAL*8 F7_0(0:Gamma_Mesh)
      REAL*8 F7_1(0:Gamma_Mesh)
      REAL*8 F7_2(0:Gamma_Mesh)
      REAL*8 F7_3(0:Gamma_Mesh)
      REAL*8 F7_4(0:Gamma_Mesh)
      REAL*8 F7_5(0:Gamma_Mesh)
      REAL*8 F7_6(0:Gamma_Mesh)
      INCLUDE "F0_77.Inc"
      INCLUDE "F7_77.Inc"
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=TwoPi5x2/(RTE*SQRT(RPE))
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         IF(T.LT.1D0)THEN
            J=AINT(T*Gamma_Grid)
            G7=F7_0(J)+T*(F7_1(J)+T*(F7_2(J)+T*(F7_3(J)
     >        +T*(F7_4(J)+T*(F7_5(J)+T*F7_6(J))))))
            ET=DEXP(-T)
            TwoT=2.0D0*T
            G6=0.7692307692307693D-01*(TwoT*G7+ET)
            G5=0.9090909090909091D-01*(TwoT*G6+ET)
            G4=0.1111111111111111D+00*(TwoT*G5+ET)
            G3=0.1428571428571428D+00*(TwoT*G4+ET)
            G2=0.2000000000000000D+00*(TwoT*G3+ET)
            G1=0.3333333333333333D+00*(TwoT*G2+ET)
            G0=TwoT*G1+ET
         ELSEIF(T.LT.Gamma_Switch)THEN
            J=AINT(T*Gamma_Grid)
            G0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >        +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))
            ET=DEXP(-T)
            OneOvr2T=5D-1/T
            G1=OneOvr2T*(0.1000000000000000D+01*G0-ET)
            G2=OneOvr2T*(0.3000000000000000D+01*G1-ET)
            G3=OneOvr2T*(0.5000000000000000D+01*G2-ET)
            G4=OneOvr2T*(0.7000000000000000D+01*G3-ET)
            G5=OneOvr2T*(0.9000000000000000D+01*G4-ET)
            G6=OneOvr2T*(0.1100000000000000D+02*G5-ET)
            G7=OneOvr2T*(0.1300000000000000D+02*G6-ET)
         ELSE
            SqrtT=DSQRT(T)
            OneOvT=One/T
            G0=SqrtPi/(Two*SqrtT)
            G1=G0*0.5000000000000000D+00*OneOvT
            G2=G1*0.1500000000000000D+01*OneOvT
            G3=G2*0.2500000000000000D+01*OneOvT
            G4=G3*0.3500000000000000D+01*OneOvT
            G5=G4*0.4500000000000000D+01*OneOvT
            G6=G5*0.5500000000000000D+01*OneOvT
            G7=G6*0.6500000000000000D+01*OneOvT
         ENDIF
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*G0
         o1=o2*o1
         AuxR(1,I)=o1*G1
         o1=o2*o1
         AuxR(2,I)=o1*G2
         o1=o2*o1
         AuxR(3,I)=o1*G3
         o1=o2*o1
         AuxR(4,I)=o1*G4
         o1=o2*o1
         AuxR(5,I)=o1*G5
         o1=o2*o1
         AuxR(6,I)=o1*G6
         o1=o2*o1
         AuxR(7,I)=o1*G7
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen8(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K,IF0,IF8
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:8,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,TwoPi5x2
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "GammaGrid_77.Inc"
      INCLUDE "GammaDimensions_77.Inc"
      REAL*8 F0_0(0:Gamma_Mesh)
      REAL*8 F0_1(0:Gamma_Mesh)
      REAL*8 F0_2(0:Gamma_Mesh)
      REAL*8 F0_3(0:Gamma_Mesh)
      REAL*8 F0_4(0:Gamma_Mesh)
      REAL*8 F0_5(0:Gamma_Mesh)
      REAL*8 F0_6(0:Gamma_Mesh)
      REAL*8 F8_0(0:Gamma_Mesh)
      REAL*8 F8_1(0:Gamma_Mesh)
      REAL*8 F8_2(0:Gamma_Mesh)
      REAL*8 F8_3(0:Gamma_Mesh)
      REAL*8 F8_4(0:Gamma_Mesh)
      REAL*8 F8_5(0:Gamma_Mesh)
      REAL*8 F8_6(0:Gamma_Mesh)
      INCLUDE "F0_77.Inc"
      INCLUDE "F8_77.Inc"
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=TwoPi5x2/(RTE*SQRT(RPE))
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         IF(T.LT.1D0)THEN
            J=AINT(T*Gamma_Grid)
            G8=F8_0(J)+T*(F8_1(J)+T*(F8_2(J)+T*(F8_3(J)
     >        +T*(F8_4(J)+T*(F8_5(J)+T*F8_6(J))))))
            ET=DEXP(-T)
            TwoT=2.0D0*T
            G7=0.6666666666666667D-01*(TwoT*G8+ET)
            G6=0.7692307692307693D-01*(TwoT*G7+ET)
            G5=0.9090909090909091D-01*(TwoT*G6+ET)
            G4=0.1111111111111111D+00*(TwoT*G5+ET)
            G3=0.1428571428571428D+00*(TwoT*G4+ET)
            G2=0.2000000000000000D+00*(TwoT*G3+ET)
            G1=0.3333333333333333D+00*(TwoT*G2+ET)
            G0=TwoT*G1+ET
         ELSEIF(T.LT.Gamma_Switch)THEN
            J=AINT(T*Gamma_Grid)
            G0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >        +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))
            ET=DEXP(-T)
            OneOvr2T=5D-1/T
            G1=OneOvr2T*(0.1000000000000000D+01*G0-ET)
            G2=OneOvr2T*(0.3000000000000000D+01*G1-ET)
            G3=OneOvr2T*(0.5000000000000000D+01*G2-ET)
            G4=OneOvr2T*(0.7000000000000000D+01*G3-ET)
            G5=OneOvr2T*(0.9000000000000000D+01*G4-ET)
            G6=OneOvr2T*(0.1100000000000000D+02*G5-ET)
            G7=OneOvr2T*(0.1300000000000000D+02*G6-ET)
            G8=OneOvr2T*(0.1500000000000000D+02*G7-ET)
         ELSE
            SqrtT=DSQRT(T)
            OneOvT=One/T
            G0=SqrtPi/(Two*SqrtT)
            G1=G0*0.5000000000000000D+00*OneOvT
            G2=G1*0.1500000000000000D+01*OneOvT
            G3=G2*0.2500000000000000D+01*OneOvT
            G4=G3*0.3500000000000000D+01*OneOvT
            G5=G4*0.4500000000000000D+01*OneOvT
            G6=G5*0.5500000000000000D+01*OneOvT
            G7=G6*0.6500000000000000D+01*OneOvT
            G8=G7*0.7500000000000000D+01*OneOvT
         ENDIF
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*G0
         o1=o2*o1
         AuxR(1,I)=o1*G1
         o1=o2*o1
         AuxR(2,I)=o1*G2
         o1=o2*o1
         AuxR(3,I)=o1*G3
         o1=o2*o1
         AuxR(4,I)=o1*G4
         o1=o2*o1
         AuxR(5,I)=o1*G5
         o1=o2*o1
         AuxR(6,I)=o1*G6
         o1=o2*o1
         AuxR(7,I)=o1*G7
         o1=o2*o1
         AuxR(8,I)=o1*G8
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen9(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K,IF0,IF9
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:9,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,TwoPi5x2
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "GammaGrid_77.Inc"
      INCLUDE "GammaDimensions_77.Inc"
      REAL*8 F0_0(0:Gamma_Mesh)
      REAL*8 F0_1(0:Gamma_Mesh)
      REAL*8 F0_2(0:Gamma_Mesh)
      REAL*8 F0_3(0:Gamma_Mesh)
      REAL*8 F0_4(0:Gamma_Mesh)
      REAL*8 F0_5(0:Gamma_Mesh)
      REAL*8 F0_6(0:Gamma_Mesh)
      REAL*8 F9_0(0:Gamma_Mesh)
      REAL*8 F9_1(0:Gamma_Mesh)
      REAL*8 F9_2(0:Gamma_Mesh)
      REAL*8 F9_3(0:Gamma_Mesh)
      REAL*8 F9_4(0:Gamma_Mesh)
      REAL*8 F9_5(0:Gamma_Mesh)
      REAL*8 F9_6(0:Gamma_Mesh)
      INCLUDE "F0_77.Inc"
      INCLUDE "F9_77.Inc"
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=TwoPi5x2/(RTE*SQRT(RPE))
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         IF(T.LT.1D0)THEN
            J=AINT(T*Gamma_Grid)
            G9=F9_0(J)+T*(F9_1(J)+T*(F9_2(J)+T*(F9_3(J)
     >        +T*(F9_4(J)+T*(F9_5(J)+T*F9_6(J))))))
            ET=DEXP(-T)
            TwoT=2.0D0*T
            G8=0.5882352941176471D-01*(TwoT*G9+ET)
            G7=0.6666666666666667D-01*(TwoT*G8+ET)
            G6=0.7692307692307693D-01*(TwoT*G7+ET)
            G5=0.9090909090909091D-01*(TwoT*G6+ET)
            G4=0.1111111111111111D+00*(TwoT*G5+ET)
            G3=0.1428571428571428D+00*(TwoT*G4+ET)
            G2=0.2000000000000000D+00*(TwoT*G3+ET)
            G1=0.3333333333333333D+00*(TwoT*G2+ET)
            G0=TwoT*G1+ET
         ELSEIF(T.LT.Gamma_Switch)THEN
            J=AINT(T*Gamma_Grid)
            G0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >        +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))
            ET=DEXP(-T)
            OneOvr2T=5D-1/T
            G1=OneOvr2T*(0.1000000000000000D+01*G0-ET)
            G2=OneOvr2T*(0.3000000000000000D+01*G1-ET)
            G3=OneOvr2T*(0.5000000000000000D+01*G2-ET)
            G4=OneOvr2T*(0.7000000000000000D+01*G3-ET)
            G5=OneOvr2T*(0.9000000000000000D+01*G4-ET)
            G6=OneOvr2T*(0.1100000000000000D+02*G5-ET)
            G7=OneOvr2T*(0.1300000000000000D+02*G6-ET)
            G8=OneOvr2T*(0.1500000000000000D+02*G7-ET)
            G9=OneOvr2T*(0.1700000000000000D+02*G8-ET)
         ELSE
            SqrtT=DSQRT(T)
            OneOvT=One/T
            G0=SqrtPi/(Two*SqrtT)
            G1=G0*0.5000000000000000D+00*OneOvT
            G2=G1*0.1500000000000000D+01*OneOvT
            G3=G2*0.2500000000000000D+01*OneOvT
            G4=G3*0.3500000000000000D+01*OneOvT
            G5=G4*0.4500000000000000D+01*OneOvT
            G6=G5*0.5500000000000000D+01*OneOvT
            G7=G6*0.6500000000000000D+01*OneOvT
            G8=G7*0.7500000000000000D+01*OneOvT
            G9=G8*0.8500000000000000D+01*OneOvT
         ENDIF
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*G0
         o1=o2*o1
         AuxR(1,I)=o1*G1
         o1=o2*o1
         AuxR(2,I)=o1*G2
         o1=o2*o1
         AuxR(3,I)=o1*G3
         o1=o2*o1
         AuxR(4,I)=o1*G4
         o1=o2*o1
         AuxR(5,I)=o1*G5
         o1=o2*o1
         AuxR(6,I)=o1*G6
         o1=o2*o1
         AuxR(7,I)=o1*G7
         o1=o2*o1
         AuxR(8,I)=o1*G8
         o1=o2*o1
         AuxR(9,I)=o1*G9
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen10(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K,IF0,IF10
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:10,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,TwoPi5x2
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "GammaGrid_77.Inc"
      INCLUDE "GammaDimensions_77.Inc"
      REAL*8 F0_0(0:Gamma_Mesh)
      REAL*8 F0_1(0:Gamma_Mesh)
      REAL*8 F0_2(0:Gamma_Mesh)
      REAL*8 F0_3(0:Gamma_Mesh)
      REAL*8 F0_4(0:Gamma_Mesh)
      REAL*8 F0_5(0:Gamma_Mesh)
      REAL*8 F0_6(0:Gamma_Mesh)
      REAL*8 F10_0(0:Gamma_Mesh)
      REAL*8 F10_1(0:Gamma_Mesh)
      REAL*8 F10_2(0:Gamma_Mesh)
      REAL*8 F10_3(0:Gamma_Mesh)
      REAL*8 F10_4(0:Gamma_Mesh)
      REAL*8 F10_5(0:Gamma_Mesh)
      REAL*8 F10_6(0:Gamma_Mesh)
      INCLUDE "F0_77.Inc"
      INCLUDE "F10_77.Inc"
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=TwoPi5x2/(RTE*SQRT(RPE))
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         IF(T.LT.1D0)THEN
            J=AINT(T*Gamma_Grid)
            G10=F10_0(J)+T*(F10_1(J)+T*(F10_2(J)+T*(F10_3(J)
     >        +T*(F10_4(J)+T*(F10_5(J)+T*F10_6(J))))))
            ET=DEXP(-T)
            TwoT=2.0D0*T
            G9=0.5263157894736842D-01*(TwoT*G10+ET)
            G8=0.5882352941176471D-01*(TwoT*G9+ET)
            G7=0.6666666666666667D-01*(TwoT*G8+ET)
            G6=0.7692307692307693D-01*(TwoT*G7+ET)
            G5=0.9090909090909091D-01*(TwoT*G6+ET)
            G4=0.1111111111111111D+00*(TwoT*G5+ET)
            G3=0.1428571428571428D+00*(TwoT*G4+ET)
            G2=0.2000000000000000D+00*(TwoT*G3+ET)
            G1=0.3333333333333333D+00*(TwoT*G2+ET)
            G0=TwoT*G1+ET
         ELSEIF(T.LT.Gamma_Switch)THEN
            J=AINT(T*Gamma_Grid)
            G0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >        +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))
            ET=DEXP(-T)
            OneOvr2T=5D-1/T
            G1=OneOvr2T*(0.1000000000000000D+01*G0-ET)
            G2=OneOvr2T*(0.3000000000000000D+01*G1-ET)
            G3=OneOvr2T*(0.5000000000000000D+01*G2-ET)
            G4=OneOvr2T*(0.7000000000000000D+01*G3-ET)
            G5=OneOvr2T*(0.9000000000000000D+01*G4-ET)
            G6=OneOvr2T*(0.1100000000000000D+02*G5-ET)
            G7=OneOvr2T*(0.1300000000000000D+02*G6-ET)
            G8=OneOvr2T*(0.1500000000000000D+02*G7-ET)
            G9=OneOvr2T*(0.1700000000000000D+02*G8-ET)
            G10=OneOvr2T*(0.1900000000000000D+02*G9-ET)
         ELSE
            SqrtT=DSQRT(T)
            OneOvT=One/T
            G0=SqrtPi/(Two*SqrtT)
            G1=G0*0.5000000000000000D+00*OneOvT
            G2=G1*0.1500000000000000D+01*OneOvT
            G3=G2*0.2500000000000000D+01*OneOvT
            G4=G3*0.3500000000000000D+01*OneOvT
            G5=G4*0.4500000000000000D+01*OneOvT
            G6=G5*0.5500000000000000D+01*OneOvT
            G7=G6*0.6500000000000000D+01*OneOvT
            G8=G7*0.7500000000000000D+01*OneOvT
            G9=G8*0.8500000000000000D+01*OneOvT
            G10=G9*0.9500000000000000D+01*OneOvT
         ENDIF
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*G0
         o1=o2*o1
         AuxR(1,I)=o1*G1
         o1=o2*o1
         AuxR(2,I)=o1*G2
         o1=o2*o1
         AuxR(3,I)=o1*G3
         o1=o2*o1
         AuxR(4,I)=o1*G4
         o1=o2*o1
         AuxR(5,I)=o1*G5
         o1=o2*o1
         AuxR(6,I)=o1*G6
         o1=o2*o1
         AuxR(7,I)=o1*G7
         o1=o2*o1
         AuxR(8,I)=o1*G8
         o1=o2*o1
         AuxR(9,I)=o1*G9
         o1=o2*o1
         AuxR(10,I)=o1*G10
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen11(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K,IF0,IF11
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:11,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,TwoPi5x2
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "GammaGrid_77.Inc"
      INCLUDE "GammaDimensions_77.Inc"
      REAL*8 F0_0(0:Gamma_Mesh)
      REAL*8 F0_1(0:Gamma_Mesh)
      REAL*8 F0_2(0:Gamma_Mesh)
      REAL*8 F0_3(0:Gamma_Mesh)
      REAL*8 F0_4(0:Gamma_Mesh)
      REAL*8 F0_5(0:Gamma_Mesh)
      REAL*8 F0_6(0:Gamma_Mesh)
      REAL*8 F11_0(0:Gamma_Mesh)
      REAL*8 F11_1(0:Gamma_Mesh)
      REAL*8 F11_2(0:Gamma_Mesh)
      REAL*8 F11_3(0:Gamma_Mesh)
      REAL*8 F11_4(0:Gamma_Mesh)
      REAL*8 F11_5(0:Gamma_Mesh)
      REAL*8 F11_6(0:Gamma_Mesh)
      INCLUDE "F0_77.Inc"
      INCLUDE "F11_77.Inc"
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=TwoPi5x2/(RTE*SQRT(RPE))
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         IF(T.LT.1D0)THEN
            J=AINT(T*Gamma_Grid)
            G11=F11_0(J)+T*(F11_1(J)+T*(F11_2(J)+T*(F11_3(J)
     >        +T*(F11_4(J)+T*(F11_5(J)+T*F11_6(J))))))
            ET=DEXP(-T)
            TwoT=2.0D0*T
            G10=0.4761904761904762D-01*(TwoT*G11+ET)
            G9=0.5263157894736842D-01*(TwoT*G10+ET)
            G8=0.5882352941176471D-01*(TwoT*G9+ET)
            G7=0.6666666666666667D-01*(TwoT*G8+ET)
            G6=0.7692307692307693D-01*(TwoT*G7+ET)
            G5=0.9090909090909091D-01*(TwoT*G6+ET)
            G4=0.1111111111111111D+00*(TwoT*G5+ET)
            G3=0.1428571428571428D+00*(TwoT*G4+ET)
            G2=0.2000000000000000D+00*(TwoT*G3+ET)
            G1=0.3333333333333333D+00*(TwoT*G2+ET)
            G0=TwoT*G1+ET
         ELSEIF(T.LT.Gamma_Switch)THEN
            J=AINT(T*Gamma_Grid)
            G0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >        +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))
            ET=DEXP(-T)
            OneOvr2T=5D-1/T
            G1=OneOvr2T*(0.1000000000000000D+01*G0-ET)
            G2=OneOvr2T*(0.3000000000000000D+01*G1-ET)
            G3=OneOvr2T*(0.5000000000000000D+01*G2-ET)
            G4=OneOvr2T*(0.7000000000000000D+01*G3-ET)
            G5=OneOvr2T*(0.9000000000000000D+01*G4-ET)
            G6=OneOvr2T*(0.1100000000000000D+02*G5-ET)
            G7=OneOvr2T*(0.1300000000000000D+02*G6-ET)
            G8=OneOvr2T*(0.1500000000000000D+02*G7-ET)
            G9=OneOvr2T*(0.1700000000000000D+02*G8-ET)
            G10=OneOvr2T*(0.1900000000000000D+02*G9-ET)
            G11=OneOvr2T*(0.2100000000000000D+02*G10-ET)
         ELSE
            SqrtT=DSQRT(T)
            OneOvT=One/T
            G0=SqrtPi/(Two*SqrtT)
            G1=G0*0.5000000000000000D+00*OneOvT
            G2=G1*0.1500000000000000D+01*OneOvT
            G3=G2*0.2500000000000000D+01*OneOvT
            G4=G3*0.3500000000000000D+01*OneOvT
            G5=G4*0.4500000000000000D+01*OneOvT
            G6=G5*0.5500000000000000D+01*OneOvT
            G7=G6*0.6500000000000000D+01*OneOvT
            G8=G7*0.7500000000000000D+01*OneOvT
            G9=G8*0.8500000000000000D+01*OneOvT
            G10=G9*0.9500000000000000D+01*OneOvT
            G11=G10*0.1050000000000000D+02*OneOvT
         ENDIF
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*G0
         o1=o2*o1
         AuxR(1,I)=o1*G1
         o1=o2*o1
         AuxR(2,I)=o1*G2
         o1=o2*o1
         AuxR(3,I)=o1*G3
         o1=o2*o1
         AuxR(4,I)=o1*G4
         o1=o2*o1
         AuxR(5,I)=o1*G5
         o1=o2*o1
         AuxR(6,I)=o1*G6
         o1=o2*o1
         AuxR(7,I)=o1*G7
         o1=o2*o1
         AuxR(8,I)=o1*G8
         o1=o2*o1
         AuxR(9,I)=o1*G9
         o1=o2*o1
         AuxR(10,I)=o1*G10
         o1=o2*o1
         AuxR(11,I)=o1*G11
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen12(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K,IF0,IF12
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:12,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,TwoPi5x2
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "GammaGrid_77.Inc"
      INCLUDE "GammaDimensions_77.Inc"
      REAL*8 F0_0(0:Gamma_Mesh)
      REAL*8 F0_1(0:Gamma_Mesh)
      REAL*8 F0_2(0:Gamma_Mesh)
      REAL*8 F0_3(0:Gamma_Mesh)
      REAL*8 F0_4(0:Gamma_Mesh)
      REAL*8 F0_5(0:Gamma_Mesh)
      REAL*8 F0_6(0:Gamma_Mesh)
      REAL*8 F12_0(0:Gamma_Mesh)
      REAL*8 F12_1(0:Gamma_Mesh)
      REAL*8 F12_2(0:Gamma_Mesh)
      REAL*8 F12_3(0:Gamma_Mesh)
      REAL*8 F12_4(0:Gamma_Mesh)
      REAL*8 F12_5(0:Gamma_Mesh)
      REAL*8 F12_6(0:Gamma_Mesh)
      INCLUDE "F0_77.Inc"
      INCLUDE "F12_77.Inc"
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=TwoPi5x2/(RTE*SQRT(RPE))
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         IF(T.LT.1D0)THEN
            J=AINT(T*Gamma_Grid)
            G12=F12_0(J)+T*(F12_1(J)+T*(F12_2(J)+T*(F12_3(J)
     >        +T*(F12_4(J)+T*(F12_5(J)+T*F12_6(J))))))
            ET=DEXP(-T)
            TwoT=2.0D0*T
            G11=0.4347826086956522D-01*(TwoT*G12+ET)
            G10=0.4761904761904762D-01*(TwoT*G11+ET)
            G9=0.5263157894736842D-01*(TwoT*G10+ET)
            G8=0.5882352941176471D-01*(TwoT*G9+ET)
            G7=0.6666666666666667D-01*(TwoT*G8+ET)
            G6=0.7692307692307693D-01*(TwoT*G7+ET)
            G5=0.9090909090909091D-01*(TwoT*G6+ET)
            G4=0.1111111111111111D+00*(TwoT*G5+ET)
            G3=0.1428571428571428D+00*(TwoT*G4+ET)
            G2=0.2000000000000000D+00*(TwoT*G3+ET)
            G1=0.3333333333333333D+00*(TwoT*G2+ET)
            G0=TwoT*G1+ET
         ELSEIF(T.LT.Gamma_Switch)THEN
            J=AINT(T*Gamma_Grid)
            G0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >        +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))
            ET=DEXP(-T)
            OneOvr2T=5D-1/T
            G1=OneOvr2T*(0.1000000000000000D+01*G0-ET)
            G2=OneOvr2T*(0.3000000000000000D+01*G1-ET)
            G3=OneOvr2T*(0.5000000000000000D+01*G2-ET)
            G4=OneOvr2T*(0.7000000000000000D+01*G3-ET)
            G5=OneOvr2T*(0.9000000000000000D+01*G4-ET)
            G6=OneOvr2T*(0.1100000000000000D+02*G5-ET)
            G7=OneOvr2T*(0.1300000000000000D+02*G6-ET)
            G8=OneOvr2T*(0.1500000000000000D+02*G7-ET)
            G9=OneOvr2T*(0.1700000000000000D+02*G8-ET)
            G10=OneOvr2T*(0.1900000000000000D+02*G9-ET)
            G11=OneOvr2T*(0.2100000000000000D+02*G10-ET)
            G12=OneOvr2T*(0.2300000000000000D+02*G11-ET)
         ELSE
            SqrtT=DSQRT(T)
            OneOvT=One/T
            G0=SqrtPi/(Two*SqrtT)
            G1=G0*0.5000000000000000D+00*OneOvT
            G2=G1*0.1500000000000000D+01*OneOvT
            G3=G2*0.2500000000000000D+01*OneOvT
            G4=G3*0.3500000000000000D+01*OneOvT
            G5=G4*0.4500000000000000D+01*OneOvT
            G6=G5*0.5500000000000000D+01*OneOvT
            G7=G6*0.6500000000000000D+01*OneOvT
            G8=G7*0.7500000000000000D+01*OneOvT
            G9=G8*0.8500000000000000D+01*OneOvT
            G10=G9*0.9500000000000000D+01*OneOvT
            G11=G10*0.1050000000000000D+02*OneOvT
            G12=G11*0.1150000000000000D+02*OneOvT
         ENDIF
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*G0
         o1=o2*o1
         AuxR(1,I)=o1*G1
         o1=o2*o1
         AuxR(2,I)=o1*G2
         o1=o2*o1
         AuxR(3,I)=o1*G3
         o1=o2*o1
         AuxR(4,I)=o1*G4
         o1=o2*o1
         AuxR(5,I)=o1*G5
         o1=o2*o1
         AuxR(6,I)=o1*G6
         o1=o2*o1
         AuxR(7,I)=o1*G7
         o1=o2*o1
         AuxR(8,I)=o1*G8
         o1=o2*o1
         AuxR(9,I)=o1*G9
         o1=o2*o1
         AuxR(10,I)=o1*G10
         o1=o2*o1
         AuxR(11,I)=o1*G11
         o1=o2*o1
         AuxR(12,I)=o1*G12
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen0(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:0,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(Pi=3.1415926535897932385D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=(Pi/RPE)**15D-1
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         ET=EXP(-T)
         TwoT=2.0D0*T
         o1=Upq
         AuxR(0,I)=o1*ET
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen1(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:1,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(Pi=3.1415926535897932385D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=(Pi/RPE)**15D-1
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         ET=EXP(-T)
         TwoT=2.0D0*T
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*ET
         o1=o2*o1
         AuxR(1,I)=o1*ET
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen2(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:2,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(Pi=3.1415926535897932385D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=(Pi/RPE)**15D-1
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         ET=EXP(-T)
         TwoT=2.0D0*T
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*ET
         o1=o2*o1
         AuxR(1,I)=o1*ET
         o1=o2*o1
         AuxR(2,I)=o1*ET
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen3(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:3,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(Pi=3.1415926535897932385D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=(Pi/RPE)**15D-1
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         ET=EXP(-T)
         TwoT=2.0D0*T
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*ET
         o1=o2*o1
         AuxR(1,I)=o1*ET
         o1=o2*o1
         AuxR(2,I)=o1*ET
         o1=o2*o1
         AuxR(3,I)=o1*ET
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen4(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:4,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(Pi=3.1415926535897932385D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=(Pi/RPE)**15D-1
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         ET=EXP(-T)
         TwoT=2.0D0*T
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*ET
         o1=o2*o1
         AuxR(1,I)=o1*ET
         o1=o2*o1
         AuxR(2,I)=o1*ET
         o1=o2*o1
         AuxR(3,I)=o1*ET
         o1=o2*o1
         AuxR(4,I)=o1*ET
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen5(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:5,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(Pi=3.1415926535897932385D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=(Pi/RPE)**15D-1
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         ET=EXP(-T)
         TwoT=2.0D0*T
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*ET
         o1=o2*o1
         AuxR(1,I)=o1*ET
         o1=o2*o1
         AuxR(2,I)=o1*ET
         o1=o2*o1
         AuxR(3,I)=o1*ET
         o1=o2*o1
         AuxR(4,I)=o1*ET
         o1=o2*o1
         AuxR(5,I)=o1*ET
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen6(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:6,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(Pi=3.1415926535897932385D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=(Pi/RPE)**15D-1
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         ET=EXP(-T)
         TwoT=2.0D0*T
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*ET
         o1=o2*o1
         AuxR(1,I)=o1*ET
         o1=o2*o1
         AuxR(2,I)=o1*ET
         o1=o2*o1
         AuxR(3,I)=o1*ET
         o1=o2*o1
         AuxR(4,I)=o1*ET
         o1=o2*o1
         AuxR(5,I)=o1*ET
         o1=o2*o1
         AuxR(6,I)=o1*ET
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen7(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:7,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(Pi=3.1415926535897932385D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=(Pi/RPE)**15D-1
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         ET=EXP(-T)
         TwoT=2.0D0*T
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*ET
         o1=o2*o1
         AuxR(1,I)=o1*ET
         o1=o2*o1
         AuxR(2,I)=o1*ET
         o1=o2*o1
         AuxR(3,I)=o1*ET
         o1=o2*o1
         AuxR(4,I)=o1*ET
         o1=o2*o1
         AuxR(5,I)=o1*ET
         o1=o2*o1
         AuxR(6,I)=o1*ET
         o1=o2*o1
         AuxR(7,I)=o1*ET
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen8(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:8,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(Pi=3.1415926535897932385D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=(Pi/RPE)**15D-1
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         ET=EXP(-T)
         TwoT=2.0D0*T
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*ET
         o1=o2*o1
         AuxR(1,I)=o1*ET
         o1=o2*o1
         AuxR(2,I)=o1*ET
         o1=o2*o1
         AuxR(3,I)=o1*ET
         o1=o2*o1
         AuxR(4,I)=o1*ET
         o1=o2*o1
         AuxR(5,I)=o1*ET
         o1=o2*o1
         AuxR(6,I)=o1*ET
         o1=o2*o1
         AuxR(7,I)=o1*ET
         o1=o2*o1
         AuxR(8,I)=o1*ET
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen9(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:9,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(Pi=3.1415926535897932385D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=(Pi/RPE)**15D-1
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         ET=EXP(-T)
         TwoT=2.0D0*T
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*ET
         o1=o2*o1
         AuxR(1,I)=o1*ET
         o1=o2*o1
         AuxR(2,I)=o1*ET
         o1=o2*o1
         AuxR(3,I)=o1*ET
         o1=o2*o1
         AuxR(4,I)=o1*ET
         o1=o2*o1
         AuxR(5,I)=o1*ET
         o1=o2*o1
         AuxR(6,I)=o1*ET
         o1=o2*o1
         AuxR(7,I)=o1*ET
         o1=o2*o1
         AuxR(8,I)=o1*ET
         o1=o2*o1
         AuxR(9,I)=o1*ET
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen10(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:10,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(Pi=3.1415926535897932385D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=(Pi/RPE)**15D-1
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         ET=EXP(-T)
         TwoT=2.0D0*T
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*ET
         o1=o2*o1
         AuxR(1,I)=o1*ET
         o1=o2*o1
         AuxR(2,I)=o1*ET
         o1=o2*o1
         AuxR(3,I)=o1*ET
         o1=o2*o1
         AuxR(4,I)=o1*ET
         o1=o2*o1
         AuxR(5,I)=o1*ET
         o1=o2*o1
         AuxR(6,I)=o1*ET
         o1=o2*o1
         AuxR(7,I)=o1*ET
         o1=o2*o1
         AuxR(8,I)=o1*ET
         o1=o2*o1
         AuxR(9,I)=o1*ET
         o1=o2*o1
         AuxR(10,I)=o1*ET
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen11(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:11,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(Pi=3.1415926535897932385D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=(Pi/RPE)**15D-1
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         ET=EXP(-T)
         TwoT=2.0D0*T
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*ET
         o1=o2*o1
         AuxR(1,I)=o1*ET
         o1=o2*o1
         AuxR(2,I)=o1*ET
         o1=o2*o1
         AuxR(3,I)=o1*ET
         o1=o2*o1
         AuxR(4,I)=o1*ET
         o1=o2*o1
         AuxR(5,I)=o1*ET
         o1=o2*o1
         AuxR(6,I)=o1*ET
         o1=o2*o1
         AuxR(7,I)=o1*ET
         o1=o2*o1
         AuxR(8,I)=o1*ET
         o1=o2*o1
         AuxR(9,I)=o1*ET
         o1=o2*o1
         AuxR(10,I)=o1*ET
         o1=o2*o1
         AuxR(11,I)=o1*ET
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen12(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)
      IMPLICIT DOUBLE PRECISION(G,O,S,F)
      INTEGER Nc,Ns,Mesh,I,J,K
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:12,*)
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi
      REAL*8 SqrtT,SqrtPi,One,Two
      PARAMETER(Pi=3.1415926535897932385D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      DO 100 I=1,Nc
         PQx(I)=-(P(1)-Q(1,I))
         PQy(I)=-(P(2)-Q(2,I))
         PQz(I)=-(P(3)-Q(3,I))
         RTE=ZetaP*ZetaQ(I)
         RPE=ZetaP+ZetaQ(I)
         Omega=RTE/RPE
         Upq=(Pi/RPE)**15D-1
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))
         ET=EXP(-T)
         TwoT=2.0D0*T
         o1=Upq
         o2=-2D0*Omega
         AuxR(0,I)=o1*ET
         o1=o2*o1
         AuxR(1,I)=o1*ET
         o1=o2*o1
         AuxR(2,I)=o1*ET
         o1=o2*o1
         AuxR(3,I)=o1*ET
         o1=o2*o1
         AuxR(4,I)=o1*ET
         o1=o2*o1
         AuxR(5,I)=o1*ET
         o1=o2*o1
         AuxR(6,I)=o1*ET
         o1=o2*o1
         AuxR(7,I)=o1*ET
         o1=o2*o1
         AuxR(8,I)=o1*ET
         o1=o2*o1
         AuxR(9,I)=o1*ET
         o1=o2*o1
         AuxR(10,I)=o1*ET
         o1=o2*o1
         AuxR(11,I)=o1*ET
         o1=o2*o1
         AuxR(12,I)=o1*ET
 100  CONTINUE
      RETURN
      END
