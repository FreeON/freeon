        IF (T1<Gamma_Switch)THEN
          MG=AINT(T1*Gamma_Grid)
          ME=AINT(T1*Exp_Grid)
          T2=T1*T1
          T3=T2*T1
          T4=T2*T2
          TwoT=2.0D0*T1
          ET=Exp_0(ME)+T1*Exp_1(ME)+T2*Exp_2(ME)+T3*Exp_3(ME)+T4*Exp_4(ME)
          G3=F3_0(MG) +T1*F3_1(MG) +T2*F3_2(MG) +T3*F3_3(MG) +T4*F3_4(MG)
          G2=.200000000000000D+00*(TwoT*G3+ET)
          G1=.333333333333333D+00*(TwoT*G2+ET)
          G0=TwoT*G1+ET
          R1=Rkk*G0
          R2=Rkk*G1
          R3=Rkk*G2
          R4=Rkk*G3
        ELSE
          T2=1.0D0/T1
          T3=DSQRT(T2)
          R1=Rkk*GammAss(0)*T3
          T3=T3*T2
          R2=Rkk*GammAss(1)*T3
          T3=T3*T2
          R3=Rkk*GammAss(2)*T3
          T3=T3*T2
          R4=Rkk*GammAss(3)*T3
        ENDIF

