        IF (T1<Gamma_Switch)THEN
          MG=AINT(T1*Gamma_Grid)
          T2=T1*T1
          T3=T2*T1
          T4=T2*T2
          G1=F1_0(MG) +T1*F1_1(MG) +T2*F1_2(MG) +T3*F1_3(MG) +T4*F1_4(MG)
          G0=F0_0(MG) +T1*F0_1(MG) +T2*F0_2(MG) +T3*F0_3(MG) +T4*F0_4(MG)
          R1=Rkk*G0
          R2=Rkk*G1
        ELSE
          T2=1.0D0/T1
          T3=DSQRT(T2)
          R1=Rkk*GammAss(0)*T3
          T3=T3*T2
          R2=Rkk*GammAss(1)*T3
        ENDIF

