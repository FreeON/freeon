C     -----------------------------------------------------------------
C     Reza Ahmadi and Jan Almlof
C     Chemical Physics Letters; 1 Dec. 1995; vol.246, no.4-5, p.364-70
C     Buffed by M. Challacombe
C     -----------------------------------------------------------------

      SUBROUTINE RAhmadiJAlmlof95c(EllP,LenP,NukeX,Threshold,ZP,P,Ket,C)
C
      IMPLICIT NONE
C
      INTEGER I,J,K,L
      INTEGER Nc,EllP,EllQ,EllPQ,LenP,LenQ,LenPQ,PQEll,PQLen,MaxEll
      REAL*8  ZP,NukeX,Threshold,P(*),Ket(*),C(*)
C
      INTEGER iOff,iNc,iEst,iZQ,iQ,iCo
      INTEGER iPQx,iPQy,iPQz,iAux,iArr,iEnd
C
      INTEGER HGLength(0:14)
      DATA (HGLength(I),I=0,14)/1,4,10,20,35,56,84,120,165,220,286,
     >                         364,455,560,680/
C
      INTEGER iClust,MxEll,MxLen,iWork
      PARAMETER (iClust=511)
      PARAMETER (MxEll=13)
      PARAMETER (MxLen=(MxEll+1)*(MxEll+2)*(MxEll+3)/6)
      PARAMETER (iWork=3*iClust+iClust*(MxLen+MxEll+1))
      REAL*8 Wrk(iWork)
C
      MaxEll=C(1)
      iOff=2
C
      DO EllQ=0,MaxEll
C
         EllPQ=EllP+EllQ
         PQEll=EllP+EllQ
         LenQ=HGLength(EllQ)
         PQLen=HGLength(EllPQ)
         LenPQ=LenP*LenQ
C
         iNc   =iOff
         Nc    =C(iNc)
         iEst  =iNc+1
         iZQ   =iEst+Nc
         iQ    =iZQ+Nc
         iCo   =iQ+3*Nc
         iOff  =iCo+LenQ*Nc
C
c         IF(.NOT.(ZP.EQ.NukeX.AND.EllQ==1))CYCLE


         IF(Nc==0.OR.C(iEst).LT.Threshold)THEN
            CYCLE
         ELSEIF(C(iEst+Nc-1).LT.Threshold)THEN
           K=0
           L=Nc
           Nc=Nc/2
           DO WHILE(K.LT.Nc.AND.L.GT.Nc)
              IF(C(iEst+Nc-1).LE.Threshold)THEN
                  L=Nc
               ELSE
                  K=Nc
               ENDIF
               Nc=(K+L)/2
            ENDDO
         ENDIF
C
         iPQx=1
         iPQy=iPQx+Nc
         iPQz=iPQy+Nc
         iAux=iPQz+Nc
         iArr=iAux+Nc*(PQEll+1)
         iEnd=iArr+Nc*(PQLen+1)
C
         IF(iEnd>iWork)
     >        STOP 'Local Memory exceeded in AA95 '
C
         IF(ZP.EQ.NukeX.AND.EllQ==0)THEN
C           CALL SPECIAL PURPOSE CODE THAT AVOIDS THE NUCLEAR NUCLEAR SELF INTERACTION
            IF(EllP==0)THEN
               CALL ENuE0(Nc,P,C(iQ),ZP,C(iZQ),C(iCo),Ket)
            ELSEIF(EllP==1)THEN
               CALL ENuE1(Nc,P,C(iQ),ZP,C(iZQ),C(iCo),Ket)
            ELSE
               WRITE(*,*)' ZP = ',ZP,' EllP = ',EllP,' EllQ = ',EllQ
               STOP 'Unsupported nuclear derivatives in AA95c.F'
            ENDIF
         ELSE
         SELECT CASE(EllPQ)
         CASE(0)
C           CALL PAPIEX_START(10,"< CASE 0 > ")
            CALL AuxIGen0(Nc,P,C(iQ),ZP,C(iZQ),Wrk(iPQx),Wrk(iPQy),
     >           Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR0(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >           Wrk(iAux),Wrk(iArr))
            CALL KTrax_0_0(Nc,Wrk(iArr),C(iCo),Ket)
C           CALL PAPIEX_STOP(10)
         CASE(1)
C           CALL PAPIEX_START(11,"< CASE 1 >")
            CALL AuxIGen1(Nc,P,C(iQ),ZP,C(iZQ),Wrk(iPQx),Wrk(iPQy),
     >           Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR1(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >           Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               CALL KTrax_0_1(Nc,Wrk(iArr),C(iCo),Ket)
            ELSE
               CALL KTrax_1_0(Nc,Wrk(iArr),C(iCo),Ket)
            ENDIF
C           CALL PAPIEX_STOP(11)
         CASE(2)
C           CALL PAPIEX_START(12,"< CASE 2 > ")
            CALL AuxIGen2(Nc,P,C(iQ),ZP,C(iZQ),Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR2(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               CALL KTrax_0_2(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_1(Nc,Wrk(iArr),C(iCo),Ket)
            ELSE
               CALL KTrax_2_0(Nc,Wrk(iArr),C(iCo),Ket)
            ENDIF
C           CALL PAPIEX_STOP(12)
         CASE(3)
            CALL AuxIGen3(Nc,P,C(iQ),ZP,C(iZQ),Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR3(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               CALL KTrax_0_3(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_2(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_1(Nc,Wrk(iArr),C(iCo),Ket)
            ELSE
               CALL KTrax_3_0(Nc,Wrk(iArr),C(iCo),Ket)
            ENDIF
         CASE(4)
            CALL AuxIGen4(Nc,P,C(iQ),ZP,C(iZQ),Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR4(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               CALL KTrax_0_4(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_3(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_2(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_1(Nc,Wrk(iArr),C(iCo),Ket)
            ELSE
               CALL KTrax_4_0(Nc,Wrk(iArr),C(iCo),Ket)
            ENDIF
         CASE(5)
            CALL AuxIGen5(Nc,P,C(iQ),ZP,C(iZQ),Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR5(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               CALL KTrax_0_5(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_4(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_3(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_2(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==4)THEN
               CALL KTrax_4_1(Nc,Wrk(iArr),C(iCo),Ket)
            ELSE
               CALL KTrax_5_0(Nc,Wrk(iArr),C(iCo),Ket)
            ENDIF
         CASE(6)
            CALL AuxIGen6(Nc,P,C(iQ),ZP,C(iZQ),Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR6(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               CALL KTrax_0_6(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_5(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_4(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_3(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==4)THEN
               CALL KTrax_4_2(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==5)THEN
               CALL KTrax_5_1(Nc,Wrk(iArr),C(iCo),Ket)
            ELSE
               CALL KTrax_6_0(Nc,Wrk(iArr),C(iCo),Ket)
            ENDIF
         CASE(7)
            CALL AuxIGen7(Nc,P,C(iQ),ZP,C(iZQ),Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR7(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               STOP '(0|3+4) in AA95 '
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_6(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_5(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_4(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==4)THEN
               CALL KTrax_4_3(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==5)THEN
               CALL KTrax_5_2(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==6)THEN
               CALL KTrax_6_1(Nc,Wrk(iArr),C(iCo),Ket)
            ELSE
               CALL KTrax_7_0(Nc,Wrk(iArr),C(iCo),Ket)
            ENDIF
         CASE(8)
            CALL AuxIGen8(Nc,P,C(iQ),ZP,C(iZQ),Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR8(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               STOP '(0|8) in AA95 '
            ELSEIF(EllP==1)THEN
               STOP '(1|8) in AA95 '
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_6(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_5(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==4)THEN
               CALL KTrax_4_4(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==5)THEN
               CALL KTrax_5_3(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==6)THEN
               CALL KTrax_6_2(Nc,Wrk(iArr),C(iCo),Ket)
            ELSEIF(EllP==7)THEN
               CALL KTrax_7_1(Nc,Wrk(iArr),C(iCo),Ket)
            ELSE
               STOP '(8+|0) in AA95'
            ENDIF


         CASE DEFAULT
            STOP 'AA95 '
         END SELECT
C
      ENDIF

      ENDDO
C
      RETURN
      END
C
C     SPECIAL PURPOSE CODE THAT AVOIDS THE NUCLEAR NUCLEAR SELF INTERACTION
C
      SUBROUTINE ENuE0(Nc,P,Q,ZetaP,ZetaQ,QCo,Ket)
C
      INTEGER Nc,I,J,IF0
      REAL*8 P(3),Q(3,*),ZetaP,ZetaQ(*),QCo(*),Ket(*)
      REAL*8 PQx,PQy,PQz,RTE,RPE,Omega,Upq,SqrtPi,TwoPi5x2
      REAL*8 PQ2,T,T1,T2,T3

      REAL*8 Switch,F0Asymp,F1Asymp,F2Asymp,F3Asymp,F4Asymp,
     >              F5Asymp,F6Asymp,F7Asymp,F8Asymp,F9Asymp,
     >              F10Asymp,F11Asymp,F12Asymp,F13Asymp,Grid

C
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)
      PARAMETER(SqrtPi=1.7724538509055160273D0)
C
      INCLUDE "GammaDimensions_77.Inc"
      INCLUDE "GammaGrid_77.Inc"
      REAL*8 F0_0(0:Gamma_Mesh)
      REAL*8 F0_1(0:Gamma_Mesh)
      REAL*8 F0_2(0:Gamma_Mesh)
      REAL*8 F0_3(0:Gamma_Mesh)
      REAL*8 F0_4(0:Gamma_Mesh)
      REAL*8 F0_5(0:Gamma_Mesh)
      REAL*8 F0_6(0:Gamma_Mesh)
      INCLUDE "F0_77.Inc"
C
      DO 100 I=1,Nc
         PQx=-(P(1)-Q(1,I))
         PQy=-(P(2)-Q(2,I))
         PQz=-(P(3)-Q(3,I))
         PQ2=PQx*PQx+PQy*PQy+PQz*PQz
         IF(.NOT.(PQ2.LT.1D-16.AND.ZetaQ(I).EQ.ZetaP))THEN
            RTE=ZetaP*ZetaQ(I)
            RPE=ZetaP+ZetaQ(I)
            Omega=RTE/RPE
            Upq=TwoPi5x2/(RTE*SQRT(RPE))
            T=Omega*PQ2
            IF(T.LT.Gamma_Switch)THEN
               j=AINT(T*Gamma_Grid)
               Ket(1)=Ket(1)+QCo(I)*Upq*(
     >                F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >               +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J)))))))
            ELSE
               Ket(1)=Ket(1)+QCo(I)*Upq*SqrtPi/(2D0*DSQRT(T))
            ENDIF
         ENDIF
 100  CONTINUE
      RETURN
      END

      SUBROUTINE ENuE1(Nc,P,Q,ZetaP,ZetaQ,QCo,Ket)
C
      INTEGER Nc,I,J,IF0,IF1
      REAL*8 P(3),Q(3,*),ZetaP,ZetaQ(*),QCo(*),Ket(*)
      REAL*8 PQx,PQy,PQz,RTE,RPE,Omega,Upq,SqrtPi,TwoPi5x2,OneOvr2T
      REAL*8 PQ2,T,T1,T2,T3,F0,F1,Aux0,Aux1,One,Two,ET,TwoT,SQRTT,OneOvT

      REAL*8 Switch,F0Asymp,F1Asymp,F2Asymp,F3Asymp,F4Asymp,
     >              F5Asymp,F6Asymp,F7Asymp,F8Asymp,F9Asymp,
     >              F10Asymp,F11Asymp,F12Asymp,F13Asymp,Grid
C
      PARAMETER(One=1D0,Two=2D0)
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)
      PARAMETER(SqrtPi=1.7724538509055160273D0)
C
      INCLUDE "GammaDimensions_77.Inc"
      INCLUDE "GammaGrid_77.Inc"
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
c
      INCLUDE "F0_77.Inc"
      INCLUDE "F1_77.Inc"
C
      DO 100 I=1,Nc
         PQx=-(P(1)-Q(1,I))
         PQy=-(P(2)-Q(2,I))
         PQz=-(P(3)-Q(3,I))
         PQ2=PQx*PQx+PQy*PQy+PQz*PQz
         IF(.NOT.(PQ2.LT.1D-16.AND.ZetaQ(I).EQ.ZetaP))THEN
            RTE=ZetaP*ZetaQ(I)
            RPE=ZetaP+ZetaQ(I)
            Omega=RTE/RPE
            Upq=TwoPi5x2/(RTE*SQRT(RPE))
            T=Omega*PQ2
            IF(T.LT.1D0)THEN
               J=AINT(T*Gamma_Grid)
               F1=F1_0(J)+T*(F1_1(J)+T*(F1_2(J)+T*(F1_3(J)
     >              +T*(F1_4(J)+T*(F1_5(J)+T*F1_6(J))))))
               ET=DEXP(-T)
               TwoT=2.0D0*T
               F0=TwoT*F1+ET
            ELSEIF(T.LT.Gamma_Switch)THEN
               J=AINT(T*Gamma_Grid)
               F0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)
     >              +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))
               ET=DEXP(-T)
               OneOvr2T=5D-1/T
               F1=OneOvr2T*(0.1000000000000000D+01*F0-ET)
            ELSE
               SqrtT=DSQRT(T)
               OneOvT=One/T
               F0=SqrtPi/(Two*SqrtT)
               F1=F0*0.5000000000000000D+00*OneOvT
            ENDIF

            Aux0=F0*Upq
            Aux1=-2D0*Omega*F1*Upq
            Ket(1)=Ket(1)+QCo(I)*Aux0
            Ket(2)=Ket(2)+QCo(I)*PQx*Aux1
            Ket(3)=Ket(3)+QCo(I)*PQy*Aux1
            Ket(4)=Ket(4)+QCo(I)*PQz*Aux1
         ENDIF
 100  CONTINUE
      RETURN
      END
