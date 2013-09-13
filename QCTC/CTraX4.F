      SUBROUTINE CTX4_0_0_C(Cp,Cpq,Spq,Cq,Sq)
        IMPLICIT NONE
        REAL*8 Cp,Cpq(0:14),Spq(0:14),Cq(0:14),Sq(0:14)
        REAL*8 COne,CTwo,Two,W
        COMMON /CTraXScratch/ W(1000)
        PARAMETER (Two=2D0)
        W(1)=Cpq(0)*Cq(0)+Cpq(1)*Cq(1)
        W(2)=Cpq(3)*Cq(3)+Cpq(6)*Cq(6)+Cpq(10)*Cq(10)
        COne=W(1)+W(2)
        W(1)=Cpq(2)*Cq(2)+Cpq(4)*Cq(4)
        W(2)=Cpq(5)*Cq(5)+Cpq(7)*Cq(7)+Cpq(8)*Cq(8)
        W(3)=Cpq(9)*Cq(9)+Cpq(11)*Cq(11)
        W(4)=Cpq(12)*Cq(12)+Cpq(13)*Cq(13)+Cpq(14)*Cq(14)
        W(5)=Spq(2)*Sq(2)+Spq(4)*Sq(4)
        W(6)=Spq(5)*Sq(5)+Spq(7)*Sq(7)+Spq(8)*Sq(8)
        W(7)=Spq(9)*Sq(9)+Spq(11)*Sq(11)
        W(8)=Spq(12)*Sq(12)+Spq(13)*Sq(13)+Spq(14)*Sq(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        Cp=COne+Cp+CTwo*Two
      RETURN
      END
      SUBROUTINE IrRegular4(PQx,PQy,PQz,Cpq,Spq)
         IMPLICIT NONE
         REAL*8 Cpq(0:14),Spq(0:14),LegendreP(0:14)
         REAL*8 Sine(0:4),Cosine(0:4),CoFact(2:4)
         REAL*8 PQx,PQy,PQz,PQ2,pqx2,pqy2,pqxy,PQ
         REAL*8 Zero,One,Two,TT,TwoC,OneOvPQ
         REAL*8 CoTan,OneOvPQxy,RS,SQ,PQToThMnsL
         PARAMETER (Two=2D0,One=1D0,Zero=0D0)
         PQx2=PQx*PQx
         PQy2=PQy*PQy
         PQ=DSQRT(PQx2+PQy2+PQz*PQz)
         OneOvPQ=One/PQ
         CoTan=PQz*OneOvPQ
!        Compute Cosine and Sine by recursion
         Sine(0)  =Zero
         Cosine(0)=One
         PQxy=DSQRT(PQx2+PQy2)
         IF(PQxy/=Zero)THEN
            OneOvPQxy=One/PQxy
            Sine(1)  =PQx*OneOvPQxy
            Cosine(1)=PQy*OneOvPQxy
         ELSE
            Sine(1)  =0.70710678118654752D0
            Cosine(1)=0.70710678118654752D0
         ENDIF
         TwoC=Two*Cosine(1)
         Sine(2)=TwoC*Sine(1)-Sine(0)
         Cosine(2)=TwoC*Cosine(1)-Cosine(0)
         Sine(3)=TwoC*Sine(2)-Sine(1)
         Cosine(3)=TwoC*Cosine(2)-Cosine(1)
         Sine(4)=TwoC*Sine(3)-Sine(2)
         Cosine(4)=TwoC*Cosine(3)-Cosine(2)
!        ----------------------------------------------------------------
!        Associated Legendre Polynomials by recursion
!
!        M=0,Ell; Compute LegendreP(M,M)
         RS=One
         Sq=SQRT(ABS(One-CoTan**2))
         LegendreP(0)=+1.000000000000000D+00*RS
         RS=RS*Sq
         LegendreP(2)=-1.000000000000000D+00*RS
         RS=RS*Sq
         LegendreP(5)=+3.000000000000000D+00*RS
         RS=RS*Sq
         LegendreP(9)=-1.500000000000000D+01*RS
         RS=RS*Sq
         LegendreP(14)=+1.050000000000000D+02*RS
!        ----------------------------------------------------------------
!        M=0,Ell-1; Compute LegendreP(M+1,M)
         LegendreP(1)=+1.000000000000000D+00*CoTan*LegendreP(0)
         LegendreP(4)=+3.000000000000000D+00*CoTan*LegendreP(2)
         LegendreP(8)=+5.000000000000000D+00*CoTan*LegendreP(5)
         LegendreP(13)=+7.000000000000000D+00*CoTan*LegendreP(9)
!        ----------------------------------------------------------------
         CoFact(2)=+3.000000000000000D+00*CoTan
         CoFact(3)=+5.000000000000000D+00*CoTan
         CoFact(4)=+7.000000000000000D+00*CoTan
!        ----------------------------------------------------------------
!        M=0,Ell-2; L=M+2,Ell; Compute LegendreP(L,M)
         LegendreP(3)=CoFact(2)*LegendreP(1)-1.000000000000000D+00*LegendreP(0)
         LegendreP(6)=CoFact(3)*LegendreP(3)-4.000000000000000D+00*LegendreP(1)
         LegendreP(10)=CoFact(4)*LegendreP(6)-9.000000000000000D+00*LegendreP(3)
         LegendreP(7)=CoFact(3)*LegendreP(4)-3.000000000000000D+00*LegendreP(2)
         LegendreP(11)=CoFact(4)*LegendreP(7)-8.000000000000000D+00*LegendreP(4)
         LegendreP(12)=CoFact(4)*LegendreP(8)-5.000000000000000D+00*LegendreP(5)
!        ----------------------------------------------------------------
!        Compute Cpq(L,M) and Spq(L,M)
         PQToThMnsL=OneOvPQ
         TT=PQToThMnsL*LegendreP(0)
         Spq(0)=TT*Sine(0)
         Cpq(0)=TT*Cosine(0)
         PQToThMnsL=PQToThMnsL*OneOvPQ
         TT=PQToThMnsL*LegendreP(1)
         Spq(1)=TT*Sine(0)
         Cpq(1)=TT*Cosine(0)
         TT=PQToThMnsL*LegendreP(2)
         Spq(2)=TT*Sine(1)
         Cpq(2)=TT*Cosine(1)
         PQToThMnsL=PQToThMnsL*OneOvPQ
         TT=PQToThMnsL*LegendreP(3)
         Spq(3)=TT*Sine(0)
         Cpq(3)=TT*Cosine(0)
         TT=PQToThMnsL*LegendreP(4)
         Spq(4)=TT*Sine(1)
         Cpq(4)=TT*Cosine(1)
         TT=PQToThMnsL*LegendreP(5)
         Spq(5)=TT*Sine(2)
         Cpq(5)=TT*Cosine(2)
         PQToThMnsL=PQToThMnsL*OneOvPQ
         TT=PQToThMnsL*LegendreP(6)
         Spq(6)=TT*Sine(0)
         Cpq(6)=TT*Cosine(0)
         TT=PQToThMnsL*LegendreP(7)
         Spq(7)=TT*Sine(1)
         Cpq(7)=TT*Cosine(1)
         TT=PQToThMnsL*LegendreP(8)
         Spq(8)=TT*Sine(2)
         Cpq(8)=TT*Cosine(2)
         TT=PQToThMnsL*LegendreP(9)
         Spq(9)=TT*Sine(3)
         Cpq(9)=TT*Cosine(3)
         PQToThMnsL=PQToThMnsL*OneOvPQ
         TT=PQToThMnsL*LegendreP(10)
         Spq(10)=TT*Sine(0)
         Cpq(10)=TT*Cosine(0)
         TT=PQToThMnsL*LegendreP(11)
         Spq(11)=TT*Sine(1)
         Cpq(11)=TT*Cosine(1)
         TT=PQToThMnsL*LegendreP(12)
         Spq(12)=TT*Sine(2)
         Cpq(12)=TT*Cosine(2)
         TT=PQToThMnsL*LegendreP(13)
         Spq(13)=TT*Sine(3)
         Cpq(13)=TT*Cosine(3)
         TT=PQToThMnsL*LegendreP(14)
         Spq(14)=TT*Sine(4)
         Cpq(14)=TT*Cosine(4)
      RETURN
      END
      SUBROUTINE CTraX4(PQx,PQy,PQz,Cq,Sq,Cp,Sp)
         IMPLICIT NONE
         REAL*8 PQx,PQy,PQz
         REAL*8 Cq(0:14),Sq(0:14),Cp(0:0),Sp(0:0)
         REAL*8 Cpq(0:14),Spq(0:14)
         CALL IrRegular4(PQx,PQy,PQz,Cpq,Spq)
         CALL CTX4_0_0_C(Cp(0),Cpq,Spq,Cq,Sq)
      RETURN
      END
      SUBROUTINE CTX104_0_0_C(Cp,Cpq,Spq,Cq,Sq)
        IMPLICIT NONE
        REAL*8 Cp,Cpq(0:20),Spq(0:20),Cq(0:14),Sq(0:14)
        REAL*8 COne,CTwo,Two,W
        COMMON /CTraXScratch/ W(1000)
        PARAMETER (Two=2D0)
        W(1)=Cpq(0)*Cq(0)+Cpq(1)*Cq(1)
        W(2)=Cpq(3)*Cq(3)+Cpq(6)*Cq(6)+Cpq(10)*Cq(10)
        COne=W(1)+W(2)
        W(1)=Cpq(2)*Cq(2)+Cpq(4)*Cq(4)
        W(2)=Cpq(5)*Cq(5)+Cpq(7)*Cq(7)+Cpq(8)*Cq(8)
        W(3)=Cpq(9)*Cq(9)+Cpq(11)*Cq(11)
        W(4)=Cpq(12)*Cq(12)+Cpq(13)*Cq(13)+Cpq(14)*Cq(14)
        W(5)=Spq(2)*Sq(2)+Spq(4)*Sq(4)
        W(6)=Spq(5)*Sq(5)+Spq(7)*Sq(7)+Spq(8)*Sq(8)
        W(7)=Spq(9)*Sq(9)+Spq(11)*Sq(11)
        W(8)=Spq(12)*Sq(12)+Spq(13)*Sq(13)+Spq(14)*Sq(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        Cp=COne+Cp+CTwo*Two
      RETURN
      END
      SUBROUTINE CTX104_0_1_C(Cp,Cpq,Spq,Cq,Sq)
        IMPLICIT NONE
        REAL*8 Cp,Cpq(0:20),Spq(0:20),Cq(0:14),Sq(0:14)
        REAL*8 COne,CTwo,Two,W
        COMMON /CTraXScratch/ W(1000)
        PARAMETER (Two=2D0)
        W(1)=-(Cpq(1)*Cq(0))-Cpq(3)*Cq(1)
        W(2)=-(Cpq(6)*Cq(3))-Cpq(10)*Cq(6)-Cpq(15)*Cq(10)
        COne=W(1)+W(2)
        W(1)=-(Cpq(4)*Cq(2))-Cpq(7)*Cq(4)
        W(2)=-(Cpq(8)*Cq(5))-Cpq(11)*Cq(7)-Cpq(12)*Cq(8)
        W(3)=-(Cpq(13)*Cq(9))-Cpq(16)*Cq(11)
        W(4)=-(Cpq(17)*Cq(12))-Cpq(18)*Cq(13)-Cpq(19)*Cq(14)
        W(5)=-(Spq(4)*Sq(2))-Spq(7)*Sq(4)
        W(6)=-(Spq(8)*Sq(5))-Spq(11)*Sq(7)-Spq(12)*Sq(8)
        W(7)=-(Spq(13)*Sq(9))-Spq(16)*Sq(11)
        W(8)=-(Spq(17)*Sq(12))-Spq(18)*Sq(13)-Spq(19)*Sq(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        Cp=COne+Cp+CTwo*Two
      RETURN
      END
      SUBROUTINE CTX104_1_1_C(Cp,Cpq,Spq,Cq,Sq)
        IMPLICIT NONE
        REAL*8 Cp,Cpq(0:20),Spq(0:20),Cq(0:14),Sq(0:14)
        REAL*8 COne,CTwo,Two,W
        COMMON /CTraXScratch/ W(1000)
        PARAMETER (Two=2D0)
        W(1)=-(Cpq(2)*Cq(0))-Cpq(4)*Cq(1)
        W(2)=Cpq(3)*Cq(2)-Cpq(5)*Cq(2)
        W(3)=Cpq(6)*Cq(4)-Cpq(8)*Cq(4)
        W(4)=-(Cpq(9)*Cq(5))+Cpq(7)*(-Cq(3)+Cq(5))
        W(5)=Cpq(10)*Cq(7)-Cpq(13)*Cq(8)
        W(6)=Cpq(11)*(-Cq(6)+Cq(8))-Cpq(14)*Cq(9)
        W(7)=Cpq(12)*(-Cq(7)+Cq(9))+Cpq(15)*Cq(11)
        W(8)=Cpq(16)*(-Cq(10)+Cq(12))-Cpq(19)*Cq(13)
        W(9)=Cpq(17)*(-Cq(11)+Cq(13))-Cpq(20)*Cq(14)
        W(10)=Cpq(18)*(-Cq(12)+Cq(14))-Spq(5)*Sq(2)
        W(11)=-(Spq(8)*Sq(4))+Spq(7)*Sq(5)
        W(12)=-(Spq(9)*Sq(5))+Spq(11)*Sq(8)
        W(13)=-(Spq(13)*Sq(8))-Spq(14)*Sq(9)
        W(14)=Spq(12)*(-Sq(7)+Sq(9))+Spq(16)*Sq(12)
        W(15)=-(Spq(19)*Sq(13))+Spq(17)*(-Sq(11)+Sq(13))
        W(16)=-(Spq(20)*Sq(14))+Spq(18)*(-Sq(12)+Sq(14))
        W(17)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(18)=W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)+W(16)
        CTwo=W(17)+W(18)
        Cp=Cp+CTwo*Two
      RETURN
      END
      SUBROUTINE CTX104_1_1_S(Sp,Cpq,Spq,Cq,Sq)
        IMPLICIT NONE
        REAL*8 Sp,Cpq(0:20),Spq(0:20),Cq(0:14),Sq(0:14)
        REAL*8 SOne,STwo,Two,W
        COMMON /CTraXScratch/ W(1000)
        PARAMETER (Two=2D0)
        W(1)=-(Cq(0)*Spq(2))-Cq(1)*Spq(4)
        W(2)=-(Cq(2)*Spq(5))-(Cq(3)+Cq(5))*Spq(7)
        W(3)=-(Cq(4)*Spq(8))-Cq(5)*Spq(9)
        W(4)=-((Cq(6)+Cq(8))*Spq(11))-(Cq(7)+Cq(9))*Spq(12)
        W(5)=-(Cq(8)*Spq(13))-Cq(9)*Spq(14)
        W(6)=-((Cq(10)+Cq(12))*Spq(16))-(Cq(11)+Cq(13))*Spq(17)
        W(7)=-((Cq(12)+Cq(14))*Spq(18))-Cq(13)*Spq(19)
        W(8)=-(Cq(14)*Spq(20))+Cpq(3)*Sq(2)
        W(9)=Cpq(5)*Sq(2)+Cpq(6)*Sq(4)+Cpq(8)*Sq(4)+Cpq(7)*Sq(5)
        W(10)=Cpq(9)*Sq(5)+Cpq(10)*Sq(7)+Cpq(11)*Sq(8)+Cpq(13)*Sq(8)
        W(11)=Cpq(14)*Sq(9)+Cpq(12)*(Sq(7)+Sq(9))
        W(12)=Cpq(15)*Sq(11)+Cpq(16)*Sq(12)
        W(13)=Cpq(19)*Sq(13)+Cpq(17)*(Sq(11)+Sq(13))
        W(14)=Cpq(20)*Sq(14)+Cpq(18)*(Sq(12)+Sq(14))
        W(15)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        W(16)=W(8)+W(9)+W(10)+W(11)+W(12)+W(13)+W(14)
        STwo=W(15)+W(16)
        Sp=Sp+STwo*Two
      RETURN
      END
      SUBROUTINE IrRegular104(PQx,PQy,PQz,Cpq,Spq)
         IMPLICIT NONE
         REAL*8 Cpq(0:20),Spq(0:20),LegendreP(0:20)
         REAL*8 Sine(0:5),Cosine(0:5),CoFact(2:5)
         REAL*8 PQx,PQy,PQz,PQ2,pqx2,pqy2,pqxy,PQ
         REAL*8 Zero,One,Two,TT,TwoC,OneOvPQ
         REAL*8 CoTan,OneOvPQxy,RS,SQ,PQToThMnsL
         PARAMETER (Two=2D0,One=1D0,Zero=0D0)
         PQx2=PQx*PQx
         PQy2=PQy*PQy
         PQ=DSQRT(PQx2+PQy2+PQz*PQz)
         OneOvPQ=One/PQ
         CoTan=PQz*OneOvPQ
!        Compute Cosine and Sine by recursion
         Sine(0)  =Zero
         Cosine(0)=One
         PQxy=DSQRT(PQx2+PQy2)
         IF(PQxy/=Zero)THEN
            OneOvPQxy=One/PQxy
            Sine(1)  =PQx*OneOvPQxy
            Cosine(1)=PQy*OneOvPQxy
         ELSE
            Sine(1)  =0.70710678118654752D0
            Cosine(1)=0.70710678118654752D0
         ENDIF
         TwoC=Two*Cosine(1)
         Sine(2)=TwoC*Sine(1)-Sine(0)
         Cosine(2)=TwoC*Cosine(1)-Cosine(0)
         Sine(3)=TwoC*Sine(2)-Sine(1)
         Cosine(3)=TwoC*Cosine(2)-Cosine(1)
         Sine(4)=TwoC*Sine(3)-Sine(2)
         Cosine(4)=TwoC*Cosine(3)-Cosine(2)
         Sine(5)=TwoC*Sine(4)-Sine(3)
         Cosine(5)=TwoC*Cosine(4)-Cosine(3)
!        ----------------------------------------------------------------
!        Associated Legendre Polynomials by recursion
!
!        M=0,Ell; Compute LegendreP(M,M)
         RS=One
         Sq=SQRT(ABS(One-CoTan**2))
         LegendreP(0)=+1.000000000000000D+00*RS
         RS=RS*Sq
         LegendreP(2)=-1.000000000000000D+00*RS
         RS=RS*Sq
         LegendreP(5)=+3.000000000000000D+00*RS
         RS=RS*Sq
         LegendreP(9)=-1.500000000000000D+01*RS
         RS=RS*Sq
         LegendreP(14)=+1.050000000000000D+02*RS
         RS=RS*Sq
         LegendreP(20)=-9.450000000000000D+02*RS
!        ----------------------------------------------------------------
!        M=0,Ell-1; Compute LegendreP(M+1,M)
         LegendreP(1)=+1.000000000000000D+00*CoTan*LegendreP(0)
         LegendreP(4)=+3.000000000000000D+00*CoTan*LegendreP(2)
         LegendreP(8)=+5.000000000000000D+00*CoTan*LegendreP(5)
         LegendreP(13)=+7.000000000000000D+00*CoTan*LegendreP(9)
         LegendreP(19)=+9.000000000000000D+00*CoTan*LegendreP(14)
!        ----------------------------------------------------------------
         CoFact(2)=+3.000000000000000D+00*CoTan
         CoFact(3)=+5.000000000000000D+00*CoTan
         CoFact(4)=+7.000000000000000D+00*CoTan
         CoFact(5)=+9.000000000000000D+00*CoTan
!        ----------------------------------------------------------------
!        M=0,Ell-2; L=M+2,Ell; Compute LegendreP(L,M)
         LegendreP(3)=CoFact(2)*LegendreP(1)-1.000000000000000D+00*LegendreP(0)
         LegendreP(6)=CoFact(3)*LegendreP(3)-4.000000000000000D+00*LegendreP(1)
         LegendreP(10)=CoFact(4)*LegendreP(6)-9.000000000000000D+00*LegendreP(3)
         LegendreP(15)=CoFact(5)*LegendreP(10)-1.600000000000000D+01*LegendreP(6)
         LegendreP(7)=CoFact(3)*LegendreP(4)-3.000000000000000D+00*LegendreP(2)
         LegendreP(11)=CoFact(4)*LegendreP(7)-8.000000000000000D+00*LegendreP(4)
         LegendreP(16)=CoFact(5)*LegendreP(11)-1.500000000000000D+01*LegendreP(7)
         LegendreP(12)=CoFact(4)*LegendreP(8)-5.000000000000000D+00*LegendreP(5)
         LegendreP(17)=CoFact(5)*LegendreP(12)-1.200000000000000D+01*LegendreP(8)
         LegendreP(18)=CoFact(5)*LegendreP(13)-7.000000000000000D+00*LegendreP(9)
!        ----------------------------------------------------------------
!        Compute Cpq(L,M) and Spq(L,M)
         PQToThMnsL=OneOvPQ
         TT=PQToThMnsL*LegendreP(0)
         Spq(0)=TT*Sine(0)
         Cpq(0)=TT*Cosine(0)
         PQToThMnsL=PQToThMnsL*OneOvPQ
         TT=PQToThMnsL*LegendreP(1)
         Spq(1)=TT*Sine(0)
         Cpq(1)=TT*Cosine(0)
         TT=PQToThMnsL*LegendreP(2)
         Spq(2)=TT*Sine(1)
         Cpq(2)=TT*Cosine(1)
         PQToThMnsL=PQToThMnsL*OneOvPQ
         TT=PQToThMnsL*LegendreP(3)
         Spq(3)=TT*Sine(0)
         Cpq(3)=TT*Cosine(0)
         TT=PQToThMnsL*LegendreP(4)
         Spq(4)=TT*Sine(1)
         Cpq(4)=TT*Cosine(1)
         TT=PQToThMnsL*LegendreP(5)
         Spq(5)=TT*Sine(2)
         Cpq(5)=TT*Cosine(2)
         PQToThMnsL=PQToThMnsL*OneOvPQ
         TT=PQToThMnsL*LegendreP(6)
         Spq(6)=TT*Sine(0)
         Cpq(6)=TT*Cosine(0)
         TT=PQToThMnsL*LegendreP(7)
         Spq(7)=TT*Sine(1)
         Cpq(7)=TT*Cosine(1)
         TT=PQToThMnsL*LegendreP(8)
         Spq(8)=TT*Sine(2)
         Cpq(8)=TT*Cosine(2)
         TT=PQToThMnsL*LegendreP(9)
         Spq(9)=TT*Sine(3)
         Cpq(9)=TT*Cosine(3)
         PQToThMnsL=PQToThMnsL*OneOvPQ
         TT=PQToThMnsL*LegendreP(10)
         Spq(10)=TT*Sine(0)
         Cpq(10)=TT*Cosine(0)
         TT=PQToThMnsL*LegendreP(11)
         Spq(11)=TT*Sine(1)
         Cpq(11)=TT*Cosine(1)
         TT=PQToThMnsL*LegendreP(12)
         Spq(12)=TT*Sine(2)
         Cpq(12)=TT*Cosine(2)
         TT=PQToThMnsL*LegendreP(13)
         Spq(13)=TT*Sine(3)
         Cpq(13)=TT*Cosine(3)
         TT=PQToThMnsL*LegendreP(14)
         Spq(14)=TT*Sine(4)
         Cpq(14)=TT*Cosine(4)
         PQToThMnsL=PQToThMnsL*OneOvPQ
         TT=PQToThMnsL*LegendreP(15)
         Spq(15)=TT*Sine(0)
         Cpq(15)=TT*Cosine(0)
         TT=PQToThMnsL*LegendreP(16)
         Spq(16)=TT*Sine(1)
         Cpq(16)=TT*Cosine(1)
         TT=PQToThMnsL*LegendreP(17)
         Spq(17)=TT*Sine(2)
         Cpq(17)=TT*Cosine(2)
         TT=PQToThMnsL*LegendreP(18)
         Spq(18)=TT*Sine(3)
         Cpq(18)=TT*Cosine(3)
         TT=PQToThMnsL*LegendreP(19)
         Spq(19)=TT*Sine(4)
         Cpq(19)=TT*Cosine(4)
         TT=PQToThMnsL*LegendreP(20)
         Spq(20)=TT*Sine(5)
         Cpq(20)=TT*Cosine(5)
      RETURN
      END
      SUBROUTINE CTraX104(PQx,PQy,PQz,Cq,Sq,Cp,Sp)
         IMPLICIT NONE
         REAL*8 PQx,PQy,PQz
         REAL*8 Cq(0:14),Sq(0:14),Cp(0:2),Sp(0:2)
         REAL*8 Cpq(0:20),Spq(0:20)
         CALL IrRegular104(PQx,PQy,PQz,Cpq,Spq)
         CALL CTX104_0_0_C(Cp(0),Cpq,Spq,Cq,Sq)
         CALL CTX104_0_1_C(Cp(1),Cpq,Spq,Cq,Sq)
         CALL CTX104_1_1_C(Cp(2),Cpq,Spq,Cq,Sq)
         CALL CTX104_1_1_S(Sp(2),Cpq,Spq,Cq,Sq)
      RETURN
      END
      SUBROUTINE CTX204_0_0_C(Cp,Cpq,Spq,Cq,Sq)
        IMPLICIT NONE
        REAL*8 Cp,Cpq(0:27),Spq(0:27),Cq(0:14),Sq(0:14)
        REAL*8 COne,CTwo,Two,W
        COMMON /CTraXScratch/ W(1000)
        PARAMETER (Two=2D0)
        W(1)=Cpq(0)*Cq(0)+Cpq(1)*Cq(1)
        W(2)=Cpq(3)*Cq(3)+Cpq(6)*Cq(6)+Cpq(10)*Cq(10)
        COne=W(1)+W(2)
        W(1)=Cpq(2)*Cq(2)+Cpq(4)*Cq(4)
        W(2)=Cpq(5)*Cq(5)+Cpq(7)*Cq(7)+Cpq(8)*Cq(8)
        W(3)=Cpq(9)*Cq(9)+Cpq(11)*Cq(11)
        W(4)=Cpq(12)*Cq(12)+Cpq(13)*Cq(13)+Cpq(14)*Cq(14)
        W(5)=Spq(2)*Sq(2)+Spq(4)*Sq(4)
        W(6)=Spq(5)*Sq(5)+Spq(7)*Sq(7)+Spq(8)*Sq(8)
        W(7)=Spq(9)*Sq(9)+Spq(11)*Sq(11)
        W(8)=Spq(12)*Sq(12)+Spq(13)*Sq(13)+Spq(14)*Sq(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        Cp=COne+Cp+CTwo*Two
      RETURN
      END
      SUBROUTINE CTX204_0_1_C(Cp,Cpq,Spq,Cq,Sq)
        IMPLICIT NONE
        REAL*8 Cp,Cpq(0:27),Spq(0:27),Cq(0:14),Sq(0:14)
        REAL*8 COne,CTwo,Two,W
        COMMON /CTraXScratch/ W(1000)
        PARAMETER (Two=2D0)
        W(1)=-(Cpq(1)*Cq(0))-Cpq(3)*Cq(1)
        W(2)=-(Cpq(6)*Cq(3))-Cpq(10)*Cq(6)-Cpq(15)*Cq(10)
        COne=W(1)+W(2)
        W(1)=-(Cpq(4)*Cq(2))-Cpq(7)*Cq(4)
        W(2)=-(Cpq(8)*Cq(5))-Cpq(11)*Cq(7)-Cpq(12)*Cq(8)
        W(3)=-(Cpq(13)*Cq(9))-Cpq(16)*Cq(11)
        W(4)=-(Cpq(17)*Cq(12))-Cpq(18)*Cq(13)-Cpq(19)*Cq(14)
        W(5)=-(Spq(4)*Sq(2))-Spq(7)*Sq(4)
        W(6)=-(Spq(8)*Sq(5))-Spq(11)*Sq(7)-Spq(12)*Sq(8)
        W(7)=-(Spq(13)*Sq(9))-Spq(16)*Sq(11)
        W(8)=-(Spq(17)*Sq(12))-Spq(18)*Sq(13)-Spq(19)*Sq(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        Cp=COne+Cp+CTwo*Two
      RETURN
      END
      SUBROUTINE CTX204_1_1_C(Cp,Cpq,Spq,Cq,Sq)
        IMPLICIT NONE
        REAL*8 Cp,Cpq(0:27),Spq(0:27),Cq(0:14),Sq(0:14)
        REAL*8 COne,CTwo,Two,W
        COMMON /CTraXScratch/ W(1000)
        PARAMETER (Two=2D0)
        W(1)=-(Cpq(2)*Cq(0))-Cpq(4)*Cq(1)
        W(2)=Cpq(3)*Cq(2)-Cpq(5)*Cq(2)
        W(3)=Cpq(6)*Cq(4)-Cpq(8)*Cq(4)
        W(4)=-(Cpq(9)*Cq(5))+Cpq(7)*(-Cq(3)+Cq(5))
        W(5)=Cpq(10)*Cq(7)-Cpq(13)*Cq(8)
        W(6)=Cpq(11)*(-Cq(6)+Cq(8))-Cpq(14)*Cq(9)
        W(7)=Cpq(12)*(-Cq(7)+Cq(9))+Cpq(15)*Cq(11)
        W(8)=Cpq(16)*(-Cq(10)+Cq(12))-Cpq(19)*Cq(13)
        W(9)=Cpq(17)*(-Cq(11)+Cq(13))-Cpq(20)*Cq(14)
        W(10)=Cpq(18)*(-Cq(12)+Cq(14))-Spq(5)*Sq(2)
        W(11)=-(Spq(8)*Sq(4))+Spq(7)*Sq(5)
        W(12)=-(Spq(9)*Sq(5))+Spq(11)*Sq(8)
        W(13)=-(Spq(13)*Sq(8))-Spq(14)*Sq(9)
        W(14)=Spq(12)*(-Sq(7)+Sq(9))+Spq(16)*Sq(12)
        W(15)=-(Spq(19)*Sq(13))+Spq(17)*(-Sq(11)+Sq(13))
        W(16)=-(Spq(20)*Sq(14))+Spq(18)*(-Sq(12)+Sq(14))
        W(17)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(18)=W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)+W(16)
        CTwo=W(17)+W(18)
        Cp=Cp+CTwo*Two
      RETURN
      END
      SUBROUTINE CTX204_1_1_S(Sp,Cpq,Spq,Cq,Sq)
        IMPLICIT NONE
        REAL*8 Sp,Cpq(0:27),Spq(0:27),Cq(0:14),Sq(0:14)
        REAL*8 SOne,STwo,Two,W
        COMMON /CTraXScratch/ W(1000)
        PARAMETER (Two=2D0)
        W(1)=-(Cq(0)*Spq(2))-Cq(1)*Spq(4)
        W(2)=-(Cq(2)*Spq(5))-(Cq(3)+Cq(5))*Spq(7)
        W(3)=-(Cq(4)*Spq(8))-Cq(5)*Spq(9)
        W(4)=-((Cq(6)+Cq(8))*Spq(11))-(Cq(7)+Cq(9))*Spq(12)
        W(5)=-(Cq(8)*Spq(13))-Cq(9)*Spq(14)
        W(6)=-((Cq(10)+Cq(12))*Spq(16))-(Cq(11)+Cq(13))*Spq(17)
        W(7)=-((Cq(12)+Cq(14))*Spq(18))-Cq(13)*Spq(19)
        W(8)=-(Cq(14)*Spq(20))+Cpq(3)*Sq(2)
        W(9)=Cpq(5)*Sq(2)+Cpq(6)*Sq(4)+Cpq(8)*Sq(4)+Cpq(7)*Sq(5)
        W(10)=Cpq(9)*Sq(5)+Cpq(10)*Sq(7)+Cpq(11)*Sq(8)+Cpq(13)*Sq(8)
        W(11)=Cpq(14)*Sq(9)+Cpq(12)*(Sq(7)+Sq(9))
        W(12)=Cpq(15)*Sq(11)+Cpq(16)*Sq(12)
        W(13)=Cpq(19)*Sq(13)+Cpq(17)*(Sq(11)+Sq(13))
        W(14)=Cpq(20)*Sq(14)+Cpq(18)*(Sq(12)+Sq(14))
        W(15)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        W(16)=W(8)+W(9)+W(10)+W(11)+W(12)+W(13)+W(14)
        STwo=W(15)+W(16)
        Sp=Sp+STwo*Two
      RETURN
      END
      SUBROUTINE CTX204_0_2_C(Cp,Cpq,Spq,Cq,Sq)
        IMPLICIT NONE
        REAL*8 Cp,Cpq(0:27),Spq(0:27),Cq(0:14),Sq(0:14)
        REAL*8 COne,CTwo,Two,W
        COMMON /CTraXScratch/ W(1000)
        PARAMETER (Two=2D0)
        W(1)=Cpq(3)*Cq(0)+Cpq(6)*Cq(1)
        W(2)=Cpq(10)*Cq(3)+Cpq(15)*Cq(6)+Cpq(21)*Cq(10)
        COne=W(1)+W(2)
        W(1)=Cpq(7)*Cq(2)+Cpq(11)*Cq(4)
        W(2)=Cpq(12)*Cq(5)+Cpq(16)*Cq(7)+Cpq(17)*Cq(8)
        W(3)=Cpq(18)*Cq(9)+Cpq(22)*Cq(11)
        W(4)=Cpq(23)*Cq(12)+Cpq(24)*Cq(13)+Cpq(25)*Cq(14)
        W(5)=Spq(7)*Sq(2)+Spq(11)*Sq(4)
        W(6)=Spq(12)*Sq(5)+Spq(16)*Sq(7)+Spq(17)*Sq(8)
        W(7)=Spq(18)*Sq(9)+Spq(22)*Sq(11)
        W(8)=Spq(23)*Sq(12)+Spq(24)*Sq(13)+Spq(25)*Sq(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        Cp=COne+Cp+CTwo*Two
      RETURN
      END
      SUBROUTINE CTX204_1_2_C(Cp,Cpq,Spq,Cq,Sq)
        IMPLICIT NONE
        REAL*8 Cp,Cpq(0:27),Spq(0:27),Cq(0:14),Sq(0:14)
        REAL*8 COne,CTwo,Two,W
        COMMON /CTraXScratch/ W(1000)
        PARAMETER (Two=2D0)
        W(1)=Cpq(4)*Cq(0)+Cpq(7)*Cq(1)
        W(2)=-(Cpq(6)*Cq(2))+Cpq(8)*Cq(2)
        W(3)=-(Cpq(10)*Cq(4))+Cpq(12)*Cq(4)
        W(4)=Cpq(11)*(Cq(3)-Cq(5))+Cpq(13)*Cq(5)
        W(5)=-(Cpq(15)*Cq(7))+Cpq(16)*(Cq(6)-Cq(8))
        W(6)=Cpq(18)*Cq(8)+Cpq(17)*(Cq(7)-Cq(9))
        W(7)=Cpq(19)*Cq(9)-Cpq(21)*Cq(11)
        W(8)=Cpq(22)*(Cq(10)-Cq(12))+Cpq(23)*(Cq(11)-Cq(13))
        W(9)=Cpq(25)*Cq(13)+Cpq(24)*(Cq(12)-Cq(14))
        W(10)=Cpq(26)*Cq(14)+Spq(8)*Sq(2)
        W(11)=Spq(12)*Sq(4)-Spq(11)*Sq(5)
        W(12)=Spq(13)*Sq(5)-Spq(16)*Sq(8)
        W(13)=Spq(18)*Sq(8)+Spq(17)*(Sq(7)-Sq(9))
        W(14)=Spq(19)*Sq(9)-Spq(22)*Sq(12)
        W(15)=Spq(23)*(Sq(11)-Sq(13))+Spq(25)*Sq(13)
        W(16)=Spq(24)*(Sq(12)-Sq(14))+Spq(26)*Sq(14)
        W(17)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(18)=W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)+W(16)
        CTwo=W(17)+W(18)
        Cp=Cp+CTwo*Two
      RETURN
      END
      SUBROUTINE CTX204_1_2_S(Sp,Cpq,Spq,Cq,Sq)
        IMPLICIT NONE
        REAL*8 Sp,Cpq(0:27),Spq(0:27),Cq(0:14),Sq(0:14)
        REAL*8 SOne,STwo,Two,W
        COMMON /CTraXScratch/ W(1000)
        PARAMETER (Two=2D0)
        W(1)=Cq(0)*Spq(4)+Cq(1)*Spq(7)
        W(2)=Cq(2)*Spq(8)+(Cq(3)+Cq(5))*Spq(11)
        W(3)=Cq(4)*Spq(12)+Cq(5)*Spq(13)
        W(4)=(Cq(6)+Cq(8))*Spq(16)+(Cq(7)+Cq(9))*Spq(17)
        W(5)=Cq(8)*Spq(18)+Cq(9)*Spq(19)
        W(6)=(Cq(10)+Cq(12))*Spq(22)+(Cq(11)+Cq(13))*Spq(23)
        W(7)=(Cq(12)+Cq(14))*Spq(24)+Cq(13)*Spq(25)
        W(8)=Cq(14)*Spq(26)-Cpq(6)*Sq(2)
        W(9)=-(Cpq(8)*Sq(2))-Cpq(10)*Sq(4)
        W(10)=-(Cpq(12)*Sq(4))-Cpq(11)*Sq(5)
        W(11)=-(Cpq(13)*Sq(5))-Cpq(15)*Sq(7)
        W(12)=-(Cpq(16)*Sq(8))-Cpq(18)*Sq(8)
        W(13)=-(Cpq(19)*Sq(9))-Cpq(17)*(Sq(7)+Sq(9))
        W(14)=-(Cpq(21)*Sq(11))-Cpq(22)*Sq(12)
        W(15)=-(Cpq(25)*Sq(13))-Cpq(23)*(Sq(11)+Sq(13))
        W(16)=-(Cpq(26)*Sq(14))-Cpq(24)*(Sq(12)+Sq(14))
        W(17)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(18)=W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)+W(16)
        STwo=W(17)+W(18)
        Sp=Sp+STwo*Two
      RETURN
      END
      SUBROUTINE CTX204_2_2_C(Cp,Cpq,Spq,Cq,Sq)
        IMPLICIT NONE
        REAL*8 Cp,Cpq(0:27),Spq(0:27),Cq(0:14),Sq(0:14)
        REAL*8 COne,CTwo,Two,W
        COMMON /CTraXScratch/ W(1000)
        PARAMETER (Two=2D0)
        W(1)=Cpq(5)*Cq(0)+Cpq(8)*Cq(1)
        W(2)=-(Cpq(7)*Cq(2))+Cpq(9)*Cq(2)
        W(3)=Cpq(12)*Cq(3)-Cpq(11)*Cq(4)
        W(4)=Cpq(13)*Cq(4)+Cpq(10)*Cq(5)+Cpq(14)*Cq(5)
        W(5)=Cpq(17)*Cq(6)+Cpq(18)*Cq(7)+Cpq(15)*Cq(8)+Cpq(19)*Cq(8)
        W(6)=Cpq(20)*Cq(9)+Cpq(16)*(-Cq(7)+Cq(9))
        W(7)=Cpq(24)*Cq(11)+Cpq(21)*Cq(12)+Cpq(25)*Cq(12)
        W(8)=Cpq(26)*Cq(13)+Cpq(22)*(-Cq(11)+Cq(13))
        W(9)=Cpq(27)*Cq(14)+Cpq(23)*(Cq(10)+Cq(14))
        W(10)=Spq(7)*Sq(2)+Spq(9)*Sq(2)
        W(11)=Spq(11)*Sq(4)+Spq(13)*Sq(4)+Spq(14)*Sq(5)
        W(12)=Spq(18)*Sq(7)+Spq(19)*Sq(8)
        W(13)=Spq(20)*Sq(9)+Spq(16)*(Sq(7)+Sq(9))+Spq(24)*Sq(11)
        W(14)=Spq(25)*Sq(12)+Spq(26)*Sq(13)
        W(15)=Spq(22)*(Sq(11)+Sq(13))+Spq(23)*Sq(14)+Spq(27)*Sq(14)
        W(16)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        W(17)=W(8)+W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)
        CTwo=W(16)+W(17)
        Cp=Cp+CTwo*Two
      RETURN
      END
      SUBROUTINE CTX204_2_2_S(Sp,Cpq,Spq,Cq,Sq)
        IMPLICIT NONE
        REAL*8 Sp,Cpq(0:27),Spq(0:27),Cq(0:14),Sq(0:14)
        REAL*8 SOne,STwo,Two,W
        COMMON /CTraXScratch/ W(1000)
        PARAMETER (Two=2D0)
        W(1)=Cq(0)*Spq(5)-Cq(2)*Spq(7)
        W(2)=Cq(1)*Spq(8)+Cq(2)*Spq(9)
        W(3)=-(Cq(4)*Spq(11))+Cq(3)*Spq(12)
        W(4)=Cq(4)*Spq(13)+Cq(5)*Spq(14)-(Cq(7)+Cq(9))*Spq(16)
        W(5)=Cq(6)*Spq(17)+Cq(7)*Spq(18)+Cq(8)*Spq(19)+Cq(9)*Spq(20)
        W(6)=-((Cq(11)+Cq(13))*Spq(22))+(Cq(10)-Cq(14))*Spq(23)
        W(7)=Cq(11)*Spq(24)+Cq(12)*Spq(25)+Cq(13)*Spq(26)
        W(8)=Cq(14)*Spq(27)-Cpq(7)*Sq(2)
        W(9)=-(Cpq(9)*Sq(2))-Cpq(11)*Sq(4)
        W(10)=-(Cpq(13)*Sq(4))+Cpq(10)*Sq(5)
        W(11)=-(Cpq(14)*Sq(5))-Cpq(18)*Sq(7)+Cpq(15)*Sq(8)
        W(12)=-(Cpq(19)*Sq(8))-Cpq(20)*Sq(9)
        W(13)=Cpq(16)*(-Sq(7)+Sq(9))
        W(14)=-(Cpq(24)*Sq(11))+Cpq(21)*Sq(12)
        W(15)=-(Cpq(25)*Sq(12))-Cpq(26)*Sq(13)
        W(16)=Cpq(22)*(-Sq(11)+Sq(13))
        W(17)=Cpq(23)*Sq(14)-Cpq(27)*Sq(14)
        W(18)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(19)=W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)+W(16)+W(17)
        STwo=W(18)+W(19)
        Sp=Sp+STwo*Two
      RETURN
      END
      SUBROUTINE IrRegular204(PQx,PQy,PQz,Cpq,Spq)
         IMPLICIT NONE
         REAL*8 Cpq(0:27),Spq(0:27),LegendreP(0:27)
         REAL*8 Sine(0:6),Cosine(0:6),CoFact(2:6)
         REAL*8 PQx,PQy,PQz,PQ2,pqx2,pqy2,pqxy,PQ
         REAL*8 Zero,One,Two,TT,TwoC,OneOvPQ
         REAL*8 CoTan,OneOvPQxy,RS,SQ,PQToThMnsL
         PARAMETER (Two=2D0,One=1D0,Zero=0D0)
         PQx2=PQx*PQx
         PQy2=PQy*PQy
         PQ=DSQRT(PQx2+PQy2+PQz*PQz)
         OneOvPQ=One/PQ
         CoTan=PQz*OneOvPQ
!        Compute Cosine and Sine by recursion
         Sine(0)  =Zero
         Cosine(0)=One
         PQxy=DSQRT(PQx2+PQy2)
         IF(PQxy/=Zero)THEN
            OneOvPQxy=One/PQxy
            Sine(1)  =PQx*OneOvPQxy
            Cosine(1)=PQy*OneOvPQxy
         ELSE
            Sine(1)  =0.70710678118654752D0
            Cosine(1)=0.70710678118654752D0
         ENDIF
         TwoC=Two*Cosine(1)
         Sine(2)=TwoC*Sine(1)-Sine(0)
         Cosine(2)=TwoC*Cosine(1)-Cosine(0)
         Sine(3)=TwoC*Sine(2)-Sine(1)
         Cosine(3)=TwoC*Cosine(2)-Cosine(1)
         Sine(4)=TwoC*Sine(3)-Sine(2)
         Cosine(4)=TwoC*Cosine(3)-Cosine(2)
         Sine(5)=TwoC*Sine(4)-Sine(3)
         Cosine(5)=TwoC*Cosine(4)-Cosine(3)
         Sine(6)=TwoC*Sine(5)-Sine(4)
         Cosine(6)=TwoC*Cosine(5)-Cosine(4)
!        ----------------------------------------------------------------
!        Associated Legendre Polynomials by recursion
!
!        M=0,Ell; Compute LegendreP(M,M)
         RS=One
         Sq=SQRT(ABS(One-CoTan**2))
         LegendreP(0)=+1.000000000000000D+00*RS
         RS=RS*Sq
         LegendreP(2)=-1.000000000000000D+00*RS
         RS=RS*Sq
         LegendreP(5)=+3.000000000000000D+00*RS
         RS=RS*Sq
         LegendreP(9)=-1.500000000000000D+01*RS
         RS=RS*Sq
         LegendreP(14)=+1.050000000000000D+02*RS
         RS=RS*Sq
         LegendreP(20)=-9.450000000000000D+02*RS
         RS=RS*Sq
         LegendreP(27)=+1.039500000000000D+04*RS
!        ----------------------------------------------------------------
!        M=0,Ell-1; Compute LegendreP(M+1,M)
         LegendreP(1)=+1.000000000000000D+00*CoTan*LegendreP(0)
         LegendreP(4)=+3.000000000000000D+00*CoTan*LegendreP(2)
         LegendreP(8)=+5.000000000000000D+00*CoTan*LegendreP(5)
         LegendreP(13)=+7.000000000000000D+00*CoTan*LegendreP(9)
         LegendreP(19)=+9.000000000000000D+00*CoTan*LegendreP(14)
         LegendreP(26)=+1.100000000000000D+01*CoTan*LegendreP(20)
!        ----------------------------------------------------------------
         CoFact(2)=+3.000000000000000D+00*CoTan
         CoFact(3)=+5.000000000000000D+00*CoTan
         CoFact(4)=+7.000000000000000D+00*CoTan
         CoFact(5)=+9.000000000000000D+00*CoTan
         CoFact(6)=+1.100000000000000D+01*CoTan
!        ----------------------------------------------------------------
!        M=0,Ell-2; L=M+2,Ell; Compute LegendreP(L,M)
         LegendreP(3)=CoFact(2)*LegendreP(1)-1.000000000000000D+00*LegendreP(0)
         LegendreP(6)=CoFact(3)*LegendreP(3)-4.000000000000000D+00*LegendreP(1)
         LegendreP(10)=CoFact(4)*LegendreP(6)-9.000000000000000D+00*LegendreP(3)
         LegendreP(15)=CoFact(5)*LegendreP(10)-1.600000000000000D+01*LegendreP(6)
         LegendreP(21)=CoFact(6)*LegendreP(15)-2.500000000000000D+01*LegendreP(10)
         LegendreP(7)=CoFact(3)*LegendreP(4)-3.000000000000000D+00*LegendreP(2)
         LegendreP(11)=CoFact(4)*LegendreP(7)-8.000000000000000D+00*LegendreP(4)
         LegendreP(16)=CoFact(5)*LegendreP(11)-1.500000000000000D+01*LegendreP(7)
         LegendreP(22)=CoFact(6)*LegendreP(16)-2.400000000000000D+01*LegendreP(11)
         LegendreP(12)=CoFact(4)*LegendreP(8)-5.000000000000000D+00*LegendreP(5)
         LegendreP(17)=CoFact(5)*LegendreP(12)-1.200000000000000D+01*LegendreP(8)
         LegendreP(23)=CoFact(6)*LegendreP(17)-2.100000000000000D+01*LegendreP(12)
         LegendreP(18)=CoFact(5)*LegendreP(13)-7.000000000000000D+00*LegendreP(9)
         LegendreP(24)=CoFact(6)*LegendreP(18)-1.600000000000000D+01*LegendreP(13)
         LegendreP(25)=CoFact(6)*LegendreP(19)-9.000000000000000D+00*LegendreP(14)
!        ----------------------------------------------------------------
!        Compute Cpq(L,M) and Spq(L,M)
         PQToThMnsL=OneOvPQ
         TT=PQToThMnsL*LegendreP(0)
         Spq(0)=TT*Sine(0)
         Cpq(0)=TT*Cosine(0)
         PQToThMnsL=PQToThMnsL*OneOvPQ
         TT=PQToThMnsL*LegendreP(1)
         Spq(1)=TT*Sine(0)
         Cpq(1)=TT*Cosine(0)
         TT=PQToThMnsL*LegendreP(2)
         Spq(2)=TT*Sine(1)
         Cpq(2)=TT*Cosine(1)
         PQToThMnsL=PQToThMnsL*OneOvPQ
         TT=PQToThMnsL*LegendreP(3)
         Spq(3)=TT*Sine(0)
         Cpq(3)=TT*Cosine(0)
         TT=PQToThMnsL*LegendreP(4)
         Spq(4)=TT*Sine(1)
         Cpq(4)=TT*Cosine(1)
         TT=PQToThMnsL*LegendreP(5)
         Spq(5)=TT*Sine(2)
         Cpq(5)=TT*Cosine(2)
         PQToThMnsL=PQToThMnsL*OneOvPQ
         TT=PQToThMnsL*LegendreP(6)
         Spq(6)=TT*Sine(0)
         Cpq(6)=TT*Cosine(0)
         TT=PQToThMnsL*LegendreP(7)
         Spq(7)=TT*Sine(1)
         Cpq(7)=TT*Cosine(1)
         TT=PQToThMnsL*LegendreP(8)
         Spq(8)=TT*Sine(2)
         Cpq(8)=TT*Cosine(2)
         TT=PQToThMnsL*LegendreP(9)
         Spq(9)=TT*Sine(3)
         Cpq(9)=TT*Cosine(3)
         PQToThMnsL=PQToThMnsL*OneOvPQ
         TT=PQToThMnsL*LegendreP(10)
         Spq(10)=TT*Sine(0)
         Cpq(10)=TT*Cosine(0)
         TT=PQToThMnsL*LegendreP(11)
         Spq(11)=TT*Sine(1)
         Cpq(11)=TT*Cosine(1)
         TT=PQToThMnsL*LegendreP(12)
         Spq(12)=TT*Sine(2)
         Cpq(12)=TT*Cosine(2)
         TT=PQToThMnsL*LegendreP(13)
         Spq(13)=TT*Sine(3)
         Cpq(13)=TT*Cosine(3)
         TT=PQToThMnsL*LegendreP(14)
         Spq(14)=TT*Sine(4)
         Cpq(14)=TT*Cosine(4)
         PQToThMnsL=PQToThMnsL*OneOvPQ
         TT=PQToThMnsL*LegendreP(15)
         Spq(15)=TT*Sine(0)
         Cpq(15)=TT*Cosine(0)
         TT=PQToThMnsL*LegendreP(16)
         Spq(16)=TT*Sine(1)
         Cpq(16)=TT*Cosine(1)
         TT=PQToThMnsL*LegendreP(17)
         Spq(17)=TT*Sine(2)
         Cpq(17)=TT*Cosine(2)
         TT=PQToThMnsL*LegendreP(18)
         Spq(18)=TT*Sine(3)
         Cpq(18)=TT*Cosine(3)
         TT=PQToThMnsL*LegendreP(19)
         Spq(19)=TT*Sine(4)
         Cpq(19)=TT*Cosine(4)
         TT=PQToThMnsL*LegendreP(20)
         Spq(20)=TT*Sine(5)
         Cpq(20)=TT*Cosine(5)
         PQToThMnsL=PQToThMnsL*OneOvPQ
         TT=PQToThMnsL*LegendreP(21)
         Spq(21)=TT*Sine(0)
         Cpq(21)=TT*Cosine(0)
         TT=PQToThMnsL*LegendreP(22)
         Spq(22)=TT*Sine(1)
         Cpq(22)=TT*Cosine(1)
         TT=PQToThMnsL*LegendreP(23)
         Spq(23)=TT*Sine(2)
         Cpq(23)=TT*Cosine(2)
         TT=PQToThMnsL*LegendreP(24)
         Spq(24)=TT*Sine(3)
         Cpq(24)=TT*Cosine(3)
         TT=PQToThMnsL*LegendreP(25)
         Spq(25)=TT*Sine(4)
         Cpq(25)=TT*Cosine(4)
         TT=PQToThMnsL*LegendreP(26)
         Spq(26)=TT*Sine(5)
         Cpq(26)=TT*Cosine(5)
         TT=PQToThMnsL*LegendreP(27)
         Spq(27)=TT*Sine(6)
         Cpq(27)=TT*Cosine(6)
      RETURN
      END
      SUBROUTINE CTraX204(PQx,PQy,PQz,Cq,Sq,Cp,Sp)
         IMPLICIT NONE
         REAL*8 PQx,PQy,PQz
         REAL*8 Cq(0:14),Sq(0:14),Cp(0:5),Sp(0:5)
         REAL*8 Cpq(0:27),Spq(0:27)
         CALL IrRegular204(PQx,PQy,PQz,Cpq,Spq)
         CALL CTX204_0_0_C(Cp(0),Cpq,Spq,Cq,Sq)
         CALL CTX204_0_1_C(Cp(1),Cpq,Spq,Cq,Sq)
         CALL CTX204_1_1_C(Cp(2),Cpq,Spq,Cq,Sq)
         CALL CTX204_1_1_S(Sp(2),Cpq,Spq,Cq,Sq)
         CALL CTX204_0_2_C(Cp(3),Cpq,Spq,Cq,Sq)
         CALL CTX204_1_2_C(Cp(4),Cpq,Spq,Cq,Sq)
         CALL CTX204_1_2_S(Sp(4),Cpq,Spq,Cq,Sq)
         CALL CTX204_2_2_C(Cp(5),Cpq,Spq,Cq,Sq)
         CALL CTX204_2_2_S(Sp(5),Cpq,Spq,Cq,Sq)
      RETURN
      END
