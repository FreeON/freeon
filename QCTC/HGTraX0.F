      SUBROUTINE MD3TRR__0_0(PQx,PQy,PQz,AuxR,R)
      IMPLICIT NONE
      REAL*8 AuxR(0:0)
      REAL*8 R(1)
      REAL*8 PQx,PQy,PQz
      R(1)=AuxR(0)
      RETURN
      END
      SUBROUTINE KTrax__0_0(R,Co,Ket)
      IMPLICIT NONE
      REAL*8 R(1)
      REAL*8 Co(1)
      REAL*8 Ket(1)
      Ket(1)=Ket(1)+R(1)*Co(1)
      RETURN
      END
      SUBROUTINE MD3TRR__0_1(PQx,PQy,PQz,AuxR,R)
      IMPLICIT NONE
      REAL*8 AuxR(0:1)
      REAL*8 R(4)
      REAL*8 PQx,PQy,PQz
      R(1)=AuxR(1)
      R(2)=PQx*R(1)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(0)
      RETURN
      END
      SUBROUTINE KTrax__0_1(R,Co,Ket)
      IMPLICIT NONE
      REAL*8 R(4)
      REAL*8 Co(4)
      REAL*8 Ket(1)
      Ket(1)=Ket(1)+R(1)*Co(1)
      Ket(1)=Ket(1)+R(2)*Co(2)
      Ket(1)=Ket(1)+R(3)*Co(3)
      Ket(1)=Ket(1)+R(4)*Co(4)
      RETURN
      END
      SUBROUTINE MD3TRR__0_2(PQx,PQy,PQz,AuxR,R)
      IMPLICIT NONE
      REAL*8 AuxR(0:2)
      REAL*8 R(10)
      REAL*8 PQx,PQy,PQz
      R(1)=AuxR(2)
      R(2)=PQx*R(1)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(1)
      R(5)=PQx*R(2)+R(1)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(6)=PQx*R(3)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(0)
      RETURN
      END
      SUBROUTINE KTrax__0_2(R,Co,Ket)
      IMPLICIT NONE
      REAL*8 R(10)
      REAL*8 Co(10)
      REAL*8 Ket(1)
      Ket(1)=Ket(1)+R(1)*Co(1)
      Ket(1)=Ket(1)+R(2)*Co(2)
      Ket(1)=Ket(1)+R(3)*Co(3)
      Ket(1)=Ket(1)+R(4)*Co(4)
      Ket(1)=Ket(1)+R(5)*Co(5)
      Ket(1)=Ket(1)+R(6)*Co(6)
      Ket(1)=Ket(1)+R(7)*Co(7)
      Ket(1)=Ket(1)+R(8)*Co(8)
      Ket(1)=Ket(1)+R(9)*Co(9)
      Ket(1)=Ket(1)+R(10)*Co(10)
      RETURN
      END
      SUBROUTINE MD3TRR__0_3(PQx,PQy,PQz,AuxR,R)
      IMPLICIT NONE
      REAL*8 AuxR(0:3)
      REAL*8 R(20)
      REAL*8 PQx,PQy,PQz
      R(1)=AuxR(3)
      R(2)=PQx*R(1)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(2)
      R(5)=PQx*R(2)+R(1)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(6)=PQx*R(3)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(1)
      R(11)=PQx*R(5)+0.20D+01*R(2)
      R(14)=PQy*R(7)+0.20D+01*R(3)
      R(20)=PQz*R(10)+0.20D+01*R(4)
      R(5)=PQx*R(2)+R(1)
      R(13)=PQx*R(7)
      R(19)=PQy*R(10)
      R(18)=PQx*R(10)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(12)=PQx*R(6)+R(3)
      R(17)=PQy*R(9)+R(4)
      R(15)=PQx*R(8)+R(4)
      R(6)=PQx*R(3)
      R(16)=PQx*R(9)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(0)
      RETURN
      END
      SUBROUTINE KTrax__0_3(R,Co,Ket)
      IMPLICIT NONE
      REAL*8 R(20)
      REAL*8 Co(20)
      REAL*8 Ket(1)
      Ket(1)=Ket(1)+R(1)*Co(1)
      Ket(1)=Ket(1)+R(2)*Co(2)
      Ket(1)=Ket(1)+R(3)*Co(3)
      Ket(1)=Ket(1)+R(4)*Co(4)
      Ket(1)=Ket(1)+R(5)*Co(5)
      Ket(1)=Ket(1)+R(6)*Co(6)
      Ket(1)=Ket(1)+R(7)*Co(7)
      Ket(1)=Ket(1)+R(8)*Co(8)
      Ket(1)=Ket(1)+R(9)*Co(9)
      Ket(1)=Ket(1)+R(10)*Co(10)
      Ket(1)=Ket(1)+R(11)*Co(11)
      Ket(1)=Ket(1)+R(12)*Co(12)
      Ket(1)=Ket(1)+R(13)*Co(13)
      Ket(1)=Ket(1)+R(14)*Co(14)
      Ket(1)=Ket(1)+R(15)*Co(15)
      Ket(1)=Ket(1)+R(16)*Co(16)
      Ket(1)=Ket(1)+R(17)*Co(17)
      Ket(1)=Ket(1)+R(18)*Co(18)
      Ket(1)=Ket(1)+R(19)*Co(19)
      Ket(1)=Ket(1)+R(20)*Co(20)
      RETURN
      END
      SUBROUTINE MD3TRR__0_4(PQx,PQy,PQz,AuxR,R)
      IMPLICIT NONE
      REAL*8 AuxR(0:4)
      REAL*8 R(35)
      REAL*8 PQx,PQy,PQz
      R(1)=AuxR(4)
      R(2)=PQx*R(1)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(3)
      R(5)=PQx*R(2)+R(1)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(6)=PQx*R(3)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(2)
      R(11)=PQx*R(5)+0.20D+01*R(2)
      R(14)=PQy*R(7)+0.20D+01*R(3)
      R(20)=PQz*R(10)+0.20D+01*R(4)
      R(5)=PQx*R(2)+R(1)
      R(13)=PQx*R(7)
      R(19)=PQy*R(10)
      R(18)=PQx*R(10)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(12)=PQx*R(6)+R(3)
      R(17)=PQy*R(9)+R(4)
      R(15)=PQx*R(8)+R(4)
      R(6)=PQx*R(3)
      R(16)=PQx*R(9)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(1)
      R(21)=PQx*R(11)+0.30D+01*R(5)
      R(25)=PQy*R(14)+0.30D+01*R(7)
      R(35)=PQz*R(20)+0.30D+01*R(10)
      R(11)=PQx*R(5)+0.20D+01*R(2)
      R(24)=PQx*R(14)
      R(34)=PQy*R(20)
      R(33)=PQx*R(20)
      R(14)=PQy*R(7)+0.20D+01*R(3)
      R(20)=PQz*R(10)+0.20D+01*R(4)
      R(5)=PQx*R(2)+R(1)
      R(23)=PQx*R(13)+R(7)
      R(32)=PQy*R(19)+R(10)
      R(30)=PQx*R(18)+R(10)
      R(13)=PQx*R(7)
      R(31)=PQx*R(19)
      R(19)=PQy*R(10)
      R(18)=PQx*R(10)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(22)=PQx*R(12)+0.20D+01*R(6)
      R(29)=PQy*R(17)+0.20D+01*R(9)
      R(26)=PQx*R(15)+0.20D+01*R(8)
      R(12)=PQx*R(6)+R(3)
      R(28)=PQx*R(17)
      R(17)=PQy*R(9)+R(4)
      R(15)=PQx*R(8)+R(4)
      R(6)=PQx*R(3)
      R(27)=PQx*R(16)+R(9)
      R(16)=PQx*R(9)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(0)
      RETURN
      END
      SUBROUTINE KTrax__0_4(R,Co,Ket)
      IMPLICIT NONE
      REAL*8 R(35)
      REAL*8 Co(35)
      REAL*8 Ket(1)
      Ket(1)=Ket(1)+R(1)*Co(1)
      Ket(1)=Ket(1)+R(2)*Co(2)
      Ket(1)=Ket(1)+R(3)*Co(3)
      Ket(1)=Ket(1)+R(4)*Co(4)
      Ket(1)=Ket(1)+R(5)*Co(5)
      Ket(1)=Ket(1)+R(6)*Co(6)
      Ket(1)=Ket(1)+R(7)*Co(7)
      Ket(1)=Ket(1)+R(8)*Co(8)
      Ket(1)=Ket(1)+R(9)*Co(9)
      Ket(1)=Ket(1)+R(10)*Co(10)
      Ket(1)=Ket(1)+R(11)*Co(11)
      Ket(1)=Ket(1)+R(12)*Co(12)
      Ket(1)=Ket(1)+R(13)*Co(13)
      Ket(1)=Ket(1)+R(14)*Co(14)
      Ket(1)=Ket(1)+R(15)*Co(15)
      Ket(1)=Ket(1)+R(16)*Co(16)
      Ket(1)=Ket(1)+R(17)*Co(17)
      Ket(1)=Ket(1)+R(18)*Co(18)
      Ket(1)=Ket(1)+R(19)*Co(19)
      Ket(1)=Ket(1)+R(20)*Co(20)
      Ket(1)=Ket(1)+R(21)*Co(21)
      Ket(1)=Ket(1)+R(22)*Co(22)
      Ket(1)=Ket(1)+R(23)*Co(23)
      Ket(1)=Ket(1)+R(24)*Co(24)
      Ket(1)=Ket(1)+R(25)*Co(25)
      Ket(1)=Ket(1)+R(26)*Co(26)
      Ket(1)=Ket(1)+R(27)*Co(27)
      Ket(1)=Ket(1)+R(28)*Co(28)
      Ket(1)=Ket(1)+R(29)*Co(29)
      Ket(1)=Ket(1)+R(30)*Co(30)
      Ket(1)=Ket(1)+R(31)*Co(31)
      Ket(1)=Ket(1)+R(32)*Co(32)
      Ket(1)=Ket(1)+R(33)*Co(33)
      Ket(1)=Ket(1)+R(34)*Co(34)
      Ket(1)=Ket(1)+R(35)*Co(35)
      RETURN
      END
      SUBROUTINE HGTraX00(PQx,PQy,PQz,U,O,Co,Ket)
      IMPLICIT NONE
      REAL*8 PQx,PQy,PQz,U,O,o1,o2,OneOvT
      REAL*8 G(0:0)
      REAL*8 AuxR(0:0)
      REAL*8 R(1)
      REAL*8 Co(1)
      REAL*8 Ket(1)
      INTEGER Mesh,I,J,K
      REAL*8 ET,TwoT,T,T2,T3,T4
      REAL*8 SqrtT,SqrtPi,One,Two
      REAL*8 Switch,Grid,F0Asymp,F1Asymp,F2Asymp,F3Asymp
      REAL*8 F4Asymp,F5Asymp,F6Asymp,F7Asymp,F8Asymp
      REAL*8 F9Asymp,F10Asymp,F11Asymp,F12Asymp
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "Gamma_Asymptotics.Inc"
      INCLUDE "Mesh.Inc"
      REAL*8 F0_0(0:Mesh)
      REAL*8 F0_1(0:Mesh)
      REAL*8 F0_2(0:Mesh)
      REAL*8 F0_3(0:Mesh)
      INCLUDE "Gamma_0.Inc"
      T=O*(PQx*PQx+PQy*PQy+PQz*PQz)
      IF(T.LT.Switch)THEN
         T2=T*T
         T3=T*T2
         J=AINT(T*Grid)
         G(0)=(F0_0(J)+T*F0_1(J)+T2*F0_2(J)+T3*F0_3(J))
         o1=U
         AuxR(0)=o1*G(0)
      ELSE
         SqrtT=DSQRT(T)
         OneOvT=One/T
         G(0)=SqrtPi/(Two*SqrtT)
         o1=U
         AuxR(0)=o1*G(0)
      ENDIF
      CALL MD3TRR__0_0(PQx,PQy,PQz,AuxR,R)
      CALL KTraX__0_0(R,Co,Ket)
      RETURN
      END

      SUBROUTINE HGTraX01(PQx,PQy,PQz,U,O,Co,Ket)
      IMPLICIT NONE
      REAL*8 PQx,PQy,PQz,U,O,o1,o2,OneOvT
      REAL*8 G(0:1)
      REAL*8 AuxR(0:1)
      REAL*8 R(4)
      REAL*8 Co(4)
      REAL*8 Ket(1)
      INTEGER Mesh,I,J,K
      REAL*8 ET,TwoT,T,T2,T3,T4
      REAL*8 SqrtT,SqrtPi,One,Two
      REAL*8 Switch,Grid,F0Asymp,F1Asymp,F2Asymp,F3Asymp
      REAL*8 F4Asymp,F5Asymp,F6Asymp,F7Asymp,F8Asymp
      REAL*8 F9Asymp,F10Asymp,F11Asymp,F12Asymp
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "Gamma_Asymptotics.Inc"
      INCLUDE "Mesh.Inc"
      REAL*8 F1_0(0:Mesh)
      REAL*8 F1_1(0:Mesh)
      REAL*8 F1_2(0:Mesh)
      REAL*8 F1_3(0:Mesh)
      INCLUDE "Gamma_1.Inc"
      T=O*(PQx*PQx+PQy*PQy+PQz*PQz)
      IF(T.LT.Switch)THEN
         T2=T*T
         T3=T*T2
         J=AINT(T*Grid)
         G(1)=(F1_0(J)+T*F1_1(J)+T2*F1_2(J)+T3*F1_3(J))
         ET=DEXP(-T)
         TwoT=2.0D0*T
         G(0)=TwoT*G(1)+ET
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         AuxR(1)=o1*G(1)
      ELSE
         SqrtT=DSQRT(T)
         OneOvT=One/T
         G(0)=SqrtPi/(Two*SqrtT)
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         G(1)=G(0)*0.5000000000000000D+00*OneOvT
         AuxR(1)=o1*G(1)
      ENDIF
      CALL MD3TRR__0_1(PQx,PQy,PQz,AuxR,R)
      CALL KTraX__0_1(R,Co,Ket)
      RETURN
      END

      SUBROUTINE HGTraX02(PQx,PQy,PQz,U,O,Co,Ket)
      IMPLICIT NONE
      REAL*8 PQx,PQy,PQz,U,O,o1,o2,OneOvT
      REAL*8 G(0:2)
      REAL*8 AuxR(0:2)
      REAL*8 R(10)
      REAL*8 Co(10)
      REAL*8 Ket(1)
      INTEGER Mesh,I,J,K
      REAL*8 ET,TwoT,T,T2,T3,T4
      REAL*8 SqrtT,SqrtPi,One,Two
      REAL*8 Switch,Grid,F0Asymp,F1Asymp,F2Asymp,F3Asymp
      REAL*8 F4Asymp,F5Asymp,F6Asymp,F7Asymp,F8Asymp
      REAL*8 F9Asymp,F10Asymp,F11Asymp,F12Asymp
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "Gamma_Asymptotics.Inc"
      INCLUDE "Mesh.Inc"
      REAL*8 F2_0(0:Mesh)
      REAL*8 F2_1(0:Mesh)
      REAL*8 F2_2(0:Mesh)
      REAL*8 F2_3(0:Mesh)
      INCLUDE "Gamma_2.Inc"
      T=O*(PQx*PQx+PQy*PQy+PQz*PQz)
      IF(T.LT.Switch)THEN
         T2=T*T
         T3=T*T2
         J=AINT(T*Grid)
         G(2)=(F2_0(J)+T*F2_1(J)+T2*F2_2(J)+T3*F2_3(J))
         ET=DEXP(-T)
         TwoT=2.0D0*T
         G(1)=0.3333333333333333D+00*(TwoT*G(2)+ET)
         G(0)=TwoT*G(1)+ET
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         AuxR(1)=o1*G(1)
         o1=o2*o1
         AuxR(2)=o1*G(2)
      ELSE
         SqrtT=DSQRT(T)
         OneOvT=One/T
         G(0)=SqrtPi/(Two*SqrtT)
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         G(1)=G(0)*0.5000000000000000D+00*OneOvT
         AuxR(1)=o1*G(1)
         o1=o2*o1
         G(2)=G(1)*0.1500000000000000D+01*OneOvT
         AuxR(2)=o1*G(2)
      ENDIF
      CALL MD3TRR__0_2(PQx,PQy,PQz,AuxR,R)
      CALL KTraX__0_2(R,Co,Ket)
      RETURN
      END

      SUBROUTINE HGTraX03(PQx,PQy,PQz,U,O,Co,Ket)
      IMPLICIT NONE
      REAL*8 PQx,PQy,PQz,U,O,o1,o2,OneOvT
      REAL*8 G(0:3)
      REAL*8 AuxR(0:3)
      REAL*8 R(20)
      REAL*8 Co(20)
      REAL*8 Ket(1)
      INTEGER Mesh,I,J,K
      REAL*8 ET,TwoT,T,T2,T3,T4
      REAL*8 SqrtT,SqrtPi,One,Two
      REAL*8 Switch,Grid,F0Asymp,F1Asymp,F2Asymp,F3Asymp
      REAL*8 F4Asymp,F5Asymp,F6Asymp,F7Asymp,F8Asymp
      REAL*8 F9Asymp,F10Asymp,F11Asymp,F12Asymp
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "Gamma_Asymptotics.Inc"
      INCLUDE "Mesh.Inc"
      REAL*8 F3_0(0:Mesh)
      REAL*8 F3_1(0:Mesh)
      REAL*8 F3_2(0:Mesh)
      REAL*8 F3_3(0:Mesh)
      INCLUDE "Gamma_3.Inc"
      T=O*(PQx*PQx+PQy*PQy+PQz*PQz)
      IF(T.LT.Switch)THEN
         T2=T*T
         T3=T*T2
         J=AINT(T*Grid)
         G(3)=(F3_0(J)+T*F3_1(J)+T2*F3_2(J)+T3*F3_3(J))
         ET=DEXP(-T)
         TwoT=2.0D0*T
         G(2)=0.2000000000000000D+00*(TwoT*G(3)+ET)
         G(1)=0.3333333333333333D+00*(TwoT*G(2)+ET)
         G(0)=TwoT*G(1)+ET
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         AuxR(1)=o1*G(1)
         o1=o2*o1
         AuxR(2)=o1*G(2)
         o1=o2*o1
         AuxR(3)=o1*G(3)
      ELSE
         SqrtT=DSQRT(T)
         OneOvT=One/T
         G(0)=SqrtPi/(Two*SqrtT)
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         G(1)=G(0)*0.5000000000000000D+00*OneOvT
         AuxR(1)=o1*G(1)
         o1=o2*o1
         G(2)=G(1)*0.1500000000000000D+01*OneOvT
         AuxR(2)=o1*G(2)
         o1=o2*o1
         G(3)=G(2)*0.2500000000000000D+01*OneOvT
         AuxR(3)=o1*G(3)
      ENDIF
      CALL MD3TRR__0_3(PQx,PQy,PQz,AuxR,R)
      CALL KTraX__0_3(R,Co,Ket)
      RETURN
      END

      SUBROUTINE HGTraX04(PQx,PQy,PQz,U,O,Co,Ket)
      IMPLICIT NONE
      REAL*8 PQx,PQy,PQz,U,O,o1,o2,OneOvT
      REAL*8 G(0:4)
      REAL*8 AuxR(0:4)
      REAL*8 R(35)
      REAL*8 Co(35)
      REAL*8 Ket(1)
      INTEGER Mesh,I,J,K
      REAL*8 ET,TwoT,T,T2,T3,T4
      REAL*8 SqrtT,SqrtPi,One,Two
      REAL*8 Switch,Grid,F0Asymp,F1Asymp,F2Asymp,F3Asymp
      REAL*8 F4Asymp,F5Asymp,F6Asymp,F7Asymp,F8Asymp
      REAL*8 F9Asymp,F10Asymp,F11Asymp,F12Asymp
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "Gamma_Asymptotics.Inc"
      INCLUDE "Mesh.Inc"
      REAL*8 F4_0(0:Mesh)
      REAL*8 F4_1(0:Mesh)
      REAL*8 F4_2(0:Mesh)
      REAL*8 F4_3(0:Mesh)
      INCLUDE "Gamma_4.Inc"
      T=O*(PQx*PQx+PQy*PQy+PQz*PQz)
      IF(T.LT.Switch)THEN
         T2=T*T
         T3=T*T2
         J=AINT(T*Grid)
         G(4)=(F4_0(J)+T*F4_1(J)+T2*F4_2(J)+T3*F4_3(J))
         ET=DEXP(-T)
         TwoT=2.0D0*T
         G(3)=0.1428571428571428D+00*(TwoT*G(4)+ET)
         G(2)=0.2000000000000000D+00*(TwoT*G(3)+ET)
         G(1)=0.3333333333333333D+00*(TwoT*G(2)+ET)
         G(0)=TwoT*G(1)+ET
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         AuxR(1)=o1*G(1)
         o1=o2*o1
         AuxR(2)=o1*G(2)
         o1=o2*o1
         AuxR(3)=o1*G(3)
         o1=o2*o1
         AuxR(4)=o1*G(4)
      ELSE
         SqrtT=DSQRT(T)
         OneOvT=One/T
         G(0)=SqrtPi/(Two*SqrtT)
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         G(1)=G(0)*0.5000000000000000D+00*OneOvT
         AuxR(1)=o1*G(1)
         o1=o2*o1
         G(2)=G(1)*0.1500000000000000D+01*OneOvT
         AuxR(2)=o1*G(2)
         o1=o2*o1
         G(3)=G(2)*0.2500000000000000D+01*OneOvT
         AuxR(3)=o1*G(3)
         o1=o2*o1
         G(4)=G(3)*0.3500000000000000D+01*OneOvT
         AuxR(4)=o1*G(4)
      ENDIF
      CALL MD3TRR__0_4(PQx,PQy,PQz,AuxR,R)
      CALL KTraX__0_4(R,Co,Ket)
      RETURN
      END
