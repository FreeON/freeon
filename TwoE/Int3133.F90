! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (P S|P P) 
! ---------------------------------------------------------- 
   SUBROUTINE Int3133(Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz, &  
                      Dx,Dy,Dz,ShlPrAC2,ShlPrBD2,I) 
      USE DerivedTypes
      USE GlobalScalars
      USE GammaF3
      USE ShellPairStruct
      IMPLICIT REAL(DOUBLE) (V,W)
      TYPE(ShellPair), POINTER :: ShlPrAC2,ShlPrBD2 
      REAL(DOUBLE),DIMENSION(3) :: AuxR
      REAL(DOUBLE),DIMENSION(4,10) :: MBarN=0D0
      REAL(DOUBLE),DIMENSION(4,1,4,4) :: I
      REAL(DOUBLE)  :: Zeta,Eta,r1xZpE,HfxZpE,r1x2E,r1x2Z,ExZpE,ZxZpE,Omega,Up,Uq,Upq
      REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz,Wx,Wy,Wz
      REAL(DOUBLE)  :: QCx,QCy,QCz,PAx,PAy,PAz,PQx,PQy,PQz,WPx,WPy,WPz,WQx,WQy,WQz   
      REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT
      INTEGER       :: J,K,L
      DO J=1,ShlPrBD2%L ! K^2 VRR |N0) loop 
         Eta=ShlPrBD2%SP(1,J)
         Qx =ShlPrBD2%SP(2,J)
         Qy =ShlPrBD2%SP(3,J)
         Qz =ShlPrBD2%SP(4,J)
         Uq =ShlPrBD2%SP(5,J)
         QCx=Qx-Cx
         QCy=Qy-Cy
         QCz=Qz-Cz
         DO K=1,ShlPrAC2%L ! K^2 VRR (M0| loop 
            Zeta=ShlPrAC2%SP(1,K)
            Px  =ShlPrAC2%SP(2,K)
            Py  =ShlPrAC2%SP(3,K)
            Pz  =ShlPrAC2%SP(4,K)
            Up  =ShlPrAC2%SP(5,K)
            r1xZpE=One/(Zeta+Eta)
            Upq=SQRT(r1xZpE)*Up*Uq
            HfxZpE=Half/(Zeta+Eta)
            r1x2E=Half/Eta
            r1x2Z=Half/Zeta
            ExZpE=Eta*r1xZpE
            ZxZpE=Zeta*r1xZpE
            Omega=ExZpe+ZxZpE
            Wx=(Zeta*Px+Eta*Qx)*r1xZpE
            Wy=(Zeta*Py+Eta*Qy)*r1xZpE
            Wz=(Zeta*Pz+Eta*Qz)*r1xZpE
            PAx=Px-Ax
            PAy=Py-Ay
            PAz=Pz-Az
            PQx=Px-Qx
            PQy=Py-Qy
            PQz=Pz-Qz
            WPx=Wx-Px
            WPy=Wy-Py
            WPz=Wz-Pz
            WQx=Wx-Qx
            WQy=Wy-Qy
            WQz=Wz-Qz
            T=Omega*(PQx*PQx+PQy*PQy+PQz*PQz)
            IF(T<Gamma_Switch)THEN
              L=AINT(T*Gamma_Grid)
              V1=Upq
              V2=-Two*Omega
              ET=EXP(-T)
              TwoT=Two*T
              W3=(F3_0(L)+T*(F3_1(L)+T*(F3_2(L)+T*(F3_3(L)+T*F3_4(L)))))
              W2=+2.000000000000000D-01*(TwoT*W3+ET)
              W1=+3.333333333333333D-01*(TwoT*W2+ET)
              W0=TwoT*W1+ET
              AuxR(0)=V1*W0
              V1=V2*V1
              AuxR(1)=V1*W1
              V1=V2*V1
              AuxR(2)=V1*W2
              V1=V2*V1
              AuxR(3)=V1*W3
            ELSE
              InvT=One/T
              SqInvT=Upq*DSQRT(InvT)
              AuxR(0)=+8.862269254527580D-01*SqInvT
              SqInvT=SqInvT*InvT
              AuxR(1)=+4.431134627263790D-01*SqInvT
              SqInvT=SqInvT*InvT
              AuxR(2)=+6.646701940895685D-01*SqInvT
              SqInvT=SqInvT*InvT
              AuxR(3)=+1.661675485223921D+00*SqInvT
            ENDIF
            V1=AuxR(0)
            V2=PAx*V1
            V3=AuxR(1)
            V4=V3*WPx
            V5=PAy*V1
            V6=V3*WPy
            V7=PAz*V1
            V8=V3*WPz
            V9=QCx*V1
            V10=V3*WQx
            V11=HfxZpE*V3
            V12=V2+V4
            V13=QCx*V12
            V14=PAx*V3
            V15=AuxR(2)
            V16=V15*WPx
            V17=V14+V16
            V18=V17*WQx
            V19=V5+V6
            V20=QCx*V19
            V21=PAy*V3
            V22=V15*WPy
            V23=V21+V22
            V24=V23*WQx
            V25=V7+V8
            V26=QCx*V25
            V27=PAz*V3
            V28=V15*WPz
            V29=V27+V28
            V30=V29*WQx
            V31=QCy*V12
            V32=V17*WQy
            V33=QCy*V19
            V34=V23*WQy
            V35=QCy*V25
            V36=V29*WQy
            V37=QCz*V12
            V38=V17*WQz
            V39=QCz*V19
            V40=V23*WQz
            V41=QCz*V25
            V42=V29*WQz
            V43=QCx*V3
            V44=V15*WQx
            V45=V43+V44
            V46=-(V17*ZxZpE)
            V47=V2+V4+V46
            V48=r1x2E*V47
            V49=HfxZpE*V15
            V50=PAx*V15
            V51=AuxR(3)
            V52=V51*WPx
            V53=V50+V52
            V54=-(V23*ZxZpE)
            V55=V5+V54+V6
            V56=r1x2E*V55
            V57=PAy*V15
            V58=V51*WPy
            V59=V57+V58
            V60=-(V29*ZxZpE)
            V61=V60+V7+V8
            V62=r1x2E*V61
            V63=PAz*V15
            V64=V51*WPz
            V65=V63+V64
            V66=V31+V32
            V67=QCy*V17
            V68=V53*WQy
            V69=V67+V68
            V70=V11+V33+V34
            V71=QCy*V23
            V72=V59*WQy
            V73=V49+V71+V72
            V74=V35+V36
            V75=QCy*V29
            V76=V65*WQy
            V77=V75+V76
            V78=V37+V38
            V79=QCz*V17
            V80=V53*WQz
            V81=V79+V80
            V82=V39+V40
            V83=QCz*V23
            V84=V59*WQz
            V85=V83+V84
            V86=V11+V41+V42
            V87=QCz*V29
            V88=V65*WQz
            V89=V49+V87+V88
            MBarN(1,1)=V1+MBarN(1,1)
            MBarN(2,1)=V2+V4+MBarN(2,1)
            MBarN(3,1)=V5+V6+MBarN(3,1)
            MBarN(4,1)=V7+V8+MBarN(4,1)
            MBarN(1,2)=V10+V9+MBarN(1,2)
            MBarN(2,2)=V11+V13+V18+MBarN(2,2)
            MBarN(3,2)=V20+V24+MBarN(3,2)
            MBarN(4,2)=V26+V30+MBarN(4,2)
            MBarN(1,3)=V1+MBarN(1,3)
            MBarN(2,3)=V31+V32+MBarN(2,3)
            MBarN(3,3)=V11+V33+V34+MBarN(3,3)
            MBarN(4,3)=V35+V36+MBarN(4,3)
            MBarN(1,4)=V1+MBarN(1,4)
            MBarN(2,4)=V37+V38+MBarN(2,4)
            MBarN(3,4)=V39+V40+MBarN(3,4)
            MBarN(4,4)=V11+V41+V42+MBarN(4,4)
            W1=QCx*(V10+V9)+V45*WQx
            W2=r1x2E*(V1-V3*ZxZpE)+MBarN(1,5)
            MBarN(1,5)=W1+W2
            W1=QCx*(V11+V13+V18)+HfxZpE*V45
            W2=V48+WQx*(QCx*V17+V49+V53*WQx)+MBarN(2,5)
            MBarN(2,5)=W1+W2
            W1=QCx*(V20+V24)+V56
            W2=WQx*(QCx*V23+V59*WQx)+MBarN(3,5)
            MBarN(3,5)=W1+W2
            W1=QCx*(V26+V30)+V62
            W2=WQx*(QCx*V29+V65*WQx)+MBarN(4,5)
            MBarN(4,5)=W1+W2
            MBarN(1,6)=V10+V9+MBarN(1,6)
            MBarN(2,6)=V11+QCx*V66+V69*WQx+MBarN(2,6)
            MBarN(3,6)=QCx*V70+V73*WQx+MBarN(3,6)
            MBarN(4,6)=QCx*V74+V77*WQx+MBarN(4,6)
            MBarN(1,7)=V1+MBarN(1,7)
            MBarN(2,7)=V48+QCy*V66+V69*WQy+MBarN(2,7)
            MBarN(3,7)=V11+V56+QCy*V70+V73*WQy+MBarN(3,7)
            MBarN(4,7)=V62+QCy*V74+V77*WQy+MBarN(4,7)
            MBarN(1,8)=V10+V9+MBarN(1,8)
            MBarN(2,8)=V11+QCx*V78+V81*WQx+MBarN(2,8)
            MBarN(3,8)=QCx*V82+V85*WQx+MBarN(3,8)
            MBarN(4,8)=QCx*V86+V89*WQx+MBarN(4,8)
            MBarN(1,9)=V1+MBarN(1,9)
            MBarN(2,9)=QCy*V78+V81*WQy+MBarN(2,9)
            MBarN(3,9)=V11+QCy*V82+V85*WQy+MBarN(3,9)
            MBarN(4,9)=QCy*V86+V89*WQy+MBarN(4,9)
            MBarN(1,10)=V1+MBarN(1,10)
            MBarN(2,10)=V48+QCz*V78+V81*WQz+MBarN(2,10)
            MBarN(3,10)=V56+QCz*V82+V85*WQz+MBarN(3,10)
            MBarN(4,10)=V11+V62+QCz*V86+V89*WQz+MBarN(4,10)
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      I(2,1,2,2)=MBarN(5,5)+MBarN(5,2)*CDx+(MBarN(2,5)+MBarN(2,2)*CDx)*ABx
      I(3,1,2,2)=MBarN(6,5)+MBarN(6,2)*CDx+(MBarN(3,5)+MBarN(3,2)*CDx)*ABx
      I(4,1,2,2)=MBarN(8,5)+MBarN(8,2)*CDx+(MBarN(4,5)+MBarN(4,2)*CDx)*ABx
      I(2,1,3,2)=MBarN(5,6)+MBarN(5,3)*CDx+(MBarN(2,6)+MBarN(2,3)*CDx)*ABx
      I(3,1,3,2)=MBarN(6,6)+MBarN(6,3)*CDx+(MBarN(3,6)+MBarN(3,3)*CDx)*ABx
      I(4,1,3,2)=MBarN(8,6)+MBarN(8,3)*CDx+(MBarN(4,6)+MBarN(4,3)*CDx)*ABx
      I(2,1,4,2)=MBarN(5,8)+MBarN(5,4)*CDx+(MBarN(2,8)+MBarN(2,4)*CDx)*ABx
      I(3,1,4,2)=MBarN(6,8)+MBarN(6,4)*CDx+(MBarN(3,8)+MBarN(3,4)*CDx)*ABx
      I(4,1,4,2)=MBarN(8,8)+MBarN(8,4)*CDx+(MBarN(4,8)+MBarN(4,4)*CDx)*ABx
      I(2,1,2,3)=MBarN(5,6)+MBarN(5,2)*CDy+ABx*(MBarN(2,6)+MBarN(2,2)*CDy)
      I(3,1,2,3)=MBarN(6,6)+MBarN(6,2)*CDy+ABx*(MBarN(3,6)+MBarN(3,2)*CDy)
      I(4,1,2,3)=MBarN(8,6)+MBarN(8,2)*CDy+ABx*(MBarN(4,6)+MBarN(4,2)*CDy)
      I(2,1,3,3)=MBarN(5,7)+MBarN(5,3)*CDy+ABx*(MBarN(2,7)+MBarN(2,3)*CDy)
      I(3,1,3,3)=MBarN(6,7)+MBarN(6,3)*CDy+ABx*(MBarN(3,7)+MBarN(3,3)*CDy)
      I(4,1,3,3)=MBarN(8,7)+MBarN(8,3)*CDy+ABx*(MBarN(4,7)+MBarN(4,3)*CDy)
      I(2,1,4,3)=MBarN(5,9)+MBarN(5,4)*CDy+ABx*(MBarN(2,9)+MBarN(2,4)*CDy)
      I(3,1,4,3)=MBarN(6,9)+MBarN(6,4)*CDy+ABx*(MBarN(3,9)+MBarN(3,4)*CDy)
      I(4,1,4,3)=MBarN(8,9)+MBarN(8,4)*CDy+ABx*(MBarN(4,9)+MBarN(4,4)*CDy)
      I(2,1,2,4)=MBarN(5,8)+MBarN(5,2)*CDz+ABx*(MBarN(2,8)+MBarN(2,2)*CDz)
      I(3,1,2,4)=MBarN(6,8)+MBarN(6,2)*CDz+ABx*(MBarN(3,8)+MBarN(3,2)*CDz)
      I(4,1,2,4)=MBarN(8,8)+MBarN(8,2)*CDz+ABx*(MBarN(4,8)+MBarN(4,2)*CDz)
      I(2,1,3,4)=MBarN(5,9)+MBarN(5,3)*CDz+ABx*(MBarN(2,9)+MBarN(2,3)*CDz)
      I(3,1,3,4)=MBarN(6,9)+MBarN(6,3)*CDz+ABx*(MBarN(3,9)+MBarN(3,3)*CDz)
      I(4,1,3,4)=MBarN(8,9)+MBarN(8,3)*CDz+ABx*(MBarN(4,9)+MBarN(4,3)*CDz)
      I(2,1,4,4)=MBarN(5,10)+MBarN(5,4)*CDz+ABx*(MBarN(2,10)+MBarN(2,4)*CDz)
      I(3,1,4,4)=MBarN(6,10)+MBarN(6,4)*CDz+ABx*(MBarN(3,10)+MBarN(3,4)*CDz)
      I(4,1,4,4)=MBarN(8,10)+MBarN(8,4)*CDz+ABx*(MBarN(4,10)+MBarN(4,4)*CDz)
   END SUBROUTINE Int3133
