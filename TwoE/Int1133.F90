! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (S S|P P) 
! ---------------------------------------------------------- 
   SUBROUTINE Int1133(Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz, &  
                      Dx,Dy,Dz,ShlPrAC2,ShlPrBD2,I) 
      USE DerivedTypes
      USE GlobalScalars
      USE GammaF2
      USE ShellPairStruct
      IMPLICIT REAL(DOUBLE) (V,W)
      TYPE(ShellPair), POINTER :: ShlPrAC2,ShlPrBD2 
      REAL(DOUBLE),DIMENSION(2) :: AuxR
      REAL(DOUBLE),DIMENSION(1,10) :: MBarN=0D0
      REAL(DOUBLE),DIMENSION(1,1,4,4) :: I
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
              W2=(F2_0(L)+T*(F2_1(L)+T*(F2_2(L)+T*(F2_3(L)+T*F2_4(L)))))
              W1=+3.333333333333333D-01*(TwoT*W2+ET)
              W0=TwoT*W1+ET
              AuxR(0)=V1*W0
              V1=V2*V1
              AuxR(1)=V1*W1
              V1=V2*V1
              AuxR(2)=V1*W2
            ELSE
              InvT=One/T
              SqInvT=Upq*DSQRT(InvT)
              AuxR(0)=+8.862269254527580D-01*SqInvT
              SqInvT=SqInvT*InvT
              AuxR(1)=+4.431134627263790D-01*SqInvT
              SqInvT=SqInvT*InvT
              AuxR(2)=+6.646701940895685D-01*SqInvT
            ENDIF
            V1=AuxR(0)
            V2=QCx*V1
            V3=AuxR(1)
            V4=V3*WQx
            MBarN(1,1)=V1+MBarN(1,1)
            MBarN(1,2)=V2+V4+MBarN(1,2)
            MBarN(1,3)=V1+MBarN(1,3)
            MBarN(1,4)=V1+MBarN(1,4)
            W1=QCx*(V2+V4)+r1x2E*(V1-V3*ZxZpE)
            W2=MBarN(1,5)+WQx*(QCx*V3+WQx*AuxR(2))
            MBarN(1,5)=W1+W2
            MBarN(1,6)=V2+V4+MBarN(1,6)
            MBarN(1,7)=V1+MBarN(1,7)
            MBarN(1,8)=V2+V4+MBarN(1,8)
            MBarN(1,9)=V1+MBarN(1,9)
            MBarN(1,10)=V1+MBarN(1,10)
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      I(1,1,2,2)=MBarN(2,5)+MBarN(2,2)*CDx+(MBarN(1,5)+MBarN(1,2)*CDx)*ABx+I(1,1,2,2)
      I(1,1,3,2)=MBarN(2,6)+MBarN(2,3)*CDx+(MBarN(1,6)+MBarN(1,3)*CDx)*ABx+I(1,1,3,2)
      I(1,1,4,2)=MBarN(2,8)+MBarN(2,4)*CDx+(MBarN(1,8)+MBarN(1,4)*CDx)*ABx+I(1,1,4,2)
      I(1,1,2,3)=MBarN(2,6)+MBarN(2,2)*CDy+ABx*(MBarN(1,6)+MBarN(1,2)*CDy)+I(1,1,2,3)
      I(1,1,3,3)=MBarN(2,7)+MBarN(2,3)*CDy+ABx*(MBarN(1,7)+MBarN(1,3)*CDy)+I(1,1,3,3)
      I(1,1,4,3)=MBarN(2,9)+MBarN(2,4)*CDy+ABx*(MBarN(1,9)+MBarN(1,4)*CDy)+I(1,1,4,3)
      I(1,1,2,4)=MBarN(2,8)+MBarN(2,2)*CDz+ABx*(MBarN(1,8)+MBarN(1,2)*CDz)+I(1,1,2,4)
      I(1,1,3,4)=MBarN(2,9)+MBarN(2,3)*CDz+ABx*(MBarN(1,9)+MBarN(1,3)*CDz)+I(1,1,3,4)
      I(1,1,4,4)=MBarN(2,10)+MBarN(2,4)*CDz+ABx*(MBarN(1,10)+MBarN(1,4)*CDz)+I(1,1,4,4)
   END SUBROUTINE Int1133
