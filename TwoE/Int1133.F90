! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (S S|P P) 
! ---------------------------------------------------------- 
   SUBROUTINE Int1133(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo,I) 
      USE DerivedTypes
      USE GlobalScalars
      USE GammaF2
      USE ONX2DataType
      IMPLICIT REAL(DOUBLE) (V,W)
      INTEGER        :: LBra,LKet
      REAL(DOUBLE)   :: PrmBufB(5,LBra),PrmBufK(5,LKet)
      TYPE(AtomInfo) :: ACInfo,BDInfo
      REAL(DOUBLE),DIMENSION(0:2) :: AuxR
      REAL(DOUBLE),DIMENSION(1,10) :: MBarN
      REAL(DOUBLE),DIMENSION(1,1,4,4) :: I
      REAL(DOUBLE)  :: Zeta,Eta,r1xZpE,HfxZpE,r1x2E,r1x2Z,ExZpE,ZxZpE,Omega,Up,Uq,Upq
      REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz,Wx,Wy,Wz
      REAL(DOUBLE)  :: QCx,QCy,QCz,PAx,PAy,PAz,PQx,PQy,PQz,WPx,WPy,WPz,WQx,WQy,WQz   
      REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT
      INTEGER       :: J,K,L
      MBarN=0.0d0
      Ax=ACInfo%Atm1X
      Ay=ACInfo%Atm1Y
      Az=ACInfo%Atm1Z
      Bx=BDInfo%Atm1X
      By=BDInfo%Atm1Y
      Bz=BDInfo%Atm1Z
      Cx=ACInfo%Atm2X
      Cy=ACInfo%Atm2Y
      Cz=ACInfo%Atm2Z
      Dx=BDInfo%Atm2X
      Dy=BDInfo%Atm2Y
      Dz=BDInfo%Atm2Z
      ABx=Ax-Bx
      ABy=Ay-By
      ABz=Az-Bz
      CDx=Cx-Dx
      CDy=Cy-Dy
      CDz=Cz-Dz
      DO J=1,LKet ! K^2 VRR |N0) loop 
         Eta=PrmBufK(1,J)
         Qx =PrmBufK(2,J)
         Qy =PrmBufK(3,J)
         Qz =PrmBufK(4,J)
         Uq =PrmBufK(5,J)
         QCx=Qx-Cx
         QCy=Qy-Cy
         QCz=Qz-Cz
         DO K=1,LBra ! K^2 VRR (M0| loop 
            Zeta=PrmBufB(1,K)
            Px  =PrmBufB(2,K)
            Py  =PrmBufB(3,K)
            Pz  =PrmBufB(4,K)
            Up  =PrmBufB(5,K)
            r1xZpE=One/(Zeta+Eta)
            Upq=SQRT(r1xZpE)*Up*Uq
            HfxZpE=Half/(Zeta+Eta)
            r1x2E=Half/Eta
            r1x2Z=Half/Zeta
            ExZpE=Eta*r1xZpE
            ZxZpE=Zeta*r1xZpE
            Omega=Eta*Zeta*r1xZpE
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
              V1=One
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
              SqInvT=DSQRT(InvT)
              AuxR(0)=+8.862269254527580D-01*SqInvT
              SqInvT=SqInvT*InvT
              AuxR(1)=+4.431134627263790D-01*SqInvT
              SqInvT=SqInvT*InvT
              AuxR(2)=+6.646701940895685D-01*SqInvT
            ENDIF
            V1=Upq*AuxR(0)
            V2=QCx*V1
            V3=Upq*AuxR(1)
            V4=V3*WQx
            MBarN(1,1)=V1+MBarN(1,1)
            MBarN(1,2)=V2+V4+MBarN(1,2)
            MBarN(1,3)=V1+MBarN(1,3)
            MBarN(1,4)=V1+MBarN(1,4)
            W1=QCx*(V2+V4)+r1x2E*(V1-V3*ZxZpE)
            W2=MBarN(1,5)+WQx*(QCx*V3+WQx*Upq*AuxR(2))
            MBarN(1,5)=W1+W2
            MBarN(1,6)=V2+V4+MBarN(1,6)
            MBarN(1,7)=V1+MBarN(1,7)
            MBarN(1,8)=V2+V4+MBarN(1,8)
            MBarN(1,9)=V1+MBarN(1,9)
            MBarN(1,10)=V1+MBarN(1,10)
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      V1=MBarN(1,2)
      V2=MBarN(1,3)
      V3=MBarN(1,6)
      V4=MBarN(1,4)
      V5=MBarN(1,8)
      V6=MBarN(1,9)
      I(1,1,2,2)=CDx*V1+MBarN(1,5)
      I(1,1,3,2)=CDx*V2+V3
      I(1,1,4,2)=CDx*V4+V5
      I(1,1,2,3)=CDy*V1+V3
      I(1,1,3,3)=CDy*V2+MBarN(1,7)
      I(1,1,4,3)=CDy*V4+V6
      I(1,1,2,4)=CDz*V1+V5
      I(1,1,3,4)=CDz*V2+V6
      I(1,1,4,4)=CDz*V4+MBarN(1,10)
   END SUBROUTINE Int1133
