! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (S S|S P) 
! ---------------------------------------------------------- 
   SUBROUTINE Int1113(Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz, &  
                      Dx,Dy,Dz,ShlPrAC2,ShlPrBD2,I) 
      USE DerivedTypes
      USE GlobalScalars
      USE GammaF0
      USE GammaF1
      USE ShellPairStruct
      IMPLICIT REAL(DOUBLE) (V,W)
      TYPE(ShellPair), POINTER :: ShlPrAC2,ShlPrBD2 
      REAL(DOUBLE),DIMENSION(0:1) :: AuxR
      REAL(DOUBLE),DIMENSION(1,4) :: MBarN=0D0
      REAL(DOUBLE),DIMENSION(1,1,1,4) :: I
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
            Omega=ExZpE*ZxZpE
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
              AuxR(0)=(F0_0(L)+T*(F0_1(L)+T*(F0_2(L)+T*(F0_3(L)+T*F0_4(L)))))
              AuxR(1)=(F1_0(L)+T*(F1_1(L)+T*(F1_2(L)+T*(F1_3(L)+T*F1_4(L)))))
            ELSE
              InvT=One/T
              SqInvT=Upq*DSQRT(InvT)
              AuxR(0)=+8.862269254527580D-01*SqInvT
              SqInvT=SqInvT*InvT
              AuxR(1)=+4.431134627263790D-01*SqInvT
            ENDIF
            V1=Upq*AuxR(0)
            MBarN(1,1)=V1+MBarN(1,1)
            MBarN(1,2)=QCx*V1+MBarN(1,2)+WQx*Upq*AuxR(1)
            MBarN(1,3)=V1+MBarN(1,3)
            MBarN(1,4)=V1+MBarN(1,4)
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      V1=MBarN(1,1)
      I(1,1,1,2)=CDx*V1+MBarN(1,2)
      I(1,1,1,3)=CDy*V1+MBarN(1,3)
      I(1,1,1,4)=CDz*V1+MBarN(1,4)
   END SUBROUTINE Int1113
