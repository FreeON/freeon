! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (S P|P S) 
! ---------------------------------------------------------- 
   SUBROUTINE Int1331(Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz, &  
                      Dx,Dy,Dz,ShlPrAC2,ShlPrBD2,I) 
      USE DerivedTypes
      USE GlobalScalars
      USE GammaF2
      USE ShellPairStruct
      IMPLICIT REAL(DOUBLE) (V,W)
      TYPE(ShellPair), POINTER :: ShlPrAC2,ShlPrBD2 
      REAL(DOUBLE),DIMENSION(0:2) :: AuxR
      REAL(DOUBLE),DIMENSION(4,4) :: MBarN=0D0
      REAL(DOUBLE),DIMENSION(1,4,4,1) :: I
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
              SqInvT=Upq*DSQRT(InvT)
              AuxR(0)=+8.862269254527580D-01*SqInvT
              SqInvT=SqInvT*InvT
              AuxR(1)=+4.431134627263790D-01*SqInvT
              SqInvT=SqInvT*InvT
              AuxR(2)=+6.646701940895685D-01*SqInvT
            ENDIF
            V1=Upq*AuxR(0)
            V2=PAx*V1
            V3=Upq*AuxR(1)
            V4=V3*WPx
            V5=PAy*V1
            V6=V3*WPy
            V7=PAz*V1
            V8=V3*WPz
            V9=HfxZpE*V3
            V10=V2+V4
            V11=PAx*V3
            V12=Upq*AuxR(2)
            V13=V12*WPx
            V14=V11+V13
            V15=V5+V6
            V16=PAy*V3
            V17=V12*WPy
            V18=V16+V17
            V19=V7+V8
            V20=PAz*V3
            V21=V12*WPz
            V22=V20+V21
            MBarN(1,1)=V1+MBarN(1,1)
            MBarN(2,1)=V2+V4+MBarN(2,1)
            MBarN(3,1)=V5+V6+MBarN(3,1)
            MBarN(4,1)=V7+V8+MBarN(4,1)
            MBarN(1,2)=QCx*V1+V3*WQx+MBarN(1,2)
            MBarN(2,2)=QCx*V10+V9+V14*WQx+MBarN(2,2)
            MBarN(3,2)=QCx*V15+V18*WQx+MBarN(3,2)
            MBarN(4,2)=QCx*V19+V22*WQx+MBarN(4,2)
            MBarN(1,3)=V1+MBarN(1,3)
            MBarN(2,3)=QCy*V10+V14*WQy+MBarN(2,3)
            MBarN(3,3)=QCy*V15+V9+V18*WQy+MBarN(3,3)
            MBarN(4,3)=QCy*V19+V22*WQy+MBarN(4,3)
            MBarN(1,4)=V1+MBarN(1,4)
            MBarN(2,4)=QCz*V10+V14*WQz+MBarN(2,4)
            MBarN(3,4)=QCz*V15+V18*WQz+MBarN(3,4)
            MBarN(4,4)=QCz*V19+V9+V22*WQz+MBarN(4,4)
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      V1=MBarN(1,2)
      V2=MBarN(1,3)
      V3=MBarN(1,4)
      I(1,2,2,1)=ABx*V1+MBarN(2,2)
      I(1,3,2,1)=ABy*V1+MBarN(3,2)
      I(1,4,2,1)=ABz*V1+MBarN(4,2)
      I(1,2,3,1)=ABx*V2+MBarN(2,3)
      I(1,3,3,1)=ABy*V2+MBarN(3,3)
      I(1,4,3,1)=ABz*V2+MBarN(4,3)
      I(1,2,4,1)=ABx*V3+MBarN(2,4)
      I(1,3,4,1)=ABy*V3+MBarN(3,4)
      I(1,4,4,1)=ABz*V3+MBarN(4,4)
   END SUBROUTINE Int1331
