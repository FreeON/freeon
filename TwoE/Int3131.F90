! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (P S|P S) 
! ---------------------------------------------------------- 
   SUBROUTINE Int3131(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo,I) 
      USE DerivedTypes
      USE GlobalScalars
      USE GammaF2
      USE ONX2DataType
      IMPLICIT REAL(DOUBLE) (A,I,V,W)
      INTEGER        :: LBra,LKet
      REAL(DOUBLE)   :: PrmBufB(5,LBra),PrmBufK(5,LKet)
      TYPE(AtomInfo) :: ACInfo,BDInfo
      REAL(DOUBLE),DIMENSION(4,1,4,1) :: I
      REAL(DOUBLE)  :: Zeta,Eta,r1xZpE,HfxZpE,r1x2E,r1x2Z,ExZpE,ZxZpE,Omega,Up,Uq,Upq
      REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz,Wx,Wy,Wz
      REAL(DOUBLE)  :: QCx,QCy,QCz,PAx,PAy,PAz,PQx,PQy,PQz,WPx,WPy,WPz,WQx,WQy,WQz   
      REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT
      INTEGER       :: J,K,L
      I1Bar1=Zero
      I2Bar1=Zero
      I3Bar1=Zero
      I4Bar1=Zero
      I1Bar2=Zero
      I2Bar2=Zero
      I3Bar2=Zero
      I4Bar2=Zero
      I1Bar3=Zero
      I2Bar3=Zero
      I3Bar3=Zero
      I4Bar3=Zero
      I1Bar4=Zero
      I2Bar4=Zero
      I3Bar4=Zero
      I4Bar4=Zero
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
              ET=EXP(-T)
              TwoT=Two*T
              W2=(F2_0(L)+T*(F2_1(L)+T*(F2_2(L)+T*(F2_3(L)+T*F2_4(L)))))
              W1=+3.333333333333333D-01*(TwoT*W2+ET)
              W0=TwoT*W1+ET
              AuxR0=Upq*W0
              AuxR1=Upq*W1
              AuxR2=Upq*W2
            ELSE
              InvT=One/T
              SqInvT=DSQRT(InvT)
              AuxR0=+8.862269254527580D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              AuxR1=+4.431134627263790D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              AuxR2=+6.646701940895685D-01*Upq*SqInvT
            ENDIF
            V1=AuxR0*PAx
            V2=AuxR1*WPx
            V3=AuxR0*PAy
            V4=AuxR1*WPy
            V5=AuxR0*PAz
            V6=AuxR1*WPz
            V7=AuxR1*HfxZpE
            V8=V1+V2
            V9=AuxR1*PAx
            V10=AuxR2*WPx
            V11=V10+V9
            V12=V3+V4
            V13=AuxR1*PAy
            V14=AuxR2*WPy
            V15=V13+V14
            V16=V5+V6
            V17=AuxR1*PAz
            V18=AuxR2*WPz
            V19=V17+V18
            I1Bar1=AuxR0+I1Bar1
            I2Bar1=V1+V2+I2Bar1
            I3Bar1=V3+V4+I3Bar1
            I4Bar1=V5+V6+I4Bar1
            I1Bar2=AuxR0*QCx+AuxR1*WQx+I1Bar2
            I2Bar2=V7+QCx*V8+V11*WQx+I2Bar2
            I3Bar2=QCx*V12+V15*WQx+I3Bar2
            I4Bar2=QCx*V16+V19*WQx+I4Bar2
            I1Bar3=AuxR0*QCy+AuxR1*WQy+I1Bar3
            I2Bar3=QCy*V8+V11*WQy+I2Bar3
            I3Bar3=QCy*V12+V7+V15*WQy+I3Bar3
            I4Bar3=QCy*V16+V19*WQy+I4Bar3
            I1Bar4=AuxR0*QCz+AuxR1*WQz+I1Bar4
            I2Bar4=QCz*V8+V11*WQz+I2Bar4
            I3Bar4=QCz*V12+V15*WQz+I3Bar4
            I4Bar4=QCz*V16+V7+V19*WQz+I4Bar4
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      I(2,1,2,1)=I2Bar2
      I(3,1,2,1)=I3Bar2
      I(4,1,2,1)=I4Bar2
      I(2,1,3,1)=I2Bar3
      I(3,1,3,1)=I3Bar3
      I(4,1,3,1)=I4Bar3
      I(2,1,4,1)=I2Bar4
      I(3,1,4,1)=I3Bar4
      I(4,1,4,1)=I4Bar4
   END SUBROUTINE Int3131
