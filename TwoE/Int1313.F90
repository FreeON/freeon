! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (S P|S P) 
! ---------------------------------------------------------- 
   SUBROUTINE Int1313(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, & 
                              OA,LDA,OB,LDB,OC,LDC,OD,LDD,PBC,I) 
      USE DerivedTypes
      USE VScratch
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF2
      IMPLICIT REAL(DOUBLE) (A,I,W)
      INTEGER        :: LBra,LKet
      REAL(DOUBLE)   :: PrmBufB(7,LBra),PrmBufK(7,LKet)
      TYPE(SmallAtomInfo) :: ACInfo,BDInfo
      TYPE(PBCInfo) :: PBC
      REAL(DOUBLE) :: I(*)
      REAL(DOUBLE)  :: Zeta,Eta,r1xZpE,HfxZpE,r1x2E,r1x2Z,ExZpE,ZxZpE,Omega,Up,Uq,Upq
      REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz
      REAL(DOUBLE)  :: QCx,QCy,QCz,PAx,PAy,PAz,PQx,PQy,PQz,WPx,WPy,WPz,WQx,WQy,WQz   
      REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT,ABx,ABy,ABz,CDx,CDy,CDz
      INTEGER       :: OffSet,OA,LDA,OB,LDB,OC,LDC,OD,LDD,J,K,L
      REAL(DOUBLE)  :: FPQx,FPQy,FPQz
      I1Bar1=0.0d0
      I2Bar1=0.0d0
      I3Bar1=0.0d0
      I4Bar1=0.0d0
      I1Bar2=0.0d0
      I2Bar2=0.0d0
      I3Bar2=0.0d0
      I4Bar2=0.0d0
      I1Bar3=0.0d0
      I2Bar3=0.0d0
      I3Bar3=0.0d0
      I4Bar3=0.0d0
      I1Bar4=0.0d0
      I2Bar4=0.0d0
      I3Bar4=0.0d0
      I4Bar4=0.0d0
      Ax=ACInfo%Atm1X
      Ay=ACInfo%Atm1Y
      Az=ACInfo%Atm1Z
      Bx=ACInfo%Atm2X
      By=ACInfo%Atm2Y
      Bz=ACInfo%Atm2Z
      Cx=BDInfo%Atm1X
      Cy=BDInfo%Atm1Y
      Cz=BDInfo%Atm1Z
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
            PAx=Px-Ax
            PAy=Py-Ay
            PAz=Pz-Az
            PQx=Px-Qx
            PQy=Py-Qy
            PQz=Pz-Qz
      ! Need to be improve...
            FPQx = PQx*PBC%InvBoxSh%D(1,1)+PQy*PBC%InvBoxSh%D(1,2)+PQz*PBC%InvBoxSh%D(1,3)
            FPQy = PQy*PBC%InvBoxSh%D(2,2)+PQz*PBC%InvBoxSh%D(2,3)
            FPQz = PQz*PBC%InvBoxSh%D(3,3)
            IF(PBC%AutoW%I(1)==1) FPQx = FPQx-ANINT(ANINT(FPQx*1d9)*1d-9)
            IF(PBC%AutoW%I(2)==1) FPQy = FPQy-ANINT(ANINT(FPQy*1d9)*1d-9)
            IF(PBC%AutoW%I(3)==1) FPQz = FPQz-ANINT(ANINT(FPQz*1d9)*1d-9)
            PQx  = FPQx*PBC%BoxShape%D(1,1)+FPQy*PBC%BoxShape%D(1,2)+FPQz*PBC%BoxShape%D(1,3)
            PQy  = FPQy*PBC%BoxShape%D(2,2)+FPQz*PBC%BoxShape%D(2,3)
            PQz  = FPQz*PBC%BoxShape%D(3,3)
      !
            WPx = -Eta*PQx*r1xZpE
            WPy = -Eta*PQy*r1xZpE
            WPz = -Eta*PQz*r1xZpE
            WQx = Zeta*PQx*r1xZpE
            WQy = Zeta*PQy*r1xZpE
            WQz = Zeta*PQz*r1xZpE
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
      V(1)=AuxR0*PAx
      V(2)=AuxR1*WPx
      V(3)=AuxR0*PAy
      V(4)=AuxR1*WPy
      V(5)=AuxR0*PAz
      V(6)=AuxR1*WPz
      V(7)=AuxR1*HfxZpE
      V(8)=V(1)+V(2)
      V(9)=AuxR1*PAx
      V(10)=AuxR2*WPx
      V(11)=V(9)+V(10)
      V(12)=V(3)+V(4)
      V(13)=AuxR1*PAy
      V(14)=AuxR2*WPy
      V(15)=V(13)+V(14)
      V(16)=V(5)+V(6)
      V(17)=AuxR1*PAz
      V(18)=AuxR2*WPz
      V(19)=V(17)+V(18)
      I1Bar1=AuxR0+I1Bar1
      I2Bar1=I2Bar1+V(1)+V(2)
      I3Bar1=I3Bar1+V(3)+V(4)
      I4Bar1=I4Bar1+V(5)+V(6)
      I1Bar2=AuxR0*QCx+AuxR1*WQx+I1Bar2
      I2Bar2=I2Bar2+V(7)+QCx*V(8)+WQx*V(11)
      I3Bar2=I3Bar2+QCx*V(12)+WQx*V(15)
      I4Bar2=I4Bar2+QCx*V(16)+WQx*V(19)
      I1Bar3=AuxR0*QCy+AuxR1*WQy+I1Bar3
      I2Bar3=I2Bar3+QCy*V(8)+WQy*V(11)
      I3Bar3=I3Bar3+V(7)+QCy*V(12)+WQy*V(15)
      I4Bar3=I4Bar3+QCy*V(16)+WQy*V(19)
      I1Bar4=AuxR0*QCz+AuxR1*WQz+I1Bar4
      I2Bar4=I2Bar4+QCz*V(8)+WQz*V(11)
      I3Bar4=I3Bar4+QCz*V(12)+WQz*V(15)
      I4Bar4=I4Bar4+V(7)+QCz*V(16)+WQz*V(19)
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      V(1)=CDx*I1Bar1
      V(2)=I1Bar2+V(1)
      V(3)=CDy*I1Bar1
      V(4)=I1Bar3+V(3)
      V(5)=CDz*I1Bar1
      V(6)=I1Bar4+V(5)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=CDx*I2Bar1+I2Bar2+I(OffSet)+ABx*V(2)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=CDx*I3Bar1+I3Bar2+I(OffSet)+ABy*V(2)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=CDx*I4Bar1+I4Bar2+I(OffSet)+ABz*V(2)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+1)*LDD 
      I(OffSet)=CDy*I2Bar1+I2Bar3+I(OffSet)+ABx*V(4)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+(OC+0)*LDC+(OD+1)*LDD 
      I(OffSet)=CDy*I3Bar1+I3Bar3+I(OffSet)+ABy*V(4)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+(OC+0)*LDC+(OD+1)*LDD 
      I(OffSet)=CDy*I4Bar1+I4Bar3+I(OffSet)+ABz*V(4)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+2)*LDD 
      I(OffSet)=CDz*I2Bar1+I2Bar4+I(OffSet)+ABx*V(6)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+(OC+0)*LDC+(OD+2)*LDD 
      I(OffSet)=CDz*I3Bar1+I3Bar4+I(OffSet)+ABy*V(6)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+(OC+0)*LDC+(OD+2)*LDD 
      I(OffSet)=CDz*I4Bar1+I4Bar4+I(OffSet)+ABz*V(6)
   END SUBROUTINE Int1313
