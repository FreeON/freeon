! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (P S|S S) 
! ---------------------------------------------------------- 
   SUBROUTINE dInt3111(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, & 
                              OA,LDA,OB,LDB,OC,LDC,OD,LDD,PBC,I) 
      USE DerivedTypes
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF0
      USE GammaF1
      IMPLICIT REAL(DOUBLE) (A,I,V,W)
      INTEGER        :: LBra,LKet
      REAL(DOUBLE)   :: PrmBufB(5,LBra),PrmBufK(5,LKet)
      TYPE(SmallAtomInfo) :: ACInfo,BDInfo
      TYPE(PBCInfo) :: PBC
      REAL(DOUBLE) :: I(*)
      REAL(DOUBLE)  :: Zeta,Eta,r1xZpE,HfxZpE,r1x2E,r1x2Z,ExZpE,ZxZpE,Omega,Up,Uq,Upq
      REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz
      REAL(DOUBLE)  :: QCx,QCy,QCz,PAx,PAy,PAz,PQx,PQy,PQz,WPx,WPy,WPz,WQx,WQy,WQz   
      REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT,ABx,ABy,ABz,CDx,CDy,CDz
      INTEGER       :: OA,LDA,OB,LDB,OC,LDC,OD,LDD,J,K,L
      REAL(DOUBLE)  :: FPQx,FPQy,FPQz
      I1Bar1=0.0d0
      I2Bar1=0.0d0
      I3Bar1=0.0d0
      I4Bar1=0.0d0
      I5Bar1=0.0d0
      I6Bar1=0.0d0
      I7Bar1=0.0d0
      I8Bar1=0.0d0
      I9Bar1=0.0d0
      I10Bar1=0.0d0
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
               Ia10Bar1=0D0
               Ia9Bar1=0D0
               Ia7Bar1=0D0
               Ia8Bar1=0D0
               Ia6Bar1=0D0
               Ia5Bar1=0D0
               Ib4Bar1=0D0
               Ib10Bar1=0D0
               Ib3Bar1=0D0
               Ib9Bar1=0D0
               Ib7Bar1=0D0
               Ib2Bar1=0D0
               Ib8Bar1=0D0
               Ib6Bar1=0D0
               Ib5Bar1=0D0
               Ic4Bar4=0D0
               Ic4Bar3=0D0
               Ic4Bar2=0D0
               Ic3Bar4=0D0
               Ic3Bar3=0D0
               Ic3Bar2=0D0
               Ic2Bar4=0D0
               Ic2Bar3=0D0
               Ic2Bar2=0D0
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
            V1=AuxR0*PAx
            V2=AuxR1*WPx
            V3=AuxR0*PAy
            V4=AuxR1*WPy
            V5=AuxR0*PAz
            V6=AuxR1*WPz
            V7=-(AuxR1*ExZpE)
            V8=AuxR0+V7
            V9=r1x2Z*V8
            V10=V1+V2
            V11=AuxR1*PAx
            V12=AuxR2*WPx
            V13=V11+V12
            V14=V3+V4
            V15=AuxR1*PAy
            V16=AuxR2*WPy
            V17=V15+V16
            V18=V5+V6
            V19=AuxR1*PAz
            V20=AuxR2*WPz
            V21=V19+V20
            V22=AuxR1*HfxZpE
            I1Bar1=AuxR0+I1Bar1
            I2Bar1=V1+V2+I2Bar1
            I3Bar1=V3+V4+I3Bar1
            I4Bar1=V5+V6+I4Bar1
            I5Bar1=PAx*V10+V9+V13*WPx+I5Bar1
            I6Bar1=PAx*V14+V17*WPx+I6Bar1
            I7Bar1=PAy*V14+V9+V17*WPy+I7Bar1
            I8Bar1=PAx*V18+V21*WPx+I8Bar1
            I9Bar1=PAy*V18+V21*WPy+I9Bar1
            I10Bar1=PAz*V18+V9+V21*WPz+I10Bar1
            I1Bar2=AuxR0*QCx+AuxR1*WQx+I1Bar2
            I2Bar2=QCx*V10+V22+V13*WQx+I2Bar2
            I3Bar2=QCx*V14+V17*WQx+I3Bar2
            I4Bar2=QCx*V18+V21*WQx+I4Bar2
            I1Bar3=AuxR0*QCy+AuxR1*WQy+I1Bar3
            I2Bar3=QCy*V10+V13*WQy+I2Bar3
            I3Bar3=QCy*V14+V22+V17*WQy+I3Bar3
            I4Bar3=QCy*V18+V21*WQy+I4Bar3
            I1Bar4=AuxR0*QCz+AuxR1*WQz+I1Bar4
            I2Bar4=QCz*V10+V13*WQz+I2Bar4
            I3Bar4=QCz*V14+V17*WQz+I3Bar4
            I4Bar4=QCz*V18+V22+V21*WQz+I4Bar4
            Ia10Bar1=Ia10Bar1+Alpha*I10Bar1
            Ia9Bar1=Ia9Bar1+Alpha*I9Bar1
            Ia7Bar1=Ia7Bar1+Alpha*I7Bar1
            Ia8Bar1=Ia8Bar1+Alpha*I8Bar1
            Ia6Bar1=Ia6Bar1+Alpha*I6Bar1
            Ia5Bar1=Ia5Bar1+Alpha*I5Bar1
            Ib4Bar1=Ib4Bar1+Beta*I4Bar1
            Ib10Bar1=Ib10Bar1+Beta*I10Bar1
            Ib3Bar1=Ib3Bar1+Beta*I3Bar1
            Ib9Bar1=Ib9Bar1+Beta*I9Bar1
            Ib7Bar1=Ib7Bar1+Beta*I7Bar1
            Ib2Bar1=Ib2Bar1+Beta*I2Bar1
            Ib8Bar1=Ib8Bar1+Beta*I8Bar1
            Ib6Bar1=Ib6Bar1+Beta*I6Bar1
            Ib5Bar1=Ib5Bar1+Beta*I5Bar1
            Ic4Bar4=Ic4Bar4+Gamma*I4Bar4
            Ic4Bar3=Ic4Bar3+Gamma*I4Bar3
            Ic4Bar2=Ic4Bar2+Gamma*I4Bar2
            Ic3Bar4=Ic3Bar4+Gamma*I3Bar4
            Ic3Bar3=Ic3Bar3+Gamma*I3Bar3
            Ic3Bar2=Ic3Bar2+Gamma*I3Bar2
            Ic2Bar4=Ic2Bar4+Gamma*I2Bar4
            Ic2Bar3=Ic2Bar3+Gamma*I2Bar3
            Ic2Bar2=Ic2Bar2+Gamma*I2Bar2
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      V1=LDA*OA
      V2=LDB*OB
      V3=LDC*OC
      V4=LDD*OD
      V5=V1+V2+V3+V4
      V6=-I1Bar1
      V7=1+OA
      V8=LDA*V7
      V9=V2+V3+V4+V8
      V10=2+OA
      V11=LDA*V10
      V12=V11+V2+V3+V4
      OffSet=V5
      dI(1,OffSet)=Ia5Bar1+V6+dI(1,OffSet)
      dI(4,OffSet)=ABx*Ib2Bar1+Ib5Bar1+dI(4,OffSet)
      dI(7,OffSet)=Ic2Bar2+dI(7,OffSet)
      W1=-dI(1,V5)-dI(4,V5)
      W2=-dI(7,V5)+dI(10,OffSet)
      dI(10,OffSet)=W1+W2
      dI(2,OffSet)=Ia6Bar1+dI(2,OffSet)
      dI(5,OffSet)=ABy*Ib2Bar1+Ib6Bar1+dI(5,OffSet)
      dI(8,OffSet)=Ic2Bar3+dI(8,OffSet)
      W1=-dI(2,V5)-dI(5,V5)
      W2=-dI(8,V5)+dI(11,OffSet)
      dI(11,OffSet)=W1+W2
      dI(3,OffSet)=Ia8Bar1+dI(3,OffSet)
      dI(6,OffSet)=ABz*Ib2Bar1+Ib8Bar1+dI(6,OffSet)
      dI(9,OffSet)=Ic2Bar4+dI(9,OffSet)
      W1=-dI(3,V5)-dI(6,V5)
      W2=-dI(9,V5)+dI(12,OffSet)
      dI(12,OffSet)=W1+W2
      OffSet=LDA*(1.D0+OA)+V2+V3+V4
      dI(1,OffSet)=Ia6Bar1+dI(1,OffSet)
      dI(4,OffSet)=ABx*Ib3Bar1+Ib6Bar1+dI(4,OffSet)
      dI(7,OffSet)=Ic3Bar2+dI(7,OffSet)
      W1=-dI(1,V9)-dI(4,V9)
      W2=-dI(7,V9)+dI(10,OffSet)
      dI(10,OffSet)=W1+W2
      dI(2,OffSet)=Ia7Bar1+V6+dI(2,OffSet)
      dI(5,OffSet)=ABy*Ib3Bar1+Ib7Bar1+dI(5,OffSet)
      dI(8,OffSet)=Ic3Bar3+dI(8,OffSet)
      W1=-dI(2,V9)-dI(5,V9)
      W2=-dI(8,V9)+dI(11,OffSet)
      dI(11,OffSet)=W1+W2
      dI(3,OffSet)=Ia9Bar1+dI(3,OffSet)
      dI(6,OffSet)=ABz*Ib3Bar1+Ib9Bar1+dI(6,OffSet)
      dI(9,OffSet)=Ic3Bar4+dI(9,OffSet)
      W1=-dI(3,V9)-dI(6,V9)
      W2=-dI(9,V9)+dI(12,OffSet)
      dI(12,OffSet)=W1+W2
      OffSet=LDA*(2.D0+OA)+V2+V3+V4
      dI(1,OffSet)=Ia8Bar1+dI(1,OffSet)
      dI(4,OffSet)=ABx*Ib4Bar1+Ib8Bar1+dI(4,OffSet)
      dI(7,OffSet)=Ic4Bar2+dI(7,OffSet)
      W1=-dI(1,V12)-dI(4,V12)
      W2=-dI(7,V12)+dI(10,OffSet)
      dI(10,OffSet)=W1+W2
      dI(2,OffSet)=Ia9Bar1+dI(2,OffSet)
      dI(5,OffSet)=ABy*Ib4Bar1+Ib9Bar1+dI(5,OffSet)
      dI(8,OffSet)=Ic4Bar3+dI(8,OffSet)
      W1=-dI(2,V12)-dI(5,V12)
      W2=-dI(8,V12)+dI(11,OffSet)
      dI(11,OffSet)=W1+W2
      dI(3,OffSet)=Ia10Bar1+V6+dI(3,OffSet)
      dI(6,OffSet)=Ib10Bar1+ABz*Ib4Bar1+dI(6,OffSet)
      dI(9,OffSet)=Ic4Bar4+dI(9,OffSet)
      W1=-dI(3,V12)-dI(6,V12)
      W2=-dI(9,V12)+dI(12,OffSet)
      dI(12,OffSet)=W1+W2
   END SUBROUTINE Int3111
