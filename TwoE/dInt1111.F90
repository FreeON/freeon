! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (S S|S S) 
! ---------------------------------------------------------- 
   SUBROUTINE dInt1111(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, & 
                              OA,LDA,OB,LDB,OC,LDC,OD,LDD,PBC,I) 
      USE DerivedTypes
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF0
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
      I1Bar2=0.0d0
      I1Bar3=0.0d0
      I1Bar4=0.0d0
               Ia4Bar1=0D0
               Ia3Bar1=0D0
               Ia2Bar1=0D0
               Ib1Bar1=0D0
               Ib4Bar1=0D0
               Ib3Bar1=0D0
               Ib2Bar1=0D0
               Ic1Bar4=0D0
               Ic1Bar3=0D0
               Ic1Bar2=0D0
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
              AuxR0=Upq*(F0_0(L)+T*(F0_1(L)+T*(F0_2(L)+T*(F0_3(L)+T*F0_4(L)))))
              AuxR1=Upq*(F1_0(L)+T*(F1_1(L)+T*(F1_2(L)+T*(F1_3(L)+T*F1_4(L)))))
            ELSE
              InvT=One/T
              SqInvT=DSQRT(InvT)
              AuxR0=+8.862269254527580D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              AuxR1=+4.431134627263790D-01*Upq*SqInvT
            ENDIF
            I1Bar1=AuxR0+I1Bar1
            I2Bar1=AuxR0*PAx+AuxR1*WPx+I2Bar1
            I3Bar1=AuxR0*PAy+AuxR1*WPy+I3Bar1
            I4Bar1=AuxR0*PAz+AuxR1*WPz+I4Bar1
            I1Bar2=AuxR0*QCx+AuxR1*WQx+I1Bar2
            I1Bar3=AuxR0*QCy+AuxR1*WQy+I1Bar3
            I1Bar4=AuxR0*QCz+AuxR1*WQz+I1Bar4
            Ia4Bar1=Ia4Bar1+Alpha*I4Bar1
            Ia3Bar1=Ia3Bar1+Alpha*I3Bar1
            Ia2Bar1=Ia2Bar1+Alpha*I2Bar1
            Ib1Bar1=Ib1Bar1+Beta*I1Bar1
            Ib4Bar1=Ib4Bar1+Beta*I4Bar1
            Ib3Bar1=Ib3Bar1+Beta*I3Bar1
            Ib2Bar1=Ib2Bar1+Beta*I2Bar1
            Ic1Bar4=Ic1Bar4+Gamma*I1Bar4
            Ic1Bar3=Ic1Bar3+Gamma*I1Bar3
            Ic1Bar2=Ic1Bar2+Gamma*I1Bar2
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      V1=LDA*OA
      V2=LDB*OB
      V3=LDC*OC
      V4=LDD*OD
      V5=V1+V2+V3+V4
      OffSet=V5
      dI(1,OffSet)=Ia2Bar1+dI(1,OffSet)
      dI(4,OffSet)=ABx*Ib1Bar1+Ib2Bar1+dI(4,OffSet)
      dI(7,OffSet)=Ic1Bar2+dI(7,OffSet)
      W1=-dI(1,V5)-dI(4,V5)
      W2=-dI(7,V5)+dI(10,OffSet)
      dI(10,OffSet)=W1+W2
      dI(2,OffSet)=Ia3Bar1+dI(2,OffSet)
      dI(5,OffSet)=ABy*Ib1Bar1+Ib3Bar1+dI(5,OffSet)
      dI(8,OffSet)=Ic1Bar3+dI(8,OffSet)
      W1=-dI(2,V5)-dI(5,V5)
      W2=-dI(8,V5)+dI(11,OffSet)
      dI(11,OffSet)=W1+W2
      dI(3,OffSet)=Ia4Bar1+dI(3,OffSet)
      dI(6,OffSet)=ABz*Ib1Bar1+Ib4Bar1+dI(6,OffSet)
      dI(9,OffSet)=Ic1Bar4+dI(9,OffSet)
      W1=-dI(3,V5)-dI(6,V5)
      W2=-dI(9,V5)+dI(12,OffSet)
      dI(12,OffSet)=W1+W2
   END SUBROUTINE Int1111
