! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (S S|S S) 
! ---------------------------------------------------------- 
   SUBROUTINE dInt1111(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, & 
                              OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NINT,PBC,dI) 
      USE DerivedTypes
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF0
      USE GammaF1
      IMPLICIT REAL(DOUBLE) (A,I,V,W)
      INTEGER        :: LBra,LKet,NINT
      REAL(DOUBLE)   :: PrmBufB(7,LBra),PrmBufK(7,LKet)
      TYPE(SmallAtomInfo) :: ACInfo,BDInfo
      TYPE(PBCInfo) :: PBC
      REAL(DOUBLE) :: dI(NINT,12)
      REAL(DOUBLE)  :: Zeta,Eta,r1xZpE,HfxZpE,r1x2E,r1x2Z,ExZpE,ZxZpE,Omega,Up,Uq,Upq
      REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz
      REAL(DOUBLE)  :: QCx,QCy,QCz,PAx,PAy,PAz,PQx,PQy,PQz,WPx,WPy,WPz,WQx,WQy,WQz   
      REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT,ABx,ABy,ABz,CDx,CDy,CDz
      REAL(DOUBLE)  :: Alpha,Beta,Gamma
      INTEGER       :: OffSet,CrtSet
      INTEGER       :: OA,LDA,OB,LDB,OC,LDC,OD,LDD,J,K,L
      INTEGER       :: GOA,GOB,GOC,GOD
      REAL(DOUBLE)  :: FPQx,FPQy,FPQz
      CrtSet1=GOA
      CrtSet2=GOA+1
      CrtSet3=GOA+2
      CrtSet4=GOB+3
      CrtSet5=GOB+4
      CrtSet6=GOB+5
      CrtSet7=GOC+6
      CrtSet8=GOC+7
      CrtSet9=GOC+8
      CrtSet10=GOD+9
      CrtSet11=GOD+10
      CrtSet12=GOD+11
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
         Gamma =PrmBufK(6,J)
         QCx=Qx-Cx
         QCy=Qy-Cy
         QCz=Qz-Cz
         DO K=1,LBra ! K^2 VRR (M0| loop 
            Zeta=PrmBufB(1,K)
            Px  =PrmBufB(2,K)
            Py  =PrmBufB(3,K)
            Pz  =PrmBufB(4,K)
            Up  =PrmBufB(5,K)
            Alpha =PrmBufB(6,K)
            Beta  =PrmBufB(7,K)
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
      ! Need to be improved...
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
            RawI1Bar1=AuxR0
            I1Bar1=RawI1Bar1+I1Bar1
            RawI2Bar1=AuxR0*PAx+AuxR1*WPx
            I2Bar1=RawI2Bar1+I2Bar1
            RawI3Bar1=AuxR0*PAy+AuxR1*WPy
            I3Bar1=RawI3Bar1+I3Bar1
            RawI4Bar1=AuxR0*PAz+AuxR1*WPz
            I4Bar1=RawI4Bar1+I4Bar1
            RawI1Bar2=AuxR0*QCx+AuxR1*WQx
            I1Bar2=RawI1Bar2+I1Bar2
            RawI1Bar3=AuxR0*QCy+AuxR1*WQy
            I1Bar3=RawI1Bar3+I1Bar3
            RawI1Bar4=AuxR0*QCz+AuxR1*WQz
            I1Bar4=RawI1Bar4+I1Bar4
            Ia4Bar1=Ia4Bar1+Alpha*RawI4Bar1
            Ia3Bar1=Ia3Bar1+Alpha*RawI3Bar1
            Ia2Bar1=Ia2Bar1+Alpha*RawI2Bar1
            Ib1Bar1=Ib1Bar1+Beta*RawI1Bar1
            Ib4Bar1=Ib4Bar1+Beta*RawI4Bar1
            Ib3Bar1=Ib3Bar1+Beta*RawI3Bar1
            Ib2Bar1=Ib2Bar1+Beta*RawI2Bar1
            Ic1Bar4=Ic1Bar4+Gamma*RawI1Bar4
            Ic1Bar3=Ic1Bar3+Gamma*RawI1Bar3
            Ic1Bar2=Ic1Bar2+Gamma*RawI1Bar2
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD
      dI(OffSet,CrtSet1)=Ia2Bar1+dI(OffSet,CrtSet1)
      dI(OffSet,CrtSet4)=ABx*Ib1Bar1+Ib2Bar1+dI(OffSet,CrtSet4)
      dI(OffSet,CrtSet7)=Ic1Bar2+dI(OffSet,CrtSet7)
      Tmp=-(dI(OffSet,CrtSet1)+dI(OffSet,CrtSet4)+dI(OffSet,CrtSet7))+dI(OffSet,CrtSet10)
      dI(OffSet,CrtSet10)=Tmp
            dI(OffSet,CrtSet2)=Ia3Bar1+dI(OffSet,CrtSet2)
      dI(OffSet,CrtSet5)=ABy*Ib1Bar1+Ib3Bar1+dI(OffSet,CrtSet5)
      dI(OffSet,CrtSet8)=Ic1Bar3+dI(OffSet,CrtSet8)
      Tmp=-(dI(OffSet,CrtSet2)+dI(OffSet,CrtSet5)+dI(OffSet,CrtSet8))+dI(OffSet,CrtSet11)
      dI(OffSet,CrtSet11)=Tmp
            dI(OffSet,CrtSet3)=Ia4Bar1+dI(OffSet,CrtSet3)
      dI(OffSet,CrtSet6)=ABz*Ib1Bar1+Ib4Bar1+dI(OffSet,CrtSet6)
      dI(OffSet,CrtSet9)=Ic1Bar4+dI(OffSet,CrtSet9)
      Tmp=-(dI(OffSet,CrtSet3)+dI(OffSet,CrtSet6)+dI(OffSet,CrtSet9))+dI(OffSet,CrtSet12)
      dI(OffSet,CrtSet12)=Tmp
   END SUBROUTINE dInt1111
