! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (S S|P S) 
! ---------------------------------------------------------- 
   SUBROUTINE dInt1131(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, & 
                              OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NINT,PBC,dI) 
      USE DerivedTypes
      USE VScratch
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF1
      USE GammaF2
      IMPLICIT REAL(DOUBLE) (A,I,W,R)
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
      INTEGER       :: CrtSet1,CrtSet2,CrtSet3,CrtSet4 ,CrtSet5 ,CrtSet6
      INTEGER       :: CrtSet7,CrtSet8,CrtSet9,CrtSet10,CrtSet11,CrtSet12
      INTEGER       :: OffSet,GOA,GOB,GOC,GOD
      INTEGER       :: OA,LDA,OB,LDB,OC,LDC,OD,LDD,J,K,L
      REAL(DOUBLE)  :: FPQx,FPQy,FPQz
      CrtSet1=GOA
      CrtSet2=GOA+1
      CrtSet3=GOA+2
      CrtSet4=GOB
      CrtSet5=GOB+1
      CrtSet6=GOB+2
      CrtSet7=GOC
      CrtSet8=GOC+1
      CrtSet9=GOC+2
      CrtSet10=GOD
      CrtSet11=GOD+1
      CrtSet12=GOD+2
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
      I1Bar5=0.0d0
      I1Bar6=0.0d0
      I1Bar7=0.0d0
      I1Bar8=0.0d0
      I1Bar9=0.0d0
      I1Bar10=0.0d0
      Ia4Bar4=0D0
      Ia4Bar3=0D0
      Ia4Bar2=0D0
      Ia3Bar4=0D0
      Ia3Bar3=0D0
      Ia3Bar2=0D0
      Ia2Bar4=0D0
      Ia2Bar3=0D0
      Ia2Bar2=0D0
      Ib1Bar4=0D0
      Ib1Bar3=0D0
      Ib1Bar2=0D0
      Ib4Bar4=0D0
      Ib4Bar3=0D0
      Ib4Bar2=0D0
      Ib3Bar4=0D0
      Ib3Bar3=0D0
      Ib3Bar2=0D0
      Ib2Bar4=0D0
      Ib2Bar3=0D0
      Ib2Bar2=0D0
      Ic1Bar10=0D0
      Ic1Bar9=0D0
      Ic1Bar7=0D0
      Ic1Bar8=0D0
      Ic1Bar6=0D0
      Ic1Bar5=0D0
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
      V(3)=V(1)+V(2)
      V(4)=AuxR0*PAy
      V(5)=AuxR1*WPy
      V(6)=V(4)+V(5)
      V(7)=AuxR0*PAz
      V(8)=AuxR1*WPz
      V(9)=V(7)+V(8)
      V(10)=AuxR0*QCx
      V(11)=AuxR1*WQx
      V(12)=V(10)+V(11)
      V(13)=AuxR1*HfxZpE
      V(14)=AuxR1*PAx
      V(15)=AuxR2*WPx
      V(16)=V(14)+V(15)
      V(17)=AuxR1*PAy
      V(18)=AuxR2*WPy
      V(19)=V(17)+V(18)
      V(20)=AuxR1*PAz
      V(21)=AuxR2*WPz
      V(22)=V(20)+V(21)
      V(23)=AuxR0*QCy
      V(24)=AuxR1*WQy
      V(25)=V(23)+V(24)
      V(26)=AuxR0*QCz
      V(27)=AuxR1*WQz
      V(28)=V(26)+V(27)
      V(29)=AuxR1*ZxZpE
      V(30)=-V(29)
      V(31)=AuxR0+V(30)
      V(32)=r1x2E*V(31)
      V(33)=AuxR1*QCy
      V(34)=AuxR2*WQy
      V(35)=V(33)+V(34)
      V(36)=AuxR1*QCz
      V(37)=AuxR2*WQz
      V(38)=V(36)+V(37)
      RawI1Bar1=AuxR0
      I1Bar1=RawI1Bar1+I1Bar1
      RawI2Bar1=V(3)
      I2Bar1=RawI2Bar1+I2Bar1
      RawI3Bar1=V(6)
      I3Bar1=RawI3Bar1+I3Bar1
      RawI4Bar1=V(9)
      I4Bar1=RawI4Bar1+I4Bar1
      RawI1Bar2=V(12)
      I1Bar2=RawI1Bar2+I1Bar2
      RawI2Bar2=QCx*V(3)+V(13)+WQx*V(16)
      I2Bar2=RawI2Bar2+I2Bar2
      RawI3Bar2=QCx*V(6)+WQx*V(19)
      I3Bar2=RawI3Bar2+I3Bar2
      RawI4Bar2=QCx*V(9)+WQx*V(22)
      I4Bar2=RawI4Bar2+I4Bar2
      RawI1Bar3=V(25)
      I1Bar3=RawI1Bar3+I1Bar3
      RawI2Bar3=QCy*V(3)+WQy*V(16)
      I2Bar3=RawI2Bar3+I2Bar3
      RawI3Bar3=QCy*V(6)+V(13)+WQy*V(19)
      I3Bar3=RawI3Bar3+I3Bar3
      RawI4Bar3=QCy*V(9)+WQy*V(22)
      I4Bar3=RawI4Bar3+I4Bar3
      RawI1Bar4=V(28)
      I1Bar4=RawI1Bar4+I1Bar4
      RawI2Bar4=QCz*V(3)+WQz*V(16)
      I2Bar4=RawI2Bar4+I2Bar4
      RawI3Bar4=QCz*V(6)+WQz*V(19)
      I3Bar4=RawI3Bar4+I3Bar4
      RawI4Bar4=QCz*V(9)+V(13)+WQz*V(22)
      I4Bar4=RawI4Bar4+I4Bar4
      W1=WQx*(AuxR1*QCx+AuxR2*WQx)
      W2=QCx*V(12)+V(32)
      RawI1Bar5=W1+W2
      I1Bar5=RawI1Bar5+I1Bar5
      RawI1Bar6=QCx*V(25)+WQx*V(35)
      I1Bar6=RawI1Bar6+I1Bar6
      RawI1Bar7=QCy*V(25)+V(32)+WQy*V(35)
      I1Bar7=RawI1Bar7+I1Bar7
      RawI1Bar8=QCx*V(28)+WQx*V(38)
      I1Bar8=RawI1Bar8+I1Bar8
      RawI1Bar9=QCy*V(28)+WQy*V(38)
      I1Bar9=RawI1Bar9+I1Bar9
      RawI1Bar10=QCz*V(28)+V(32)+WQz*V(38)
      I1Bar10=RawI1Bar10+I1Bar10
            Ia4Bar4=Ia4Bar4+Alpha*RawI4Bar4
            Ia4Bar3=Ia4Bar3+Alpha*RawI4Bar3
            Ia4Bar2=Ia4Bar2+Alpha*RawI4Bar2
            Ia3Bar4=Ia3Bar4+Alpha*RawI3Bar4
            Ia3Bar3=Ia3Bar3+Alpha*RawI3Bar3
            Ia3Bar2=Ia3Bar2+Alpha*RawI3Bar2
            Ia2Bar4=Ia2Bar4+Alpha*RawI2Bar4
            Ia2Bar3=Ia2Bar3+Alpha*RawI2Bar3
            Ia2Bar2=Ia2Bar2+Alpha*RawI2Bar2
            Ib1Bar4=Ib1Bar4+Beta*RawI1Bar4
            Ib1Bar3=Ib1Bar3+Beta*RawI1Bar3
            Ib1Bar2=Ib1Bar2+Beta*RawI1Bar2
            Ib4Bar4=Ib4Bar4+Beta*RawI4Bar4
            Ib4Bar3=Ib4Bar3+Beta*RawI4Bar3
            Ib4Bar2=Ib4Bar2+Beta*RawI4Bar2
            Ib3Bar4=Ib3Bar4+Beta*RawI3Bar4
            Ib3Bar3=Ib3Bar3+Beta*RawI3Bar3
            Ib3Bar2=Ib3Bar2+Beta*RawI3Bar2
            Ib2Bar4=Ib2Bar4+Beta*RawI2Bar4
            Ib2Bar3=Ib2Bar3+Beta*RawI2Bar3
            Ib2Bar2=Ib2Bar2+Beta*RawI2Bar2
            Ic1Bar10=Ic1Bar10+Gamma*RawI1Bar10
            Ic1Bar9=Ic1Bar9+Gamma*RawI1Bar9
            Ic1Bar7=Ic1Bar7+Gamma*RawI1Bar7
            Ic1Bar8=Ic1Bar8+Gamma*RawI1Bar8
            Ic1Bar6=Ic1Bar6+Gamma*RawI1Bar6
            Ic1Bar5=Ic1Bar5+Gamma*RawI1Bar5
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      V(1)=ABx*Ib1Bar2
      V(2)=-I1Bar1
      V(3)=ABy*Ib1Bar2
      V(4)=-Ic1Bar6
      V(5)=ABz*Ib1Bar2
      V(6)=-Ic1Bar8
      V(7)=ABx*Ib1Bar3
      V(8)=ABy*Ib1Bar3
      V(9)=ABz*Ib1Bar3
      V(10)=-Ic1Bar9
      V(11)=ABx*Ib1Bar4
      V(12)=ABy*Ib1Bar4
      V(13)=ABz*Ib1Bar4
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD !=ZippyForPres
      dI(OffSet,CrtSet1)=Ia2Bar2+dI(OffSet,CrtSet1)
      dI(OffSet,CrtSet4)=Ib2Bar2+dI(OffSet,CrtSet4)+V(1)
      dI(OffSet,CrtSet7)=Ic1Bar5+dI(OffSet,CrtSet7)+V(2)
      W1=I1Bar1-Ia2Bar2-Ib2Bar2
      W2=-Ic1Bar5+dI(OffSet,CrtSet10)-V(1)
      dI(OffSet,CrtSet10)=W1+W2
      dI(OffSet,CrtSet2)=Ia3Bar2+dI(OffSet,CrtSet2)
      dI(OffSet,CrtSet5)=Ib3Bar2+dI(OffSet,CrtSet5)+V(3)
      dI(OffSet,CrtSet8)=Ic1Bar6+dI(OffSet,CrtSet8)
      W1=-Ia3Bar2-Ib3Bar2
      W2=dI(OffSet,CrtSet11)-V(3)+V(4)
      dI(OffSet,CrtSet11)=W1+W2
      dI(OffSet,CrtSet3)=Ia4Bar2+dI(OffSet,CrtSet3)
      dI(OffSet,CrtSet6)=Ib4Bar2+dI(OffSet,CrtSet6)+V(5)
      dI(OffSet,CrtSet9)=Ic1Bar8+dI(OffSet,CrtSet9)
      W1=-Ia4Bar2-Ib4Bar2
      W2=dI(OffSet,CrtSet12)-V(5)+V(6)
      dI(OffSet,CrtSet12)=W1+W2
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+1)*LDC+(OD+0)*LDD !=ZippyForPres
      dI(OffSet,CrtSet1)=Ia2Bar3+dI(OffSet,CrtSet1)
      dI(OffSet,CrtSet4)=Ib2Bar3+dI(OffSet,CrtSet4)+V(7)
      dI(OffSet,CrtSet7)=Ic1Bar6+dI(OffSet,CrtSet7)
      W1=-Ia2Bar3-Ib2Bar3
      W2=dI(OffSet,CrtSet10)+V(4)-V(7)
      dI(OffSet,CrtSet10)=W1+W2
      dI(OffSet,CrtSet2)=Ia3Bar3+dI(OffSet,CrtSet2)
      dI(OffSet,CrtSet5)=Ib3Bar3+dI(OffSet,CrtSet5)+V(8)
      dI(OffSet,CrtSet8)=Ic1Bar7+dI(OffSet,CrtSet8)+V(2)
      W1=I1Bar1-Ia3Bar3-Ib3Bar3
      W2=-Ic1Bar7+dI(OffSet,CrtSet11)-V(8)
      dI(OffSet,CrtSet11)=W1+W2
      dI(OffSet,CrtSet3)=Ia4Bar3+dI(OffSet,CrtSet3)
      dI(OffSet,CrtSet6)=Ib4Bar3+dI(OffSet,CrtSet6)+V(9)
      dI(OffSet,CrtSet9)=Ic1Bar9+dI(OffSet,CrtSet9)
      W1=-Ia4Bar3-Ib4Bar3
      W2=dI(OffSet,CrtSet12)-V(9)+V(10)
      dI(OffSet,CrtSet12)=W1+W2
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+2)*LDC+(OD+0)*LDD !=ZippyForPres
      dI(OffSet,CrtSet1)=Ia2Bar4+dI(OffSet,CrtSet1)
      dI(OffSet,CrtSet4)=Ib2Bar4+dI(OffSet,CrtSet4)+V(11)
      dI(OffSet,CrtSet7)=Ic1Bar8+dI(OffSet,CrtSet7)
      W1=-Ia2Bar4-Ib2Bar4
      W2=dI(OffSet,CrtSet10)+V(6)-V(11)
      dI(OffSet,CrtSet10)=W1+W2
      dI(OffSet,CrtSet2)=Ia3Bar4+dI(OffSet,CrtSet2)
      dI(OffSet,CrtSet5)=Ib3Bar4+dI(OffSet,CrtSet5)+V(12)
      dI(OffSet,CrtSet8)=Ic1Bar9+dI(OffSet,CrtSet8)
      W1=-Ia3Bar4-Ib3Bar4
      W2=dI(OffSet,CrtSet11)+V(10)-V(12)
      dI(OffSet,CrtSet11)=W1+W2
      dI(OffSet,CrtSet3)=Ia4Bar4+dI(OffSet,CrtSet3)
      dI(OffSet,CrtSet6)=Ib4Bar4+dI(OffSet,CrtSet6)+V(13)
      dI(OffSet,CrtSet9)=Ic1Bar10+dI(OffSet,CrtSet9)+V(2)
      W1=I1Bar1-Ia4Bar4-Ib4Bar4
      W2=-Ic1Bar10+dI(OffSet,CrtSet12)-V(13)
      dI(OffSet,CrtSet12)=W1+W2
   END SUBROUTINE dInt1131
