! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (P S|S S) 
! ---------------------------------------------------------- 
   SUBROUTINE dInt3111(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, & 
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
      REAL(DOUBLE)  :: FPQx,FPQy,FPQz,Dum
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
            IF(PBC%AutoW%I(1)==1) FPQx=FPQx-DNINT(FPQx-SIGN(1.D0,FPQx)*1.D-14)
            IF(PBC%AutoW%I(2)==1) FPQy=FPQy-DNINT(FPQy-SIGN(1.D0,FPQy)*1.D-14)
            IF(PBC%AutoW%I(3)==1) FPQz=FPQz-DNINT(FPQz-SIGN(1.D0,FPQz)*1.D-14)
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
      V(10)=AuxR1*ExZpE
      V(11)=-V(10)
      V(12)=AuxR0+V(11)
      V(13)=r1x2Z*V(12)
      V(14)=AuxR1*PAx
      V(15)=AuxR2*WPx
      V(16)=V(14)+V(15)
      V(17)=AuxR1*PAy
      V(18)=AuxR2*WPy
      V(19)=V(17)+V(18)
      V(20)=AuxR1*PAz
      V(21)=AuxR2*WPz
      V(22)=V(20)+V(21)
      V(23)=AuxR1*HfxZpE
      RawI1Bar1=AuxR0
      I1Bar1=RawI1Bar1+I1Bar1
      RawI2Bar1=V(3)
      I2Bar1=RawI2Bar1+I2Bar1
      RawI3Bar1=V(6)
      I3Bar1=RawI3Bar1+I3Bar1
      RawI4Bar1=V(9)
      I4Bar1=RawI4Bar1+I4Bar1
      RawI5Bar1=PAx*V(3)+V(13)+WPx*V(16)
      I5Bar1=RawI5Bar1+I5Bar1
      RawI6Bar1=PAx*V(6)+WPx*V(19)
      I6Bar1=RawI6Bar1+I6Bar1
      RawI7Bar1=PAy*V(6)+V(13)+WPy*V(19)
      I7Bar1=RawI7Bar1+I7Bar1
      RawI8Bar1=PAx*V(9)+WPx*V(22)
      I8Bar1=RawI8Bar1+I8Bar1
      RawI9Bar1=PAy*V(9)+WPy*V(22)
      I9Bar1=RawI9Bar1+I9Bar1
      RawI10Bar1=PAz*V(9)+V(13)+WPz*V(22)
      I10Bar1=RawI10Bar1+I10Bar1
      RawI1Bar2=AuxR0*QCx+AuxR1*WQx
      I1Bar2=RawI1Bar2+I1Bar2
      RawI2Bar2=QCx*V(3)+WQx*V(16)+V(23)
      I2Bar2=RawI2Bar2+I2Bar2
      RawI3Bar2=QCx*V(6)+WQx*V(19)
      I3Bar2=RawI3Bar2+I3Bar2
      RawI4Bar2=QCx*V(9)+WQx*V(22)
      I4Bar2=RawI4Bar2+I4Bar2
      RawI1Bar3=AuxR0*QCy+AuxR1*WQy
      I1Bar3=RawI1Bar3+I1Bar3
      RawI2Bar3=QCy*V(3)+WQy*V(16)
      I2Bar3=RawI2Bar3+I2Bar3
      RawI3Bar3=QCy*V(6)+WQy*V(19)+V(23)
      I3Bar3=RawI3Bar3+I3Bar3
      RawI4Bar3=QCy*V(9)+WQy*V(22)
      I4Bar3=RawI4Bar3+I4Bar3
      RawI1Bar4=AuxR0*QCz+AuxR1*WQz
      I1Bar4=RawI1Bar4+I1Bar4
      RawI2Bar4=QCz*V(3)+WQz*V(16)
      I2Bar4=RawI2Bar4+I2Bar4
      RawI3Bar4=QCz*V(6)+WQz*V(19)
      I3Bar4=RawI3Bar4+I3Bar4
      RawI4Bar4=QCz*V(9)+WQz*V(22)+V(23)
      I4Bar4=RawI4Bar4+I4Bar4
            Ia10Bar1=Ia10Bar1+Alpha*RawI10Bar1
            Ia9Bar1=Ia9Bar1+Alpha*RawI9Bar1
            Ia7Bar1=Ia7Bar1+Alpha*RawI7Bar1
            Ia8Bar1=Ia8Bar1+Alpha*RawI8Bar1
            Ia6Bar1=Ia6Bar1+Alpha*RawI6Bar1
            Ia5Bar1=Ia5Bar1+Alpha*RawI5Bar1
            Ib4Bar1=Ib4Bar1+Beta*RawI4Bar1
            Ib10Bar1=Ib10Bar1+Beta*RawI10Bar1
            Ib3Bar1=Ib3Bar1+Beta*RawI3Bar1
            Ib9Bar1=Ib9Bar1+Beta*RawI9Bar1
            Ib7Bar1=Ib7Bar1+Beta*RawI7Bar1
            Ib2Bar1=Ib2Bar1+Beta*RawI2Bar1
            Ib8Bar1=Ib8Bar1+Beta*RawI8Bar1
            Ib6Bar1=Ib6Bar1+Beta*RawI6Bar1
            Ib5Bar1=Ib5Bar1+Beta*RawI5Bar1
            Ic4Bar4=Ic4Bar4+Gamma*RawI4Bar4
            Ic4Bar3=Ic4Bar3+Gamma*RawI4Bar3
            Ic4Bar2=Ic4Bar2+Gamma*RawI4Bar2
            Ic3Bar4=Ic3Bar4+Gamma*RawI3Bar4
            Ic3Bar3=Ic3Bar3+Gamma*RawI3Bar3
            Ic3Bar2=Ic3Bar2+Gamma*RawI3Bar2
            Ic2Bar4=Ic2Bar4+Gamma*RawI2Bar4
            Ic2Bar3=Ic2Bar3+Gamma*RawI2Bar3
            Ic2Bar2=Ic2Bar2+Gamma*RawI2Bar2
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      V(1)=-I1Bar1
      V(2)=ABx*Ib2Bar1
      V(3)=ABy*Ib2Bar1
      V(4)=-Ia6Bar1
      V(5)=-Ib6Bar1
      V(6)=ABz*Ib2Bar1
      V(7)=-Ia8Bar1
      V(8)=-Ib8Bar1
      V(9)=ABx*Ib3Bar1
      V(10)=ABy*Ib3Bar1
      V(11)=ABz*Ib3Bar1
      V(12)=-Ia9Bar1
      V(13)=-Ib9Bar1
      V(14)=ABx*Ib4Bar1
      V(15)=ABy*Ib4Bar1
      V(16)=ABz*Ib4Bar1
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD !=ZippyForPres
      dI(OffSet,CrtSet1)=Ia5Bar1+dI(OffSet,CrtSet1)+V(1)
      dI(OffSet,CrtSet4)=Ib5Bar1+dI(OffSet,CrtSet4)+V(2)
      dI(OffSet,CrtSet7)=Ic2Bar2+dI(OffSet,CrtSet7)
      W1=I1Bar1-Ia5Bar1-Ib5Bar1
      W2=-Ic2Bar2+dI(OffSet,CrtSet10)-V(2)
      dI(OffSet,CrtSet10)=W1+W2
      dI(OffSet,CrtSet2)=Ia6Bar1+dI(OffSet,CrtSet2)
      dI(OffSet,CrtSet5)=Ib6Bar1+dI(OffSet,CrtSet5)+V(3)
      dI(OffSet,CrtSet8)=Ic2Bar3+dI(OffSet,CrtSet8)
      W1=-Ic2Bar3+dI(OffSet,CrtSet11)
      W2=-V(3)+V(4)+V(5)
      dI(OffSet,CrtSet11)=W1+W2
      dI(OffSet,CrtSet3)=Ia8Bar1+dI(OffSet,CrtSet3)
      dI(OffSet,CrtSet6)=Ib8Bar1+dI(OffSet,CrtSet6)+V(6)
      dI(OffSet,CrtSet9)=Ic2Bar4+dI(OffSet,CrtSet9)
      W1=-Ic2Bar4+dI(OffSet,CrtSet12)
      W2=-V(6)+V(7)+V(8)
      dI(OffSet,CrtSet12)=W1+W2
      OffSet=(OA+1)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD !=ZippyForPres
      dI(OffSet,CrtSet1)=Ia6Bar1+dI(OffSet,CrtSet1)
      dI(OffSet,CrtSet4)=Ib6Bar1+dI(OffSet,CrtSet4)+V(9)
      dI(OffSet,CrtSet7)=Ic3Bar2+dI(OffSet,CrtSet7)
      W1=-Ic3Bar2+dI(OffSet,CrtSet10)
      W2=V(4)+V(5)-V(9)
      dI(OffSet,CrtSet10)=W1+W2
      dI(OffSet,CrtSet2)=Ia7Bar1+dI(OffSet,CrtSet2)+V(1)
      dI(OffSet,CrtSet5)=Ib7Bar1+dI(OffSet,CrtSet5)+V(10)
      dI(OffSet,CrtSet8)=Ic3Bar3+dI(OffSet,CrtSet8)
      W1=I1Bar1-Ia7Bar1-Ib7Bar1
      W2=-Ic3Bar3+dI(OffSet,CrtSet11)-V(10)
      dI(OffSet,CrtSet11)=W1+W2
      dI(OffSet,CrtSet3)=Ia9Bar1+dI(OffSet,CrtSet3)
      dI(OffSet,CrtSet6)=Ib9Bar1+dI(OffSet,CrtSet6)+V(11)
      dI(OffSet,CrtSet9)=Ic3Bar4+dI(OffSet,CrtSet9)
      W1=-Ic3Bar4+dI(OffSet,CrtSet12)
      W2=-V(11)+V(12)+V(13)
      dI(OffSet,CrtSet12)=W1+W2
      OffSet=(OA+2)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD !=ZippyForPres
      dI(OffSet,CrtSet1)=Ia8Bar1+dI(OffSet,CrtSet1)
      dI(OffSet,CrtSet4)=Ib8Bar1+dI(OffSet,CrtSet4)+V(14)
      dI(OffSet,CrtSet7)=Ic4Bar2+dI(OffSet,CrtSet7)
      W1=-Ic4Bar2+dI(OffSet,CrtSet10)
      W2=V(7)+V(8)-V(14)
      dI(OffSet,CrtSet10)=W1+W2
      dI(OffSet,CrtSet2)=Ia9Bar1+dI(OffSet,CrtSet2)
      dI(OffSet,CrtSet5)=Ib9Bar1+dI(OffSet,CrtSet5)+V(15)
      dI(OffSet,CrtSet8)=Ic4Bar3+dI(OffSet,CrtSet8)
      W1=-Ic4Bar3+dI(OffSet,CrtSet11)
      W2=V(12)+V(13)-V(15)
      dI(OffSet,CrtSet11)=W1+W2
      dI(OffSet,CrtSet3)=Ia10Bar1+dI(OffSet,CrtSet3)+V(1)
      dI(OffSet,CrtSet6)=Ib10Bar1+dI(OffSet,CrtSet6)+V(16)
      dI(OffSet,CrtSet9)=Ic4Bar4+dI(OffSet,CrtSet9)
      W1=I1Bar1-Ia10Bar1-Ib10Bar1
      W2=-Ic4Bar4+dI(OffSet,CrtSet12)-V(16)
      dI(OffSet,CrtSet12)=W1+W2
   END SUBROUTINE dInt3111
