! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (D P|S S) 
! ---------------------------------------------------------- 
   SUBROUTINE Int6311(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, & 
                              OA,LDA,OB,LDB,OC,LDC,OD,LDD,PBC,I) 
      USE DerivedTypes
      USE VScratch
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF3
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
      I5Bar1=0.0d0
      I6Bar1=0.0d0
      I7Bar1=0.0d0
      I8Bar1=0.0d0
      I9Bar1=0.0d0
      I10Bar1=0.0d0
      I11Bar1=0.0d0
      I12Bar1=0.0d0
      I13Bar1=0.0d0
      I14Bar1=0.0d0
      I15Bar1=0.0d0
      I16Bar1=0.0d0
      I17Bar1=0.0d0
      I18Bar1=0.0d0
      I19Bar1=0.0d0
      I20Bar1=0.0d0
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
              W3=(F3_0(L)+T*(F3_1(L)+T*(F3_2(L)+T*(F3_3(L)+T*F3_4(L)))))
              W2=+2.000000000000000D-01*(TwoT*W3+ET)
              W1=+3.333333333333333D-01*(TwoT*W2+ET)
              W0=TwoT*W1+ET
              AuxR0=Upq*W0
              AuxR1=Upq*W1
              AuxR2=Upq*W2
              AuxR3=Upq*W3
            ELSE
              InvT=One/T
              SqInvT=DSQRT(InvT)
              AuxR0=+8.862269254527580D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              AuxR1=+4.431134627263790D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              AuxR2=+6.646701940895685D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              AuxR3=+1.661675485223921D+00*Upq*SqInvT
            ENDIF
      V(1)=AuxR0*PAx
      V(2)=AuxR1*WPx
      V(3)=AuxR0*PAy
      V(4)=AuxR1*WPy
      V(5)=AuxR0*PAz
      V(6)=AuxR1*WPz
      V(7)=AuxR1*ExZpE
      V(8)=-V(7)
      V(9)=AuxR0+V(8)
      V(10)=r1x2Z*V(9)
      V(11)=V(1)+V(2)
      V(12)=PAx*V(11)
      V(13)=AuxR1*PAx
      V(14)=AuxR2*WPx
      V(15)=V(13)+V(14)
      V(16)=WPx*V(15)
      V(17)=V(3)+V(4)
      V(18)=PAx*V(17)
      V(19)=AuxR1*PAy
      V(20)=AuxR2*WPy
      V(21)=V(19)+V(20)
      V(22)=WPx*V(21)
      V(23)=PAy*V(17)
      V(24)=WPy*V(21)
      V(25)=V(5)+V(6)
      V(26)=PAx*V(25)
      V(27)=AuxR1*PAz
      V(28)=AuxR2*WPz
      V(29)=V(27)+V(28)
      V(30)=WPx*V(29)
      V(31)=PAy*V(25)
      V(32)=WPy*V(29)
      V(33)=PAz*V(25)
      V(34)=WPz*V(29)
      V(35)=ExZpE*V(15)
      V(36)=-V(35)
      V(37)=V(1)+V(2)+V(36)
      V(38)=r1x2Z*V(37)
      V(39)=AuxR2*ExZpE
      V(40)=-V(39)
      V(41)=AuxR1+V(40)
      V(42)=r1x2Z*V(41)
      V(43)=ExZpE*V(21)
      V(44)=-V(43)
      V(45)=V(3)+V(4)+V(44)
      V(46)=r1x2Z*V(45)
      V(47)=V(18)+V(22)
      V(48)=PAx*V(21)
      V(49)=AuxR2*PAy
      V(50)=AuxR3*WPy
      V(51)=V(49)+V(50)
      V(52)=WPx*V(51)
      V(53)=V(48)+V(52)
      V(54)=ExZpE*V(29)
      V(55)=-V(54)
      V(56)=V(5)+V(6)+V(55)
      V(57)=r1x2Z*V(56)
      V(58)=V(26)+V(30)
      V(59)=PAx*V(29)
      V(60)=AuxR2*PAz
      V(61)=AuxR3*WPz
      V(62)=V(60)+V(61)
      V(63)=WPx*V(62)
      V(64)=V(59)+V(63)
      V(65)=V(31)+V(32)
      V(66)=PAy*V(29)
      V(67)=WPy*V(62)
      V(68)=V(66)+V(67)
      I1Bar1=AuxR0+I1Bar1
      I2Bar1=I2Bar1+V(1)+V(2)
      I3Bar1=I3Bar1+V(3)+V(4)
      I4Bar1=I4Bar1+V(5)+V(6)
      I5Bar1=I5Bar1+V(10)+V(12)+V(16)
      I6Bar1=I6Bar1+V(18)+V(22)
      I7Bar1=I7Bar1+V(10)+V(23)+V(24)
      I8Bar1=I8Bar1+V(26)+V(30)
      I9Bar1=I9Bar1+V(31)+V(32)
      I10Bar1=I10Bar1+V(10)+V(33)+V(34)
      W1=I11Bar1+PAx*(V(10)+V(12)+V(16))
      W2=2.D0*V(38)
      W3=WPx*(WPx*(AuxR2*PAx+AuxR3*WPx)+PAx*V(15)+V(42))
      I11Bar1=W1+W2+W3
      I12Bar1=I12Bar1+V(46)+PAx*V(47)+WPx*V(53)
      I13Bar1=I13Bar1+V(38)+PAy*V(47)+WPy*V(53)
      W1=I14Bar1+PAy*(V(10)+V(23)+V(24))
      W2=2.D0*V(46)
      W3=WPy*(PAy*V(21)+V(42)+WPy*V(51))
      I14Bar1=W1+W2+W3
      I15Bar1=I15Bar1+V(57)+PAx*V(58)+WPx*V(64)
      I16Bar1=I16Bar1+PAx*V(65)+WPx*V(68)
      I17Bar1=I17Bar1+V(57)+PAy*V(65)+WPy*V(68)
      I18Bar1=I18Bar1+V(38)+PAz*V(58)+WPz*V(64)
      I19Bar1=I19Bar1+V(46)+PAz*V(65)+WPz*V(68)
      W1=I20Bar1+PAz*(V(10)+V(33)+V(34))
      W2=2.D0*V(57)
      W3=WPz*(PAz*V(29)+V(42)+WPz*V(62))
      I20Bar1=W1+W2+W3
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I11Bar1+ABx*I5Bar1+I(OffSet)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I12Bar1+ABx*I6Bar1+I(OffSet)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I13Bar1+ABx*I7Bar1+I(OffSet)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I15Bar1+ABx*I8Bar1+I(OffSet)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I16Bar1+ABx*I9Bar1+I(OffSet)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=ABx*I10Bar1+I18Bar1+I(OffSet)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I12Bar1+ABy*I5Bar1+I(OffSet)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I13Bar1+ABy*I6Bar1+I(OffSet)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I14Bar1+ABy*I7Bar1+I(OffSet)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I16Bar1+ABy*I8Bar1+I(OffSet)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I17Bar1+ABy*I9Bar1+I(OffSet)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=ABy*I10Bar1+I19Bar1+I(OffSet)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I15Bar1+ABz*I5Bar1+I(OffSet)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I16Bar1+ABz*I6Bar1+I(OffSet)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I17Bar1+ABz*I7Bar1+I(OffSet)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I18Bar1+ABz*I8Bar1+I(OffSet)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I19Bar1+ABz*I9Bar1+I(OffSet)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=ABz*I10Bar1+I20Bar1+I(OffSet)
   END SUBROUTINE Int6311
