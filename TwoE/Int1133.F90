! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (S S|P P) 
! ---------------------------------------------------------- 
   SUBROUTINE Int1133(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, & 
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
      I1Bar2=0.0d0
      I1Bar3=0.0d0
      I1Bar4=0.0d0
      I1Bar5=0.0d0
      I1Bar6=0.0d0
      I1Bar7=0.0d0
      I1Bar8=0.0d0
      I1Bar9=0.0d0
      I1Bar10=0.0d0
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
      V(1)=AuxR0*QCx
      V(2)=AuxR1*WQx
      V(3)=AuxR0*QCy
      V(4)=AuxR1*WQy
      V(5)=AuxR0*QCz
      V(6)=AuxR1*WQz
      V(7)=AuxR1*ZxZpE
      V(8)=-V(7)
      V(9)=AuxR0+V(8)
      V(10)=r1x2E*V(9)
      V(11)=V(3)+V(4)
      V(12)=AuxR1*QCy
      V(13)=AuxR2*WQy
      V(14)=V(12)+V(13)
      V(15)=V(5)+V(6)
      V(16)=AuxR1*QCz
      V(17)=AuxR2*WQz
      V(18)=V(16)+V(17)
      I1Bar1=AuxR0+I1Bar1
      I1Bar2=I1Bar2+V(1)+V(2)
      I1Bar3=I1Bar3+V(3)+V(4)
      I1Bar4=I1Bar4+V(5)+V(6)
      W1=WQx*(AuxR1*QCx+AuxR2*WQx)+I1Bar5
      W2=QCx*(V(1)+V(2))+V(10)
      I1Bar5=W1+W2
      I1Bar6=I1Bar6+QCx*V(11)+WQx*V(14)
      I1Bar7=I1Bar7+V(10)+QCy*V(11)+WQy*V(14)
      I1Bar8=I1Bar8+QCx*V(15)+WQx*V(18)
      I1Bar9=I1Bar9+QCy*V(15)+WQy*V(18)
      I1Bar10=I1Bar10+V(10)+QCz*V(15)+WQz*V(18)
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=CDx*I1Bar2+I1Bar5+I(OffSet)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+1)*LDC+(OD+0)*LDD 
      I(OffSet)=CDx*I1Bar3+I1Bar6+I(OffSet)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+2)*LDC+(OD+0)*LDD 
      I(OffSet)=CDx*I1Bar4+I1Bar8+I(OffSet)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+1)*LDD 
      I(OffSet)=CDy*I1Bar2+I1Bar6+I(OffSet)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+1)*LDC+(OD+1)*LDD 
      I(OffSet)=CDy*I1Bar3+I1Bar7+I(OffSet)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+2)*LDC+(OD+1)*LDD 
      I(OffSet)=CDy*I1Bar4+I1Bar9+I(OffSet)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+2)*LDD 
      I(OffSet)=CDz*I1Bar2+I1Bar8+I(OffSet)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+1)*LDC+(OD+2)*LDD 
      I(OffSet)=CDz*I1Bar3+I1Bar9+I(OffSet)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+2)*LDC+(OD+2)*LDD 
      I(OffSet)=I1Bar10+CDz*I1Bar4+I(OffSet)
   END SUBROUTINE Int1133
