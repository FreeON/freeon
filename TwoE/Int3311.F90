! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (P P|S S) 
! ---------------------------------------------------------- 
   SUBROUTINE Int3311(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, & 
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
      I5Bar1=0.0d0
      I6Bar1=0.0d0
      I7Bar1=0.0d0
      I8Bar1=0.0d0
      I9Bar1=0.0d0
      I10Bar1=0.0d0
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
            INCLUDE 'ERIMIC.Inc'
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
      V(7)=AuxR1*ExZpE
      V(8)=-V(7)
      V(9)=AuxR0+V(8)
      V(10)=r1x2Z*V(9)
      V(11)=V(3)+V(4)
      V(12)=AuxR1*PAy
      V(13)=AuxR2*WPy
      V(14)=V(12)+V(13)
      V(15)=V(5)+V(6)
      V(16)=AuxR1*PAz
      V(17)=AuxR2*WPz
      V(18)=V(16)+V(17)
      I1Bar1=AuxR0+I1Bar1
      I2Bar1=I2Bar1+V(1)+V(2)
      I3Bar1=I3Bar1+V(3)+V(4)
      I4Bar1=I4Bar1+V(5)+V(6)
      W1=WPx*(AuxR1*PAx+AuxR2*WPx)+I5Bar1
      W2=PAx*(V(1)+V(2))+V(10)
      I5Bar1=W1+W2
      I6Bar1=I6Bar1+PAx*V(11)+WPx*V(14)
      I7Bar1=I7Bar1+V(10)+PAy*V(11)+WPy*V(14)
      I8Bar1=I8Bar1+PAx*V(15)+WPx*V(18)
      I9Bar1=I9Bar1+PAy*V(15)+WPy*V(18)
      I10Bar1=I10Bar1+V(10)+PAz*V(15)+WPz*V(18)
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      OffSet=(OA+0)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=ABx*I2Bar1+I5Bar1+I(OffSet)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=ABx*I3Bar1+I6Bar1+I(OffSet)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=ABx*I4Bar1+I8Bar1+I(OffSet)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=ABy*I2Bar1+I6Bar1+I(OffSet)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=ABy*I3Bar1+I7Bar1+I(OffSet)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=ABy*I4Bar1+I9Bar1+I(OffSet)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=ABz*I2Bar1+I8Bar1+I(OffSet)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=ABz*I3Bar1+I9Bar1+I(OffSet)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+(OC+0)*LDC+(OD+0)*LDD 
      I(OffSet)=I10Bar1+ABz*I4Bar1+I(OffSet)
   END SUBROUTINE Int3311
