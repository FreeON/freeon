! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (d p|d p) 
! ---------------------------------------------------------- 
SUBROUTINE dIntB6030603(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, &
 OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NINT,PBC,GRADIENTS)
       USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF7
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER        :: LBra,LKet,NINT,CDOffSet
      REAL(DOUBLE)   :: PrmBufB(10,LBra),PrmBufK(10,LKet)
      TYPE(SmallAtomInfo) :: ACInfo,BDInfo
      TYPE(PBCInfo) :: PBC
      REAL(DOUBLE)  :: GRADIENTS(NINT,12)
      REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz
      REAL(DOUBLE)  :: Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz
      REAL(DOUBLE)  :: PQx,PQy,PQz,FPQx,FPQy,FPQz
      REAL(DOUBLE)  :: Zeta,Eta,Omega,Up,Uq,Upq
      REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT
      REAL(DOUBLE)  :: Alpha,Beta,Gamma
      REAL(DOUBLE), DIMENSION(35) :: HRRTmp 
      REAL(DOUBLE), DIMENSION(20,20,4) :: HRR 
      REAL(DOUBLE), DIMENSION(35,20,4) :: HRRA,HRRB 
      REAL(DOUBLE), DIMENSION(20,35,4) :: HRRC 
      REAL(DOUBLE)  :: VRR(35,35,0:7)
      INTEGER       :: OffSet,OA,LDA,GOA,OB,LDB,GOB,OC,LDC,GOC,OD,LDD,GOD,I,J,K,L
      EXTERNAL InitDbl
      CALL InitDbl(20*20,HRR(1,1,1))
      CALL InitDbl(35*20,HRRA(1,1,1))
      CALL InitDbl(35*20,HRRB(1,1,1))
      CALL InitDbl(20*35,HRRC(1,1,1))
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
         Qx=PrmBufK(2,J)
         Qy=PrmBufK(3,J)
         Qz=PrmBufK(4,J)
         Uq=PrmBufK(5,J)
         Gamma =PrmBufK(9,J)
         QCx=Qx-Cx
         QCy=Qy-Cy
         QCz=Qz-Cz
         DO K=1,LBra ! K^2 VRR (M0| loop 
            Zeta=PrmBufB(1,K)
            Px=PrmBufB(2,K)
            Py=PrmBufB(3,K)
            Pz=PrmBufB(4,K)
            Up=PrmBufB(5,K)
            Alpha =PrmBufB(9,K)
            Beta  =PrmBufB(10,K)
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
            ! Begin Minimum Image Convention 
            FPQx = PQx*PBC%InvBoxSh%D(1,1)+PQy*PBC%InvBoxSh%D(1,2)+PQz*PBC%InvBoxSh%D(1,3)
            FPQy = PQy*PBC%InvBoxSh%D(2,2)+PQz*PBC%InvBoxSh%D(2,3)
            FPQz = PQz*PBC%InvBoxSh%D(3,3)
            IF(PBC%AutoW%I(1)==1)FPQx=FPQx-ANINT(FPQx-SIGN(1D-15,FPQx))
            IF(PBC%AutoW%I(2)==1)FPQy=FPQy-ANINT(FPQy-SIGN(1D-15,FPQy))
            IF(PBC%AutoW%I(3)==1)FPQz=FPQz-ANINT(FPQz-SIGN(1D-15,FPQz))
            PQx=FPQx*PBC%BoxShape%D(1,1)+FPQy*PBC%BoxShape%D(1,2)+FPQz*PBC%BoxShape%D(1,3)
            PQy=FPQy*PBC%BoxShape%D(2,2)+FPQz*PBC%BoxShape%D(2,3)
            PQz=FPQz*PBC%BoxShape%D(3,3)
            ! End MIC
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
              W7=(F7_0(L)+T*(F7_1(L)+T*(F7_2(L)+T*(F7_3(L)+T*F7_4(L)))))
              W6=+7.692307692307693D-02*(TwoT*W7+ET)
              W5=+9.090909090909090D-02*(TwoT*W6+ET)
              W4=+1.111111111111111D-01*(TwoT*W5+ET)
              W3=+1.428571428571428D-01*(TwoT*W4+ET)
              W2=+2.000000000000000D-01*(TwoT*W3+ET)
              W1=+3.333333333333333D-01*(TwoT*W2+ET)
              W0=TwoT*W1+ET
              VRR(1,1,0)=Upq*W0
              VRR(1,1,1)=Upq*W1
              VRR(1,1,2)=Upq*W2
              VRR(1,1,3)=Upq*W3
              VRR(1,1,4)=Upq*W4
              VRR(1,1,5)=Upq*W5
              VRR(1,1,6)=Upq*W6
              VRR(1,1,7)=Upq*W7
            ELSE
              InvT=One/T
              SqInvT=DSQRT(InvT)
              VRR(1,1,0)=+8.862269254527580D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              VRR(1,1,1)=+4.431134627263790D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              VRR(1,1,2)=+6.646701940895685D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              VRR(1,1,3)=+1.661675485223921D+00*Upq*SqInvT
              SqInvT=SqInvT*InvT
              VRR(1,1,4)=+5.815864198283724D+00*Upq*SqInvT
              SqInvT=SqInvT*InvT
              VRR(1,1,5)=+2.617138889227676D+01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              VRR(1,1,6)=+1.439426389075222D+02*Upq*SqInvT
              SqInvT=SqInvT*InvT
              VRR(1,1,7)=+9.356271528988940D+02*Upq*SqInvT
            ENDIF
            ! Generating (p0|s0)^(6)
            VRR(2,1,6)=PAx*VRR(1,1,6)+WPx*VRR(1,1,7) 
            VRR(3,1,6)=PAy*VRR(1,1,6)+WPy*VRR(1,1,7) 
            VRR(4,1,6)=PAz*VRR(1,1,6)+WPz*VRR(1,1,7) 
            ! Generating (p0|s0)^(5)
            VRR(2,1,5)=PAx*VRR(1,1,5)+WPx*VRR(1,1,6) 
            VRR(3,1,5)=PAy*VRR(1,1,5)+WPy*VRR(1,1,6) 
            VRR(4,1,5)=PAz*VRR(1,1,5)+WPz*VRR(1,1,6) 
            ! Generating (p0|s0)^(4)
            VRR(2,1,4)=PAx*VRR(1,1,4)+WPx*VRR(1,1,5) 
            VRR(3,1,4)=PAy*VRR(1,1,4)+WPy*VRR(1,1,5) 
            VRR(4,1,4)=PAz*VRR(1,1,4)+WPz*VRR(1,1,5) 
            ! Generating (p0|s0)^(3)
            VRR(2,1,3)=PAx*VRR(1,1,3)+WPx*VRR(1,1,4) 
            VRR(3,1,3)=PAy*VRR(1,1,3)+WPy*VRR(1,1,4) 
            VRR(4,1,3)=PAz*VRR(1,1,3)+WPz*VRR(1,1,4) 
            ! Generating (p0|s0)^(2)
            VRR(2,1,2)=PAx*VRR(1,1,2)+WPx*VRR(1,1,3) 
            VRR(3,1,2)=PAy*VRR(1,1,2)+WPy*VRR(1,1,3) 
            VRR(4,1,2)=PAz*VRR(1,1,2)+WPz*VRR(1,1,3) 
            ! Generating (p0|s0)^(1)
            VRR(2,1,1)=PAx*VRR(1,1,1)+WPx*VRR(1,1,2) 
            VRR(3,1,1)=PAy*VRR(1,1,1)+WPy*VRR(1,1,2) 
            VRR(4,1,1)=PAz*VRR(1,1,1)+WPz*VRR(1,1,2) 
            ! Generating (p0|s0)^(0)
            VRR(2,1,0)=PAx*VRR(1,1,0)+WPx*VRR(1,1,1) 
            VRR(3,1,0)=PAy*VRR(1,1,0)+WPy*VRR(1,1,1) 
            VRR(4,1,0)=PAz*VRR(1,1,0)+WPz*VRR(1,1,1) 
            ! Generating (d0|s0)^(5)
            VRR(5,1,5)=PAx*VRR(2,1,5)+r1x2Z*(VRR(1,1,5)-ExZpE*VRR(1,1,6))+WPx*VRR(2,1,6)
            VRR(6,1,5)=PAx*VRR(3,1,5)+WPx*VRR(3,1,6)
            VRR(7,1,5)=PAy*VRR(3,1,5)+r1x2Z*(VRR(1,1,5)-ExZpE*VRR(1,1,6))+WPy*VRR(3,1,6)
            VRR(8,1,5)=PAx*VRR(4,1,5)+WPx*VRR(4,1,6)
            VRR(9,1,5)=PAy*VRR(4,1,5)+WPy*VRR(4,1,6)
            VRR(10,1,5)=PAz*VRR(4,1,5)+r1x2Z*(VRR(1,1,5)-ExZpE*VRR(1,1,6))+WPz*VRR(4,1,6)
            ! Generating (d0|s0)^(4)
            VRR(5,1,4)=PAx*VRR(2,1,4)+r1x2Z*(VRR(1,1,4)-ExZpE*VRR(1,1,5))+WPx*VRR(2,1,5)
            VRR(6,1,4)=PAx*VRR(3,1,4)+WPx*VRR(3,1,5)
            VRR(7,1,4)=PAy*VRR(3,1,4)+r1x2Z*(VRR(1,1,4)-ExZpE*VRR(1,1,5))+WPy*VRR(3,1,5)
            VRR(8,1,4)=PAx*VRR(4,1,4)+WPx*VRR(4,1,5)
            VRR(9,1,4)=PAy*VRR(4,1,4)+WPy*VRR(4,1,5)
            VRR(10,1,4)=PAz*VRR(4,1,4)+r1x2Z*(VRR(1,1,4)-ExZpE*VRR(1,1,5))+WPz*VRR(4,1,5)
            ! Generating (d0|s0)^(3)
            VRR(5,1,3)=PAx*VRR(2,1,3)+r1x2Z*(VRR(1,1,3)-ExZpE*VRR(1,1,4))+WPx*VRR(2,1,4)
            VRR(6,1,3)=PAx*VRR(3,1,3)+WPx*VRR(3,1,4)
            VRR(7,1,3)=PAy*VRR(3,1,3)+r1x2Z*(VRR(1,1,3)-ExZpE*VRR(1,1,4))+WPy*VRR(3,1,4)
            VRR(8,1,3)=PAx*VRR(4,1,3)+WPx*VRR(4,1,4)
            VRR(9,1,3)=PAy*VRR(4,1,3)+WPy*VRR(4,1,4)
            VRR(10,1,3)=PAz*VRR(4,1,3)+r1x2Z*(VRR(1,1,3)-ExZpE*VRR(1,1,4))+WPz*VRR(4,1,4)
            ! Generating (d0|s0)^(2)
            VRR(5,1,2)=PAx*VRR(2,1,2)+r1x2Z*(VRR(1,1,2)-ExZpE*VRR(1,1,3))+WPx*VRR(2,1,3)
            VRR(6,1,2)=PAx*VRR(3,1,2)+WPx*VRR(3,1,3)
            VRR(7,1,2)=PAy*VRR(3,1,2)+r1x2Z*(VRR(1,1,2)-ExZpE*VRR(1,1,3))+WPy*VRR(3,1,3)
            VRR(8,1,2)=PAx*VRR(4,1,2)+WPx*VRR(4,1,3)
            VRR(9,1,2)=PAy*VRR(4,1,2)+WPy*VRR(4,1,3)
            VRR(10,1,2)=PAz*VRR(4,1,2)+r1x2Z*(VRR(1,1,2)-ExZpE*VRR(1,1,3))+WPz*VRR(4,1,3)
            ! Generating (d0|s0)^(1)
            VRR(5,1,1)=PAx*VRR(2,1,1)+r1x2Z*(VRR(1,1,1)-ExZpE*VRR(1,1,2))+WPx*VRR(2,1,2)
            VRR(6,1,1)=PAx*VRR(3,1,1)+WPx*VRR(3,1,2)
            VRR(7,1,1)=PAy*VRR(3,1,1)+r1x2Z*(VRR(1,1,1)-ExZpE*VRR(1,1,2))+WPy*VRR(3,1,2)
            VRR(8,1,1)=PAx*VRR(4,1,1)+WPx*VRR(4,1,2)
            VRR(9,1,1)=PAy*VRR(4,1,1)+WPy*VRR(4,1,2)
            VRR(10,1,1)=PAz*VRR(4,1,1)+r1x2Z*(VRR(1,1,1)-ExZpE*VRR(1,1,2))+WPz*VRR(4,1,2)
            ! Generating (d0|s0)^(0)
            VRR(5,1,0)=PAx*VRR(2,1,0)+r1x2Z*(VRR(1,1,0)-ExZpE*VRR(1,1,1))+WPx*VRR(2,1,1)
            VRR(6,1,0)=PAx*VRR(3,1,0)+WPx*VRR(3,1,1)
            VRR(7,1,0)=PAy*VRR(3,1,0)+r1x2Z*(VRR(1,1,0)-ExZpE*VRR(1,1,1))+WPy*VRR(3,1,1)
            VRR(8,1,0)=PAx*VRR(4,1,0)+WPx*VRR(4,1,1)
            VRR(9,1,0)=PAy*VRR(4,1,0)+WPy*VRR(4,1,1)
            VRR(10,1,0)=PAz*VRR(4,1,0)+r1x2Z*(VRR(1,1,0)-ExZpE*VRR(1,1,1))+WPz*VRR(4,1,1)
            ! Generating (f0|s0)^(4)
            CALL VRRf0s0(35,35,VRR(1,1,4),VRR(1,1,5))
            ! Generating (f0|s0)^(3)
            CALL VRRf0s0(35,35,VRR(1,1,3),VRR(1,1,4))
            ! Generating (f0|s0)^(2)
            CALL VRRf0s0(35,35,VRR(1,1,2),VRR(1,1,3))
            ! Generating (f0|s0)^(1)
            CALL VRRf0s0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (f0|s0)^(0)
            CALL VRRf0s0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (g0|s0)^(3)
            CALL VRRg0s0(35,35,VRR(1,1,3),VRR(1,1,4))
            ! Generating (g0|s0)^(2)
            CALL VRRg0s0(35,35,VRR(1,1,2),VRR(1,1,3))
            ! Generating (g0|s0)^(1)
            CALL VRRg0s0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (g0|s0)^(0)
            CALL VRRg0s0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (s0|p0)^(6)
            VRR(1,2,6)=QCx*VRR(1,1,6)+WQx*VRR(1,1,7)
            VRR(1,3,6)=QCy*VRR(1,1,6)+WQy*VRR(1,1,7)
            VRR(1,4,6)=QCz*VRR(1,1,6)+WQz*VRR(1,1,7)
            ! Generating (s0|p0)^(5)
            VRR(1,2,5)=QCx*VRR(1,1,5)+WQx*VRR(1,1,6)
            VRR(1,3,5)=QCy*VRR(1,1,5)+WQy*VRR(1,1,6)
            VRR(1,4,5)=QCz*VRR(1,1,5)+WQz*VRR(1,1,6)
            ! Generating (s0|p0)^(4)
            VRR(1,2,4)=QCx*VRR(1,1,4)+WQx*VRR(1,1,5)
            VRR(1,3,4)=QCy*VRR(1,1,4)+WQy*VRR(1,1,5)
            VRR(1,4,4)=QCz*VRR(1,1,4)+WQz*VRR(1,1,5)
            ! Generating (s0|p0)^(3)
            VRR(1,2,3)=QCx*VRR(1,1,3)+WQx*VRR(1,1,4)
            VRR(1,3,3)=QCy*VRR(1,1,3)+WQy*VRR(1,1,4)
            VRR(1,4,3)=QCz*VRR(1,1,3)+WQz*VRR(1,1,4)
            ! Generating (s0|p0)^(2)
            VRR(1,2,2)=QCx*VRR(1,1,2)+WQx*VRR(1,1,3)
            VRR(1,3,2)=QCy*VRR(1,1,2)+WQy*VRR(1,1,3)
            VRR(1,4,2)=QCz*VRR(1,1,2)+WQz*VRR(1,1,3)
            ! Generating (s0|p0)^(1)
            VRR(1,2,1)=QCx*VRR(1,1,1)+WQx*VRR(1,1,2)
            VRR(1,3,1)=QCy*VRR(1,1,1)+WQy*VRR(1,1,2)
            VRR(1,4,1)=QCz*VRR(1,1,1)+WQz*VRR(1,1,2)
            ! Generating (s0|p0)^(0)
            VRR(1,2,0)=QCx*VRR(1,1,0)+WQx*VRR(1,1,1)
            VRR(1,3,0)=QCy*VRR(1,1,0)+WQy*VRR(1,1,1)
            VRR(1,4,0)=QCz*VRR(1,1,0)+WQz*VRR(1,1,1)
            ! Generating (p0|p0)^(5)
            VRR(2,2,5)=QCx*VRR(2,1,5)+HfxZpE*VRR(1,1,6)+WQx*VRR(2,1,6) 
            VRR(2,3,5)=QCy*VRR(2,1,5)+WQy*VRR(2,1,6) 
            VRR(2,4,5)=QCz*VRR(2,1,5)+WQz*VRR(2,1,6) 
            VRR(3,2,5)=QCx*VRR(3,1,5)+WQx*VRR(3,1,6) 
            VRR(3,3,5)=QCy*VRR(3,1,5)+HfxZpE*VRR(1,1,6)+WQy*VRR(3,1,6) 
            VRR(3,4,5)=QCz*VRR(3,1,5)+WQz*VRR(3,1,6) 
            VRR(4,2,5)=QCx*VRR(4,1,5)+WQx*VRR(4,1,6) 
            VRR(4,3,5)=QCy*VRR(4,1,5)+WQy*VRR(4,1,6) 
            VRR(4,4,5)=QCz*VRR(4,1,5)+HfxZpE*VRR(1,1,6)+WQz*VRR(4,1,6) 
            ! Generating (p0|p0)^(4)
            VRR(2,2,4)=QCx*VRR(2,1,4)+HfxZpE*VRR(1,1,5)+WQx*VRR(2,1,5) 
            VRR(2,3,4)=QCy*VRR(2,1,4)+WQy*VRR(2,1,5) 
            VRR(2,4,4)=QCz*VRR(2,1,4)+WQz*VRR(2,1,5) 
            VRR(3,2,4)=QCx*VRR(3,1,4)+WQx*VRR(3,1,5) 
            VRR(3,3,4)=QCy*VRR(3,1,4)+HfxZpE*VRR(1,1,5)+WQy*VRR(3,1,5) 
            VRR(3,4,4)=QCz*VRR(3,1,4)+WQz*VRR(3,1,5) 
            VRR(4,2,4)=QCx*VRR(4,1,4)+WQx*VRR(4,1,5) 
            VRR(4,3,4)=QCy*VRR(4,1,4)+WQy*VRR(4,1,5) 
            VRR(4,4,4)=QCz*VRR(4,1,4)+HfxZpE*VRR(1,1,5)+WQz*VRR(4,1,5) 
            ! Generating (p0|p0)^(3)
            VRR(2,2,3)=QCx*VRR(2,1,3)+HfxZpE*VRR(1,1,4)+WQx*VRR(2,1,4) 
            VRR(2,3,3)=QCy*VRR(2,1,3)+WQy*VRR(2,1,4) 
            VRR(2,4,3)=QCz*VRR(2,1,3)+WQz*VRR(2,1,4) 
            VRR(3,2,3)=QCx*VRR(3,1,3)+WQx*VRR(3,1,4) 
            VRR(3,3,3)=QCy*VRR(3,1,3)+HfxZpE*VRR(1,1,4)+WQy*VRR(3,1,4) 
            VRR(3,4,3)=QCz*VRR(3,1,3)+WQz*VRR(3,1,4) 
            VRR(4,2,3)=QCx*VRR(4,1,3)+WQx*VRR(4,1,4) 
            VRR(4,3,3)=QCy*VRR(4,1,3)+WQy*VRR(4,1,4) 
            VRR(4,4,3)=QCz*VRR(4,1,3)+HfxZpE*VRR(1,1,4)+WQz*VRR(4,1,4) 
            ! Generating (p0|p0)^(2)
            VRR(2,2,2)=QCx*VRR(2,1,2)+HfxZpE*VRR(1,1,3)+WQx*VRR(2,1,3) 
            VRR(2,3,2)=QCy*VRR(2,1,2)+WQy*VRR(2,1,3) 
            VRR(2,4,2)=QCz*VRR(2,1,2)+WQz*VRR(2,1,3) 
            VRR(3,2,2)=QCx*VRR(3,1,2)+WQx*VRR(3,1,3) 
            VRR(3,3,2)=QCy*VRR(3,1,2)+HfxZpE*VRR(1,1,3)+WQy*VRR(3,1,3) 
            VRR(3,4,2)=QCz*VRR(3,1,2)+WQz*VRR(3,1,3) 
            VRR(4,2,2)=QCx*VRR(4,1,2)+WQx*VRR(4,1,3) 
            VRR(4,3,2)=QCy*VRR(4,1,2)+WQy*VRR(4,1,3) 
            VRR(4,4,2)=QCz*VRR(4,1,2)+HfxZpE*VRR(1,1,3)+WQz*VRR(4,1,3) 
            ! Generating (p0|p0)^(1)
            VRR(2,2,1)=QCx*VRR(2,1,1)+HfxZpE*VRR(1,1,2)+WQx*VRR(2,1,2) 
            VRR(2,3,1)=QCy*VRR(2,1,1)+WQy*VRR(2,1,2) 
            VRR(2,4,1)=QCz*VRR(2,1,1)+WQz*VRR(2,1,2) 
            VRR(3,2,1)=QCx*VRR(3,1,1)+WQx*VRR(3,1,2) 
            VRR(3,3,1)=QCy*VRR(3,1,1)+HfxZpE*VRR(1,1,2)+WQy*VRR(3,1,2) 
            VRR(3,4,1)=QCz*VRR(3,1,1)+WQz*VRR(3,1,2) 
            VRR(4,2,1)=QCx*VRR(4,1,1)+WQx*VRR(4,1,2) 
            VRR(4,3,1)=QCy*VRR(4,1,1)+WQy*VRR(4,1,2) 
            VRR(4,4,1)=QCz*VRR(4,1,1)+HfxZpE*VRR(1,1,2)+WQz*VRR(4,1,2) 
            ! Generating (p0|p0)^(0)
            VRR(2,2,0)=QCx*VRR(2,1,0)+HfxZpE*VRR(1,1,1)+WQx*VRR(2,1,1) 
            VRR(2,3,0)=QCy*VRR(2,1,0)+WQy*VRR(2,1,1) 
            VRR(2,4,0)=QCz*VRR(2,1,0)+WQz*VRR(2,1,1) 
            VRR(3,2,0)=QCx*VRR(3,1,0)+WQx*VRR(3,1,1) 
            VRR(3,3,0)=QCy*VRR(3,1,0)+HfxZpE*VRR(1,1,1)+WQy*VRR(3,1,1) 
            VRR(3,4,0)=QCz*VRR(3,1,0)+WQz*VRR(3,1,1) 
            VRR(4,2,0)=QCx*VRR(4,1,0)+WQx*VRR(4,1,1) 
            VRR(4,3,0)=QCy*VRR(4,1,0)+WQy*VRR(4,1,1) 
            VRR(4,4,0)=QCz*VRR(4,1,0)+HfxZpE*VRR(1,1,1)+WQz*VRR(4,1,1) 
            ! Generating (d0|p0)^(4)
            CALL VRRd0p0(35,35,VRR(1,1,4),VRR(1,1,5))
            ! Generating (d0|p0)^(3)
            CALL VRRd0p0(35,35,VRR(1,1,3),VRR(1,1,4))
            ! Generating (d0|p0)^(2)
            CALL VRRd0p0(35,35,VRR(1,1,2),VRR(1,1,3))
            ! Generating (d0|p0)^(1)
            CALL VRRd0p0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (d0|p0)^(0)
            CALL VRRd0p0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (f0|p0)^(3)
            CALL VRRf0p0(35,35,VRR(1,1,3),VRR(1,1,4))
            ! Generating (f0|p0)^(2)
            CALL VRRf0p0(35,35,VRR(1,1,2),VRR(1,1,3))
            ! Generating (f0|p0)^(1)
            CALL VRRf0p0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (f0|p0)^(0)
            CALL VRRf0p0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (g0|p0)^(2)
            CALL VRRg0p0(35,35,VRR(1,1,2),VRR(1,1,3))
            ! Generating (g0|p0)^(1)
            CALL VRRg0p0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (g0|p0)^(0)
            CALL VRRg0p0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (s0|d0)^(5)
            VRR(1,5,5)=r1x2E*VRR(1,1,5)+QCx*VRR(1,2,5)-r1x2E*ZxZpE*VRR(1,1,6)+WQx*VRR(1,2,6)
            VRR(1,6,5)=QCx*VRR(1,3,5)+WQx*VRR(1,3,6)
            VRR(1,7,5)=r1x2E*VRR(1,1,5)+QCy*VRR(1,3,5)-r1x2E*ZxZpE*VRR(1,1,6)+WQy*VRR(1,3,6)
            VRR(1,8,5)=QCx*VRR(1,4,5)+WQx*VRR(1,4,6)
            VRR(1,9,5)=QCy*VRR(1,4,5)+WQy*VRR(1,4,6)
            VRR(1,10,5)=r1x2E*VRR(1,1,5)+QCz*VRR(1,4,5)-r1x2E*ZxZpE*VRR(1,1,6)+WQz*VRR(1,4,6)
            ! Generating (s0|d0)^(4)
            VRR(1,5,4)=r1x2E*VRR(1,1,4)+QCx*VRR(1,2,4)-r1x2E*ZxZpE*VRR(1,1,5)+WQx*VRR(1,2,5)
            VRR(1,6,4)=QCx*VRR(1,3,4)+WQx*VRR(1,3,5)
            VRR(1,7,4)=r1x2E*VRR(1,1,4)+QCy*VRR(1,3,4)-r1x2E*ZxZpE*VRR(1,1,5)+WQy*VRR(1,3,5)
            VRR(1,8,4)=QCx*VRR(1,4,4)+WQx*VRR(1,4,5)
            VRR(1,9,4)=QCy*VRR(1,4,4)+WQy*VRR(1,4,5)
            VRR(1,10,4)=r1x2E*VRR(1,1,4)+QCz*VRR(1,4,4)-r1x2E*ZxZpE*VRR(1,1,5)+WQz*VRR(1,4,5)
            ! Generating (s0|d0)^(3)
            VRR(1,5,3)=r1x2E*VRR(1,1,3)+QCx*VRR(1,2,3)-r1x2E*ZxZpE*VRR(1,1,4)+WQx*VRR(1,2,4)
            VRR(1,6,3)=QCx*VRR(1,3,3)+WQx*VRR(1,3,4)
            VRR(1,7,3)=r1x2E*VRR(1,1,3)+QCy*VRR(1,3,3)-r1x2E*ZxZpE*VRR(1,1,4)+WQy*VRR(1,3,4)
            VRR(1,8,3)=QCx*VRR(1,4,3)+WQx*VRR(1,4,4)
            VRR(1,9,3)=QCy*VRR(1,4,3)+WQy*VRR(1,4,4)
            VRR(1,10,3)=r1x2E*VRR(1,1,3)+QCz*VRR(1,4,3)-r1x2E*ZxZpE*VRR(1,1,4)+WQz*VRR(1,4,4)
            ! Generating (s0|d0)^(2)
            VRR(1,5,2)=r1x2E*VRR(1,1,2)+QCx*VRR(1,2,2)-r1x2E*ZxZpE*VRR(1,1,3)+WQx*VRR(1,2,3)
            VRR(1,6,2)=QCx*VRR(1,3,2)+WQx*VRR(1,3,3)
            VRR(1,7,2)=r1x2E*VRR(1,1,2)+QCy*VRR(1,3,2)-r1x2E*ZxZpE*VRR(1,1,3)+WQy*VRR(1,3,3)
            VRR(1,8,2)=QCx*VRR(1,4,2)+WQx*VRR(1,4,3)
            VRR(1,9,2)=QCy*VRR(1,4,2)+WQy*VRR(1,4,3)
            VRR(1,10,2)=r1x2E*VRR(1,1,2)+QCz*VRR(1,4,2)-r1x2E*ZxZpE*VRR(1,1,3)+WQz*VRR(1,4,3)
            ! Generating (s0|d0)^(1)
            VRR(1,5,1)=r1x2E*VRR(1,1,1)+QCx*VRR(1,2,1)-r1x2E*ZxZpE*VRR(1,1,2)+WQx*VRR(1,2,2)
            VRR(1,6,1)=QCx*VRR(1,3,1)+WQx*VRR(1,3,2)
            VRR(1,7,1)=r1x2E*VRR(1,1,1)+QCy*VRR(1,3,1)-r1x2E*ZxZpE*VRR(1,1,2)+WQy*VRR(1,3,2)
            VRR(1,8,1)=QCx*VRR(1,4,1)+WQx*VRR(1,4,2)
            VRR(1,9,1)=QCy*VRR(1,4,1)+WQy*VRR(1,4,2)
            VRR(1,10,1)=r1x2E*VRR(1,1,1)+QCz*VRR(1,4,1)-r1x2E*ZxZpE*VRR(1,1,2)+WQz*VRR(1,4,2)
            ! Generating (s0|d0)^(0)
            VRR(1,5,0)=r1x2E*VRR(1,1,0)+QCx*VRR(1,2,0)-r1x2E*ZxZpE*VRR(1,1,1)+WQx*VRR(1,2,1)
            VRR(1,6,0)=QCx*VRR(1,3,0)+WQx*VRR(1,3,1)
            VRR(1,7,0)=r1x2E*VRR(1,1,0)+QCy*VRR(1,3,0)-r1x2E*ZxZpE*VRR(1,1,1)+WQy*VRR(1,3,1)
            VRR(1,8,0)=QCx*VRR(1,4,0)+WQx*VRR(1,4,1)
            VRR(1,9,0)=QCy*VRR(1,4,0)+WQy*VRR(1,4,1)
            VRR(1,10,0)=r1x2E*VRR(1,1,0)+QCz*VRR(1,4,0)-r1x2E*ZxZpE*VRR(1,1,1)+WQz*VRR(1,4,1)
            ! Generating (p0|d0)^(4)
            CALL VRRp0d0(35,35,VRR(1,1,4),VRR(1,1,5))
            ! Generating (p0|d0)^(3)
            CALL VRRp0d0(35,35,VRR(1,1,3),VRR(1,1,4))
            ! Generating (p0|d0)^(2)
            CALL VRRp0d0(35,35,VRR(1,1,2),VRR(1,1,3))
            ! Generating (p0|d0)^(1)
            CALL VRRp0d0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (p0|d0)^(0)
            CALL VRRp0d0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (d0|d0)^(3)
            CALL VRRd0d0(35,35,VRR(1,1,3),VRR(1,1,4))
            ! Generating (d0|d0)^(2)
            CALL VRRd0d0(35,35,VRR(1,1,2),VRR(1,1,3))
            ! Generating (d0|d0)^(1)
            CALL VRRd0d0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (d0|d0)^(0)
            CALL VRRd0d0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (f0|d0)^(2)
            CALL VRRf0d0(35,35,VRR(1,1,2),VRR(1,1,3))
            ! Generating (f0|d0)^(1)
            CALL VRRf0d0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (f0|d0)^(0)
            CALL VRRf0d0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (g0|d0)^(1)
            CALL VRRg0d0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (g0|d0)^(0)
            CALL VRRg0d0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (s0|f0)^(4)
            CALL VRRs0f0(35,35,VRR(1,1,4),VRR(1,1,5))
            ! Generating (s0|f0)^(3)
            CALL VRRs0f0(35,35,VRR(1,1,3),VRR(1,1,4))
            ! Generating (s0|f0)^(2)
            CALL VRRs0f0(35,35,VRR(1,1,2),VRR(1,1,3))
            ! Generating (s0|f0)^(1)
            CALL VRRs0f0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (s0|f0)^(0)
            CALL VRRs0f0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (p0|f0)^(3)
            CALL VRRp0f0(35,35,VRR(1,1,3),VRR(1,1,4))
            ! Generating (p0|f0)^(2)
            CALL VRRp0f0(35,35,VRR(1,1,2),VRR(1,1,3))
            ! Generating (p0|f0)^(1)
            CALL VRRp0f0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (p0|f0)^(0)
            CALL VRRp0f0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (d0|f0)^(2)
            CALL VRRd0f0(35,35,VRR(1,1,2),VRR(1,1,3))
            ! Generating (d0|f0)^(1)
            CALL VRRd0f0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (d0|f0)^(0)
            CALL VRRd0f0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (f0|f0)^(1)
            CALL VRRf0f0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (f0|f0)^(0)
            CALL VRRf0f0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (g0|f0)^(0)
            CALL VRRg0f0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (s0|g0)^(3)
            CALL VRRs0g0(35,35,VRR(1,1,3),VRR(1,1,4))
            ! Generating (s0|g0)^(2)
            CALL VRRs0g0(35,35,VRR(1,1,2),VRR(1,1,3))
            ! Generating (s0|g0)^(1)
            CALL VRRs0g0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (s0|g0)^(0)
            CALL VRRs0g0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (p0|g0)^(2)
            CALL VRRp0g0(35,35,VRR(1,1,2),VRR(1,1,3))
            ! Generating (p0|g0)^(1)
            CALL VRRp0g0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (p0|g0)^(0)
            CALL VRRp0g0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (d0|g0)^(1)
            CALL VRRd0g0(35,35,VRR(1,1,1),VRR(1,1,2))
            ! Generating (d0|g0)^(0)
            CALL VRRd0g0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Generating (f0|g0)^(0)
            CALL VRRf0g0(35,35,VRR(1,1,0),VRR(1,1,1))
            ! Contracting ... 
            CALL CNTRCTG6363(VRR,HRR,Alpha,HRRA,Beta,HRRB,Gamma,HRRC)
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! Generating (f,0|d,p)
      CALL KetHRR63(20,HRR) 
      ! Generating (g,0|d,p)^a
      CALL KetHRR63(35,HRRA) 
      ! Generating (g,0|d,p)^b
      CALL KetHRR63(35,HRRB) 
      ! Generating (f,0|f,p)^c
      CALL KetHRR103(20,HRRC) 
      DO L=2,4
      
         !K = 5
         CDOffSet=(OC+5-5)*LDC+(OD+L-2)*LDD
         ! Generating (d',p|5,L)  and (d,p'|5,L)
         CALL BraHRR63ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,5,L),&
                          HRRA(1,5,L),HRRB(1,5,L),GRADIENTS(1,1))
         ! Generating (d,p|5_x,L)  and (d,p|5,L_x)
         HRRTmp(1:20)=HRRC(1:20,11,L)-2D0*HRR(1:20,2,L)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,HRRTmp,GRADIENTS(1,1))
         ! Generating (d,p|5_y,L)  and (d,p|5,L_y)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRC(1,12,L),GRADIENTS(1,1))
         ! Generating (d,p|5_z,L)  and (d,p|5,L_z)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRC(1,15,L),GRADIENTS(1,1))
      
         !K = 6
         CDOffSet=(OC+6-5)*LDC+(OD+L-2)*LDD
         ! Generating (d',p|6,L)  and (d,p'|6,L)
         CALL BraHRR63ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,6,L),&
                          HRRA(1,6,L),HRRB(1,6,L),GRADIENTS(1,1))
         ! Generating (d,p|6_x,L)  and (d,p|6,L_x)
         HRRTmp(1:20)=HRRC(1:20,12,L)-1D0*HRR(1:20,3,L)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,HRRTmp,GRADIENTS(1,1))
         ! Generating (d,p|6_y,L)  and (d,p|6,L_y)
         HRRTmp(1:20)=HRRC(1:20,13,L)-1D0*HRR(1:20,2,L)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRTmp,GRADIENTS(1,1))
         ! Generating (d,p|6_z,L)  and (d,p|6,L_z)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRC(1,16,L),GRADIENTS(1,1))
      
         !K = 7
         CDOffSet=(OC+7-5)*LDC+(OD+L-2)*LDD
         ! Generating (d',p|7,L)  and (d,p'|7,L)
         CALL BraHRR63ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,7,L),&
                          HRRA(1,7,L),HRRB(1,7,L),GRADIENTS(1,1))
         ! Generating (d,p|7_x,L)  and (d,p|7,L_x)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,HRRC(1,13,L),GRADIENTS(1,1))
         ! Generating (d,p|7_y,L)  and (d,p|7,L_y)
         HRRTmp(1:20)=HRRC(1:20,14,L)-2D0*HRR(1:20,3,L)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRTmp,GRADIENTS(1,1))
         ! Generating (d,p|7_z,L)  and (d,p|7,L_z)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRC(1,17,L),GRADIENTS(1,1))
      
         !K = 8
         CDOffSet=(OC+8-5)*LDC+(OD+L-2)*LDD
         ! Generating (d',p|8,L)  and (d,p'|8,L)
         CALL BraHRR63ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,8,L),&
                          HRRA(1,8,L),HRRB(1,8,L),GRADIENTS(1,1))
         ! Generating (d,p|8_x,L)  and (d,p|8,L_x)
         HRRTmp(1:20)=HRRC(1:20,15,L)-1D0*HRR(1:20,4,L)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,HRRTmp,GRADIENTS(1,1))
         ! Generating (d,p|8_y,L)  and (d,p|8,L_y)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRC(1,16,L),GRADIENTS(1,1))
         ! Generating (d,p|8_z,L)  and (d,p|8,L_z)
         HRRTmp(1:20)=HRRC(1:20,18,L)-1D0*HRR(1:20,2,L)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRTmp,GRADIENTS(1,1))
      
         !K = 9
         CDOffSet=(OC+9-5)*LDC+(OD+L-2)*LDD
         ! Generating (d',p|9,L)  and (d,p'|9,L)
         CALL BraHRR63ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,9,L),&
                          HRRA(1,9,L),HRRB(1,9,L),GRADIENTS(1,1))
         ! Generating (d,p|9_x,L)  and (d,p|9,L_x)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,HRRC(1,16,L),GRADIENTS(1,1))
         ! Generating (d,p|9_y,L)  and (d,p|9,L_y)
         HRRTmp(1:20)=HRRC(1:20,17,L)-1D0*HRR(1:20,4,L)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRTmp,GRADIENTS(1,1))
         ! Generating (d,p|9_z,L)  and (d,p|9,L_z)
         HRRTmp(1:20)=HRRC(1:20,19,L)-1D0*HRR(1:20,3,L)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRTmp,GRADIENTS(1,1))
      
         !K = 10
         CDOffSet=(OC+10-5)*LDC+(OD+L-2)*LDD
         ! Generating (d',p|10,L)  and (d,p'|10,L)
         CALL BraHRR63ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,10,L),&
                          HRRA(1,10,L),HRRB(1,10,L),GRADIENTS(1,1))
         ! Generating (d,p|10_x,L)  and (d,p|10,L_x)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,HRRC(1,18,L),GRADIENTS(1,1))
         ! Generating (d,p|10_y,L)  and (d,p|10,L_y)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRC(1,19,L),GRADIENTS(1,1))
         ! Generating (d,p|10_z,L)  and (d,p|10,L_z)
         HRRTmp(1:20)=HRRC(1:20,20,L)-2D0*HRR(1:20,4,L)
         CALL BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRTmp,GRADIENTS(1,1))
      ENDDO 
    END SUBROUTINE dIntB6030603
    SUBROUTINE CNTRCTG6363(VRR,HRR,Alpha,HRRA,Beta,HRRB,Gamma,HRRC)
      USE DerivedTypes
      USE VScratchB
      INTEGER :: K
      REAL(DOUBLE)  :: Alpha,Beta,Gamma
      REAL(DOUBLE), DIMENSION(20,20,4) :: HRR 
      REAL(DOUBLE), DIMENSION(35,20,4) :: HRRA,HRRB 
      REAL(DOUBLE), DIMENSION(20,35,4) :: HRRC 
      REAL(DOUBLE)  :: VRR(35,35,0:7)
      DO K=1,20
         HRR(1,K,1)=HRR(1,K,1)+VRR(1,K,0)
         HRRC(1,K,1)=HRRC(1,K,1)+Gamma*VRR(1,K,0)
         HRRA(1,K,1)=HRRA(1,K,1)+Alpha*VRR(1,K,0)
         HRRB(1,K,1)=HRRB(1,K,1)+Beta*VRR(1,K,0)
         HRR(2,K,1)=HRR(2,K,1)+VRR(2,K,0)
         HRRC(2,K,1)=HRRC(2,K,1)+Gamma*VRR(2,K,0)
         HRRA(2,K,1)=HRRA(2,K,1)+Alpha*VRR(2,K,0)
         HRRB(2,K,1)=HRRB(2,K,1)+Beta*VRR(2,K,0)
         HRR(3,K,1)=HRR(3,K,1)+VRR(3,K,0)
         HRRC(3,K,1)=HRRC(3,K,1)+Gamma*VRR(3,K,0)
         HRRA(3,K,1)=HRRA(3,K,1)+Alpha*VRR(3,K,0)
         HRRB(3,K,1)=HRRB(3,K,1)+Beta*VRR(3,K,0)
         HRR(4,K,1)=HRR(4,K,1)+VRR(4,K,0)
         HRRC(4,K,1)=HRRC(4,K,1)+Gamma*VRR(4,K,0)
         HRRA(4,K,1)=HRRA(4,K,1)+Alpha*VRR(4,K,0)
         HRRB(4,K,1)=HRRB(4,K,1)+Beta*VRR(4,K,0)
         HRR(5,K,1)=HRR(5,K,1)+VRR(5,K,0)
         HRRC(5,K,1)=HRRC(5,K,1)+Gamma*VRR(5,K,0)
         HRRA(5,K,1)=HRRA(5,K,1)+Alpha*VRR(5,K,0)
         HRRB(5,K,1)=HRRB(5,K,1)+Beta*VRR(5,K,0)
         HRR(6,K,1)=HRR(6,K,1)+VRR(6,K,0)
         HRRC(6,K,1)=HRRC(6,K,1)+Gamma*VRR(6,K,0)
         HRRA(6,K,1)=HRRA(6,K,1)+Alpha*VRR(6,K,0)
         HRRB(6,K,1)=HRRB(6,K,1)+Beta*VRR(6,K,0)
         HRR(7,K,1)=HRR(7,K,1)+VRR(7,K,0)
         HRRC(7,K,1)=HRRC(7,K,1)+Gamma*VRR(7,K,0)
         HRRA(7,K,1)=HRRA(7,K,1)+Alpha*VRR(7,K,0)
         HRRB(7,K,1)=HRRB(7,K,1)+Beta*VRR(7,K,0)
         HRR(8,K,1)=HRR(8,K,1)+VRR(8,K,0)
         HRRC(8,K,1)=HRRC(8,K,1)+Gamma*VRR(8,K,0)
         HRRA(8,K,1)=HRRA(8,K,1)+Alpha*VRR(8,K,0)
         HRRB(8,K,1)=HRRB(8,K,1)+Beta*VRR(8,K,0)
         HRR(9,K,1)=HRR(9,K,1)+VRR(9,K,0)
         HRRC(9,K,1)=HRRC(9,K,1)+Gamma*VRR(9,K,0)
         HRRA(9,K,1)=HRRA(9,K,1)+Alpha*VRR(9,K,0)
         HRRB(9,K,1)=HRRB(9,K,1)+Beta*VRR(9,K,0)
         HRR(10,K,1)=HRR(10,K,1)+VRR(10,K,0)
         HRRC(10,K,1)=HRRC(10,K,1)+Gamma*VRR(10,K,0)
         HRRA(10,K,1)=HRRA(10,K,1)+Alpha*VRR(10,K,0)
         HRRB(10,K,1)=HRRB(10,K,1)+Beta*VRR(10,K,0)
         HRR(11,K,1)=HRR(11,K,1)+VRR(11,K,0)
         HRRC(11,K,1)=HRRC(11,K,1)+Gamma*VRR(11,K,0)
         HRRA(11,K,1)=HRRA(11,K,1)+Alpha*VRR(11,K,0)
         HRRB(11,K,1)=HRRB(11,K,1)+Beta*VRR(11,K,0)
         HRR(12,K,1)=HRR(12,K,1)+VRR(12,K,0)
         HRRC(12,K,1)=HRRC(12,K,1)+Gamma*VRR(12,K,0)
         HRRA(12,K,1)=HRRA(12,K,1)+Alpha*VRR(12,K,0)
         HRRB(12,K,1)=HRRB(12,K,1)+Beta*VRR(12,K,0)
         HRR(13,K,1)=HRR(13,K,1)+VRR(13,K,0)
         HRRC(13,K,1)=HRRC(13,K,1)+Gamma*VRR(13,K,0)
         HRRA(13,K,1)=HRRA(13,K,1)+Alpha*VRR(13,K,0)
         HRRB(13,K,1)=HRRB(13,K,1)+Beta*VRR(13,K,0)
         HRR(14,K,1)=HRR(14,K,1)+VRR(14,K,0)
         HRRC(14,K,1)=HRRC(14,K,1)+Gamma*VRR(14,K,0)
         HRRA(14,K,1)=HRRA(14,K,1)+Alpha*VRR(14,K,0)
         HRRB(14,K,1)=HRRB(14,K,1)+Beta*VRR(14,K,0)
         HRR(15,K,1)=HRR(15,K,1)+VRR(15,K,0)
         HRRC(15,K,1)=HRRC(15,K,1)+Gamma*VRR(15,K,0)
         HRRA(15,K,1)=HRRA(15,K,1)+Alpha*VRR(15,K,0)
         HRRB(15,K,1)=HRRB(15,K,1)+Beta*VRR(15,K,0)
         HRR(16,K,1)=HRR(16,K,1)+VRR(16,K,0)
         HRRC(16,K,1)=HRRC(16,K,1)+Gamma*VRR(16,K,0)
         HRRA(16,K,1)=HRRA(16,K,1)+Alpha*VRR(16,K,0)
         HRRB(16,K,1)=HRRB(16,K,1)+Beta*VRR(16,K,0)
         HRR(17,K,1)=HRR(17,K,1)+VRR(17,K,0)
         HRRC(17,K,1)=HRRC(17,K,1)+Gamma*VRR(17,K,0)
         HRRA(17,K,1)=HRRA(17,K,1)+Alpha*VRR(17,K,0)
         HRRB(17,K,1)=HRRB(17,K,1)+Beta*VRR(17,K,0)
         HRR(18,K,1)=HRR(18,K,1)+VRR(18,K,0)
         HRRC(18,K,1)=HRRC(18,K,1)+Gamma*VRR(18,K,0)
         HRRA(18,K,1)=HRRA(18,K,1)+Alpha*VRR(18,K,0)
         HRRB(18,K,1)=HRRB(18,K,1)+Beta*VRR(18,K,0)
         HRR(19,K,1)=HRR(19,K,1)+VRR(19,K,0)
         HRRC(19,K,1)=HRRC(19,K,1)+Gamma*VRR(19,K,0)
         HRRA(19,K,1)=HRRA(19,K,1)+Alpha*VRR(19,K,0)
         HRRB(19,K,1)=HRRB(19,K,1)+Beta*VRR(19,K,0)
         HRR(20,K,1)=HRR(20,K,1)+VRR(20,K,0)
         HRRC(20,K,1)=HRRC(20,K,1)+Gamma*VRR(20,K,0)
         HRRA(20,K,1)=HRRA(20,K,1)+Alpha*VRR(20,K,0)
         HRRB(20,K,1)=HRRB(20,K,1)+Beta*VRR(20,K,0)
         HRRA(21,K,1)=HRRA(21,K,1)+Alpha*VRR(21,K,0)
         HRRB(21,K,1)=HRRB(21,K,1)+Beta*VRR(21,K,0)
         HRRA(22,K,1)=HRRA(22,K,1)+Alpha*VRR(22,K,0)
         HRRB(22,K,1)=HRRB(22,K,1)+Beta*VRR(22,K,0)
         HRRA(23,K,1)=HRRA(23,K,1)+Alpha*VRR(23,K,0)
         HRRB(23,K,1)=HRRB(23,K,1)+Beta*VRR(23,K,0)
         HRRA(24,K,1)=HRRA(24,K,1)+Alpha*VRR(24,K,0)
         HRRB(24,K,1)=HRRB(24,K,1)+Beta*VRR(24,K,0)
         HRRA(25,K,1)=HRRA(25,K,1)+Alpha*VRR(25,K,0)
         HRRB(25,K,1)=HRRB(25,K,1)+Beta*VRR(25,K,0)
         HRRA(26,K,1)=HRRA(26,K,1)+Alpha*VRR(26,K,0)
         HRRB(26,K,1)=HRRB(26,K,1)+Beta*VRR(26,K,0)
         HRRA(27,K,1)=HRRA(27,K,1)+Alpha*VRR(27,K,0)
         HRRB(27,K,1)=HRRB(27,K,1)+Beta*VRR(27,K,0)
         HRRA(28,K,1)=HRRA(28,K,1)+Alpha*VRR(28,K,0)
         HRRB(28,K,1)=HRRB(28,K,1)+Beta*VRR(28,K,0)
         HRRA(29,K,1)=HRRA(29,K,1)+Alpha*VRR(29,K,0)
         HRRB(29,K,1)=HRRB(29,K,1)+Beta*VRR(29,K,0)
         HRRA(30,K,1)=HRRA(30,K,1)+Alpha*VRR(30,K,0)
         HRRB(30,K,1)=HRRB(30,K,1)+Beta*VRR(30,K,0)
         HRRA(31,K,1)=HRRA(31,K,1)+Alpha*VRR(31,K,0)
         HRRB(31,K,1)=HRRB(31,K,1)+Beta*VRR(31,K,0)
         HRRA(32,K,1)=HRRA(32,K,1)+Alpha*VRR(32,K,0)
         HRRB(32,K,1)=HRRB(32,K,1)+Beta*VRR(32,K,0)
         HRRA(33,K,1)=HRRA(33,K,1)+Alpha*VRR(33,K,0)
         HRRB(33,K,1)=HRRB(33,K,1)+Beta*VRR(33,K,0)
         HRRA(34,K,1)=HRRA(34,K,1)+Alpha*VRR(34,K,0)
         HRRB(34,K,1)=HRRB(34,K,1)+Beta*VRR(34,K,0)
         HRRA(35,K,1)=HRRA(35,K,1)+Alpha*VRR(35,K,0)
         HRRB(35,K,1)=HRRB(35,K,1)+Beta*VRR(35,K,0)
      ENDDO
      DO K=21,35
         HRRC(1,K,1)=HRRC(1,K,1)+Gamma*VRR(1,K,0)
         HRRC(2,K,1)=HRRC(2,K,1)+Gamma*VRR(2,K,0)
         HRRC(3,K,1)=HRRC(3,K,1)+Gamma*VRR(3,K,0)
         HRRC(4,K,1)=HRRC(4,K,1)+Gamma*VRR(4,K,0)
         HRRC(5,K,1)=HRRC(5,K,1)+Gamma*VRR(5,K,0)
         HRRC(6,K,1)=HRRC(6,K,1)+Gamma*VRR(6,K,0)
         HRRC(7,K,1)=HRRC(7,K,1)+Gamma*VRR(7,K,0)
         HRRC(8,K,1)=HRRC(8,K,1)+Gamma*VRR(8,K,0)
         HRRC(9,K,1)=HRRC(9,K,1)+Gamma*VRR(9,K,0)
         HRRC(10,K,1)=HRRC(10,K,1)+Gamma*VRR(10,K,0)
         HRRC(11,K,1)=HRRC(11,K,1)+Gamma*VRR(11,K,0)
         HRRC(12,K,1)=HRRC(12,K,1)+Gamma*VRR(12,K,0)
         HRRC(13,K,1)=HRRC(13,K,1)+Gamma*VRR(13,K,0)
         HRRC(14,K,1)=HRRC(14,K,1)+Gamma*VRR(14,K,0)
         HRRC(15,K,1)=HRRC(15,K,1)+Gamma*VRR(15,K,0)
         HRRC(16,K,1)=HRRC(16,K,1)+Gamma*VRR(16,K,0)
         HRRC(17,K,1)=HRRC(17,K,1)+Gamma*VRR(17,K,0)
         HRRC(18,K,1)=HRRC(18,K,1)+Gamma*VRR(18,K,0)
         HRRC(19,K,1)=HRRC(19,K,1)+Gamma*VRR(19,K,0)
         HRRC(20,K,1)=HRRC(20,K,1)+Gamma*VRR(20,K,0)
      ENDDO
    END SUBROUTINE CNTRCTG6363
