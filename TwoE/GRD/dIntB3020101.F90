! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (p sp|s s) 
! ---------------------------------------------------------- 
SUBROUTINE dIntB3020101(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, &
 OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NINT,PBC,GRADIENTS)
       USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF3
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
      REAL(DOUBLE), DIMENSION(20) :: HRRTmp 
      REAL(DOUBLE), DIMENSION(13,1,1) :: HRR 
      REAL(DOUBLE), DIMENSION(20,1,1) :: HRRA,HRRB 
      REAL(DOUBLE), DIMENSION(13,4,1) :: HRRC 
      REAL(DOUBLE)  :: VRR(20,4,0:3)
      INTEGER       :: OffSet,OA,LDA,GOA,OB,LDB,GOB,OC,LDC,GOC,OD,LDD,GOD,I,J,K,L
      EXTERNAL InitDbl
      CALL InitDbl(13*1,HRR(1,1,1))
      CALL InitDbl(20*1,HRRA(1,1,1))
      CALL InitDbl(20*1,HRRB(1,1,1))
      CALL InitDbl(13*4,HRRC(1,1,1))
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
            FnSpB=PrmBufB(6,K)
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
              W3=(F3_0(L)+T*(F3_1(L)+T*(F3_2(L)+T*(F3_3(L)+T*F3_4(L)))))
              W2=+2.000000000000000D-01*(TwoT*W3+ET)
              W1=+3.333333333333333D-01*(TwoT*W2+ET)
              W0=TwoT*W1+ET
              VRR(1,1,0)=Upq*W0
              VRR(1,1,1)=Upq*W1
              VRR(1,1,2)=Upq*W2
              VRR(1,1,3)=Upq*W3
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
            ENDIF
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
            ! Generating (f0|s0)^(0)
            CALL VRRf0s0(20,4,VRR(1,1,0),VRR(1,1,1))
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
            ! Generating (d0|p0)^(0)
            CALL VRRd0p0(20,4,VRR(1,1,0),VRR(1,1,1))
            ! Contracting ... 
            CALL CNTRCTG3211(VRR,HRR,Alpha,HRRA,Beta,HRRB,Gamma,HRRC)
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! Dont need to generate (spd,0|s,s)
      ! Dont need to generate (spdf,0|s,s)^a
      ! Dont need to generate (spdf,0|s,s)^b
      ! Dont need to generate (spd,0|p,s)^c
      DO L=1,1
      
         !K = 1
         CDOffSet=(OC+1-1)*LDC+(OD+L-1)*LDD
         ! Generating (p',sp|1,L)  and (p,sp'|1,L)
         CALL BraHRR32ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,1,L),&
                          HRRA(1,1,L),HRRB(1,1,L),GRADIENTS(1,1))
         ! Generating (p,sp|1_x,L)  and (p,sp|1,L_x)
         CALL BraHRR32cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,HRRC(1,2,L),GRADIENTS(1,1))
         ! Generating (p,sp|1_y,L)  and (p,sp|1,L_y)
         CALL BraHRR32cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRC(1,3,L),GRADIENTS(1,1))
         ! Generating (p,sp|1_z,L)  and (p,sp|1,L_z)
         CALL BraHRR32cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRC(1,4,L),GRADIENTS(1,1))
      ENDDO 
    END SUBROUTINE dIntB3020101
    SUBROUTINE CNTRCTG3211(VRR,HRR,Alpha,HRRA,Beta,HRRB,Gamma,HRRC)
      USE DerivedTypes
      USE VScratchB
      REAL(DOUBLE)  :: Alpha,Beta,Gamma
      REAL(DOUBLE), DIMENSION(13,1,1) :: HRR 
      REAL(DOUBLE), DIMENSION(20,1,1) :: HRRA,HRRB 
      REAL(DOUBLE), DIMENSION(13,4,1) :: HRRC 
      REAL(DOUBLE)  :: VRR(20,4,0:3)
      HRR(1,1,1)=HRR(1,1,1)+VRR(1,1,0)
      HRRA(1,1,1)=HRRA(1,1,1)+Alpha*VRR(1,1,0)
      HRRB(1,1,1)=HRRB(1,1,1)+Beta*VRR(1,1,0)
      HRRC(1,1,1)=HRRC(1,1,1)+Gamma*VRR(1,1,0)
      HRRC(1,2,1)=HRRC(1,2,1)+Gamma*VRR(1,2,0)
      HRRC(1,3,1)=HRRC(1,3,1)+Gamma*VRR(1,3,0)
      HRRC(1,4,1)=HRRC(1,4,1)+Gamma*VRR(1,4,0)
      HRR(2,1,1)=HRR(2,1,1)+VRR(2,1,0)
      HRRA(2,1,1)=HRRA(2,1,1)+Alpha*VRR(2,1,0)
      HRRB(2,1,1)=HRRB(2,1,1)+Beta*VRR(2,1,0)
      HRRC(2,1,1)=HRRC(2,1,1)+Gamma*VRR(2,1,0)
      HRR(11,1,1)=HRR(11,1,1)+FnSpB*VRR(2,1,0)
      HRRC(2,2,1)=HRRC(2,2,1)+Gamma*VRR(2,2,0)
      HRR(11,2,1)=HRR(11,2,1)+FnSpB*VRR(2,2,0)
      HRRC(2,3,1)=HRRC(2,3,1)+Gamma*VRR(2,3,0)
      HRR(11,3,1)=HRR(11,3,1)+FnSpB*VRR(2,3,0)
      HRRC(2,4,1)=HRRC(2,4,1)+Gamma*VRR(2,4,0)
      HRR(11,4,1)=HRR(11,4,1)+FnSpB*VRR(2,4,0)
      HRR(3,1,1)=HRR(3,1,1)+VRR(3,1,0)
      HRRA(3,1,1)=HRRA(3,1,1)+Alpha*VRR(3,1,0)
      HRRB(3,1,1)=HRRB(3,1,1)+Beta*VRR(3,1,0)
      HRRC(3,1,1)=HRRC(3,1,1)+Gamma*VRR(3,1,0)
      HRR(12,1,1)=HRR(12,1,1)+FnSpB*VRR(3,1,0)
      HRRC(3,2,1)=HRRC(3,2,1)+Gamma*VRR(3,2,0)
      HRR(12,2,1)=HRR(12,2,1)+FnSpB*VRR(3,2,0)
      HRRC(3,3,1)=HRRC(3,3,1)+Gamma*VRR(3,3,0)
      HRR(12,3,1)=HRR(12,3,1)+FnSpB*VRR(3,3,0)
      HRRC(3,4,1)=HRRC(3,4,1)+Gamma*VRR(3,4,0)
      HRR(12,4,1)=HRR(12,4,1)+FnSpB*VRR(3,4,0)
      HRR(4,1,1)=HRR(4,1,1)+VRR(4,1,0)
      HRRA(4,1,1)=HRRA(4,1,1)+Alpha*VRR(4,1,0)
      HRRB(4,1,1)=HRRB(4,1,1)+Beta*VRR(4,1,0)
      HRRC(4,1,1)=HRRC(4,1,1)+Gamma*VRR(4,1,0)
      HRR(13,1,1)=HRR(13,1,1)+FnSpB*VRR(4,1,0)
      HRRC(4,2,1)=HRRC(4,2,1)+Gamma*VRR(4,2,0)
      HRR(13,2,1)=HRR(13,2,1)+FnSpB*VRR(4,2,0)
      HRRC(4,3,1)=HRRC(4,3,1)+Gamma*VRR(4,3,0)
      HRR(13,3,1)=HRR(13,3,1)+FnSpB*VRR(4,3,0)
      HRRC(4,4,1)=HRRC(4,4,1)+Gamma*VRR(4,4,0)
      HRR(13,4,1)=HRR(13,4,1)+FnSpB*VRR(4,4,0)
      HRR(5,1,1)=HRR(5,1,1)+VRR(5,1,0)
      HRRA(5,1,1)=HRRA(5,1,1)+Alpha*VRR(5,1,0)
      HRRB(5,1,1)=HRRB(5,1,1)+Beta*VRR(5,1,0)
      HRRC(5,1,1)=HRRC(5,1,1)+Gamma*VRR(5,1,0)
      HRRC(5,2,1)=HRRC(5,2,1)+Gamma*VRR(5,2,0)
      HRRC(5,3,1)=HRRC(5,3,1)+Gamma*VRR(5,3,0)
      HRRC(5,4,1)=HRRC(5,4,1)+Gamma*VRR(5,4,0)
      HRR(6,1,1)=HRR(6,1,1)+VRR(6,1,0)
      HRRA(6,1,1)=HRRA(6,1,1)+Alpha*VRR(6,1,0)
      HRRB(6,1,1)=HRRB(6,1,1)+Beta*VRR(6,1,0)
      HRRC(6,1,1)=HRRC(6,1,1)+Gamma*VRR(6,1,0)
      HRRC(6,2,1)=HRRC(6,2,1)+Gamma*VRR(6,2,0)
      HRRC(6,3,1)=HRRC(6,3,1)+Gamma*VRR(6,3,0)
      HRRC(6,4,1)=HRRC(6,4,1)+Gamma*VRR(6,4,0)
      HRR(7,1,1)=HRR(7,1,1)+VRR(7,1,0)
      HRRA(7,1,1)=HRRA(7,1,1)+Alpha*VRR(7,1,0)
      HRRB(7,1,1)=HRRB(7,1,1)+Beta*VRR(7,1,0)
      HRRC(7,1,1)=HRRC(7,1,1)+Gamma*VRR(7,1,0)
      HRRC(7,2,1)=HRRC(7,2,1)+Gamma*VRR(7,2,0)
      HRRC(7,3,1)=HRRC(7,3,1)+Gamma*VRR(7,3,0)
      HRRC(7,4,1)=HRRC(7,4,1)+Gamma*VRR(7,4,0)
      HRR(8,1,1)=HRR(8,1,1)+VRR(8,1,0)
      HRRA(8,1,1)=HRRA(8,1,1)+Alpha*VRR(8,1,0)
      HRRB(8,1,1)=HRRB(8,1,1)+Beta*VRR(8,1,0)
      HRRC(8,1,1)=HRRC(8,1,1)+Gamma*VRR(8,1,0)
      HRRC(8,2,1)=HRRC(8,2,1)+Gamma*VRR(8,2,0)
      HRRC(8,3,1)=HRRC(8,3,1)+Gamma*VRR(8,3,0)
      HRRC(8,4,1)=HRRC(8,4,1)+Gamma*VRR(8,4,0)
      HRR(9,1,1)=HRR(9,1,1)+VRR(9,1,0)
      HRRA(9,1,1)=HRRA(9,1,1)+Alpha*VRR(9,1,0)
      HRRB(9,1,1)=HRRB(9,1,1)+Beta*VRR(9,1,0)
      HRRC(9,1,1)=HRRC(9,1,1)+Gamma*VRR(9,1,0)
      HRRC(9,2,1)=HRRC(9,2,1)+Gamma*VRR(9,2,0)
      HRRC(9,3,1)=HRRC(9,3,1)+Gamma*VRR(9,3,0)
      HRRC(9,4,1)=HRRC(9,4,1)+Gamma*VRR(9,4,0)
      HRR(10,1,1)=HRR(10,1,1)+VRR(10,1,0)
      HRRA(10,1,1)=HRRA(10,1,1)+Alpha*VRR(10,1,0)
      HRRB(10,1,1)=HRRB(10,1,1)+Beta*VRR(10,1,0)
      HRRC(10,1,1)=HRRC(10,1,1)+Gamma*VRR(10,1,0)
      HRRC(10,2,1)=HRRC(10,2,1)+Gamma*VRR(10,2,0)
      HRRC(10,3,1)=HRRC(10,3,1)+Gamma*VRR(10,3,0)
      HRRC(10,4,1)=HRRC(10,4,1)+Gamma*VRR(10,4,0)
      HRR(11,1,1)=HRR(11,1,1)+VRR(11,1,0)
      HRRA(11,1,1)=HRRA(11,1,1)+Alpha*VRR(11,1,0)
      HRRB(11,1,1)=HRRB(11,1,1)+Beta*VRR(11,1,0)
      HRRC(11,1,1)=HRRC(11,1,1)+Gamma*VRR(11,1,0)
      HRRC(11,2,1)=HRRC(11,2,1)+Gamma*VRR(11,2,0)
      HRRC(11,3,1)=HRRC(11,3,1)+Gamma*VRR(11,3,0)
      HRRC(11,4,1)=HRRC(11,4,1)+Gamma*VRR(11,4,0)
      HRR(12,1,1)=HRR(12,1,1)+VRR(12,1,0)
      HRRA(12,1,1)=HRRA(12,1,1)+Alpha*VRR(12,1,0)
      HRRB(12,1,1)=HRRB(12,1,1)+Beta*VRR(12,1,0)
      HRRC(12,1,1)=HRRC(12,1,1)+Gamma*VRR(12,1,0)
      HRRC(12,2,1)=HRRC(12,2,1)+Gamma*VRR(12,2,0)
      HRRC(12,3,1)=HRRC(12,3,1)+Gamma*VRR(12,3,0)
      HRRC(12,4,1)=HRRC(12,4,1)+Gamma*VRR(12,4,0)
      HRR(13,1,1)=HRR(13,1,1)+VRR(13,1,0)
      HRRA(13,1,1)=HRRA(13,1,1)+Alpha*VRR(13,1,0)
      HRRB(13,1,1)=HRRB(13,1,1)+Beta*VRR(13,1,0)
      HRRC(13,1,1)=HRRC(13,1,1)+Gamma*VRR(13,1,0)
      HRRC(13,2,1)=HRRC(13,2,1)+Gamma*VRR(13,2,0)
      HRRC(13,3,1)=HRRC(13,3,1)+Gamma*VRR(13,3,0)
      HRRC(13,4,1)=HRRC(13,4,1)+Gamma*VRR(13,4,0)
      HRRA(14,1,1)=HRRA(14,1,1)+Alpha*VRR(14,1,0)
      HRRB(14,1,1)=HRRB(14,1,1)+Beta*VRR(14,1,0)
      HRRA(15,1,1)=HRRA(15,1,1)+Alpha*VRR(15,1,0)
      HRRB(15,1,1)=HRRB(15,1,1)+Beta*VRR(15,1,0)
      HRRA(16,1,1)=HRRA(16,1,1)+Alpha*VRR(16,1,0)
      HRRB(16,1,1)=HRRB(16,1,1)+Beta*VRR(16,1,0)
      HRRA(17,1,1)=HRRA(17,1,1)+Alpha*VRR(17,1,0)
      HRRB(17,1,1)=HRRB(17,1,1)+Beta*VRR(17,1,0)
      HRRA(18,1,1)=HRRA(18,1,1)+Alpha*VRR(18,1,0)
      HRRB(18,1,1)=HRRB(18,1,1)+Beta*VRR(18,1,0)
      HRRA(19,1,1)=HRRA(19,1,1)+Alpha*VRR(19,1,0)
      HRRB(19,1,1)=HRRB(19,1,1)+Beta*VRR(19,1,0)
      HRRA(20,1,1)=HRRA(20,1,1)+Alpha*VRR(20,1,0)
      HRRB(20,1,1)=HRRB(20,1,1)+Beta*VRR(20,1,0)
    END SUBROUTINE CNTRCTG3211
