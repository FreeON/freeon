! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (p s|s s) 
! ---------------------------------------------------------- 
SUBROUTINE dIntB3010101(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, &
 OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NINT,PBC,GRADIENTS)
       USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF2
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
      REAL(DOUBLE), DIMENSION(10) :: HRRTmp 
      REAL(DOUBLE), DIMENSION(4,1,1) :: HRR 
      REAL(DOUBLE), DIMENSION(10,1,1) :: HRRA,HRRB 
      REAL(DOUBLE), DIMENSION(4,4,1) :: HRRC 
      REAL(DOUBLE)  :: VRR(10,4,0:2)
      INTEGER       :: OffSet,OA,LDA,GOA,OB,LDB,GOB,OC,LDC,GOC,OD,LDD,GOD,I,J,K,L
      EXTERNAL InitDbl
      CALL InitDbl(4*1,HRR(1,1,1))
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
              W2=(F2_0(L)+T*(F2_1(L)+T*(F2_2(L)+T*(F2_3(L)+T*F2_4(L)))))
              W1=+3.333333333333333D-01*(TwoT*W2+ET)
              W0=TwoT*W1+ET
              VRR(1,1,0)=Upq*W0
              VRR(1,1,1)=Upq*W1
              VRR(1,1,2)=Upq*W2
            ELSE
              InvT=One/T
              SqInvT=DSQRT(InvT)
              VRR(1,1,0)=+8.862269254527580D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              VRR(1,1,1)=+4.431134627263790D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              VRR(1,1,2)=+6.646701940895685D-01*Upq*SqInvT
            ENDIF
            ! Generating (p0|s0)^(1)
            CALL VRRp0s0(10,4,VRR(1,1,1),VRR(1,1,2))
            ! Generating (p0|s0)^(0)
            CALL VRRp0s0(10,4,VRR(1,1,0),VRR(1,1,1))
            ! Generating (d0|s0)^(0)
            CALL VRRd0s0(10,4,VRR(1,1,0),VRR(1,1,1))
            ! Generating (s0|p0)^(1)
            CALL VRRs0p0(10,4,VRR(1,1,1),VRR(1,1,2))
            ! Generating (s0|p0)^(0)
            CALL VRRs0p0(10,4,VRR(1,1,0),VRR(1,1,1))
            ! Generating (p0|p0)^(0)
            CALL VRRp0p0(10,4,VRR(1,1,0),VRR(1,1,1))
            ! Contracting ... 
            CALL DBLAXPY(4,HRR(1,1,1),VRR(1,1,0)) 
            CALL DBLAXPZY(10,HRRA(1,1,1),Alpha,VRR(1,1,0)) 
            CALL DBLAXPZY(10,HRRB(1,1,1),Beta, VRR(1,1,0)) 
            CALL DBLAXPZY(40,HRRC(1,1,1),Gamma,VRR(1,1,0)) 
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! Dont need to generate (p,0|s,s)
      ! Dont need to generate (d,0|s,s)^a
      ! Dont need to generate (d,0|s,s)^b
      ! Dont need to generate (p,0|p,s)^c
      CALL KetHRR31(4,HRRC) 
      DO L=1,1
      
         !K = 1
         CDOffSet=(OC+1-1)*LDC+(OD+L-1)*LDD
         ! Generating (p',s|1,L)  and (p,s'|1,L)
         CALL BraHRR31ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,1,L),&
                          HRRA(1,1,L),HRRB(1,1,L),GRADIENTS(1,1))
         ! Generating (p,s|1_x,L)  and (p,s|1,L_x)
         CALL BraHRR31cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRC(1,2,L),GRADIENTS(1,1))
         ! Generating (p,s|1_y,L)  and (p,s|1,L_y)
         CALL BraHRR31cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRC(1,3,L),GRADIENTS(1,1))
         ! Generating (p,s|1_z,L)  and (p,s|1,L_z)
         CALL BraHRR31cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,3,HRRC(1,4,L),GRADIENTS(1,1))
      ENDDO 
    END SUBROUTINE dIntB3010101
