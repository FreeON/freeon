! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (p p|p s) 
! ---------------------------------------------------------- 
SUBROUTINE dIntB3030301(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, &
 OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NINT,PBC,GRADIENTS)
       USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF4
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
      REAL(DOUBLE), DIMENSION(20,10,4) :: HRR,HRRA,HRRB,HRRC 
      REAL(DOUBLE)  :: VRR(20,10,0:4)
      INTEGER       :: OffSet,OA,LDA,GOA,OB,LDB,GOB,OC,LDC,GOC,OD,LDD,GOD,I,J,K,L
      EXTERNAL InitDbl
      CALL InitDbl(20*10,HRR(1,1,1))
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
              W4=(F4_0(L)+T*(F4_1(L)+T*(F4_2(L)+T*(F4_3(L)+T*F4_4(L)))))
              W3=+1.428571428571428D-01*(TwoT*W4+ET)
              W2=+2.000000000000000D-01*(TwoT*W3+ET)
              W1=+3.333333333333333D-01*(TwoT*W2+ET)
              W0=TwoT*W1+ET
              VRR(1,1,0)=Upq*W0
              VRR(1,1,1)=Upq*W1
              VRR(1,1,2)=Upq*W2
              VRR(1,1,3)=Upq*W3
              VRR(1,1,4)=Upq*W4
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
            ENDIF
            ! Generating (p0|s0)^(3)
            CALL VRRp0s0(20,10,VRR(1,1,3),VRR(1,1,4))
            ! Generating (p0|s0)^(2)
            CALL VRRp0s0(20,10,VRR(1,1,2),VRR(1,1,3))
            ! Generating (p0|s0)^(1)
            CALL VRRp0s0(20,10,VRR(1,1,1),VRR(1,1,2))
            ! Generating (p0|s0)^(0)
            CALL VRRp0s0(20,10,VRR(1,1,0),VRR(1,1,1))
            ! Generating (d0|s0)^(2)
            CALL VRRd0s0(20,10,VRR(1,1,2),VRR(1,1,3))
            ! Generating (d0|s0)^(1)
            CALL VRRd0s0(20,10,VRR(1,1,1),VRR(1,1,2))
            ! Generating (d0|s0)^(0)
            CALL VRRd0s0(20,10,VRR(1,1,0),VRR(1,1,1))
            ! Generating (f0|s0)^(1)
            CALL VRRf0s0(20,10,VRR(1,1,1),VRR(1,1,2))
            ! Generating (f0|s0)^(0)
            CALL VRRf0s0(20,10,VRR(1,1,0),VRR(1,1,1))
            ! Generating (s0|p0)^(3)
            CALL VRRs0p0(20,10,VRR(1,1,3),VRR(1,1,4))
            ! Generating (s0|p0)^(2)
            CALL VRRs0p0(20,10,VRR(1,1,2),VRR(1,1,3))
            ! Generating (s0|p0)^(1)
            CALL VRRs0p0(20,10,VRR(1,1,1),VRR(1,1,2))
            ! Generating (s0|p0)^(0)
            CALL VRRs0p0(20,10,VRR(1,1,0),VRR(1,1,1))
            ! Generating (p0|p0)^(2)
            CALL VRRp0p0(20,10,VRR(1,1,2),VRR(1,1,3))
            ! Generating (p0|p0)^(1)
            CALL VRRp0p0(20,10,VRR(1,1,1),VRR(1,1,2))
            ! Generating (p0|p0)^(0)
            CALL VRRp0p0(20,10,VRR(1,1,0),VRR(1,1,1))
            ! Generating (d0|p0)^(1)
            CALL VRRd0p0(20,10,VRR(1,1,1),VRR(1,1,2))
            ! Generating (d0|p0)^(0)
            CALL VRRd0p0(20,10,VRR(1,1,0),VRR(1,1,1))
            ! Generating (f0|p0)^(0)
            CALL VRRf0p0(20,10,VRR(1,1,0),VRR(1,1,1))
            ! Generating (s0|d0)^(2)
            CALL VRRs0d0(20,10,VRR(1,1,2),VRR(1,1,3))
            ! Generating (s0|d0)^(1)
            CALL VRRs0d0(20,10,VRR(1,1,1),VRR(1,1,2))
            ! Generating (s0|d0)^(0)
            CALL VRRs0d0(20,10,VRR(1,1,0),VRR(1,1,1))
            ! Generating (p0|d0)^(1)
            CALL VRRp0d0(20,10,VRR(1,1,1),VRR(1,1,2))
            ! Generating (p0|d0)^(0)
            CALL VRRp0d0(20,10,VRR(1,1,0),VRR(1,1,1))
            ! Generating (d0|d0)^(0)
            CALL VRRd0d0(20,10,VRR(1,1,0),VRR(1,1,1))
            ! Contracting ... 
            CALL DBLAXPY( 200,HRR(1,1,1),       VRR(1,1,0)) 
            CALL DBLAXPZY(200,HRRA(1,1,1),Alpha,VRR(1,1,0)) 
            CALL DBLAXPZY(200,HRRB(1,1,1),Beta, VRR(1,1,0)) 
            CALL DBLAXPZY(200,HRRC(1,1,1),Gamma,VRR(1,1,0)) 
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! Generating (d,0|p,s)
      CALL KetHRR31(10,HRR) 
      ! Generating (f,0|p,s)^a
      CALL KetHRR31(10,HRRA) 
      ! Generating (f,0|p,s)^b
      CALL KetHRR31(10,HRRB) 
      ! Generating (d,0|d,s)^c
      CALL KetHRR61(10,HRRC) 
      DO L=1,1
      
         !K = 2
         CDOffSet=(OC+2-2)*LDC+(OD+L-1)*LDD
         ! Generating (p',p|2,L)  and (p,p'|2,L)
         CALL BraHRR33ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,2,L),&
                          HRRA(1,2,L),HRRB(1,2,L),GRADIENTS(1,1))
         ! Generating (p,p|2_x,L)  and (p,p|2,L_x)
         HRRTmp(1:10)=HRRC(1:10,5,L)-1D0*HRR(1:10,1,L)
         CALL BraHRR33cd(OA,OB,LDA,LDB,CDOffSet,1,HRRTmp,GRADIENTS(1,1))
         ! Generating (p,p|2_y,L)  and (p,p|2,L_y)
         CALL BraHRR33cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRC(1,6,L),GRADIENTS(1,1))
         ! Generating (p,p|2_z,L)  and (p,p|2,L_z)
         CALL BraHRR33cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRC(1,8,L),GRADIENTS(1,1))
      
         !K = 3
         CDOffSet=(OC+3-2)*LDC+(OD+L-1)*LDD
         ! Generating (p',p|3,L)  and (p,p'|3,L)
         CALL BraHRR33ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,3,L),&
                          HRRA(1,3,L),HRRB(1,3,L),GRADIENTS(1,1))
         ! Generating (p,p|3_x,L)  and (p,p|3,L_x)
         CALL BraHRR33cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRC(1,6,L),GRADIENTS(1,1))
         ! Generating (p,p|3_y,L)  and (p,p|3,L_y)
         HRRTmp(1:10)=HRRC(1:10,7,L)-1D0*HRR(1:10,1,L)
         CALL BraHRR33cd(OA,OB,LDA,LDB,CDOffSet,2,HRRTmp,GRADIENTS(1,1))
         ! Generating (p,p|3_z,L)  and (p,p|3,L_z)
         CALL BraHRR33cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRC(1,9,L),GRADIENTS(1,1))
      
         !K = 4
         CDOffSet=(OC+4-2)*LDC+(OD+L-1)*LDD
         ! Generating (p',p|4,L)  and (p,p'|4,L)
         CALL BraHRR33ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,4,L),&
                          HRRA(1,4,L),HRRB(1,4,L),GRADIENTS(1,1))
         ! Generating (p,p|4_x,L)  and (p,p|4,L_x)
         CALL BraHRR33cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRC(1,8,L),GRADIENTS(1,1))
         ! Generating (p,p|4_y,L)  and (p,p|4,L_y)
         CALL BraHRR33cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRC(1,9,L),GRADIENTS(1,1))
         ! Generating (p,p|4_z,L)  and (p,p|4,L_z)
         HRRTmp(1:10)=HRRC(1:10,10,L)-1D0*HRR(1:10,1,L)
         CALL BraHRR33cd(OA,OB,LDA,LDB,CDOffSet,2,HRRTmp,GRADIENTS(1,1))
      ENDDO 
    END SUBROUTINE dIntB3030301
