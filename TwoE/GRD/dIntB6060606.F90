! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (d d|d d) 
! ---------------------------------------------------------- 
SUBROUTINE dIntB6060606(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, &
 OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NINT,PBC,GRADIENTS)
       USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF9
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
      REAL(DOUBLE), DIMENSION(56) :: HRRTmp 
      REAL(DOUBLE), DIMENSION(35,35,10) :: HRR 
      REAL(DOUBLE), DIMENSION(56,35,10) :: HRRA,HRRB 
      REAL(DOUBLE), DIMENSION(35,56,10) :: HRRC 
      REAL(DOUBLE)  :: VRR(56,56,0:9)
      INTEGER       :: OffSet,OA,LDA,GOA,OB,LDB,GOB,OC,LDC,GOC,OD,LDD,GOD,I,J,K,L
      EXTERNAL InitDbl
      CALL InitDbl(35*35,HRR(1,1,1))
      CALL InitDbl(56*35,HRRA(1,1,1))
      CALL InitDbl(56*35,HRRB(1,1,1))
      CALL InitDbl(35*56,HRRC(1,1,1))
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
              W9=(F9_0(L)+T*(F9_1(L)+T*(F9_2(L)+T*(F9_3(L)+T*F9_4(L)))))
              W8=+5.882352941176471D-02*(TwoT*W9+ET)
              W7=+6.666666666666666D-02*(TwoT*W8+ET)
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
              VRR(1,1,8)=Upq*W8
              VRR(1,1,9)=Upq*W9
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
              SqInvT=SqInvT*InvT
              VRR(1,1,8)=+7.017203646741708D+03*Upq*SqInvT
              SqInvT=SqInvT*InvT
              VRR(1,1,9)=+5.964623099730450D+04*Upq*SqInvT
            ENDIF
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! Generating (g,0|d,d)
      CALL KetHRR66(35,HRR) 
      ! Generating (h,0|d,d)^a
      CALL KetHRR66(56,HRRA) 
      ! Generating (h,0|d,d)^b
      CALL KetHRR66(56,HRRB) 
      ! Generating (g,0|f,d)^c
      CALL KetHRR106(35,HRRC) 
      DO L=5,10
      
         !K = 5
         CDOffSet=(OC+5-5)*LDC+(OD+L-5)*LDD
         ! Generating (d',d|5,L)  and (d,d'|5,L)
         CALL BraHRR66ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,5,L),&
                          HRRA(1,5,L),HRRB(1,5,L),GRADIENTS(1,1))
         ! Generating (d,d|5_x,L)  and (d,d|5,L_x)
         HRRTmp(1:35)=HRRC(1:35,11,L)-2D0*HRR(1:35,2,L)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,HRRTmp,GRADIENTS(1,1))
         ! Generating (d,d|5_y,L)  and (d,d|5,L_y)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRC(1,12,L),GRADIENTS(1,1))
         ! Generating (d,d|5_z,L)  and (d,d|5,L_z)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRC(1,15,L),GRADIENTS(1,1))
      
         !K = 6
         CDOffSet=(OC+6-5)*LDC+(OD+L-5)*LDD
         ! Generating (d',d|6,L)  and (d,d'|6,L)
         CALL BraHRR66ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,6,L),&
                          HRRA(1,6,L),HRRB(1,6,L),GRADIENTS(1,1))
         ! Generating (d,d|6_x,L)  and (d,d|6,L_x)
         HRRTmp(1:35)=HRRC(1:35,12,L)-1D0*HRR(1:35,3,L)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,HRRTmp,GRADIENTS(1,1))
         ! Generating (d,d|6_y,L)  and (d,d|6,L_y)
         HRRTmp(1:35)=HRRC(1:35,13,L)-1D0*HRR(1:35,2,L)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRTmp,GRADIENTS(1,1))
         ! Generating (d,d|6_z,L)  and (d,d|6,L_z)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRC(1,16,L),GRADIENTS(1,1))
      
         !K = 7
         CDOffSet=(OC+7-5)*LDC+(OD+L-5)*LDD
         ! Generating (d',d|7,L)  and (d,d'|7,L)
         CALL BraHRR66ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,7,L),&
                          HRRA(1,7,L),HRRB(1,7,L),GRADIENTS(1,1))
         ! Generating (d,d|7_x,L)  and (d,d|7,L_x)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,HRRC(1,13,L),GRADIENTS(1,1))
         ! Generating (d,d|7_y,L)  and (d,d|7,L_y)
         HRRTmp(1:35)=HRRC(1:35,14,L)-2D0*HRR(1:35,3,L)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRTmp,GRADIENTS(1,1))
         ! Generating (d,d|7_z,L)  and (d,d|7,L_z)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRC(1,17,L),GRADIENTS(1,1))
      
         !K = 8
         CDOffSet=(OC+8-5)*LDC+(OD+L-5)*LDD
         ! Generating (d',d|8,L)  and (d,d'|8,L)
         CALL BraHRR66ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,8,L),&
                          HRRA(1,8,L),HRRB(1,8,L),GRADIENTS(1,1))
         ! Generating (d,d|8_x,L)  and (d,d|8,L_x)
         HRRTmp(1:35)=HRRC(1:35,15,L)-1D0*HRR(1:35,4,L)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,HRRTmp,GRADIENTS(1,1))
         ! Generating (d,d|8_y,L)  and (d,d|8,L_y)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRC(1,16,L),GRADIENTS(1,1))
         ! Generating (d,d|8_z,L)  and (d,d|8,L_z)
         HRRTmp(1:35)=HRRC(1:35,18,L)-1D0*HRR(1:35,2,L)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRTmp,GRADIENTS(1,1))
      
         !K = 9
         CDOffSet=(OC+9-5)*LDC+(OD+L-5)*LDD
         ! Generating (d',d|9,L)  and (d,d'|9,L)
         CALL BraHRR66ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,9,L),&
                          HRRA(1,9,L),HRRB(1,9,L),GRADIENTS(1,1))
         ! Generating (d,d|9_x,L)  and (d,d|9,L_x)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,HRRC(1,16,L),GRADIENTS(1,1))
         ! Generating (d,d|9_y,L)  and (d,d|9,L_y)
         HRRTmp(1:35)=HRRC(1:35,17,L)-1D0*HRR(1:35,4,L)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRTmp,GRADIENTS(1,1))
         ! Generating (d,d|9_z,L)  and (d,d|9,L_z)
         HRRTmp(1:35)=HRRC(1:35,19,L)-1D0*HRR(1:35,3,L)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRTmp,GRADIENTS(1,1))
      
         !K = 10
         CDOffSet=(OC+10-5)*LDC+(OD+L-5)*LDD
         ! Generating (d',d|10,L)  and (d,d'|10,L)
         CALL BraHRR66ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,10,L),&
                          HRRA(1,10,L),HRRB(1,10,L),GRADIENTS(1,1))
         ! Generating (d,d|10_x,L)  and (d,d|10,L_x)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,HRRC(1,18,L),GRADIENTS(1,1))
         ! Generating (d,d|10_y,L)  and (d,d|10,L_y)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRC(1,19,L),GRADIENTS(1,1))
         ! Generating (d,d|10_z,L)  and (d,d|10,L_z)
         HRRTmp(1:35)=HRRC(1:35,20,L)-2D0*HRR(1:35,4,L)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRTmp,GRADIENTS(1,1))
      ENDDO 
    END SUBROUTINE dIntB6060606
