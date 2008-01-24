!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (d d|s s) 
! ---------------------------------------------------------- 
SUBROUTINE dIntB6060101(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, &
 OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NINT,PBC,GRADIENTS,STRESS)
       USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF5
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER        :: LBra,LKet,NINT,CDOffSet
      REAL(DOUBLE)   :: PrmBufB(10,LBra),PrmBufK(10,LKet)
      TYPE(SmallAtomInfo) :: ACInfo,BDInfo
      TYPE(PBCInfo) :: PBC
      REAL(DOUBLE)  :: GRADIENTS(NINT,12)
      REAL(DOUBLE)  :: STRESS(NINT,9)
      REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz
      REAL(DOUBLE)  :: Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz
      REAL(DOUBLE)  :: PQx,PQy,PQz,FPQx,FPQy,FPQz
      REAL(DOUBLE)  :: Zeta,Eta,Omega,Up,Uq,Upq
      REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT
      REAL(DOUBLE)  :: Alpha,Beta,Gamma
      REAL(DOUBLE), DIMENSION(56) :: HRRTmp 
      REAL(DOUBLE), DIMENSION(35,1,1) :: HRR 
      REAL(DOUBLE), DIMENSION(56,1,1) :: HRRA,HRRB 
      REAL(DOUBLE), DIMENSION(35,4,1) :: HRRC 
      REAL(DOUBLE)  :: VRR(56,4,0:5)
      REAL(DOUBLE)  :: VRRS(35,1,0:4,3)
      REAL(DOUBLE)  :: HRRS(35,1,1,9)
      REAL(DOUBLE)  :: TOm,PQJ(3),FP(9)
      INTEGER       :: OffSet,OA,LDA,GOA,OB,LDB,GOB,OC,LDC,GOC,OD,LDD,GOD,I,J,K,L,IJ
      EXTERNAL InitDbl
      CALL InitDbl(35*1,HRR(1,1,1))
      CALL InitDbl(56*1,HRRA(1,1,1))
      CALL InitDbl(56*1,HRRB(1,1,1))
      CALL InitDbl(35*4,HRRC(1,1,1))
      CALL InitDbl(9*1*35*1,HRRS(1,1,1,1))
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
      !
      !This will feel better above!
      FP(1)=PBC%InvBoxSh%D(1,1)*(Ax-Dx)+PBC%InvBoxSh%D(1,2)*(Ay-Dy)+PBC%InvBoxSh%D(1,3)*(Az-Dz)
      FP(2)=                            PBC%InvBoxSh%D(2,2)*(Ay-Dy)+PBC%InvBoxSh%D(2,3)*(Az-Dz)
      FP(3)=                                                        PBC%InvBoxSh%D(3,3)*(Az-Dz)
      FP(4)=PBC%InvBoxSh%D(1,1)*(Cx-Dx)+PBC%InvBoxSh%D(1,2)*(Cy-Dy)+PBC%InvBoxSh%D(1,3)*(Cz-Dz)
      FP(5)=                            PBC%InvBoxSh%D(2,2)*(Cy-Dy)+PBC%InvBoxSh%D(2,3)*(Cz-Dz)
      FP(6)=                                                        PBC%InvBoxSh%D(3,3)*(Cz-Dz)
      FP(7)=PBC%InvBoxSh%D(1,1)*(Bx-Dx)+PBC%InvBoxSh%D(1,2)*(By-Dy)+PBC%InvBoxSh%D(1,3)*(Bz-Dz)
      FP(8)=                            PBC%InvBoxSh%D(2,2)*(By-Dy)+PBC%InvBoxSh%D(2,3)*(Bz-Dz)
      FP(9)=                                                        PBC%InvBoxSh%D(3,3)*(Bz-Dz)
      !
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
            TOm=2.0d0*Omega
            IF(PBC%AutoW%I(1)==1) THEN
              PQJ(1)=ANINT(FPQx-SIGN(1D-15,FPQx));FPQx=FPQx-PQJ(1)
              PQJ(1)=PQJ(1)*TOm
            ELSE
              PQJ(1)=0.0D0
            ENDIF
            IF(PBC%AutoW%I(2)==1) THEN
              PQJ(2)=ANINT(FPQy-SIGN(1D-15,FPQy));FPQy=FPQy-PQJ(2)
              PQJ(2)=PQJ(2)*TOm
            ELSE
              PQJ(2)=0.0D0
            ENDIF
            IF(PBC%AutoW%I(3)==1) THEN
              PQJ(3)=ANINT(FPQz-SIGN(1D-15,FPQz));FPQz=FPQz-PQJ(3)
              PQJ(3)=PQJ(3)*TOm
            ELSE
              PQJ(3)=0.0D0
            ENDIF
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
              W5=(F5_0(L)+T*(F5_1(L)+T*(F5_2(L)+T*(F5_3(L)+T*F5_4(L)))))
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
            ENDIF
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
            ! Generating (f0|s0)^(2)
            CALL VRRf0s0(56,4,VRR(1,1,2),VRR(1,1,3))
            ! Generating (f0|s0)^(1)
            CALL VRRf0s0(56,4,VRR(1,1,1),VRR(1,1,2))
            ! Generating (f0|s0)^(0)
            CALL VRRf0s0(56,4,VRR(1,1,0),VRR(1,1,1))
            ! Generating (g0|s0)^(1)
            CALL VRRg0s0(56,4,VRR(1,1,1),VRR(1,1,2))
            ! Generating (g0|s0)^(0)
            CALL VRRg0s0(56,4,VRR(1,1,0),VRR(1,1,1))
            ! Generating (h0|s0)^(0)
            CALL VRRh0s0(56,4,VRR(1,1,0),VRR(1,1,1))
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
            ! Generating (d0|p0)^(2)
            CALL VRRd0p0(56,4,VRR(1,1,2),VRR(1,1,3))
            ! Generating (d0|p0)^(1)
            CALL VRRd0p0(56,4,VRR(1,1,1),VRR(1,1,2))
            ! Generating (d0|p0)^(0)
            CALL VRRd0p0(56,4,VRR(1,1,0),VRR(1,1,1))
            ! Generating (f0|p0)^(1)
            CALL VRRf0p0(56,4,VRR(1,1,1),VRR(1,1,2))
            ! Generating (f0|p0)^(0)
            CALL VRRf0p0(56,4,VRR(1,1,0),VRR(1,1,1))
            ! Generating (g0|p0)^(0)
            CALL VRRg0p0(56,4,VRR(1,1,0),VRR(1,1,1))
            !MAY BE BETTER TO PUT WHAT FOLLOWS IN A DO LOOP!
            IF(PBC%AutoW%I(1).EQ.1) THEN
            VRRS(1,1,0,1)=PQx*VRR(1,1,1)
            VRRS(1,1,1,1)=PQx*VRR(1,1,2)
            VRRS(1,1,2,1)=PQx*VRR(1,1,3)
            VRRS(1,1,3,1)=PQx*VRR(1,1,4)
            VRRS(1,1,4,1)=PQx*VRR(1,1,5)
            ! MIC-VRR: Generating [p0|s0]^(3)
            VRRS(2,1,3,1)=PAx*VRRS(1,1,3,1)+WPx*VRRS(1,1,4,1)+r1x2Z*VRR(1,1,4)
            VRRS(3,1,3,1)=PAy*VRRS(1,1,3,1)+WPy*VRRS(1,1,4,1) 
            VRRS(4,1,3,1)=PAz*VRRS(1,1,3,1)+WPz*VRRS(1,1,4,1) 
            ! MIC-VRR: Generating [p0|s0]^(2)
            VRRS(2,1,2,1)=PAx*VRRS(1,1,2,1)+WPx*VRRS(1,1,3,1)+r1x2Z*VRR(1,1,3)
            VRRS(3,1,2,1)=PAy*VRRS(1,1,2,1)+WPy*VRRS(1,1,3,1) 
            VRRS(4,1,2,1)=PAz*VRRS(1,1,2,1)+WPz*VRRS(1,1,3,1) 
            ! MIC-VRR: Generating [p0|s0]^(1)
            VRRS(2,1,1,1)=PAx*VRRS(1,1,1,1)+WPx*VRRS(1,1,2,1)+r1x2Z*VRR(1,1,2)
            VRRS(3,1,1,1)=PAy*VRRS(1,1,1,1)+WPy*VRRS(1,1,2,1) 
            VRRS(4,1,1,1)=PAz*VRRS(1,1,1,1)+WPz*VRRS(1,1,2,1) 
            ! MIC-VRR: Generating [p0|s0]^(0)
            VRRS(2,1,0,1)=PAx*VRRS(1,1,0,1)+WPx*VRRS(1,1,1,1)+r1x2Z*VRR(1,1,1)
            VRRS(3,1,0,1)=PAy*VRRS(1,1,0,1)+WPy*VRRS(1,1,1,1) 
            VRRS(4,1,0,1)=PAz*VRRS(1,1,0,1)+WPz*VRRS(1,1,1,1) 
            ! MIC-VRR: Generating [d0|s0]^(2)
            VRRS( 5,1,2,1)=PAx*VRRS(2,1,2,1)+r1x2Z*(VRRS(1,1,2,1)-ExZpE*VRRS(1,1,3,1))+WPx*VRRS(2,1,3,1)+r1x2Z*VRR(2,1,3)
            VRRS( 6,1,2,1)=PAx*VRRS(3,1,2,1)+WPx*VRRS(3,1,3,1)+r1x2Z*VRR(3,1,3)
            VRRS( 7,1,2,1)=PAy*VRRS(3,1,2,1)+r1x2Z*(VRRS(1,1,2,1)-ExZpE*VRRS(1,1,3,1))+WPy*VRRS(3,1,3,1)
            VRRS( 8,1,2,1)=PAx*VRRS(4,1,2,1)+WPx*VRRS(4,1,3,1)+r1x2Z*VRR(4,1,3)
            VRRS( 9,1,2,1)=PAy*VRRS(4,1,2,1)+WPy*VRRS(4,1,3,1)
            VRRS(10,1,2,1)=PAz*VRRS(4,1,2,1)+r1x2Z*(VRRS(1,1,2,1)-ExZpE*VRRS(1,1,3,1))+WPz*VRRS(4,1,3,1)
            ! MIC-VRR: Generating [d0|s0]^(1)
            VRRS( 5,1,1,1)=PAx*VRRS(2,1,1,1)+r1x2Z*(VRRS(1,1,1,1)-ExZpE*VRRS(1,1,2,1))+WPx*VRRS(2,1,2,1)+r1x2Z*VRR(2,1,2)
            VRRS( 6,1,1,1)=PAx*VRRS(3,1,1,1)+WPx*VRRS(3,1,2,1)+r1x2Z*VRR(3,1,2)
            VRRS( 7,1,1,1)=PAy*VRRS(3,1,1,1)+r1x2Z*(VRRS(1,1,1,1)-ExZpE*VRRS(1,1,2,1))+WPy*VRRS(3,1,2,1)
            VRRS( 8,1,1,1)=PAx*VRRS(4,1,1,1)+WPx*VRRS(4,1,2,1)+r1x2Z*VRR(4,1,2)
            VRRS( 9,1,1,1)=PAy*VRRS(4,1,1,1)+WPy*VRRS(4,1,2,1)
            VRRS(10,1,1,1)=PAz*VRRS(4,1,1,1)+r1x2Z*(VRRS(1,1,1,1)-ExZpE*VRRS(1,1,2,1))+WPz*VRRS(4,1,2,1)
            ! MIC-VRR: Generating [d0|s0]^(0)
            VRRS( 5,1,0,1)=PAx*VRRS(2,1,0,1)+r1x2Z*(VRRS(1,1,0,1)-ExZpE*VRRS(1,1,1,1))+WPx*VRRS(2,1,1,1)+r1x2Z*VRR(2,1,1)
            VRRS( 6,1,0,1)=PAx*VRRS(3,1,0,1)+WPx*VRRS(3,1,1,1)+r1x2Z*VRR(3,1,1)
            VRRS( 7,1,0,1)=PAy*VRRS(3,1,0,1)+r1x2Z*(VRRS(1,1,0,1)-ExZpE*VRRS(1,1,1,1))+WPy*VRRS(3,1,1,1)
            VRRS( 8,1,0,1)=PAx*VRRS(4,1,0,1)+WPx*VRRS(4,1,1,1)+r1x2Z*VRR(4,1,1)
            VRRS( 9,1,0,1)=PAy*VRRS(4,1,0,1)+WPy*VRRS(4,1,1,1)
            VRRS(10,1,0,1)=PAz*VRRS(4,1,0,1)+r1x2Z*(VRRS(1,1,0,1)-ExZpE*VRRS(1,1,1,1))+WPz*VRRS(4,1,1,1)
            ! MIC-VRR: Generating [f0|s0]^(1)
            CALL MVRRf0s0(1,35,1,VRRS(1,1,1,1),VRRS(1,1,2,1),56,4,VRR(1,1,2))
            ! MIC-VRR: Generating [f0|s0]^(0)
            CALL MVRRf0s0(1,35,1,VRRS(1,1,0,1),VRRS(1,1,1,1),56,4,VRR(1,1,1))
            ! MIC-VRR: Generating [g0|s0]^(0)
            CALL MVRRg0s0(1,35,1,VRRS(1,1,0,1),VRRS(1,1,1,1),56,4,VRR(1,1,1))
            ENDIF
            IF(PBC%AutoW%I(2).EQ.1) THEN
            VRRS(1,1,0,2)=PQy*VRR(1,1,1)
            VRRS(1,1,1,2)=PQy*VRR(1,1,2)
            VRRS(1,1,2,2)=PQy*VRR(1,1,3)
            VRRS(1,1,3,2)=PQy*VRR(1,1,4)
            VRRS(1,1,4,2)=PQy*VRR(1,1,5)
            ! MIC-VRR: Generating [p0|s0]^(3)
            VRRS(2,1,3,2)=PAx*VRRS(1,1,3,2)+WPx*VRRS(1,1,4,2)
            VRRS(3,1,3,2)=PAy*VRRS(1,1,3,2)+WPy*VRRS(1,1,4,2)+r1x2Z*VRR(1,1,4)
            VRRS(4,1,3,2)=PAz*VRRS(1,1,3,2)+WPz*VRRS(1,1,4,2)
            ! MIC-VRR: Generating [p0|s0]^(2)
            VRRS(2,1,2,2)=PAx*VRRS(1,1,2,2)+WPx*VRRS(1,1,3,2)
            VRRS(3,1,2,2)=PAy*VRRS(1,1,2,2)+WPy*VRRS(1,1,3,2)+r1x2Z*VRR(1,1,3)
            VRRS(4,1,2,2)=PAz*VRRS(1,1,2,2)+WPz*VRRS(1,1,3,2)
            ! MIC-VRR: Generating [p0|s0]^(1)
            VRRS(2,1,1,2)=PAx*VRRS(1,1,1,2)+WPx*VRRS(1,1,2,2)
            VRRS(3,1,1,2)=PAy*VRRS(1,1,1,2)+WPy*VRRS(1,1,2,2)+r1x2Z*VRR(1,1,2)
            VRRS(4,1,1,2)=PAz*VRRS(1,1,1,2)+WPz*VRRS(1,1,2,2)
            ! MIC-VRR: Generating [p0|s0]^(0)
            VRRS(2,1,0,2)=PAx*VRRS(1,1,0,2)+WPx*VRRS(1,1,1,2)
            VRRS(3,1,0,2)=PAy*VRRS(1,1,0,2)+WPy*VRRS(1,1,1,2)+r1x2Z*VRR(1,1,1)
            VRRS(4,1,0,2)=PAz*VRRS(1,1,0,2)+WPz*VRRS(1,1,1,2)
            ! MIC-VRR: Generating [d0|s0]^(2)
            VRRS( 5,1,2,2)=PAx*VRRS(2,1,2,2)+r1x2Z*(VRRS(1,1,2,2)-ExZpE*VRRS(1,1,3,2))+WPx*VRRS(2,1,3,2)
            VRRS( 6,1,2,2)=PAx*VRRS(3,1,2,2)+WPx*VRRS(3,1,3,2)
            VRRS( 7,1,2,2)=PAy*VRRS(3,1,2,2)+r1x2Z*(VRRS(1,1,2,2)-ExZpE*VRRS(1,1,3,2))+WPy*VRRS(3,1,3,2)+r1x2Z*VRR(3,1,3)
            VRRS( 8,1,2,2)=PAx*VRRS(4,1,2,2)+WPx*VRRS(4,1,3,2)
            VRRS( 9,1,2,2)=PAy*VRRS(4,1,2,2)+WPy*VRRS(4,1,3,2)+r1x2Z*VRR(4,1,3)
            VRRS(10,1,2,2)=PAz*VRRS(4,1,2,2)+r1x2Z*(VRRS(1,1,2,2)-ExZpE*VRRS(1,1,3,2))+WPz*VRRS(4,1,3,2)
            ! MIC-VRR: Generating [d0|s0]^(1)
            VRRS( 5,1,1,2)=PAx*VRRS(2,1,1,2)+r1x2Z*(VRRS(1,1,1,2)-ExZpE*VRRS(1,1,2,2))+WPx*VRRS(2,1,2,2)
            VRRS( 6,1,1,2)=PAx*VRRS(3,1,1,2)+WPx*VRRS(3,1,2,2)
            VRRS( 7,1,1,2)=PAy*VRRS(3,1,1,2)+r1x2Z*(VRRS(1,1,1,2)-ExZpE*VRRS(1,1,2,2))+WPy*VRRS(3,1,2,2)+r1x2Z*VRR(3,1,2)
            VRRS( 8,1,1,2)=PAx*VRRS(4,1,1,2)+WPx*VRRS(4,1,2,2)
            VRRS( 9,1,1,2)=PAy*VRRS(4,1,1,2)+WPy*VRRS(4,1,2,2)+r1x2Z*VRR(4,1,2)
            VRRS(10,1,1,2)=PAz*VRRS(4,1,1,2)+r1x2Z*(VRRS(1,1,1,2)-ExZpE*VRRS(1,1,2,2))+WPz*VRRS(4,1,2,2)
            ! MIC-VRR: Generating [d0|s0]^(0)
            VRRS( 5,1,0,2)=PAx*VRRS(2,1,0,2)+r1x2Z*(VRRS(1,1,0,2)-ExZpE*VRRS(1,1,1,2))+WPx*VRRS(2,1,1,2)
            VRRS( 6,1,0,2)=PAx*VRRS(3,1,0,2)+WPx*VRRS(3,1,1,2)
            VRRS( 7,1,0,2)=PAy*VRRS(3,1,0,2)+r1x2Z*(VRRS(1,1,0,2)-ExZpE*VRRS(1,1,1,2))+WPy*VRRS(3,1,1,2)+r1x2Z*VRR(3,1,1)
            VRRS( 8,1,0,2)=PAx*VRRS(4,1,0,2)+WPx*VRRS(4,1,1,2)
            VRRS( 9,1,0,2)=PAy*VRRS(4,1,0,2)+WPy*VRRS(4,1,1,2)+r1x2Z*VRR(4,1,1)
            VRRS(10,1,0,2)=PAz*VRRS(4,1,0,2)+r1x2Z*(VRRS(1,1,0,2)-ExZpE*VRRS(1,1,1,2))+WPz*VRRS(4,1,1,2)
            ! MIC-VRR: Generating [f0|s0]^(1)
            CALL MVRRf0s0(2,35,1,VRRS(1,1,1,2),VRRS(1,1,2,2),56,4,VRR(1,1,2))
            ! MIC-VRR: Generating [f0|s0]^(0)
            CALL MVRRf0s0(2,35,1,VRRS(1,1,0,2),VRRS(1,1,1,2),56,4,VRR(1,1,1))
            ! MIC-VRR: Generating [g0|s0]^(0)
            CALL MVRRg0s0(2,35,1,VRRS(1,1,0,2),VRRS(1,1,1,2),56,4,VRR(1,1,1))
            ENDIF
            IF(PBC%AutoW%I(3).EQ.1) THEN
            VRRS(1,1,0,3)=PQz*VRR(1,1,1)
            VRRS(1,1,1,3)=PQz*VRR(1,1,2)
            VRRS(1,1,2,3)=PQz*VRR(1,1,3)
            VRRS(1,1,3,3)=PQz*VRR(1,1,4)
            VRRS(1,1,4,3)=PQz*VRR(1,1,5)
            ! MIC-VRR: Generating [p0|s0]^(3)
            VRRS(2,1,3,3)=PAx*VRRS(1,1,3,3)+WPx*VRRS(1,1,4,3)
            VRRS(3,1,3,3)=PAy*VRRS(1,1,3,3)+WPy*VRRS(1,1,4,3)
            VRRS(4,1,3,3)=PAz*VRRS(1,1,3,3)+WPz*VRRS(1,1,4,3)+r1x2Z*VRR(1,1,4)
            ! MIC-VRR: Generating [p0|s0]^(2)
            VRRS(2,1,2,3)=PAx*VRRS(1,1,2,3)+WPx*VRRS(1,1,3,3)
            VRRS(3,1,2,3)=PAy*VRRS(1,1,2,3)+WPy*VRRS(1,1,3,3)
            VRRS(4,1,2,3)=PAz*VRRS(1,1,2,3)+WPz*VRRS(1,1,3,3)+r1x2Z*VRR(1,1,3)
            ! MIC-VRR: Generating [p0|s0]^(1)
            VRRS(2,1,1,3)=PAx*VRRS(1,1,1,3)+WPx*VRRS(1,1,2,3)
            VRRS(3,1,1,3)=PAy*VRRS(1,1,1,3)+WPy*VRRS(1,1,2,3)
            VRRS(4,1,1,3)=PAz*VRRS(1,1,1,3)+WPz*VRRS(1,1,2,3)+r1x2Z*VRR(1,1,2)
            ! MIC-VRR: Generating [p0|s0]^(0)
            VRRS(2,1,0,3)=PAx*VRRS(1,1,0,3)+WPx*VRRS(1,1,1,3)
            VRRS(3,1,0,3)=PAy*VRRS(1,1,0,3)+WPy*VRRS(1,1,1,3)
            VRRS(4,1,0,3)=PAz*VRRS(1,1,0,3)+WPz*VRRS(1,1,1,3)+r1x2Z*VRR(1,1,1)
            ! MIC-VRR: Generating [d0|s0]^(2)
            VRRS( 5,1,2,3)=PAx*VRRS(2,1,2,3)+r1x2Z*(VRRS(1,1,2,3)-ExZpE*VRRS(1,1,3,3))+WPx*VRRS(2,1,3,3)
            VRRS( 6,1,2,3)=PAx*VRRS(3,1,2,3)+WPx*VRRS(3,1,3,3)
            VRRS( 7,1,2,3)=PAy*VRRS(3,1,2,3)+r1x2Z*(VRRS(1,1,2,3)-ExZpE*VRRS(1,1,3,3))+WPy*VRRS(3,1,3,3)
            VRRS( 8,1,2,3)=PAx*VRRS(4,1,2,3)+WPx*VRRS(4,1,3,3)
            VRRS( 9,1,2,3)=PAy*VRRS(4,1,2,3)+WPy*VRRS(4,1,3,3)
            VRRS(10,1,2,3)=PAz*VRRS(4,1,2,3)+r1x2Z*(VRRS(1,1,2,3)-ExZpE*VRRS(1,1,3,3))+WPz*VRRS(4,1,3,3)+r1x2Z*VRR(4,1,3)
            ! MIC-VRR: Generating [d0|s0]^(1)
            VRRS( 5,1,1,3)=PAx*VRRS(2,1,1,3)+r1x2Z*(VRRS(1,1,1,3)-ExZpE*VRRS(1,1,2,3))+WPx*VRRS(2,1,2,3)
            VRRS( 6,1,1,3)=PAx*VRRS(3,1,1,3)+WPx*VRRS(3,1,2,3)
            VRRS( 7,1,1,3)=PAy*VRRS(3,1,1,3)+r1x2Z*(VRRS(1,1,1,3)-ExZpE*VRRS(1,1,2,3))+WPy*VRRS(3,1,2,3)
            VRRS( 8,1,1,3)=PAx*VRRS(4,1,1,3)+WPx*VRRS(4,1,2,3)
            VRRS( 9,1,1,3)=PAy*VRRS(4,1,1,3)+WPy*VRRS(4,1,2,3)
            VRRS(10,1,1,3)=PAz*VRRS(4,1,1,3)+r1x2Z*(VRRS(1,1,1,3)-ExZpE*VRRS(1,1,2,3))+WPz*VRRS(4,1,2,3)+r1x2Z*VRR(4,1,2)
            ! MIC-VRR: Generating [d0|s0]^(0)
            VRRS( 5,1,0,3)=PAx*VRRS(2,1,0,3)+r1x2Z*(VRRS(1,1,0,3)-ExZpE*VRRS(1,1,1,3))+WPx*VRRS(2,1,1,3)
            VRRS( 6,1,0,3)=PAx*VRRS(3,1,0,3)+WPx*VRRS(3,1,1,3)
            VRRS( 7,1,0,3)=PAy*VRRS(3,1,0,3)+r1x2Z*(VRRS(1,1,0,3)-ExZpE*VRRS(1,1,1,3))+WPy*VRRS(3,1,1,3)
            VRRS( 8,1,0,3)=PAx*VRRS(4,1,0,3)+WPx*VRRS(4,1,1,3)
            VRRS( 9,1,0,3)=PAy*VRRS(4,1,0,3)+WPy*VRRS(4,1,1,3)
            VRRS(10,1,0,3)=PAz*VRRS(4,1,0,3)+r1x2Z*(VRRS(1,1,0,3)-ExZpE*VRRS(1,1,1,3))+WPz*VRRS(4,1,1,3)+r1x2Z*VRR(4,1,1)
            ! MIC-VRR: Generating [f0|s0]^(1)
            CALL MVRRf0s0(3,35,1,VRRS(1,1,1,3),VRRS(1,1,2,3),56,4,VRR(1,1,2))
            ! MIC-VRR: Generating [f0|s0]^(0)
            CALL MVRRf0s0(3,35,1,VRRS(1,1,0,3),VRRS(1,1,1,3),56,4,VRR(1,1,1))
            ! MIC-VRR: Generating [g0|s0]^(0)
            CALL MVRRg0s0(3,35,1,VRRS(1,1,0,3),VRRS(1,1,1,3),56,4,VRR(1,1,1))
            ENDIF
            ! Contracting ... 
            CALL CNTRCTG6611(VRR,HRR,Alpha,HRRA,Beta,HRRB,Gamma,HRRC, &
                       VRRS,HRRS(1,1,1,1),PQJ(1),PBC%AutoW%I(1))
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! Dont need to generate (g,0|s,s)
      ! Dont need to generate (h,0|s,s)^a
      ! Dont need to generate (h,0|s,s)^b
      ! Dont need to generate (g,0|p,s)^c
      ! Stress: No need to generate [d,0|s,s] 
      DO L=1,1
      
         !K = 1
         CDOffSet=(OC+1-1)*LDC+(OD+L-1)*LDD
         ! Generating (d',d|1,L)  and (d,d'|1,L)
         CALL BraHRR66ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(1,1,L),&
                      HRRA(1,1,L),HRRB(1,1,L),GRADIENTS(1,1),FP(1),&
                      STRESS(1,1))
         ! Generating (d,d|1_x,L)  and (d,d|1,L_x)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,&
                      HRRC(1,2,L),GRADIENTS(1,1),FP(1),STRESS(1,1))
         ! Generating (d,d|1_y,L)  and (d,d|1,L_y)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,&
                      HRRC(1,3,L),GRADIENTS(1,1),FP(1),STRESS(1,1))
         ! Generating (d,d|1_z,L)  and (d,d|1,L_z)
         CALL BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,&
                      HRRC(1,4,L),GRADIENTS(1,1),FP(1),STRESS(1,1))
      ENDDO 
      ! Stress: Generating (d,d|s,s)^(1) 
      DO J=1,3
      DO I=1,3
      IJ=3*(J-1)+I
        DO L=1,1
          DO K=1,1
            CDOffSet=(OC+K-1)*LDC+(OD+L-1)*LDD 
            CALL BraHRR66(OA,OB,LDA,LDB,CDOffSet,HRRS(1,K,L,IJ),STRESS(1,IJ))
          ENDDO 
        ENDDO 
      ENDDO 
      ENDDO 
   END SUBROUTINE dIntB6060101
    SUBROUTINE CNTRCTG6611(VRR,HRR,Alpha,HRRA,Beta,HRRB,Gamma,HRRC,VRRS,HRRS,PQJ,IW)
      USE DerivedTypes
      USE VScratchB
      IMPLICIT NONE
      INTEGER :: K
      REAL(DOUBLE)  :: Alpha,Beta,Gamma
      REAL(DOUBLE), DIMENSION(35,1,1) :: HRR 
      REAL(DOUBLE), DIMENSION(56,1,1) :: HRRA,HRRB 
      REAL(DOUBLE), DIMENSION(35,4,1) :: HRRC 
      REAL(DOUBLE)  :: VRR(56,4,0:5)
      REAL(DOUBLE)  :: VRRS(35,1,0:4,3)
      REAL(DOUBLE)  :: HRRS(35,1,1,9),PQJ(3)
      INTEGER :: IJ,J,I,IW(3)
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
      HRRC(2,2,1)=HRRC(2,2,1)+Gamma*VRR(2,2,0)
      HRRC(2,3,1)=HRRC(2,3,1)+Gamma*VRR(2,3,0)
      HRRC(2,4,1)=HRRC(2,4,1)+Gamma*VRR(2,4,0)
      HRR(3,1,1)=HRR(3,1,1)+VRR(3,1,0)
      HRRA(3,1,1)=HRRA(3,1,1)+Alpha*VRR(3,1,0)
      HRRB(3,1,1)=HRRB(3,1,1)+Beta*VRR(3,1,0)
      HRRC(3,1,1)=HRRC(3,1,1)+Gamma*VRR(3,1,0)
      HRRC(3,2,1)=HRRC(3,2,1)+Gamma*VRR(3,2,0)
      HRRC(3,3,1)=HRRC(3,3,1)+Gamma*VRR(3,3,0)
      HRRC(3,4,1)=HRRC(3,4,1)+Gamma*VRR(3,4,0)
      HRR(4,1,1)=HRR(4,1,1)+VRR(4,1,0)
      HRRA(4,1,1)=HRRA(4,1,1)+Alpha*VRR(4,1,0)
      HRRB(4,1,1)=HRRB(4,1,1)+Beta*VRR(4,1,0)
      HRRC(4,1,1)=HRRC(4,1,1)+Gamma*VRR(4,1,0)
      HRRC(4,2,1)=HRRC(4,2,1)+Gamma*VRR(4,2,0)
      HRRC(4,3,1)=HRRC(4,3,1)+Gamma*VRR(4,3,0)
      HRRC(4,4,1)=HRRC(4,4,1)+Gamma*VRR(4,4,0)
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
      HRR(14,1,1)=HRR(14,1,1)+VRR(14,1,0)
      HRRA(14,1,1)=HRRA(14,1,1)+Alpha*VRR(14,1,0)
      HRRB(14,1,1)=HRRB(14,1,1)+Beta*VRR(14,1,0)
      HRRC(14,1,1)=HRRC(14,1,1)+Gamma*VRR(14,1,0)
      HRRC(14,2,1)=HRRC(14,2,1)+Gamma*VRR(14,2,0)
      HRRC(14,3,1)=HRRC(14,3,1)+Gamma*VRR(14,3,0)
      HRRC(14,4,1)=HRRC(14,4,1)+Gamma*VRR(14,4,0)
      HRR(15,1,1)=HRR(15,1,1)+VRR(15,1,0)
      HRRA(15,1,1)=HRRA(15,1,1)+Alpha*VRR(15,1,0)
      HRRB(15,1,1)=HRRB(15,1,1)+Beta*VRR(15,1,0)
      HRRC(15,1,1)=HRRC(15,1,1)+Gamma*VRR(15,1,0)
      HRRC(15,2,1)=HRRC(15,2,1)+Gamma*VRR(15,2,0)
      HRRC(15,3,1)=HRRC(15,3,1)+Gamma*VRR(15,3,0)
      HRRC(15,4,1)=HRRC(15,4,1)+Gamma*VRR(15,4,0)
      HRR(16,1,1)=HRR(16,1,1)+VRR(16,1,0)
      HRRA(16,1,1)=HRRA(16,1,1)+Alpha*VRR(16,1,0)
      HRRB(16,1,1)=HRRB(16,1,1)+Beta*VRR(16,1,0)
      HRRC(16,1,1)=HRRC(16,1,1)+Gamma*VRR(16,1,0)
      HRRC(16,2,1)=HRRC(16,2,1)+Gamma*VRR(16,2,0)
      HRRC(16,3,1)=HRRC(16,3,1)+Gamma*VRR(16,3,0)
      HRRC(16,4,1)=HRRC(16,4,1)+Gamma*VRR(16,4,0)
      HRR(17,1,1)=HRR(17,1,1)+VRR(17,1,0)
      HRRA(17,1,1)=HRRA(17,1,1)+Alpha*VRR(17,1,0)
      HRRB(17,1,1)=HRRB(17,1,1)+Beta*VRR(17,1,0)
      HRRC(17,1,1)=HRRC(17,1,1)+Gamma*VRR(17,1,0)
      HRRC(17,2,1)=HRRC(17,2,1)+Gamma*VRR(17,2,0)
      HRRC(17,3,1)=HRRC(17,3,1)+Gamma*VRR(17,3,0)
      HRRC(17,4,1)=HRRC(17,4,1)+Gamma*VRR(17,4,0)
      HRR(18,1,1)=HRR(18,1,1)+VRR(18,1,0)
      HRRA(18,1,1)=HRRA(18,1,1)+Alpha*VRR(18,1,0)
      HRRB(18,1,1)=HRRB(18,1,1)+Beta*VRR(18,1,0)
      HRRC(18,1,1)=HRRC(18,1,1)+Gamma*VRR(18,1,0)
      HRRC(18,2,1)=HRRC(18,2,1)+Gamma*VRR(18,2,0)
      HRRC(18,3,1)=HRRC(18,3,1)+Gamma*VRR(18,3,0)
      HRRC(18,4,1)=HRRC(18,4,1)+Gamma*VRR(18,4,0)
      HRR(19,1,1)=HRR(19,1,1)+VRR(19,1,0)
      HRRA(19,1,1)=HRRA(19,1,1)+Alpha*VRR(19,1,0)
      HRRB(19,1,1)=HRRB(19,1,1)+Beta*VRR(19,1,0)
      HRRC(19,1,1)=HRRC(19,1,1)+Gamma*VRR(19,1,0)
      HRRC(19,2,1)=HRRC(19,2,1)+Gamma*VRR(19,2,0)
      HRRC(19,3,1)=HRRC(19,3,1)+Gamma*VRR(19,3,0)
      HRRC(19,4,1)=HRRC(19,4,1)+Gamma*VRR(19,4,0)
      HRR(20,1,1)=HRR(20,1,1)+VRR(20,1,0)
      HRRA(20,1,1)=HRRA(20,1,1)+Alpha*VRR(20,1,0)
      HRRB(20,1,1)=HRRB(20,1,1)+Beta*VRR(20,1,0)
      HRRC(20,1,1)=HRRC(20,1,1)+Gamma*VRR(20,1,0)
      HRRC(20,2,1)=HRRC(20,2,1)+Gamma*VRR(20,2,0)
      HRRC(20,3,1)=HRRC(20,3,1)+Gamma*VRR(20,3,0)
      HRRC(20,4,1)=HRRC(20,4,1)+Gamma*VRR(20,4,0)
      HRR(21,1,1)=HRR(21,1,1)+VRR(21,1,0)
      HRRA(21,1,1)=HRRA(21,1,1)+Alpha*VRR(21,1,0)
      HRRB(21,1,1)=HRRB(21,1,1)+Beta*VRR(21,1,0)
      HRRC(21,1,1)=HRRC(21,1,1)+Gamma*VRR(21,1,0)
      HRRC(21,2,1)=HRRC(21,2,1)+Gamma*VRR(21,2,0)
      HRRC(21,3,1)=HRRC(21,3,1)+Gamma*VRR(21,3,0)
      HRRC(21,4,1)=HRRC(21,4,1)+Gamma*VRR(21,4,0)
      HRR(22,1,1)=HRR(22,1,1)+VRR(22,1,0)
      HRRA(22,1,1)=HRRA(22,1,1)+Alpha*VRR(22,1,0)
      HRRB(22,1,1)=HRRB(22,1,1)+Beta*VRR(22,1,0)
      HRRC(22,1,1)=HRRC(22,1,1)+Gamma*VRR(22,1,0)
      HRRC(22,2,1)=HRRC(22,2,1)+Gamma*VRR(22,2,0)
      HRRC(22,3,1)=HRRC(22,3,1)+Gamma*VRR(22,3,0)
      HRRC(22,4,1)=HRRC(22,4,1)+Gamma*VRR(22,4,0)
      HRR(23,1,1)=HRR(23,1,1)+VRR(23,1,0)
      HRRA(23,1,1)=HRRA(23,1,1)+Alpha*VRR(23,1,0)
      HRRB(23,1,1)=HRRB(23,1,1)+Beta*VRR(23,1,0)
      HRRC(23,1,1)=HRRC(23,1,1)+Gamma*VRR(23,1,0)
      HRRC(23,2,1)=HRRC(23,2,1)+Gamma*VRR(23,2,0)
      HRRC(23,3,1)=HRRC(23,3,1)+Gamma*VRR(23,3,0)
      HRRC(23,4,1)=HRRC(23,4,1)+Gamma*VRR(23,4,0)
      HRR(24,1,1)=HRR(24,1,1)+VRR(24,1,0)
      HRRA(24,1,1)=HRRA(24,1,1)+Alpha*VRR(24,1,0)
      HRRB(24,1,1)=HRRB(24,1,1)+Beta*VRR(24,1,0)
      HRRC(24,1,1)=HRRC(24,1,1)+Gamma*VRR(24,1,0)
      HRRC(24,2,1)=HRRC(24,2,1)+Gamma*VRR(24,2,0)
      HRRC(24,3,1)=HRRC(24,3,1)+Gamma*VRR(24,3,0)
      HRRC(24,4,1)=HRRC(24,4,1)+Gamma*VRR(24,4,0)
      HRR(25,1,1)=HRR(25,1,1)+VRR(25,1,0)
      HRRA(25,1,1)=HRRA(25,1,1)+Alpha*VRR(25,1,0)
      HRRB(25,1,1)=HRRB(25,1,1)+Beta*VRR(25,1,0)
      HRRC(25,1,1)=HRRC(25,1,1)+Gamma*VRR(25,1,0)
      HRRC(25,2,1)=HRRC(25,2,1)+Gamma*VRR(25,2,0)
      HRRC(25,3,1)=HRRC(25,3,1)+Gamma*VRR(25,3,0)
      HRRC(25,4,1)=HRRC(25,4,1)+Gamma*VRR(25,4,0)
      HRR(26,1,1)=HRR(26,1,1)+VRR(26,1,0)
      HRRA(26,1,1)=HRRA(26,1,1)+Alpha*VRR(26,1,0)
      HRRB(26,1,1)=HRRB(26,1,1)+Beta*VRR(26,1,0)
      HRRC(26,1,1)=HRRC(26,1,1)+Gamma*VRR(26,1,0)
      HRRC(26,2,1)=HRRC(26,2,1)+Gamma*VRR(26,2,0)
      HRRC(26,3,1)=HRRC(26,3,1)+Gamma*VRR(26,3,0)
      HRRC(26,4,1)=HRRC(26,4,1)+Gamma*VRR(26,4,0)
      HRR(27,1,1)=HRR(27,1,1)+VRR(27,1,0)
      HRRA(27,1,1)=HRRA(27,1,1)+Alpha*VRR(27,1,0)
      HRRB(27,1,1)=HRRB(27,1,1)+Beta*VRR(27,1,0)
      HRRC(27,1,1)=HRRC(27,1,1)+Gamma*VRR(27,1,0)
      HRRC(27,2,1)=HRRC(27,2,1)+Gamma*VRR(27,2,0)
      HRRC(27,3,1)=HRRC(27,3,1)+Gamma*VRR(27,3,0)
      HRRC(27,4,1)=HRRC(27,4,1)+Gamma*VRR(27,4,0)
      HRR(28,1,1)=HRR(28,1,1)+VRR(28,1,0)
      HRRA(28,1,1)=HRRA(28,1,1)+Alpha*VRR(28,1,0)
      HRRB(28,1,1)=HRRB(28,1,1)+Beta*VRR(28,1,0)
      HRRC(28,1,1)=HRRC(28,1,1)+Gamma*VRR(28,1,0)
      HRRC(28,2,1)=HRRC(28,2,1)+Gamma*VRR(28,2,0)
      HRRC(28,3,1)=HRRC(28,3,1)+Gamma*VRR(28,3,0)
      HRRC(28,4,1)=HRRC(28,4,1)+Gamma*VRR(28,4,0)
      HRR(29,1,1)=HRR(29,1,1)+VRR(29,1,0)
      HRRA(29,1,1)=HRRA(29,1,1)+Alpha*VRR(29,1,0)
      HRRB(29,1,1)=HRRB(29,1,1)+Beta*VRR(29,1,0)
      HRRC(29,1,1)=HRRC(29,1,1)+Gamma*VRR(29,1,0)
      HRRC(29,2,1)=HRRC(29,2,1)+Gamma*VRR(29,2,0)
      HRRC(29,3,1)=HRRC(29,3,1)+Gamma*VRR(29,3,0)
      HRRC(29,4,1)=HRRC(29,4,1)+Gamma*VRR(29,4,0)
      HRR(30,1,1)=HRR(30,1,1)+VRR(30,1,0)
      HRRA(30,1,1)=HRRA(30,1,1)+Alpha*VRR(30,1,0)
      HRRB(30,1,1)=HRRB(30,1,1)+Beta*VRR(30,1,0)
      HRRC(30,1,1)=HRRC(30,1,1)+Gamma*VRR(30,1,0)
      HRRC(30,2,1)=HRRC(30,2,1)+Gamma*VRR(30,2,0)
      HRRC(30,3,1)=HRRC(30,3,1)+Gamma*VRR(30,3,0)
      HRRC(30,4,1)=HRRC(30,4,1)+Gamma*VRR(30,4,0)
      HRR(31,1,1)=HRR(31,1,1)+VRR(31,1,0)
      HRRA(31,1,1)=HRRA(31,1,1)+Alpha*VRR(31,1,0)
      HRRB(31,1,1)=HRRB(31,1,1)+Beta*VRR(31,1,0)
      HRRC(31,1,1)=HRRC(31,1,1)+Gamma*VRR(31,1,0)
      HRRC(31,2,1)=HRRC(31,2,1)+Gamma*VRR(31,2,0)
      HRRC(31,3,1)=HRRC(31,3,1)+Gamma*VRR(31,3,0)
      HRRC(31,4,1)=HRRC(31,4,1)+Gamma*VRR(31,4,0)
      HRR(32,1,1)=HRR(32,1,1)+VRR(32,1,0)
      HRRA(32,1,1)=HRRA(32,1,1)+Alpha*VRR(32,1,0)
      HRRB(32,1,1)=HRRB(32,1,1)+Beta*VRR(32,1,0)
      HRRC(32,1,1)=HRRC(32,1,1)+Gamma*VRR(32,1,0)
      HRRC(32,2,1)=HRRC(32,2,1)+Gamma*VRR(32,2,0)
      HRRC(32,3,1)=HRRC(32,3,1)+Gamma*VRR(32,3,0)
      HRRC(32,4,1)=HRRC(32,4,1)+Gamma*VRR(32,4,0)
      HRR(33,1,1)=HRR(33,1,1)+VRR(33,1,0)
      HRRA(33,1,1)=HRRA(33,1,1)+Alpha*VRR(33,1,0)
      HRRB(33,1,1)=HRRB(33,1,1)+Beta*VRR(33,1,0)
      HRRC(33,1,1)=HRRC(33,1,1)+Gamma*VRR(33,1,0)
      HRRC(33,2,1)=HRRC(33,2,1)+Gamma*VRR(33,2,0)
      HRRC(33,3,1)=HRRC(33,3,1)+Gamma*VRR(33,3,0)
      HRRC(33,4,1)=HRRC(33,4,1)+Gamma*VRR(33,4,0)
      HRR(34,1,1)=HRR(34,1,1)+VRR(34,1,0)
      HRRA(34,1,1)=HRRA(34,1,1)+Alpha*VRR(34,1,0)
      HRRB(34,1,1)=HRRB(34,1,1)+Beta*VRR(34,1,0)
      HRRC(34,1,1)=HRRC(34,1,1)+Gamma*VRR(34,1,0)
      HRRC(34,2,1)=HRRC(34,2,1)+Gamma*VRR(34,2,0)
      HRRC(34,3,1)=HRRC(34,3,1)+Gamma*VRR(34,3,0)
      HRRC(34,4,1)=HRRC(34,4,1)+Gamma*VRR(34,4,0)
      HRR(35,1,1)=HRR(35,1,1)+VRR(35,1,0)
      HRRA(35,1,1)=HRRA(35,1,1)+Alpha*VRR(35,1,0)
      HRRB(35,1,1)=HRRB(35,1,1)+Beta*VRR(35,1,0)
      HRRC(35,1,1)=HRRC(35,1,1)+Gamma*VRR(35,1,0)
      HRRC(35,2,1)=HRRC(35,2,1)+Gamma*VRR(35,2,0)
      HRRC(35,3,1)=HRRC(35,3,1)+Gamma*VRR(35,3,0)
      HRRC(35,4,1)=HRRC(35,4,1)+Gamma*VRR(35,4,0)
      HRRA(36,1,1)=HRRA(36,1,1)+Alpha*VRR(36,1,0)
      HRRB(36,1,1)=HRRB(36,1,1)+Beta*VRR(36,1,0)
      HRRA(37,1,1)=HRRA(37,1,1)+Alpha*VRR(37,1,0)
      HRRB(37,1,1)=HRRB(37,1,1)+Beta*VRR(37,1,0)
      HRRA(38,1,1)=HRRA(38,1,1)+Alpha*VRR(38,1,0)
      HRRB(38,1,1)=HRRB(38,1,1)+Beta*VRR(38,1,0)
      HRRA(39,1,1)=HRRA(39,1,1)+Alpha*VRR(39,1,0)
      HRRB(39,1,1)=HRRB(39,1,1)+Beta*VRR(39,1,0)
      HRRA(40,1,1)=HRRA(40,1,1)+Alpha*VRR(40,1,0)
      HRRB(40,1,1)=HRRB(40,1,1)+Beta*VRR(40,1,0)
      HRRA(41,1,1)=HRRA(41,1,1)+Alpha*VRR(41,1,0)
      HRRB(41,1,1)=HRRB(41,1,1)+Beta*VRR(41,1,0)
      HRRA(42,1,1)=HRRA(42,1,1)+Alpha*VRR(42,1,0)
      HRRB(42,1,1)=HRRB(42,1,1)+Beta*VRR(42,1,0)
      HRRA(43,1,1)=HRRA(43,1,1)+Alpha*VRR(43,1,0)
      HRRB(43,1,1)=HRRB(43,1,1)+Beta*VRR(43,1,0)
      HRRA(44,1,1)=HRRA(44,1,1)+Alpha*VRR(44,1,0)
      HRRB(44,1,1)=HRRB(44,1,1)+Beta*VRR(44,1,0)
      HRRA(45,1,1)=HRRA(45,1,1)+Alpha*VRR(45,1,0)
      HRRB(45,1,1)=HRRB(45,1,1)+Beta*VRR(45,1,0)
      HRRA(46,1,1)=HRRA(46,1,1)+Alpha*VRR(46,1,0)
      HRRB(46,1,1)=HRRB(46,1,1)+Beta*VRR(46,1,0)
      HRRA(47,1,1)=HRRA(47,1,1)+Alpha*VRR(47,1,0)
      HRRB(47,1,1)=HRRB(47,1,1)+Beta*VRR(47,1,0)
      HRRA(48,1,1)=HRRA(48,1,1)+Alpha*VRR(48,1,0)
      HRRB(48,1,1)=HRRB(48,1,1)+Beta*VRR(48,1,0)
      HRRA(49,1,1)=HRRA(49,1,1)+Alpha*VRR(49,1,0)
      HRRB(49,1,1)=HRRB(49,1,1)+Beta*VRR(49,1,0)
      HRRA(50,1,1)=HRRA(50,1,1)+Alpha*VRR(50,1,0)
      HRRB(50,1,1)=HRRB(50,1,1)+Beta*VRR(50,1,0)
      HRRA(51,1,1)=HRRA(51,1,1)+Alpha*VRR(51,1,0)
      HRRB(51,1,1)=HRRB(51,1,1)+Beta*VRR(51,1,0)
      HRRA(52,1,1)=HRRA(52,1,1)+Alpha*VRR(52,1,0)
      HRRB(52,1,1)=HRRB(52,1,1)+Beta*VRR(52,1,0)
      HRRA(53,1,1)=HRRA(53,1,1)+Alpha*VRR(53,1,0)
      HRRB(53,1,1)=HRRB(53,1,1)+Beta*VRR(53,1,0)
      HRRA(54,1,1)=HRRA(54,1,1)+Alpha*VRR(54,1,0)
      HRRB(54,1,1)=HRRB(54,1,1)+Beta*VRR(54,1,0)
      HRRA(55,1,1)=HRRA(55,1,1)+Alpha*VRR(55,1,0)
      HRRB(55,1,1)=HRRB(55,1,1)+Beta*VRR(55,1,0)
      HRRA(56,1,1)=HRRA(56,1,1)+Alpha*VRR(56,1,0)
      HRRB(56,1,1)=HRRB(56,1,1)+Beta*VRR(56,1,0)
      IJ=1
      DO J=1,3
      IF(IW(J).EQ.1) THEN
        IF(IW(1).EQ.1) THEN
        HRRS(1,1,1,IJ)=HRRS(1,1,1,IJ)+PQJ(J)*VRRS(1,1,0,1)
        HRRS(2,1,1,IJ)=HRRS(2,1,1,IJ)+PQJ(J)*VRRS(2,1,0,1)
        HRRS(3,1,1,IJ)=HRRS(3,1,1,IJ)+PQJ(J)*VRRS(3,1,0,1)
        HRRS(4,1,1,IJ)=HRRS(4,1,1,IJ)+PQJ(J)*VRRS(4,1,0,1)
        HRRS(5,1,1,IJ)=HRRS(5,1,1,IJ)+PQJ(J)*VRRS(5,1,0,1)
        HRRS(6,1,1,IJ)=HRRS(6,1,1,IJ)+PQJ(J)*VRRS(6,1,0,1)
        HRRS(7,1,1,IJ)=HRRS(7,1,1,IJ)+PQJ(J)*VRRS(7,1,0,1)
        HRRS(8,1,1,IJ)=HRRS(8,1,1,IJ)+PQJ(J)*VRRS(8,1,0,1)
        HRRS(9,1,1,IJ)=HRRS(9,1,1,IJ)+PQJ(J)*VRRS(9,1,0,1)
        HRRS(10,1,1,IJ)=HRRS(10,1,1,IJ)+PQJ(J)*VRRS(10,1,0,1)
        HRRS(11,1,1,IJ)=HRRS(11,1,1,IJ)+PQJ(J)*VRRS(11,1,0,1)
        HRRS(12,1,1,IJ)=HRRS(12,1,1,IJ)+PQJ(J)*VRRS(12,1,0,1)
        HRRS(13,1,1,IJ)=HRRS(13,1,1,IJ)+PQJ(J)*VRRS(13,1,0,1)
        HRRS(14,1,1,IJ)=HRRS(14,1,1,IJ)+PQJ(J)*VRRS(14,1,0,1)
        HRRS(15,1,1,IJ)=HRRS(15,1,1,IJ)+PQJ(J)*VRRS(15,1,0,1)
        HRRS(16,1,1,IJ)=HRRS(16,1,1,IJ)+PQJ(J)*VRRS(16,1,0,1)
        HRRS(17,1,1,IJ)=HRRS(17,1,1,IJ)+PQJ(J)*VRRS(17,1,0,1)
        HRRS(18,1,1,IJ)=HRRS(18,1,1,IJ)+PQJ(J)*VRRS(18,1,0,1)
        HRRS(19,1,1,IJ)=HRRS(19,1,1,IJ)+PQJ(J)*VRRS(19,1,0,1)
        HRRS(20,1,1,IJ)=HRRS(20,1,1,IJ)+PQJ(J)*VRRS(20,1,0,1)
        HRRS(21,1,1,IJ)=HRRS(21,1,1,IJ)+PQJ(J)*VRRS(21,1,0,1)
        HRRS(22,1,1,IJ)=HRRS(22,1,1,IJ)+PQJ(J)*VRRS(22,1,0,1)
        HRRS(23,1,1,IJ)=HRRS(23,1,1,IJ)+PQJ(J)*VRRS(23,1,0,1)
        HRRS(24,1,1,IJ)=HRRS(24,1,1,IJ)+PQJ(J)*VRRS(24,1,0,1)
        HRRS(25,1,1,IJ)=HRRS(25,1,1,IJ)+PQJ(J)*VRRS(25,1,0,1)
        HRRS(26,1,1,IJ)=HRRS(26,1,1,IJ)+PQJ(J)*VRRS(26,1,0,1)
        HRRS(27,1,1,IJ)=HRRS(27,1,1,IJ)+PQJ(J)*VRRS(27,1,0,1)
        HRRS(28,1,1,IJ)=HRRS(28,1,1,IJ)+PQJ(J)*VRRS(28,1,0,1)
        HRRS(29,1,1,IJ)=HRRS(29,1,1,IJ)+PQJ(J)*VRRS(29,1,0,1)
        HRRS(30,1,1,IJ)=HRRS(30,1,1,IJ)+PQJ(J)*VRRS(30,1,0,1)
        HRRS(31,1,1,IJ)=HRRS(31,1,1,IJ)+PQJ(J)*VRRS(31,1,0,1)
        HRRS(32,1,1,IJ)=HRRS(32,1,1,IJ)+PQJ(J)*VRRS(32,1,0,1)
        HRRS(33,1,1,IJ)=HRRS(33,1,1,IJ)+PQJ(J)*VRRS(33,1,0,1)
        HRRS(34,1,1,IJ)=HRRS(34,1,1,IJ)+PQJ(J)*VRRS(34,1,0,1)
        HRRS(35,1,1,IJ)=HRRS(35,1,1,IJ)+PQJ(J)*VRRS(35,1,0,1)
        ENDIF
        IJ=IJ+1
        IF(IW(2).EQ.1) THEN
        HRRS(1,1,1,IJ)=HRRS(1,1,1,IJ)+PQJ(J)*VRRS(1,1,0,2)
        HRRS(2,1,1,IJ)=HRRS(2,1,1,IJ)+PQJ(J)*VRRS(2,1,0,2)
        HRRS(3,1,1,IJ)=HRRS(3,1,1,IJ)+PQJ(J)*VRRS(3,1,0,2)
        HRRS(4,1,1,IJ)=HRRS(4,1,1,IJ)+PQJ(J)*VRRS(4,1,0,2)
        HRRS(5,1,1,IJ)=HRRS(5,1,1,IJ)+PQJ(J)*VRRS(5,1,0,2)
        HRRS(6,1,1,IJ)=HRRS(6,1,1,IJ)+PQJ(J)*VRRS(6,1,0,2)
        HRRS(7,1,1,IJ)=HRRS(7,1,1,IJ)+PQJ(J)*VRRS(7,1,0,2)
        HRRS(8,1,1,IJ)=HRRS(8,1,1,IJ)+PQJ(J)*VRRS(8,1,0,2)
        HRRS(9,1,1,IJ)=HRRS(9,1,1,IJ)+PQJ(J)*VRRS(9,1,0,2)
        HRRS(10,1,1,IJ)=HRRS(10,1,1,IJ)+PQJ(J)*VRRS(10,1,0,2)
        HRRS(11,1,1,IJ)=HRRS(11,1,1,IJ)+PQJ(J)*VRRS(11,1,0,2)
        HRRS(12,1,1,IJ)=HRRS(12,1,1,IJ)+PQJ(J)*VRRS(12,1,0,2)
        HRRS(13,1,1,IJ)=HRRS(13,1,1,IJ)+PQJ(J)*VRRS(13,1,0,2)
        HRRS(14,1,1,IJ)=HRRS(14,1,1,IJ)+PQJ(J)*VRRS(14,1,0,2)
        HRRS(15,1,1,IJ)=HRRS(15,1,1,IJ)+PQJ(J)*VRRS(15,1,0,2)
        HRRS(16,1,1,IJ)=HRRS(16,1,1,IJ)+PQJ(J)*VRRS(16,1,0,2)
        HRRS(17,1,1,IJ)=HRRS(17,1,1,IJ)+PQJ(J)*VRRS(17,1,0,2)
        HRRS(18,1,1,IJ)=HRRS(18,1,1,IJ)+PQJ(J)*VRRS(18,1,0,2)
        HRRS(19,1,1,IJ)=HRRS(19,1,1,IJ)+PQJ(J)*VRRS(19,1,0,2)
        HRRS(20,1,1,IJ)=HRRS(20,1,1,IJ)+PQJ(J)*VRRS(20,1,0,2)
        HRRS(21,1,1,IJ)=HRRS(21,1,1,IJ)+PQJ(J)*VRRS(21,1,0,2)
        HRRS(22,1,1,IJ)=HRRS(22,1,1,IJ)+PQJ(J)*VRRS(22,1,0,2)
        HRRS(23,1,1,IJ)=HRRS(23,1,1,IJ)+PQJ(J)*VRRS(23,1,0,2)
        HRRS(24,1,1,IJ)=HRRS(24,1,1,IJ)+PQJ(J)*VRRS(24,1,0,2)
        HRRS(25,1,1,IJ)=HRRS(25,1,1,IJ)+PQJ(J)*VRRS(25,1,0,2)
        HRRS(26,1,1,IJ)=HRRS(26,1,1,IJ)+PQJ(J)*VRRS(26,1,0,2)
        HRRS(27,1,1,IJ)=HRRS(27,1,1,IJ)+PQJ(J)*VRRS(27,1,0,2)
        HRRS(28,1,1,IJ)=HRRS(28,1,1,IJ)+PQJ(J)*VRRS(28,1,0,2)
        HRRS(29,1,1,IJ)=HRRS(29,1,1,IJ)+PQJ(J)*VRRS(29,1,0,2)
        HRRS(30,1,1,IJ)=HRRS(30,1,1,IJ)+PQJ(J)*VRRS(30,1,0,2)
        HRRS(31,1,1,IJ)=HRRS(31,1,1,IJ)+PQJ(J)*VRRS(31,1,0,2)
        HRRS(32,1,1,IJ)=HRRS(32,1,1,IJ)+PQJ(J)*VRRS(32,1,0,2)
        HRRS(33,1,1,IJ)=HRRS(33,1,1,IJ)+PQJ(J)*VRRS(33,1,0,2)
        HRRS(34,1,1,IJ)=HRRS(34,1,1,IJ)+PQJ(J)*VRRS(34,1,0,2)
        HRRS(35,1,1,IJ)=HRRS(35,1,1,IJ)+PQJ(J)*VRRS(35,1,0,2)
        ENDIF
        IJ=IJ+1
        IF(IW(3).EQ.1) THEN
        HRRS(1,1,1,IJ)=HRRS(1,1,1,IJ)+PQJ(J)*VRRS(1,1,0,3)
        HRRS(2,1,1,IJ)=HRRS(2,1,1,IJ)+PQJ(J)*VRRS(2,1,0,3)
        HRRS(3,1,1,IJ)=HRRS(3,1,1,IJ)+PQJ(J)*VRRS(3,1,0,3)
        HRRS(4,1,1,IJ)=HRRS(4,1,1,IJ)+PQJ(J)*VRRS(4,1,0,3)
        HRRS(5,1,1,IJ)=HRRS(5,1,1,IJ)+PQJ(J)*VRRS(5,1,0,3)
        HRRS(6,1,1,IJ)=HRRS(6,1,1,IJ)+PQJ(J)*VRRS(6,1,0,3)
        HRRS(7,1,1,IJ)=HRRS(7,1,1,IJ)+PQJ(J)*VRRS(7,1,0,3)
        HRRS(8,1,1,IJ)=HRRS(8,1,1,IJ)+PQJ(J)*VRRS(8,1,0,3)
        HRRS(9,1,1,IJ)=HRRS(9,1,1,IJ)+PQJ(J)*VRRS(9,1,0,3)
        HRRS(10,1,1,IJ)=HRRS(10,1,1,IJ)+PQJ(J)*VRRS(10,1,0,3)
        HRRS(11,1,1,IJ)=HRRS(11,1,1,IJ)+PQJ(J)*VRRS(11,1,0,3)
        HRRS(12,1,1,IJ)=HRRS(12,1,1,IJ)+PQJ(J)*VRRS(12,1,0,3)
        HRRS(13,1,1,IJ)=HRRS(13,1,1,IJ)+PQJ(J)*VRRS(13,1,0,3)
        HRRS(14,1,1,IJ)=HRRS(14,1,1,IJ)+PQJ(J)*VRRS(14,1,0,3)
        HRRS(15,1,1,IJ)=HRRS(15,1,1,IJ)+PQJ(J)*VRRS(15,1,0,3)
        HRRS(16,1,1,IJ)=HRRS(16,1,1,IJ)+PQJ(J)*VRRS(16,1,0,3)
        HRRS(17,1,1,IJ)=HRRS(17,1,1,IJ)+PQJ(J)*VRRS(17,1,0,3)
        HRRS(18,1,1,IJ)=HRRS(18,1,1,IJ)+PQJ(J)*VRRS(18,1,0,3)
        HRRS(19,1,1,IJ)=HRRS(19,1,1,IJ)+PQJ(J)*VRRS(19,1,0,3)
        HRRS(20,1,1,IJ)=HRRS(20,1,1,IJ)+PQJ(J)*VRRS(20,1,0,3)
        HRRS(21,1,1,IJ)=HRRS(21,1,1,IJ)+PQJ(J)*VRRS(21,1,0,3)
        HRRS(22,1,1,IJ)=HRRS(22,1,1,IJ)+PQJ(J)*VRRS(22,1,0,3)
        HRRS(23,1,1,IJ)=HRRS(23,1,1,IJ)+PQJ(J)*VRRS(23,1,0,3)
        HRRS(24,1,1,IJ)=HRRS(24,1,1,IJ)+PQJ(J)*VRRS(24,1,0,3)
        HRRS(25,1,1,IJ)=HRRS(25,1,1,IJ)+PQJ(J)*VRRS(25,1,0,3)
        HRRS(26,1,1,IJ)=HRRS(26,1,1,IJ)+PQJ(J)*VRRS(26,1,0,3)
        HRRS(27,1,1,IJ)=HRRS(27,1,1,IJ)+PQJ(J)*VRRS(27,1,0,3)
        HRRS(28,1,1,IJ)=HRRS(28,1,1,IJ)+PQJ(J)*VRRS(28,1,0,3)
        HRRS(29,1,1,IJ)=HRRS(29,1,1,IJ)+PQJ(J)*VRRS(29,1,0,3)
        HRRS(30,1,1,IJ)=HRRS(30,1,1,IJ)+PQJ(J)*VRRS(30,1,0,3)
        HRRS(31,1,1,IJ)=HRRS(31,1,1,IJ)+PQJ(J)*VRRS(31,1,0,3)
        HRRS(32,1,1,IJ)=HRRS(32,1,1,IJ)+PQJ(J)*VRRS(32,1,0,3)
        HRRS(33,1,1,IJ)=HRRS(33,1,1,IJ)+PQJ(J)*VRRS(33,1,0,3)
        HRRS(34,1,1,IJ)=HRRS(34,1,1,IJ)+PQJ(J)*VRRS(34,1,0,3)
        HRRS(35,1,1,IJ)=HRRS(35,1,1,IJ)+PQJ(J)*VRRS(35,1,0,3)
        ENDIF
        IJ=IJ+1
      ELSE
        IJ=IJ+3
      ENDIF
      ENDDO !J
    END SUBROUTINE CNTRCTG6611
