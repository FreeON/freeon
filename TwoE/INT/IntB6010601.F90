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
! COMPUTES THE INTEGRAL CLASS (d s|d s)
! ----------------------------------------------------------
   SUBROUTINE IntB6010601(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, &
                              OA,LDA,OB,LDB,OC,LDC,OD,LDD,PBC,INTGRL)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF4
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER        :: LBra,LKet,CDOffSet
      REAL(DOUBLE)   :: PrmBufB(10,LBra),PrmBufK(10,LKet)
      TYPE(SmallAtomInfo) :: ACInfo,BDInfo
      TYPE(PBCInfo) :: PBC
      REAL(DOUBLE)  :: INTGRL(*)
      REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz
      REAL(DOUBLE)  :: Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz
      REAL(DOUBLE)  :: PQx,PQy,PQz,FPQx,FPQy,FPQz
      REAL(DOUBLE)  :: Zeta,Eta,Omega,Up,Uq,Upq
      REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT
      REAL(DOUBLE)  :: VRR(10,10,0:4)
      REAL(DOUBLE)  :: HRR(10,10,1)
      INTEGER       :: OffSet,OA,LDA,OB,LDB,OC,LDC,OD,LDD,I,J,K,L
      EXTERNAL InitDbl
      CALL InitDbl(10*10,HRR(1,1,1))
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
         QCx=Qx-Cx
         QCy=Qy-Cy
         QCz=Qz-Cz
         DO K=1,LBra ! K^2 VRR (M0| loop
            Zeta=PrmBufB(1,K)
            Px=PrmBufB(2,K)
            Py=PrmBufB(3,K)
            Pz=PrmBufB(4,K)
            Up=PrmBufB(5,K)
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
            ! Generating (d0|p0)^(1)
            CALL VRRd0p0(10,10,VRR(1,1,1),VRR(1,1,2))
            ! Generating (d0|p0)^(0)
            CALL VRRd0p0(10,10,VRR(1,1,0),VRR(1,1,1))
            ! Generating (s0|d0)^(2)
            CALL VRRs0d0(10,10,VRR(1,1,2),VRR(1,1,3))
            ! Generating (s0|d0)^(1)
            CALL VRRs0d0(10,10,VRR(1,1,1),VRR(1,1,2))
            ! Generating (s0|d0)^(0)
            CALL VRRs0d0(10,10,VRR(1,1,0),VRR(1,1,1))
            ! Generating (p0|d0)^(1)
            CALL VRRp0d0(10,10,VRR(1,1,1),VRR(1,1,2))
            ! Generating (p0|d0)^(0)
            CALL VRRp0d0(10,10,VRR(1,1,0),VRR(1,1,1))
            ! Generating (d0|d0)^(0)
            CALL VRRd0d0(10,10,VRR(1,1,0),VRR(1,1,1))
            ! Contracting ...
            CALL CNTRCT6161(VRR,HRR)
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! No need to generate (d,0|d,s)^(0)
      ! Generating (d,s|d,s)^(0)
      DO L=1,1
         DO K=5,10
            CDOffSet=(OC+K-5)*LDC+(OD+L-1)*LDD
            CALL BraHRR61(OA,OB,LDA,LDB,CDOffSet,HRR(1,K,L),INTGRL)
          ENDDO
      ENDDO
    END SUBROUTINE IntB6010601
    SUBROUTINE CNTRCT6161(VRR,HRR)
      USE DerivedTypes
      USE VScratchB
      REAL(DOUBLE)  :: VRR(10,10,0:4)
      REAL(DOUBLE)  :: HRR(10,10,1)
      HRR(1,1,1)=HRR(1,1,1)+VRR(1,1,0)
      HRR(1,2,1)=HRR(1,2,1)+VRR(1,2,0)
      HRR(1,3,1)=HRR(1,3,1)+VRR(1,3,0)
      HRR(1,4,1)=HRR(1,4,1)+VRR(1,4,0)
      HRR(1,5,1)=HRR(1,5,1)+VRR(1,5,0)
      HRR(1,6,1)=HRR(1,6,1)+VRR(1,6,0)
      HRR(1,7,1)=HRR(1,7,1)+VRR(1,7,0)
      HRR(1,8,1)=HRR(1,8,1)+VRR(1,8,0)
      HRR(1,9,1)=HRR(1,9,1)+VRR(1,9,0)
      HRR(1,10,1)=HRR(1,10,1)+VRR(1,10,0)
      HRR(2,1,1)=HRR(2,1,1)+VRR(2,1,0)
      HRR(2,2,1)=HRR(2,2,1)+VRR(2,2,0)
      HRR(2,3,1)=HRR(2,3,1)+VRR(2,3,0)
      HRR(2,4,1)=HRR(2,4,1)+VRR(2,4,0)
      HRR(2,5,1)=HRR(2,5,1)+VRR(2,5,0)
      HRR(2,6,1)=HRR(2,6,1)+VRR(2,6,0)
      HRR(2,7,1)=HRR(2,7,1)+VRR(2,7,0)
      HRR(2,8,1)=HRR(2,8,1)+VRR(2,8,0)
      HRR(2,9,1)=HRR(2,9,1)+VRR(2,9,0)
      HRR(2,10,1)=HRR(2,10,1)+VRR(2,10,0)
      HRR(3,1,1)=HRR(3,1,1)+VRR(3,1,0)
      HRR(3,2,1)=HRR(3,2,1)+VRR(3,2,0)
      HRR(3,3,1)=HRR(3,3,1)+VRR(3,3,0)
      HRR(3,4,1)=HRR(3,4,1)+VRR(3,4,0)
      HRR(3,5,1)=HRR(3,5,1)+VRR(3,5,0)
      HRR(3,6,1)=HRR(3,6,1)+VRR(3,6,0)
      HRR(3,7,1)=HRR(3,7,1)+VRR(3,7,0)
      HRR(3,8,1)=HRR(3,8,1)+VRR(3,8,0)
      HRR(3,9,1)=HRR(3,9,1)+VRR(3,9,0)
      HRR(3,10,1)=HRR(3,10,1)+VRR(3,10,0)
      HRR(4,1,1)=HRR(4,1,1)+VRR(4,1,0)
      HRR(4,2,1)=HRR(4,2,1)+VRR(4,2,0)
      HRR(4,3,1)=HRR(4,3,1)+VRR(4,3,0)
      HRR(4,4,1)=HRR(4,4,1)+VRR(4,4,0)
      HRR(4,5,1)=HRR(4,5,1)+VRR(4,5,0)
      HRR(4,6,1)=HRR(4,6,1)+VRR(4,6,0)
      HRR(4,7,1)=HRR(4,7,1)+VRR(4,7,0)
      HRR(4,8,1)=HRR(4,8,1)+VRR(4,8,0)
      HRR(4,9,1)=HRR(4,9,1)+VRR(4,9,0)
      HRR(4,10,1)=HRR(4,10,1)+VRR(4,10,0)
      HRR(5,1,1)=HRR(5,1,1)+VRR(5,1,0)
      HRR(5,2,1)=HRR(5,2,1)+VRR(5,2,0)
      HRR(5,3,1)=HRR(5,3,1)+VRR(5,3,0)
      HRR(5,4,1)=HRR(5,4,1)+VRR(5,4,0)
      HRR(5,5,1)=HRR(5,5,1)+VRR(5,5,0)
      HRR(5,6,1)=HRR(5,6,1)+VRR(5,6,0)
      HRR(5,7,1)=HRR(5,7,1)+VRR(5,7,0)
      HRR(5,8,1)=HRR(5,8,1)+VRR(5,8,0)
      HRR(5,9,1)=HRR(5,9,1)+VRR(5,9,0)
      HRR(5,10,1)=HRR(5,10,1)+VRR(5,10,0)
      HRR(6,1,1)=HRR(6,1,1)+VRR(6,1,0)
      HRR(6,2,1)=HRR(6,2,1)+VRR(6,2,0)
      HRR(6,3,1)=HRR(6,3,1)+VRR(6,3,0)
      HRR(6,4,1)=HRR(6,4,1)+VRR(6,4,0)
      HRR(6,5,1)=HRR(6,5,1)+VRR(6,5,0)
      HRR(6,6,1)=HRR(6,6,1)+VRR(6,6,0)
      HRR(6,7,1)=HRR(6,7,1)+VRR(6,7,0)
      HRR(6,8,1)=HRR(6,8,1)+VRR(6,8,0)
      HRR(6,9,1)=HRR(6,9,1)+VRR(6,9,0)
      HRR(6,10,1)=HRR(6,10,1)+VRR(6,10,0)
      HRR(7,1,1)=HRR(7,1,1)+VRR(7,1,0)
      HRR(7,2,1)=HRR(7,2,1)+VRR(7,2,0)
      HRR(7,3,1)=HRR(7,3,1)+VRR(7,3,0)
      HRR(7,4,1)=HRR(7,4,1)+VRR(7,4,0)
      HRR(7,5,1)=HRR(7,5,1)+VRR(7,5,0)
      HRR(7,6,1)=HRR(7,6,1)+VRR(7,6,0)
      HRR(7,7,1)=HRR(7,7,1)+VRR(7,7,0)
      HRR(7,8,1)=HRR(7,8,1)+VRR(7,8,0)
      HRR(7,9,1)=HRR(7,9,1)+VRR(7,9,0)
      HRR(7,10,1)=HRR(7,10,1)+VRR(7,10,0)
      HRR(8,1,1)=HRR(8,1,1)+VRR(8,1,0)
      HRR(8,2,1)=HRR(8,2,1)+VRR(8,2,0)
      HRR(8,3,1)=HRR(8,3,1)+VRR(8,3,0)
      HRR(8,4,1)=HRR(8,4,1)+VRR(8,4,0)
      HRR(8,5,1)=HRR(8,5,1)+VRR(8,5,0)
      HRR(8,6,1)=HRR(8,6,1)+VRR(8,6,0)
      HRR(8,7,1)=HRR(8,7,1)+VRR(8,7,0)
      HRR(8,8,1)=HRR(8,8,1)+VRR(8,8,0)
      HRR(8,9,1)=HRR(8,9,1)+VRR(8,9,0)
      HRR(8,10,1)=HRR(8,10,1)+VRR(8,10,0)
      HRR(9,1,1)=HRR(9,1,1)+VRR(9,1,0)
      HRR(9,2,1)=HRR(9,2,1)+VRR(9,2,0)
      HRR(9,3,1)=HRR(9,3,1)+VRR(9,3,0)
      HRR(9,4,1)=HRR(9,4,1)+VRR(9,4,0)
      HRR(9,5,1)=HRR(9,5,1)+VRR(9,5,0)
      HRR(9,6,1)=HRR(9,6,1)+VRR(9,6,0)
      HRR(9,7,1)=HRR(9,7,1)+VRR(9,7,0)
      HRR(9,8,1)=HRR(9,8,1)+VRR(9,8,0)
      HRR(9,9,1)=HRR(9,9,1)+VRR(9,9,0)
      HRR(9,10,1)=HRR(9,10,1)+VRR(9,10,0)
      HRR(10,1,1)=HRR(10,1,1)+VRR(10,1,0)
      HRR(10,2,1)=HRR(10,2,1)+VRR(10,2,0)
      HRR(10,3,1)=HRR(10,3,1)+VRR(10,3,0)
      HRR(10,4,1)=HRR(10,4,1)+VRR(10,4,0)
      HRR(10,5,1)=HRR(10,5,1)+VRR(10,5,0)
      HRR(10,6,1)=HRR(10,6,1)+VRR(10,6,0)
      HRR(10,7,1)=HRR(10,7,1)+VRR(10,7,0)
      HRR(10,8,1)=HRR(10,8,1)+VRR(10,8,0)
      HRR(10,9,1)=HRR(10,9,1)+VRR(10,9,0)
      HRR(10,10,1)=HRR(10,10,1)+VRR(10,10,0)
    END SUBROUTINE CNTRCT6161
