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
! COMPUTES THE INTEGRAL CLASS (f sp|s s) 
! ---------------------------------------------------------- 
   SUBROUTINE IntB10020101(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, & 
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
      REAL(DOUBLE)  :: VRR(35,1,0:4)
      REAL(DOUBLE)  :: HRR(45,1,1)
      INTEGER       :: OffSet,OA,LDA,OB,LDB,OC,LDC,OD,LDD,I,J,K,L
      EXTERNAL InitDbl
      CALL InitDbl(45*1,HRR(1,1,1))
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
            FnSpB=PrmBufB(6,K)
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
            ! Generating (f0|s0)^(1)
            CALL VRRf0s0(35,1,VRR(1,1,1),VRR(1,1,2))
            ! Generating (f0|s0)^(0)
            CALL VRRf0s0(35,1,VRR(1,1,0),VRR(1,1,1))
            ! Generating (g0|s0)^(0)
            CALL VRRg0s0(35,1,VRR(1,1,0),VRR(1,1,1))
            ! Contracting ... 
            CALL CNTRCT10211(VRR,HRR)
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! No need to generate (f,0|s,s)^(0) 
      ! Generating (f,sp|s,s)^(0) 
      DO L=1,1
         DO K=1,1
            CDOffSet=(OC+K-1)*LDC+(OD+L-1)*LDD 
            CALL BraHRR102(OA,OB,LDA,LDB,CDOffSet,HRR(1,K,L),INTGRL) 
          ENDDO 
      ENDDO 
    END SUBROUTINE IntB10020101
    SUBROUTINE CNTRCT10211(VRR,HRR)
      USE DerivedTypes
      USE VScratchB
      REAL(DOUBLE)  :: VRR(35,1,0:4)
      REAL(DOUBLE)  :: HRR(45,1,1)
      HRR(1,1,1)=HRR(1,1,1)+VRR(1,1,0)
      HRR(2,1,1)=HRR(2,1,1)+VRR(2,1,0)
      HRR(3,1,1)=HRR(3,1,1)+VRR(3,1,0)
      HRR(4,1,1)=HRR(4,1,1)+VRR(4,1,0)
      HRR(5,1,1)=HRR(5,1,1)+VRR(5,1,0)
      HRR(6,1,1)=HRR(6,1,1)+VRR(6,1,0)
      HRR(7,1,1)=HRR(7,1,1)+VRR(7,1,0)
      HRR(8,1,1)=HRR(8,1,1)+VRR(8,1,0)
      HRR(9,1,1)=HRR(9,1,1)+VRR(9,1,0)
      HRR(10,1,1)=HRR(10,1,1)+VRR(10,1,0)
      HRR(11,1,1)=HRR(11,1,1)+VRR(11,1,0)
      HRR(36,1,1)=HRR(36,1,1)+FnSpB*VRR(11,1,0)
      HRR(12,1,1)=HRR(12,1,1)+VRR(12,1,0)
      HRR(37,1,1)=HRR(37,1,1)+FnSpB*VRR(12,1,0)
      HRR(13,1,1)=HRR(13,1,1)+VRR(13,1,0)
      HRR(38,1,1)=HRR(38,1,1)+FnSpB*VRR(13,1,0)
      HRR(14,1,1)=HRR(14,1,1)+VRR(14,1,0)
      HRR(39,1,1)=HRR(39,1,1)+FnSpB*VRR(14,1,0)
      HRR(15,1,1)=HRR(15,1,1)+VRR(15,1,0)
      HRR(40,1,1)=HRR(40,1,1)+FnSpB*VRR(15,1,0)
      HRR(16,1,1)=HRR(16,1,1)+VRR(16,1,0)
      HRR(41,1,1)=HRR(41,1,1)+FnSpB*VRR(16,1,0)
      HRR(17,1,1)=HRR(17,1,1)+VRR(17,1,0)
      HRR(42,1,1)=HRR(42,1,1)+FnSpB*VRR(17,1,0)
      HRR(18,1,1)=HRR(18,1,1)+VRR(18,1,0)
      HRR(43,1,1)=HRR(43,1,1)+FnSpB*VRR(18,1,0)
      HRR(19,1,1)=HRR(19,1,1)+VRR(19,1,0)
      HRR(44,1,1)=HRR(44,1,1)+FnSpB*VRR(19,1,0)
      HRR(20,1,1)=HRR(20,1,1)+VRR(20,1,0)
      HRR(45,1,1)=HRR(45,1,1)+FnSpB*VRR(20,1,0)
      HRR(21,1,1)=HRR(21,1,1)+VRR(21,1,0)
      HRR(22,1,1)=HRR(22,1,1)+VRR(22,1,0)
      HRR(23,1,1)=HRR(23,1,1)+VRR(23,1,0)
      HRR(24,1,1)=HRR(24,1,1)+VRR(24,1,0)
      HRR(25,1,1)=HRR(25,1,1)+VRR(25,1,0)
      HRR(26,1,1)=HRR(26,1,1)+VRR(26,1,0)
      HRR(27,1,1)=HRR(27,1,1)+VRR(27,1,0)
      HRR(28,1,1)=HRR(28,1,1)+VRR(28,1,0)
      HRR(29,1,1)=HRR(29,1,1)+VRR(29,1,0)
      HRR(30,1,1)=HRR(30,1,1)+VRR(30,1,0)
      HRR(31,1,1)=HRR(31,1,1)+VRR(31,1,0)
      HRR(32,1,1)=HRR(32,1,1)+VRR(32,1,0)
      HRR(33,1,1)=HRR(33,1,1)+VRR(33,1,0)
      HRR(34,1,1)=HRR(34,1,1)+VRR(34,1,0)
      HRR(35,1,1)=HRR(35,1,1)+VRR(35,1,0)
    END SUBROUTINE CNTRCT10211
