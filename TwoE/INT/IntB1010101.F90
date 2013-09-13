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
! COMPUTES THE INTEGRAL CLASS (s s|s s)
! ----------------------------------------------------------
   SUBROUTINE IntB1010101(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, &
                              OA,LDA,OB,LDB,OC,LDC,OD,LDD,PBC,INTGRL)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF0
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
      REAL(DOUBLE)  :: VRR(1,1,0:0)
      REAL(DOUBLE)  :: HRR(1,1,1)
      INTEGER       :: OffSet,OA,LDA,OB,LDB,OC,LDC,OD,LDD,I,J,K,L
      EXTERNAL InitDbl
      CALL InitDbl(1*1,HRR(1,1,1))
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
              VRR(1,1,0)=Upq*(F0_0(L)+T*(F0_1(L)+T*(F0_2(L)+T*(F0_3(L)+T*F0_4(L)))))
            ELSE
              InvT=One/T
              SqInvT=DSQRT(InvT)
              VRR(1,1,0)=+8.862269254527580D-01*Upq*SqInvT
            ENDIF
            ! Contracting ...
            CALL CNTRCT1111(VRR,HRR)
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! No need to generate (s,0|s,s)^(0)
      ! Generating (s,s|s,s)^(0)
      DO L=1,1
         DO K=1,1
            CDOffSet=(OC+K-1)*LDC+(OD+L-1)*LDD
            CALL BraHRR11(OA,OB,LDA,LDB,CDOffSet,HRR(1,K,L),INTGRL)
          ENDDO
      ENDDO
    END SUBROUTINE IntB1010101
    SUBROUTINE CNTRCT1111(VRR,HRR)
      USE DerivedTypes
      USE VScratchB
      REAL(DOUBLE)  :: VRR(1,1,0:0)
      REAL(DOUBLE)  :: HRR(1,1,1)
      HRR(1,1,1)=HRR(1,1,1)+VRR(1,1,0)
    END SUBROUTINE CNTRCT1111
