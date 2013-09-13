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
   SUBROUTINE VRRg0s0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2Z*VRR0(6,1)
      V(2)=2.D0*V(1)
      V(3)=ExZpE*r1x2Z*VRR1(6,1)
      V(4)=-2.D0*V(3)
      V(5)=r1x2Z*VRR0(7,1)
      V(6)=ExZpE*r1x2Z*VRR1(7,1)
      V(7)=r1x2Z*VRR0(8,1)
      V(8)=2.D0*V(7)
      V(9)=ExZpE*r1x2Z*VRR1(8,1)
      V(10)=-2.D0*V(9)
      V(11)=r1x2Z*VRR0(9,1)
      V(12)=ExZpE*r1x2Z*VRR1(9,1)
      V(13)=2.D0*V(11)
      V(14)=-2.D0*V(12)
      V(15)=r1x2Z*VRR0(10,1)
      V(16)=ExZpE*r1x2Z*VRR1(10,1)
      V(17)=-V(16)
      VRR0(21,1)=3.D0*r1x2Z*VRR0(5,1)+PAx*VRR0(11,1)-3.D0*ExZpE*r1x2Z*VRR1(5,1)+WPx*VRR1(11,1)
      VRR0(22,1)=V(2)+V(4)+PAx*VRR0(12,1)+WPx*VRR1(12,1)
      VRR0(23,1)=V(5)-V(6)+PAx*VRR0(13,1)+WPx*VRR1(13,1)
      VRR0(24,1)=V(2)+V(4)+PAy*VRR0(13,1)+WPy*VRR1(13,1)
      VRR0(25,1)=3.D0*V(5)-3.D0*V(6)+PAy*VRR0(14,1)+WPy*VRR1(14,1)
      VRR0(26,1)=V(8)+V(10)+PAx*VRR0(15,1)+WPx*VRR1(15,1)
      VRR0(27,1)=V(11)-V(12)+PAx*VRR0(16,1)+WPx*VRR1(16,1)
      VRR0(28,1)=V(7)-V(9)+PAy*VRR0(16,1)+WPy*VRR1(16,1)
      VRR0(29,1)=V(13)+V(14)+PAy*VRR0(17,1)+WPy*VRR1(17,1)
      VRR0(30,1)=V(15)+V(17)+PAx*VRR0(18,1)+WPx*VRR1(18,1)
      VRR0(31,1)=V(1)-V(3)+PAz*VRR0(16,1)+WPz*VRR1(16,1)
      VRR0(32,1)=V(15)+V(17)+PAy*VRR0(19,1)+WPy*VRR1(19,1)
      VRR0(33,1)=V(8)+V(10)+PAz*VRR0(18,1)+WPz*VRR1(18,1)
      VRR0(34,1)=V(13)+V(14)+PAz*VRR0(19,1)+WPz*VRR1(19,1)
      VRR0(35,1)=3.D0*V(15)-3.D0*V(16)+PAz*VRR0(20,1)+WPz*VRR1(20,1)
END SUBROUTINE VRRg0s0
SUBROUTINE MVRRg0s0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(21,1)=PAx*VS0(11,1)+WPx*VS1(11,1)+r1x2Z*VR1(11,1)&
   +3D0*r1x2Z*(VS0(5,1)-ExZpE*VS1(5,1))
VS0(22,1)=PAx*VS0(12,1)+WPx*VS1(12,1)+r1x2Z*VR1(12,1)&
   +2D0*r1x2Z*(VS0(6,1)-ExZpE*VS1(6,1))
VS0(23,1)=PAx*VS0(13,1)+WPx*VS1(13,1)+r1x2Z*VR1(13,1)&
   +r1x2Z*(VS0(7,1)-ExZpE*VS1(7,1))
VS0(24,1)=PAy*VS0(13,1)+WPy*VS1(13,1)&
   +2D0*r1x2Z*(VS0(6,1)-ExZpE*VS1(6,1))
VS0(25,1)=PAy*VS0(14,1)+WPy*VS1(14,1)&
   +3D0*r1x2Z*(VS0(7,1)-ExZpE*VS1(7,1))
VS0(26,1)=PAx*VS0(15,1)+WPx*VS1(15,1)+r1x2Z*VR1(15,1)&
   +2D0*r1x2Z*(VS0(8,1)-ExZpE*VS1(8,1))
VS0(27,1)=PAx*VS0(16,1)+WPx*VS1(16,1)+r1x2Z*VR1(16,1)&
   +r1x2Z*(VS0(9,1)-ExZpE*VS1(9,1))
VS0(28,1)=PAy*VS0(16,1)+WPy*VS1(16,1)&
   +r1x2Z*(VS0(8,1)-ExZpE*VS1(8,1))
VS0(29,1)=PAy*VS0(17,1)+WPy*VS1(17,1)&
   +2D0*r1x2Z*(VS0(9,1)-ExZpE*VS1(9,1))
VS0(30,1)=PAx*VS0(18,1)+WPx*VS1(18,1)+r1x2Z*VR1(18,1)&
   +r1x2Z*(VS0(10,1)-ExZpE*VS1(10,1))
VS0(31,1)=PAz*VS0(16,1)+WPz*VS1(16,1)&
   +r1x2Z*(VS0(6,1)-ExZpE*VS1(6,1))
VS0(32,1)=PAy*VS0(19,1)+WPy*VS1(19,1)&
   +r1x2Z*(VS0(10,1)-ExZpE*VS1(10,1))
VS0(33,1)=PAz*VS0(18,1)+WPz*VS1(18,1)&
   +2D0*r1x2Z*(VS0(8,1)-ExZpE*VS1(8,1))
VS0(34,1)=PAz*VS0(19,1)+WPz*VS1(19,1)&
   +2D0*r1x2Z*(VS0(9,1)-ExZpE*VS1(9,1))
VS0(35,1)=PAz*VS0(20,1)+WPz*VS1(20,1)&
   +3D0*r1x2Z*(VS0(10,1)-ExZpE*VS1(10,1))
CASE(2)
VS0(21,1)=PAx*VS0(11,1)+WPx*VS1(11,1)&
   +3D0*r1x2Z*(VS0(5,1)-ExZpE*VS1(5,1))
VS0(22,1)=PAx*VS0(12,1)+WPx*VS1(12,1)&
   +2D0*r1x2Z*(VS0(6,1)-ExZpE*VS1(6,1))
VS0(23,1)=PAx*VS0(13,1)+WPx*VS1(13,1)&
   +r1x2Z*(VS0(7,1)-ExZpE*VS1(7,1))
VS0(24,1)=PAy*VS0(13,1)+WPy*VS1(13,1)+r1x2Z*VR1(13,1)&
   +2D0*r1x2Z*(VS0(6,1)-ExZpE*VS1(6,1))
VS0(25,1)=PAy*VS0(14,1)+WPy*VS1(14,1)+r1x2Z*VR1(14,1)&
   +3D0*r1x2Z*(VS0(7,1)-ExZpE*VS1(7,1))
VS0(26,1)=PAx*VS0(15,1)+WPx*VS1(15,1)&
   +2D0*r1x2Z*(VS0(8,1)-ExZpE*VS1(8,1))
VS0(27,1)=PAx*VS0(16,1)+WPx*VS1(16,1)&
   +r1x2Z*(VS0(9,1)-ExZpE*VS1(9,1))
VS0(28,1)=PAy*VS0(16,1)+WPy*VS1(16,1)+r1x2Z*VR1(16,1)&
   +r1x2Z*(VS0(8,1)-ExZpE*VS1(8,1))
VS0(29,1)=PAy*VS0(17,1)+WPy*VS1(17,1)+r1x2Z*VR1(17,1)&
   +2D0*r1x2Z*(VS0(9,1)-ExZpE*VS1(9,1))
VS0(30,1)=PAx*VS0(18,1)+WPx*VS1(18,1)&
   +r1x2Z*(VS0(10,1)-ExZpE*VS1(10,1))
VS0(31,1)=PAz*VS0(16,1)+WPz*VS1(16,1)&
   +r1x2Z*(VS0(6,1)-ExZpE*VS1(6,1))
VS0(32,1)=PAy*VS0(19,1)+WPy*VS1(19,1)+r1x2Z*VR1(19,1)&
   +r1x2Z*(VS0(10,1)-ExZpE*VS1(10,1))
VS0(33,1)=PAz*VS0(18,1)+WPz*VS1(18,1)&
   +2D0*r1x2Z*(VS0(8,1)-ExZpE*VS1(8,1))
VS0(34,1)=PAz*VS0(19,1)+WPz*VS1(19,1)&
   +2D0*r1x2Z*(VS0(9,1)-ExZpE*VS1(9,1))
VS0(35,1)=PAz*VS0(20,1)+WPz*VS1(20,1)&
   +3D0*r1x2Z*(VS0(10,1)-ExZpE*VS1(10,1))
CASE(3)
VS0(21,1)=PAx*VS0(11,1)+WPx*VS1(11,1)&
   +3D0*r1x2Z*(VS0(5,1)-ExZpE*VS1(5,1))
VS0(22,1)=PAx*VS0(12,1)+WPx*VS1(12,1)&
   +2D0*r1x2Z*(VS0(6,1)-ExZpE*VS1(6,1))
VS0(23,1)=PAx*VS0(13,1)+WPx*VS1(13,1)&
   +r1x2Z*(VS0(7,1)-ExZpE*VS1(7,1))
VS0(24,1)=PAy*VS0(13,1)+WPy*VS1(13,1)&
   +2D0*r1x2Z*(VS0(6,1)-ExZpE*VS1(6,1))
VS0(25,1)=PAy*VS0(14,1)+WPy*VS1(14,1)&
   +3D0*r1x2Z*(VS0(7,1)-ExZpE*VS1(7,1))
VS0(26,1)=PAx*VS0(15,1)+WPx*VS1(15,1)&
   +2D0*r1x2Z*(VS0(8,1)-ExZpE*VS1(8,1))
VS0(27,1)=PAx*VS0(16,1)+WPx*VS1(16,1)&
   +r1x2Z*(VS0(9,1)-ExZpE*VS1(9,1))
VS0(28,1)=PAy*VS0(16,1)+WPy*VS1(16,1)&
   +r1x2Z*(VS0(8,1)-ExZpE*VS1(8,1))
VS0(29,1)=PAy*VS0(17,1)+WPy*VS1(17,1)&
   +2D0*r1x2Z*(VS0(9,1)-ExZpE*VS1(9,1))
VS0(30,1)=PAx*VS0(18,1)+WPx*VS1(18,1)&
   +r1x2Z*(VS0(10,1)-ExZpE*VS1(10,1))
VS0(31,1)=PAz*VS0(16,1)+WPz*VS1(16,1)+r1x2Z*VR1(16,1)&
   +r1x2Z*(VS0(6,1)-ExZpE*VS1(6,1))
VS0(32,1)=PAy*VS0(19,1)+WPy*VS1(19,1)&
   +r1x2Z*(VS0(10,1)-ExZpE*VS1(10,1))
VS0(33,1)=PAz*VS0(18,1)+WPz*VS1(18,1)+r1x2Z*VR1(18,1)&
   +2D0*r1x2Z*(VS0(8,1)-ExZpE*VS1(8,1))
VS0(34,1)=PAz*VS0(19,1)+WPz*VS1(19,1)+r1x2Z*VR1(19,1)&
   +2D0*r1x2Z*(VS0(9,1)-ExZpE*VS1(9,1))
VS0(35,1)=PAz*VS0(20,1)+WPz*VS1(20,1)+r1x2Z*VR1(20,1)&
   +3D0*r1x2Z*(VS0(10,1)-ExZpE*VS1(10,1))
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRg0s0'
STOP
END SELECT
END SUBROUTINE MVRRg0s0
