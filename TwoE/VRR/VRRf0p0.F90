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
   SUBROUTINE VRRf0p0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2Z*VRR0(2,2)
      V(2)=ExZpE*r1x2Z*VRR1(2,2)
      V(3)=r1x2Z*VRR0(2,3)
      V(4)=ExZpE*r1x2Z*VRR1(2,3)
      V(5)=r1x2Z*VRR0(2,4)
      V(6)=ExZpE*r1x2Z*VRR1(2,4)
      V(7)=r1x2Z*VRR0(3,2)
      V(8)=ExZpE*r1x2Z*VRR1(3,2)
      V(9)=-V(8)
      V(10)=HfxZpE*VRR1(6,1)
      V(11)=r1x2Z*VRR0(3,3)
      V(12)=ExZpE*r1x2Z*VRR1(3,3)
      V(13)=-V(12)
      V(14)=r1x2Z*VRR0(3,4)
      V(15)=ExZpE*r1x2Z*VRR1(3,4)
      V(16)=-V(15)
      V(17)=-V(2)
      V(18)=-V(4)
      V(19)=-V(6)
      V(20)=r1x2Z*VRR0(4,2)
      V(21)=ExZpE*r1x2Z*VRR1(4,2)
      V(22)=-V(21)
      V(23)=HfxZpE*VRR1(8,1)
      V(24)=r1x2Z*VRR0(4,3)
      V(25)=ExZpE*r1x2Z*VRR1(4,3)
      V(26)=-V(25)
      V(27)=r1x2Z*VRR0(4,4)
      V(28)=ExZpE*r1x2Z*VRR1(4,4)
      V(29)=-V(28)
      V(30)=HfxZpE*VRR1(9,1)
      VRR0(11,2)=2.D0*V(1)-2.D0*V(2)+PAx*VRR0(5,2)+HfxZpE*VRR1(5,1)+WPx*VRR1(5,2)
      VRR0(11,3)=2.D0*V(3)-2.D0*V(4)+PAx*VRR0(5,3)+WPx*VRR1(5,3)
      VRR0(11,4)=2.D0*V(5)-2.D0*V(6)+PAx*VRR0(5,4)+WPx*VRR1(5,4)
      VRR0(12,2)=V(7)+V(9)+V(10)+PAx*VRR0(6,2)+WPx*VRR1(6,2)
      VRR0(12,3)=V(11)+V(13)+PAx*VRR0(6,3)+WPx*VRR1(6,3)
      VRR0(12,4)=V(14)+V(16)+PAx*VRR0(6,4)+WPx*VRR1(6,4)
      VRR0(13,2)=V(1)+V(17)+PAy*VRR0(6,2)+WPy*VRR1(6,2)
      VRR0(13,3)=V(3)+V(10)+V(18)+PAy*VRR0(6,3)+WPy*VRR1(6,3)
      VRR0(13,4)=V(5)+V(19)+PAy*VRR0(6,4)+WPy*VRR1(6,4)
      VRR0(14,2)=2.D0*V(7)-2.D0*V(8)+PAy*VRR0(7,2)+WPy*VRR1(7,2)
      VRR0(14,3)=2.D0*V(11)-2.D0*V(12)+PAy*VRR0(7,3)+HfxZpE*VRR1(7,1)+WPy*VRR1(7,3)
      VRR0(14,4)=2.D0*V(14)-2.D0*V(15)+PAy*VRR0(7,4)+WPy*VRR1(7,4)
      VRR0(15,2)=V(20)+V(22)+V(23)+PAx*VRR0(8,2)+WPx*VRR1(8,2)
      VRR0(15,3)=V(24)+V(26)+PAx*VRR0(8,3)+WPx*VRR1(8,3)
      VRR0(15,4)=V(27)+V(29)+PAx*VRR0(8,4)+WPx*VRR1(8,4)
      VRR0(16,2)=V(30)+QCx*VRR0(16,1)+WQx*VRR1(16,1)
      VRR0(16,3)=V(23)+QCy*VRR0(16,1)+WQy*VRR1(16,1)
      VRR0(16,4)=V(10)+QCz*VRR0(16,1)+WQz*VRR1(16,1)
      VRR0(17,2)=V(20)+V(22)+PAy*VRR0(9,2)+WPy*VRR1(9,2)
      VRR0(17,3)=V(24)+V(26)+V(30)+PAy*VRR0(9,3)+WPy*VRR1(9,3)
      VRR0(17,4)=V(27)+V(29)+PAy*VRR0(9,4)+WPy*VRR1(9,4)
      VRR0(18,2)=V(1)+V(17)+PAz*VRR0(8,2)+WPz*VRR1(8,2)
      VRR0(18,3)=V(3)+V(18)+PAz*VRR0(8,3)+WPz*VRR1(8,3)
      VRR0(18,4)=V(5)+V(19)+V(23)+PAz*VRR0(8,4)+WPz*VRR1(8,4)
      VRR0(19,2)=V(7)+V(9)+PAz*VRR0(9,2)+WPz*VRR1(9,2)
      VRR0(19,3)=V(11)+V(13)+PAz*VRR0(9,3)+WPz*VRR1(9,3)
      VRR0(19,4)=V(14)+V(16)+V(30)+PAz*VRR0(9,4)+WPz*VRR1(9,4)
      VRR0(20,2)=2.D0*V(20)-2.D0*V(21)+PAz*VRR0(10,2)+WPz*VRR1(10,2)
      VRR0(20,3)=2.D0*V(24)-2.D0*V(25)+PAz*VRR0(10,3)+WPz*VRR1(10,3)
      VRR0(20,4)=2.D0*V(27)-2.D0*V(28)+PAz*VRR0(10,4)+HfxZpE*VRR1(10,1)+WPz*VRR1(10,4)
END SUBROUTINE VRRf0p0
SUBROUTINE MVRRf0p0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(11,2)=PAx*VS0(5,2)+WPx*VS1(5,2)+r1x2Z*VR1(5,2)&
   +2D0*r1x2Z*(VS0(2,2)-ExZpE*VS1(2,2))&
   +HfxZpE*VS1(5,1)
VS0(11,3)=PAx*VS0(5,3)+WPx*VS1(5,3)+r1x2Z*VR1(5,3)&
   +2D0*r1x2Z*(VS0(2,3)-ExZpE*VS1(2,3))
VS0(11,4)=PAx*VS0(5,4)+WPx*VS1(5,4)+r1x2Z*VR1(5,4)&
   +2D0*r1x2Z*(VS0(2,4)-ExZpE*VS1(2,4))
VS0(12,2)=PAx*VS0(6,2)+WPx*VS1(6,2)+r1x2Z*VR1(6,2)&
   +r1x2Z*(VS0(3,2)-ExZpE*VS1(3,2))&
   +HfxZpE*VS1(6,1)
VS0(12,3)=PAx*VS0(6,3)+WPx*VS1(6,3)+r1x2Z*VR1(6,3)&
   +r1x2Z*(VS0(3,3)-ExZpE*VS1(3,3))
VS0(12,4)=PAx*VS0(6,4)+WPx*VS1(6,4)+r1x2Z*VR1(6,4)&
   +r1x2Z*(VS0(3,4)-ExZpE*VS1(3,4))
VS0(13,2)=PAy*VS0(6,2)+WPy*VS1(6,2)&
   +r1x2Z*(VS0(2,2)-ExZpE*VS1(2,2))
VS0(13,3)=PAy*VS0(6,3)+WPy*VS1(6,3)&
   +r1x2Z*(VS0(2,3)-ExZpE*VS1(2,3))&
   +HfxZpE*VS1(6,1)
VS0(13,4)=PAy*VS0(6,4)+WPy*VS1(6,4)&
   +r1x2Z*(VS0(2,4)-ExZpE*VS1(2,4))
VS0(14,2)=PAy*VS0(7,2)+WPy*VS1(7,2)&
   +2D0*r1x2Z*(VS0(3,2)-ExZpE*VS1(3,2))
VS0(14,3)=PAy*VS0(7,3)+WPy*VS1(7,3)&
   +2D0*r1x2Z*(VS0(3,3)-ExZpE*VS1(3,3))&
   +HfxZpE*VS1(7,1)
VS0(14,4)=PAy*VS0(7,4)+WPy*VS1(7,4)&
   +2D0*r1x2Z*(VS0(3,4)-ExZpE*VS1(3,4))
VS0(15,2)=PAx*VS0(8,2)+WPx*VS1(8,2)+r1x2Z*VR1(8,2)&
   +r1x2Z*(VS0(4,2)-ExZpE*VS1(4,2))&
   +HfxZpE*VS1(8,1)
VS0(15,3)=PAx*VS0(8,3)+WPx*VS1(8,3)+r1x2Z*VR1(8,3)&
   +r1x2Z*(VS0(4,3)-ExZpE*VS1(4,3))
VS0(15,4)=PAx*VS0(8,4)+WPx*VS1(8,4)+r1x2Z*VR1(8,4)&
   +r1x2Z*(VS0(4,4)-ExZpE*VS1(4,4))
VS0(16,2)=PAx*VS0(9,2)+WPx*VS1(9,2)+r1x2Z*VR1(9,2)&
   +HfxZpE*VS1(9,1)
VS0(16,3)=PAx*VS0(9,3)+WPx*VS1(9,3)+r1x2Z*VR1(9,3)
VS0(16,4)=PAx*VS0(9,4)+WPx*VS1(9,4)+r1x2Z*VR1(9,4)
VS0(17,2)=PAy*VS0(9,2)+WPy*VS1(9,2)&
   +r1x2Z*(VS0(4,2)-ExZpE*VS1(4,2))
VS0(17,3)=PAy*VS0(9,3)+WPy*VS1(9,3)&
   +r1x2Z*(VS0(4,3)-ExZpE*VS1(4,3))&
   +HfxZpE*VS1(9,1)
VS0(17,4)=PAy*VS0(9,4)+WPy*VS1(9,4)&
   +r1x2Z*(VS0(4,4)-ExZpE*VS1(4,4))
VS0(18,2)=PAz*VS0(8,2)+WPz*VS1(8,2)&
   +r1x2Z*(VS0(2,2)-ExZpE*VS1(2,2))
VS0(18,3)=PAz*VS0(8,3)+WPz*VS1(8,3)&
   +r1x2Z*(VS0(2,3)-ExZpE*VS1(2,3))
VS0(18,4)=PAz*VS0(8,4)+WPz*VS1(8,4)&
   +r1x2Z*(VS0(2,4)-ExZpE*VS1(2,4))&
   +HfxZpE*VS1(8,1)
VS0(19,2)=PAz*VS0(9,2)+WPz*VS1(9,2)&
   +r1x2Z*(VS0(3,2)-ExZpE*VS1(3,2))
VS0(19,3)=PAz*VS0(9,3)+WPz*VS1(9,3)&
   +r1x2Z*(VS0(3,3)-ExZpE*VS1(3,3))
VS0(19,4)=PAz*VS0(9,4)+WPz*VS1(9,4)&
   +r1x2Z*(VS0(3,4)-ExZpE*VS1(3,4))&
   +HfxZpE*VS1(9,1)
VS0(20,2)=PAz*VS0(10,2)+WPz*VS1(10,2)&
   +2D0*r1x2Z*(VS0(4,2)-ExZpE*VS1(4,2))
VS0(20,3)=PAz*VS0(10,3)+WPz*VS1(10,3)&
   +2D0*r1x2Z*(VS0(4,3)-ExZpE*VS1(4,3))
VS0(20,4)=PAz*VS0(10,4)+WPz*VS1(10,4)&
   +2D0*r1x2Z*(VS0(4,4)-ExZpE*VS1(4,4))&
   +HfxZpE*VS1(10,1)
CASE(2)
VS0(11,2)=PAx*VS0(5,2)+WPx*VS1(5,2)&
   +2D0*r1x2Z*(VS0(2,2)-ExZpE*VS1(2,2))&
   +HfxZpE*VS1(5,1)
VS0(11,3)=PAx*VS0(5,3)+WPx*VS1(5,3)&
   +2D0*r1x2Z*(VS0(2,3)-ExZpE*VS1(2,3))
VS0(11,4)=PAx*VS0(5,4)+WPx*VS1(5,4)&
   +2D0*r1x2Z*(VS0(2,4)-ExZpE*VS1(2,4))
VS0(12,2)=PAx*VS0(6,2)+WPx*VS1(6,2)&
   +r1x2Z*(VS0(3,2)-ExZpE*VS1(3,2))&
   +HfxZpE*VS1(6,1)
VS0(12,3)=PAx*VS0(6,3)+WPx*VS1(6,3)&
   +r1x2Z*(VS0(3,3)-ExZpE*VS1(3,3))
VS0(12,4)=PAx*VS0(6,4)+WPx*VS1(6,4)&
   +r1x2Z*(VS0(3,4)-ExZpE*VS1(3,4))
VS0(13,2)=PAy*VS0(6,2)+WPy*VS1(6,2)+r1x2Z*VR1(6,2)&
   +r1x2Z*(VS0(2,2)-ExZpE*VS1(2,2))
VS0(13,3)=PAy*VS0(6,3)+WPy*VS1(6,3)+r1x2Z*VR1(6,3)&
   +r1x2Z*(VS0(2,3)-ExZpE*VS1(2,3))&
   +HfxZpE*VS1(6,1)
VS0(13,4)=PAy*VS0(6,4)+WPy*VS1(6,4)+r1x2Z*VR1(6,4)&
   +r1x2Z*(VS0(2,4)-ExZpE*VS1(2,4))
VS0(14,2)=PAy*VS0(7,2)+WPy*VS1(7,2)+r1x2Z*VR1(7,2)&
   +2D0*r1x2Z*(VS0(3,2)-ExZpE*VS1(3,2))
VS0(14,3)=PAy*VS0(7,3)+WPy*VS1(7,3)+r1x2Z*VR1(7,3)&
   +2D0*r1x2Z*(VS0(3,3)-ExZpE*VS1(3,3))&
   +HfxZpE*VS1(7,1)
VS0(14,4)=PAy*VS0(7,4)+WPy*VS1(7,4)+r1x2Z*VR1(7,4)&
   +2D0*r1x2Z*(VS0(3,4)-ExZpE*VS1(3,4))
VS0(15,2)=PAx*VS0(8,2)+WPx*VS1(8,2)&
   +r1x2Z*(VS0(4,2)-ExZpE*VS1(4,2))&
   +HfxZpE*VS1(8,1)
VS0(15,3)=PAx*VS0(8,3)+WPx*VS1(8,3)&
   +r1x2Z*(VS0(4,3)-ExZpE*VS1(4,3))
VS0(15,4)=PAx*VS0(8,4)+WPx*VS1(8,4)&
   +r1x2Z*(VS0(4,4)-ExZpE*VS1(4,4))
VS0(16,2)=PAx*VS0(9,2)+WPx*VS1(9,2)&
   +HfxZpE*VS1(9,1)
VS0(16,3)=PAx*VS0(9,3)+WPx*VS1(9,3)
VS0(16,4)=PAx*VS0(9,4)+WPx*VS1(9,4)
VS0(17,2)=PAy*VS0(9,2)+WPy*VS1(9,2)+r1x2Z*VR1(9,2)&
   +r1x2Z*(VS0(4,2)-ExZpE*VS1(4,2))
VS0(17,3)=PAy*VS0(9,3)+WPy*VS1(9,3)+r1x2Z*VR1(9,3)&
   +r1x2Z*(VS0(4,3)-ExZpE*VS1(4,3))&
   +HfxZpE*VS1(9,1)
VS0(17,4)=PAy*VS0(9,4)+WPy*VS1(9,4)+r1x2Z*VR1(9,4)&
   +r1x2Z*(VS0(4,4)-ExZpE*VS1(4,4))
VS0(18,2)=PAz*VS0(8,2)+WPz*VS1(8,2)&
   +r1x2Z*(VS0(2,2)-ExZpE*VS1(2,2))
VS0(18,3)=PAz*VS0(8,3)+WPz*VS1(8,3)&
   +r1x2Z*(VS0(2,3)-ExZpE*VS1(2,3))
VS0(18,4)=PAz*VS0(8,4)+WPz*VS1(8,4)&
   +r1x2Z*(VS0(2,4)-ExZpE*VS1(2,4))&
   +HfxZpE*VS1(8,1)
VS0(19,2)=PAz*VS0(9,2)+WPz*VS1(9,2)&
   +r1x2Z*(VS0(3,2)-ExZpE*VS1(3,2))
VS0(19,3)=PAz*VS0(9,3)+WPz*VS1(9,3)&
   +r1x2Z*(VS0(3,3)-ExZpE*VS1(3,3))
VS0(19,4)=PAz*VS0(9,4)+WPz*VS1(9,4)&
   +r1x2Z*(VS0(3,4)-ExZpE*VS1(3,4))&
   +HfxZpE*VS1(9,1)
VS0(20,2)=PAz*VS0(10,2)+WPz*VS1(10,2)&
   +2D0*r1x2Z*(VS0(4,2)-ExZpE*VS1(4,2))
VS0(20,3)=PAz*VS0(10,3)+WPz*VS1(10,3)&
   +2D0*r1x2Z*(VS0(4,3)-ExZpE*VS1(4,3))
VS0(20,4)=PAz*VS0(10,4)+WPz*VS1(10,4)&
   +2D0*r1x2Z*(VS0(4,4)-ExZpE*VS1(4,4))&
   +HfxZpE*VS1(10,1)
CASE(3)
VS0(11,2)=PAx*VS0(5,2)+WPx*VS1(5,2)&
   +2D0*r1x2Z*(VS0(2,2)-ExZpE*VS1(2,2))&
   +HfxZpE*VS1(5,1)
VS0(11,3)=PAx*VS0(5,3)+WPx*VS1(5,3)&
   +2D0*r1x2Z*(VS0(2,3)-ExZpE*VS1(2,3))
VS0(11,4)=PAx*VS0(5,4)+WPx*VS1(5,4)&
   +2D0*r1x2Z*(VS0(2,4)-ExZpE*VS1(2,4))
VS0(12,2)=PAx*VS0(6,2)+WPx*VS1(6,2)&
   +r1x2Z*(VS0(3,2)-ExZpE*VS1(3,2))&
   +HfxZpE*VS1(6,1)
VS0(12,3)=PAx*VS0(6,3)+WPx*VS1(6,3)&
   +r1x2Z*(VS0(3,3)-ExZpE*VS1(3,3))
VS0(12,4)=PAx*VS0(6,4)+WPx*VS1(6,4)&
   +r1x2Z*(VS0(3,4)-ExZpE*VS1(3,4))
VS0(13,2)=PAy*VS0(6,2)+WPy*VS1(6,2)&
   +r1x2Z*(VS0(2,2)-ExZpE*VS1(2,2))
VS0(13,3)=PAy*VS0(6,3)+WPy*VS1(6,3)&
   +r1x2Z*(VS0(2,3)-ExZpE*VS1(2,3))&
   +HfxZpE*VS1(6,1)
VS0(13,4)=PAy*VS0(6,4)+WPy*VS1(6,4)&
   +r1x2Z*(VS0(2,4)-ExZpE*VS1(2,4))
VS0(14,2)=PAy*VS0(7,2)+WPy*VS1(7,2)&
   +2D0*r1x2Z*(VS0(3,2)-ExZpE*VS1(3,2))
VS0(14,3)=PAy*VS0(7,3)+WPy*VS1(7,3)&
   +2D0*r1x2Z*(VS0(3,3)-ExZpE*VS1(3,3))&
   +HfxZpE*VS1(7,1)
VS0(14,4)=PAy*VS0(7,4)+WPy*VS1(7,4)&
   +2D0*r1x2Z*(VS0(3,4)-ExZpE*VS1(3,4))
VS0(15,2)=PAx*VS0(8,2)+WPx*VS1(8,2)&
   +r1x2Z*(VS0(4,2)-ExZpE*VS1(4,2))&
   +HfxZpE*VS1(8,1)
VS0(15,3)=PAx*VS0(8,3)+WPx*VS1(8,3)&
   +r1x2Z*(VS0(4,3)-ExZpE*VS1(4,3))
VS0(15,4)=PAx*VS0(8,4)+WPx*VS1(8,4)&
   +r1x2Z*(VS0(4,4)-ExZpE*VS1(4,4))
VS0(16,2)=PAx*VS0(9,2)+WPx*VS1(9,2)&
   +HfxZpE*VS1(9,1)
VS0(16,3)=PAx*VS0(9,3)+WPx*VS1(9,3)
VS0(16,4)=PAx*VS0(9,4)+WPx*VS1(9,4)
VS0(17,2)=PAy*VS0(9,2)+WPy*VS1(9,2)&
   +r1x2Z*(VS0(4,2)-ExZpE*VS1(4,2))
VS0(17,3)=PAy*VS0(9,3)+WPy*VS1(9,3)&
   +r1x2Z*(VS0(4,3)-ExZpE*VS1(4,3))&
   +HfxZpE*VS1(9,1)
VS0(17,4)=PAy*VS0(9,4)+WPy*VS1(9,4)&
   +r1x2Z*(VS0(4,4)-ExZpE*VS1(4,4))
VS0(18,2)=PAz*VS0(8,2)+WPz*VS1(8,2)+r1x2Z*VR1(8,2)&
   +r1x2Z*(VS0(2,2)-ExZpE*VS1(2,2))
VS0(18,3)=PAz*VS0(8,3)+WPz*VS1(8,3)+r1x2Z*VR1(8,3)&
   +r1x2Z*(VS0(2,3)-ExZpE*VS1(2,3))
VS0(18,4)=PAz*VS0(8,4)+WPz*VS1(8,4)+r1x2Z*VR1(8,4)&
   +r1x2Z*(VS0(2,4)-ExZpE*VS1(2,4))&
   +HfxZpE*VS1(8,1)
VS0(19,2)=PAz*VS0(9,2)+WPz*VS1(9,2)+r1x2Z*VR1(9,2)&
   +r1x2Z*(VS0(3,2)-ExZpE*VS1(3,2))
VS0(19,3)=PAz*VS0(9,3)+WPz*VS1(9,3)+r1x2Z*VR1(9,3)&
   +r1x2Z*(VS0(3,3)-ExZpE*VS1(3,3))
VS0(19,4)=PAz*VS0(9,4)+WPz*VS1(9,4)+r1x2Z*VR1(9,4)&
   +r1x2Z*(VS0(3,4)-ExZpE*VS1(3,4))&
   +HfxZpE*VS1(9,1)
VS0(20,2)=PAz*VS0(10,2)+WPz*VS1(10,2)+r1x2Z*VR1(10,2)&
   +2D0*r1x2Z*(VS0(4,2)-ExZpE*VS1(4,2))
VS0(20,3)=PAz*VS0(10,3)+WPz*VS1(10,3)+r1x2Z*VR1(10,3)&
   +2D0*r1x2Z*(VS0(4,3)-ExZpE*VS1(4,3))
VS0(20,4)=PAz*VS0(10,4)+WPz*VS1(10,4)+r1x2Z*VR1(10,4)&
   +2D0*r1x2Z*(VS0(4,4)-ExZpE*VS1(4,4))&
   +HfxZpE*VS1(10,1)
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRf0p0'
STOP
END SELECT
END SUBROUTINE MVRRf0p0
