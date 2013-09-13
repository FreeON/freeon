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

   SUBROUTINE VRRd0d0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2E*VRR0(5,1)
      V(2)=r1x2E*ZxZpE*VRR1(5,1)
      V(3)=-V(2)
      V(4)=r1x2Z*VRR0(1,6)
      V(5)=ExZpE*r1x2Z*VRR1(1,6)
      V(6)=-V(5)
      V(7)=HfxZpE*VRR1(2,3)
      V(8)=r1x2Z*VRR0(1,8)
      V(9)=ExZpE*r1x2Z*VRR1(1,8)
      V(10)=-V(9)
      V(11)=HfxZpE*VRR1(2,4)
      V(12)=r1x2Z*VRR0(1,9)
      V(13)=ExZpE*r1x2Z*VRR1(1,9)
      V(14)=-V(13)
      V(15)=r1x2E*VRR0(6,1)
      V(16)=HfxZpE*VRR1(3,2)
      V(17)=r1x2E*ZxZpE*VRR1(6,1)
      V(18)=-V(17)
      V(19)=HfxZpE*VRR1(3,3)
      V(20)=HfxZpE*VRR1(3,4)
      V(21)=r1x2E*VRR0(7,1)
      V(22)=r1x2E*ZxZpE*VRR1(7,1)
      V(23)=-V(22)
      V(24)=r1x2E*VRR0(8,1)
      V(25)=HfxZpE*VRR1(4,2)
      V(26)=r1x2E*ZxZpE*VRR1(8,1)
      V(27)=-V(26)
      V(28)=HfxZpE*VRR1(4,3)
      V(29)=HfxZpE*VRR1(4,4)
      V(30)=r1x2E*VRR0(9,1)
      V(31)=r1x2E*ZxZpE*VRR1(9,1)
      V(32)=-V(31)
      V(33)=r1x2E*VRR0(10,1)
      V(34)=r1x2E*ZxZpE*VRR1(10,1)
      V(35)=-V(34)
      VRR0(5,5)=V(1)+V(3)+QCx*VRR0(5,2)+2.D0*HfxZpE*VRR1(2,2)+WQx*VRR1(5,2)
      VRR0(5,6)=V(4)+V(6)+V(7)+PAx*VRR0(2,6)+WPx*VRR1(2,6)
      VRR0(5,7)=V(1)+V(3)+QCy*VRR0(5,3)+WQy*VRR1(5,3)
      VRR0(5,8)=V(8)+V(10)+V(11)+PAx*VRR0(2,8)+WPx*VRR1(2,8)
      VRR0(5,9)=V(12)+V(14)+PAx*VRR0(2,9)+WPx*VRR1(2,9)
      VRR0(5,10)=V(1)+V(3)+QCz*VRR0(5,4)+WQz*VRR1(5,4)
      VRR0(6,5)=V(15)+V(16)+V(18)+QCx*VRR0(6,2)+WQx*VRR1(6,2)
      VRR0(6,6)=V(19)+QCx*VRR0(6,3)+WQx*VRR1(6,3)
      VRR0(6,7)=V(7)+V(15)+V(18)+QCy*VRR0(6,3)+WQy*VRR1(6,3)
      VRR0(6,8)=V(20)+QCx*VRR0(6,4)+WQx*VRR1(6,4)
      VRR0(6,9)=V(11)+QCy*VRR0(6,4)+WQy*VRR1(6,4)
      VRR0(6,10)=V(15)+V(18)+QCz*VRR0(6,4)+WQz*VRR1(6,4)
      VRR0(7,5)=V(21)+V(23)+QCx*VRR0(7,2)+WQx*VRR1(7,2)
      VRR0(7,6)=V(4)+V(6)+V(16)+PAy*VRR0(3,6)+WPy*VRR1(3,6)
      VRR0(7,7)=2.D0*V(19)+V(21)+V(23)+QCy*VRR0(7,3)+WQy*VRR1(7,3)
      VRR0(7,8)=V(8)+V(10)+PAy*VRR0(3,8)+WPy*VRR1(3,8)
      VRR0(7,9)=V(12)+V(14)+V(20)+PAy*VRR0(3,9)+WPy*VRR1(3,9)
      VRR0(7,10)=V(21)+V(23)+QCz*VRR0(7,4)+WQz*VRR1(7,4)
      VRR0(8,5)=V(24)+V(25)+V(27)+QCx*VRR0(8,2)+WQx*VRR1(8,2)
      VRR0(8,6)=V(28)+QCx*VRR0(8,3)+WQx*VRR1(8,3)
      VRR0(8,7)=V(24)+V(27)+QCy*VRR0(8,3)+WQy*VRR1(8,3)
      VRR0(8,8)=V(29)+QCx*VRR0(8,4)+WQx*VRR1(8,4)
      VRR0(8,9)=QCy*VRR0(8,4)+WQy*VRR1(8,4)
      VRR0(8,10)=V(11)+V(24)+V(27)+QCz*VRR0(8,4)+WQz*VRR1(8,4)
      VRR0(9,5)=V(30)+V(32)+QCx*VRR0(9,2)+WQx*VRR1(9,2)
      VRR0(9,6)=QCx*VRR0(9,3)+WQx*VRR1(9,3)
      VRR0(9,7)=V(28)+V(30)+V(32)+QCy*VRR0(9,3)+WQy*VRR1(9,3)
      VRR0(9,8)=QCx*VRR0(9,4)+WQx*VRR1(9,4)
      VRR0(9,9)=V(29)+QCy*VRR0(9,4)+WQy*VRR1(9,4)
      VRR0(9,10)=V(20)+V(30)+V(32)+QCz*VRR0(9,4)+WQz*VRR1(9,4)
      VRR0(10,5)=V(33)+V(35)+QCx*VRR0(10,2)+WQx*VRR1(10,2)
      VRR0(10,6)=V(4)+V(6)+PAz*VRR0(4,6)+WPz*VRR1(4,6)
      VRR0(10,7)=V(33)+V(35)+QCy*VRR0(10,3)+WQy*VRR1(10,3)
      VRR0(10,8)=V(8)+V(10)+V(25)+PAz*VRR0(4,8)+WPz*VRR1(4,8)
      VRR0(10,9)=V(12)+V(14)+V(28)+PAz*VRR0(4,9)+WPz*VRR1(4,9)
      VRR0(10,10)=2.D0*V(29)+V(33)+V(35)+QCz*VRR0(10,4)+WQz*VRR1(10,4)
END SUBROUTINE VRRd0d0
SUBROUTINE MVRRd0d0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(5,5)=QCx*VS0(5,2)+WQx*VS1(5,2)-r1x2E*VR1(5,2)&
   +r1x2E*(VS0(5,1)-ZxZpE*VS1(5,1))&
   +2D0*HfxZpE*VS1(2,2)
VS0(5,6)=QCx*VS0(5,3)+WQx*VS1(5,3)-r1x2E*VR1(5,3)&
   +2D0*HfxZpE*VS1(2,3)
VS0(5,7)=QCy*VS0(5,3)+WQy*VS1(5,3)&
   +r1x2E*(VS0(5,1)-ZxZpE*VS1(5,1))
VS0(5,8)=QCx*VS0(5,4)+WQx*VS1(5,4)-r1x2E*VR1(5,4)&
   +2D0*HfxZpE*VS1(2,4)
VS0(5,9)=QCy*VS0(5,4)+WQy*VS1(5,4)
VS0(5,10)=QCz*VS0(5,4)+WQz*VS1(5,4)&
   +r1x2E*(VS0(5,1)-ZxZpE*VS1(5,1))
VS0(6,5)=QCx*VS0(6,2)+WQx*VS1(6,2)-r1x2E*VR1(6,2)&
   +r1x2E*(VS0(6,1)-ZxZpE*VS1(6,1))&
   +HfxZpE*VS1(3,2)
VS0(6,6)=QCx*VS0(6,3)+WQx*VS1(6,3)-r1x2E*VR1(6,3)&
   +HfxZpE*VS1(3,3)
VS0(6,7)=QCy*VS0(6,3)+WQy*VS1(6,3)&
   +r1x2E*(VS0(6,1)-ZxZpE*VS1(6,1))&
   +HfxZpE*VS1(2,3)
VS0(6,8)=QCx*VS0(6,4)+WQx*VS1(6,4)-r1x2E*VR1(6,4)&
   +HfxZpE*VS1(3,4)
VS0(6,9)=QCy*VS0(6,4)+WQy*VS1(6,4)&
   +HfxZpE*VS1(2,4)
VS0(6,10)=QCz*VS0(6,4)+WQz*VS1(6,4)&
   +r1x2E*(VS0(6,1)-ZxZpE*VS1(6,1))
VS0(7,5)=QCx*VS0(7,2)+WQx*VS1(7,2)-r1x2E*VR1(7,2)&
   +r1x2E*(VS0(7,1)-ZxZpE*VS1(7,1))
VS0(7,6)=QCx*VS0(7,3)+WQx*VS1(7,3)-r1x2E*VR1(7,3)
VS0(7,7)=QCy*VS0(7,3)+WQy*VS1(7,3)&
   +r1x2E*(VS0(7,1)-ZxZpE*VS1(7,1))&
   +2D0*HfxZpE*VS1(3,3)
VS0(7,8)=QCx*VS0(7,4)+WQx*VS1(7,4)-r1x2E*VR1(7,4)
VS0(7,9)=QCy*VS0(7,4)+WQy*VS1(7,4)&
   +2D0*HfxZpE*VS1(3,4)
VS0(7,10)=QCz*VS0(7,4)+WQz*VS1(7,4)&
   +r1x2E*(VS0(7,1)-ZxZpE*VS1(7,1))
VS0(8,5)=QCx*VS0(8,2)+WQx*VS1(8,2)-r1x2E*VR1(8,2)&
   +r1x2E*(VS0(8,1)-ZxZpE*VS1(8,1))&
   +HfxZpE*VS1(4,2)
VS0(8,6)=QCx*VS0(8,3)+WQx*VS1(8,3)-r1x2E*VR1(8,3)&
   +HfxZpE*VS1(4,3)
VS0(8,7)=QCy*VS0(8,3)+WQy*VS1(8,3)&
   +r1x2E*(VS0(8,1)-ZxZpE*VS1(8,1))
VS0(8,8)=QCx*VS0(8,4)+WQx*VS1(8,4)-r1x2E*VR1(8,4)&
   +HfxZpE*VS1(4,4)
VS0(8,9)=QCy*VS0(8,4)+WQy*VS1(8,4)
VS0(8,10)=QCz*VS0(8,4)+WQz*VS1(8,4)&
   +r1x2E*(VS0(8,1)-ZxZpE*VS1(8,1))&
   +HfxZpE*VS1(2,4)
VS0(9,5)=QCx*VS0(9,2)+WQx*VS1(9,2)-r1x2E*VR1(9,2)&
   +r1x2E*(VS0(9,1)-ZxZpE*VS1(9,1))
VS0(9,6)=QCx*VS0(9,3)+WQx*VS1(9,3)-r1x2E*VR1(9,3)
VS0(9,7)=QCy*VS0(9,3)+WQy*VS1(9,3)&
   +r1x2E*(VS0(9,1)-ZxZpE*VS1(9,1))&
   +HfxZpE*VS1(4,3)
VS0(9,8)=QCx*VS0(9,4)+WQx*VS1(9,4)-r1x2E*VR1(9,4)
VS0(9,9)=QCy*VS0(9,4)+WQy*VS1(9,4)&
   +HfxZpE*VS1(4,4)
VS0(9,10)=QCz*VS0(9,4)+WQz*VS1(9,4)&
   +r1x2E*(VS0(9,1)-ZxZpE*VS1(9,1))&
   +HfxZpE*VS1(3,4)
VS0(10,5)=QCx*VS0(10,2)+WQx*VS1(10,2)-r1x2E*VR1(10,2)&
   +r1x2E*(VS0(10,1)-ZxZpE*VS1(10,1))
VS0(10,6)=QCx*VS0(10,3)+WQx*VS1(10,3)-r1x2E*VR1(10,3)
VS0(10,7)=QCy*VS0(10,3)+WQy*VS1(10,3)&
   +r1x2E*(VS0(10,1)-ZxZpE*VS1(10,1))
VS0(10,8)=QCx*VS0(10,4)+WQx*VS1(10,4)-r1x2E*VR1(10,4)
VS0(10,9)=QCy*VS0(10,4)+WQy*VS1(10,4)
VS0(10,10)=QCz*VS0(10,4)+WQz*VS1(10,4)&
   +r1x2E*(VS0(10,1)-ZxZpE*VS1(10,1))&
   +2D0*HfxZpE*VS1(4,4)
CASE(2)
VS0(5,5)=QCx*VS0(5,2)+WQx*VS1(5,2)&
   +r1x2E*(VS0(5,1)-ZxZpE*VS1(5,1))&
   +2D0*HfxZpE*VS1(2,2)
VS0(5,6)=QCx*VS0(5,3)+WQx*VS1(5,3)&
   +2D0*HfxZpE*VS1(2,3)
VS0(5,7)=QCy*VS0(5,3)+WQy*VS1(5,3)-r1x2E*VR1(5,3)&
   +r1x2E*(VS0(5,1)-ZxZpE*VS1(5,1))
VS0(5,8)=QCx*VS0(5,4)+WQx*VS1(5,4)&
   +2D0*HfxZpE*VS1(2,4)
VS0(5,9)=QCy*VS0(5,4)+WQy*VS1(5,4)-r1x2E*VR1(5,4)
VS0(5,10)=QCz*VS0(5,4)+WQz*VS1(5,4)&
   +r1x2E*(VS0(5,1)-ZxZpE*VS1(5,1))
VS0(6,5)=QCx*VS0(6,2)+WQx*VS1(6,2)&
   +r1x2E*(VS0(6,1)-ZxZpE*VS1(6,1))&
   +HfxZpE*VS1(3,2)
VS0(6,6)=QCx*VS0(6,3)+WQx*VS1(6,3)&
   +HfxZpE*VS1(3,3)
VS0(6,7)=QCy*VS0(6,3)+WQy*VS1(6,3)-r1x2E*VR1(6,3)&
   +r1x2E*(VS0(6,1)-ZxZpE*VS1(6,1))&
   +HfxZpE*VS1(2,3)
VS0(6,8)=QCx*VS0(6,4)+WQx*VS1(6,4)&
   +HfxZpE*VS1(3,4)
VS0(6,9)=QCy*VS0(6,4)+WQy*VS1(6,4)-r1x2E*VR1(6,4)&
   +HfxZpE*VS1(2,4)
VS0(6,10)=QCz*VS0(6,4)+WQz*VS1(6,4)&
   +r1x2E*(VS0(6,1)-ZxZpE*VS1(6,1))
VS0(7,5)=QCx*VS0(7,2)+WQx*VS1(7,2)&
   +r1x2E*(VS0(7,1)-ZxZpE*VS1(7,1))
VS0(7,6)=QCx*VS0(7,3)+WQx*VS1(7,3)
VS0(7,7)=QCy*VS0(7,3)+WQy*VS1(7,3)-r1x2E*VR1(7,3)&
   +r1x2E*(VS0(7,1)-ZxZpE*VS1(7,1))&
   +2D0*HfxZpE*VS1(3,3)
VS0(7,8)=QCx*VS0(7,4)+WQx*VS1(7,4)
VS0(7,9)=QCy*VS0(7,4)+WQy*VS1(7,4)-r1x2E*VR1(7,4)&
   +2D0*HfxZpE*VS1(3,4)
VS0(7,10)=QCz*VS0(7,4)+WQz*VS1(7,4)&
   +r1x2E*(VS0(7,1)-ZxZpE*VS1(7,1))
VS0(8,5)=QCx*VS0(8,2)+WQx*VS1(8,2)&
   +r1x2E*(VS0(8,1)-ZxZpE*VS1(8,1))&
   +HfxZpE*VS1(4,2)
VS0(8,6)=QCx*VS0(8,3)+WQx*VS1(8,3)&
   +HfxZpE*VS1(4,3)
VS0(8,7)=QCy*VS0(8,3)+WQy*VS1(8,3)-r1x2E*VR1(8,3)&
   +r1x2E*(VS0(8,1)-ZxZpE*VS1(8,1))
VS0(8,8)=QCx*VS0(8,4)+WQx*VS1(8,4)&
   +HfxZpE*VS1(4,4)
VS0(8,9)=QCy*VS0(8,4)+WQy*VS1(8,4)-r1x2E*VR1(8,4)
VS0(8,10)=QCz*VS0(8,4)+WQz*VS1(8,4)&
   +r1x2E*(VS0(8,1)-ZxZpE*VS1(8,1))&
   +HfxZpE*VS1(2,4)
VS0(9,5)=QCx*VS0(9,2)+WQx*VS1(9,2)&
   +r1x2E*(VS0(9,1)-ZxZpE*VS1(9,1))
VS0(9,6)=QCx*VS0(9,3)+WQx*VS1(9,3)
VS0(9,7)=QCy*VS0(9,3)+WQy*VS1(9,3)-r1x2E*VR1(9,3)&
   +r1x2E*(VS0(9,1)-ZxZpE*VS1(9,1))&
   +HfxZpE*VS1(4,3)
VS0(9,8)=QCx*VS0(9,4)+WQx*VS1(9,4)
VS0(9,9)=QCy*VS0(9,4)+WQy*VS1(9,4)-r1x2E*VR1(9,4)&
   +HfxZpE*VS1(4,4)
VS0(9,10)=QCz*VS0(9,4)+WQz*VS1(9,4)&
   +r1x2E*(VS0(9,1)-ZxZpE*VS1(9,1))&
   +HfxZpE*VS1(3,4)
VS0(10,5)=QCx*VS0(10,2)+WQx*VS1(10,2)&
   +r1x2E*(VS0(10,1)-ZxZpE*VS1(10,1))
VS0(10,6)=QCx*VS0(10,3)+WQx*VS1(10,3)
VS0(10,7)=QCy*VS0(10,3)+WQy*VS1(10,3)-r1x2E*VR1(10,3)&
   +r1x2E*(VS0(10,1)-ZxZpE*VS1(10,1))
VS0(10,8)=QCx*VS0(10,4)+WQx*VS1(10,4)
VS0(10,9)=QCy*VS0(10,4)+WQy*VS1(10,4)-r1x2E*VR1(10,4)
VS0(10,10)=QCz*VS0(10,4)+WQz*VS1(10,4)&
   +r1x2E*(VS0(10,1)-ZxZpE*VS1(10,1))&
   +2D0*HfxZpE*VS1(4,4)
CASE(3)
VS0(5,5)=QCx*VS0(5,2)+WQx*VS1(5,2)&
   +r1x2E*(VS0(5,1)-ZxZpE*VS1(5,1))&
   +2D0*HfxZpE*VS1(2,2)
VS0(5,6)=QCx*VS0(5,3)+WQx*VS1(5,3)&
   +2D0*HfxZpE*VS1(2,3)
VS0(5,7)=QCy*VS0(5,3)+WQy*VS1(5,3)&
   +r1x2E*(VS0(5,1)-ZxZpE*VS1(5,1))
VS0(5,8)=QCx*VS0(5,4)+WQx*VS1(5,4)&
   +2D0*HfxZpE*VS1(2,4)
VS0(5,9)=QCy*VS0(5,4)+WQy*VS1(5,4)
VS0(5,10)=QCz*VS0(5,4)+WQz*VS1(5,4)-r1x2E*VR1(5,4)&
   +r1x2E*(VS0(5,1)-ZxZpE*VS1(5,1))
VS0(6,5)=QCx*VS0(6,2)+WQx*VS1(6,2)&
   +r1x2E*(VS0(6,1)-ZxZpE*VS1(6,1))&
   +HfxZpE*VS1(3,2)
VS0(6,6)=QCx*VS0(6,3)+WQx*VS1(6,3)&
   +HfxZpE*VS1(3,3)
VS0(6,7)=QCy*VS0(6,3)+WQy*VS1(6,3)&
   +r1x2E*(VS0(6,1)-ZxZpE*VS1(6,1))&
   +HfxZpE*VS1(2,3)
VS0(6,8)=QCx*VS0(6,4)+WQx*VS1(6,4)&
   +HfxZpE*VS1(3,4)
VS0(6,9)=QCy*VS0(6,4)+WQy*VS1(6,4)&
   +HfxZpE*VS1(2,4)
VS0(6,10)=QCz*VS0(6,4)+WQz*VS1(6,4)-r1x2E*VR1(6,4)&
   +r1x2E*(VS0(6,1)-ZxZpE*VS1(6,1))
VS0(7,5)=QCx*VS0(7,2)+WQx*VS1(7,2)&
   +r1x2E*(VS0(7,1)-ZxZpE*VS1(7,1))
VS0(7,6)=QCx*VS0(7,3)+WQx*VS1(7,3)
VS0(7,7)=QCy*VS0(7,3)+WQy*VS1(7,3)&
   +r1x2E*(VS0(7,1)-ZxZpE*VS1(7,1))&
   +2D0*HfxZpE*VS1(3,3)
VS0(7,8)=QCx*VS0(7,4)+WQx*VS1(7,4)
VS0(7,9)=QCy*VS0(7,4)+WQy*VS1(7,4)&
   +2D0*HfxZpE*VS1(3,4)
VS0(7,10)=QCz*VS0(7,4)+WQz*VS1(7,4)-r1x2E*VR1(7,4)&
   +r1x2E*(VS0(7,1)-ZxZpE*VS1(7,1))
VS0(8,5)=QCx*VS0(8,2)+WQx*VS1(8,2)&
   +r1x2E*(VS0(8,1)-ZxZpE*VS1(8,1))&
   +HfxZpE*VS1(4,2)
VS0(8,6)=QCx*VS0(8,3)+WQx*VS1(8,3)&
   +HfxZpE*VS1(4,3)
VS0(8,7)=QCy*VS0(8,3)+WQy*VS1(8,3)&
   +r1x2E*(VS0(8,1)-ZxZpE*VS1(8,1))
VS0(8,8)=QCx*VS0(8,4)+WQx*VS1(8,4)&
   +HfxZpE*VS1(4,4)
VS0(8,9)=QCy*VS0(8,4)+WQy*VS1(8,4)
VS0(8,10)=QCz*VS0(8,4)+WQz*VS1(8,4)-r1x2E*VR1(8,4)&
   +r1x2E*(VS0(8,1)-ZxZpE*VS1(8,1))&
   +HfxZpE*VS1(2,4)
VS0(9,5)=QCx*VS0(9,2)+WQx*VS1(9,2)&
   +r1x2E*(VS0(9,1)-ZxZpE*VS1(9,1))
VS0(9,6)=QCx*VS0(9,3)+WQx*VS1(9,3)
VS0(9,7)=QCy*VS0(9,3)+WQy*VS1(9,3)&
   +r1x2E*(VS0(9,1)-ZxZpE*VS1(9,1))&
   +HfxZpE*VS1(4,3)
VS0(9,8)=QCx*VS0(9,4)+WQx*VS1(9,4)
VS0(9,9)=QCy*VS0(9,4)+WQy*VS1(9,4)&
   +HfxZpE*VS1(4,4)
VS0(9,10)=QCz*VS0(9,4)+WQz*VS1(9,4)-r1x2E*VR1(9,4)&
   +r1x2E*(VS0(9,1)-ZxZpE*VS1(9,1))&
   +HfxZpE*VS1(3,4)
VS0(10,5)=QCx*VS0(10,2)+WQx*VS1(10,2)&
   +r1x2E*(VS0(10,1)-ZxZpE*VS1(10,1))
VS0(10,6)=QCx*VS0(10,3)+WQx*VS1(10,3)
VS0(10,7)=QCy*VS0(10,3)+WQy*VS1(10,3)&
   +r1x2E*(VS0(10,1)-ZxZpE*VS1(10,1))
VS0(10,8)=QCx*VS0(10,4)+WQx*VS1(10,4)
VS0(10,9)=QCy*VS0(10,4)+WQy*VS1(10,4)
VS0(10,10)=QCz*VS0(10,4)+WQz*VS1(10,4)-r1x2E*VR1(10,4)&
   +r1x2E*(VS0(10,1)-ZxZpE*VS1(10,1))&
   +2D0*HfxZpE*VS1(4,4)
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRd0d0'
STOP
END SELECT
END SUBROUTINE MVRRd0d0
