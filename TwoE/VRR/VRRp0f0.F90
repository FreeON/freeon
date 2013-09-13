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
   SUBROUTINE VRRp0f0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2E*VRR0(2,2)
      V(2)=r1x2E*ZxZpE*VRR1(2,2)
      V(3)=r1x2E*VRR0(2,3)
      V(4)=HfxZpE*VRR1(1,6)
      V(5)=r1x2E*ZxZpE*VRR1(2,3)
      V(6)=-V(5)
      V(7)=-V(2)
      V(8)=r1x2E*VRR0(2,4)
      V(9)=HfxZpE*VRR1(1,8)
      V(10)=r1x2E*ZxZpE*VRR1(2,4)
      V(11)=-V(10)
      V(12)=HfxZpE*VRR1(1,9)
      V(13)=r1x2E*VRR0(3,2)
      V(14)=r1x2E*ZxZpE*VRR1(3,2)
      V(15)=r1x2E*VRR0(3,3)
      V(16)=r1x2E*ZxZpE*VRR1(3,3)
      V(17)=-V(16)
      V(18)=-V(14)
      V(19)=r1x2E*VRR0(3,4)
      V(20)=r1x2E*ZxZpE*VRR1(3,4)
      V(21)=-V(20)
      V(22)=r1x2E*VRR0(4,2)
      V(23)=r1x2E*ZxZpE*VRR1(4,2)
      V(24)=r1x2E*VRR0(4,3)
      V(25)=r1x2E*ZxZpE*VRR1(4,3)
      V(26)=-V(25)
      V(27)=-V(23)
      V(28)=r1x2E*VRR0(4,4)
      V(29)=r1x2E*ZxZpE*VRR1(4,4)
      V(30)=-V(29)
      VRR0(2,11)=2.D0*V(1)-2.D0*V(2)+QCx*VRR0(2,5)+HfxZpE*VRR1(1,5)+WQx*VRR1(2,5)
      VRR0(2,12)=V(3)+V(4)+V(6)+QCx*VRR0(2,6)+WQx*VRR1(2,6)
      VRR0(2,13)=V(1)+V(7)+QCy*VRR0(2,6)+WQy*VRR1(2,6)
      VRR0(2,14)=2.D0*V(3)-2.D0*V(5)+QCy*VRR0(2,7)+WQy*VRR1(2,7)
      VRR0(2,15)=V(8)+V(9)+V(11)+QCx*VRR0(2,8)+WQx*VRR1(2,8)
      VRR0(2,16)=V(12)+QCx*VRR0(2,9)+WQx*VRR1(2,9)
      VRR0(2,17)=V(8)+V(11)+QCy*VRR0(2,9)+WQy*VRR1(2,9)
      VRR0(2,18)=V(1)+V(7)+QCz*VRR0(2,8)+WQz*VRR1(2,8)
      VRR0(2,19)=V(3)+V(6)+QCz*VRR0(2,9)+WQz*VRR1(2,9)
      VRR0(2,20)=2.D0*V(8)-2.D0*V(10)+QCz*VRR0(2,10)+WQz*VRR1(2,10)
      VRR0(3,11)=2.D0*V(13)-2.D0*V(14)+QCx*VRR0(3,5)+WQx*VRR1(3,5)
      VRR0(3,12)=V(15)+V(17)+QCx*VRR0(3,6)+WQx*VRR1(3,6)
      VRR0(3,13)=V(4)+V(13)+V(18)+QCy*VRR0(3,6)+WQy*VRR1(3,6)
      VRR0(3,14)=2.D0*V(15)-2.D0*V(16)+QCy*VRR0(3,7)+HfxZpE*VRR1(1,7)+WQy*VRR1(3,7)
      VRR0(3,15)=V(19)+V(21)+QCx*VRR0(3,8)+WQx*VRR1(3,8)
      VRR0(3,16)=QCx*VRR0(3,9)+WQx*VRR1(3,9)
      VRR0(3,17)=V(12)+V(19)+V(21)+QCy*VRR0(3,9)+WQy*VRR1(3,9)
      VRR0(3,18)=V(13)+V(18)+QCz*VRR0(3,8)+WQz*VRR1(3,8)
      VRR0(3,19)=V(15)+V(17)+QCz*VRR0(3,9)+WQz*VRR1(3,9)
      VRR0(3,20)=2.D0*V(19)-2.D0*V(20)+QCz*VRR0(3,10)+WQz*VRR1(3,10)
      VRR0(4,11)=2.D0*V(22)-2.D0*V(23)+QCx*VRR0(4,5)+WQx*VRR1(4,5)
      VRR0(4,12)=V(24)+V(26)+QCx*VRR0(4,6)+WQx*VRR1(4,6)
      VRR0(4,13)=V(22)+V(27)+QCy*VRR0(4,6)+WQy*VRR1(4,6)
      VRR0(4,14)=2.D0*V(24)-2.D0*V(25)+QCy*VRR0(4,7)+WQy*VRR1(4,7)
      VRR0(4,15)=V(28)+V(30)+QCx*VRR0(4,8)+WQx*VRR1(4,8)
      VRR0(4,16)=QCx*VRR0(4,9)+WQx*VRR1(4,9)
      VRR0(4,17)=V(28)+V(30)+QCy*VRR0(4,9)+WQy*VRR1(4,9)
      VRR0(4,18)=V(9)+V(22)+V(27)+QCz*VRR0(4,8)+WQz*VRR1(4,8)
      VRR0(4,19)=V(12)+V(24)+V(26)+QCz*VRR0(4,9)+WQz*VRR1(4,9)
      VRR0(4,20)=2.D0*V(28)-2.D0*V(29)+QCz*VRR0(4,10)+HfxZpE*VRR1(1,10)+WQz*VRR1(4,10)
END SUBROUTINE VRRp0f0
SUBROUTINE MVRRp0f0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(2,11)=QCx*VS0(2,5)+WQx*VS1(2,5)-r1x2E*VR1(2,5)&
   +2D0*r1x2E*(VS0(2,2)-ZxZpE*VS1(2,2))&
   +HfxZpE*VS1(1,5)
VS0(2,12)=QCx*VS0(2,6)+WQx*VS1(2,6)-r1x2E*VR1(2,6)&
   +r1x2E*(VS0(2,3)-ZxZpE*VS1(2,3))&
   +HfxZpE*VS1(1,6)
VS0(2,13)=QCy*VS0(2,6)+WQy*VS1(2,6)&
   +r1x2E*(VS0(2,2)-ZxZpE*VS1(2,2))
VS0(2,14)=QCy*VS0(2,7)+WQy*VS1(2,7)&
   +2D0*r1x2E*(VS0(2,3)-ZxZpE*VS1(2,3))
VS0(2,15)=QCx*VS0(2,8)+WQx*VS1(2,8)-r1x2E*VR1(2,8)&
   +r1x2E*(VS0(2,4)-ZxZpE*VS1(2,4))&
   +HfxZpE*VS1(1,8)
VS0(2,16)=QCx*VS0(2,9)+WQx*VS1(2,9)-r1x2E*VR1(2,9)&
   +HfxZpE*VS1(1,9)
VS0(2,17)=QCy*VS0(2,9)+WQy*VS1(2,9)&
   +r1x2E*(VS0(2,4)-ZxZpE*VS1(2,4))
VS0(2,18)=QCz*VS0(2,8)+WQz*VS1(2,8)&
   +r1x2E*(VS0(2,2)-ZxZpE*VS1(2,2))
VS0(2,19)=QCz*VS0(2,9)+WQz*VS1(2,9)&
   +r1x2E*(VS0(2,3)-ZxZpE*VS1(2,3))
VS0(2,20)=QCz*VS0(2,10)+WQz*VS1(2,10)&
   +2D0*r1x2E*(VS0(2,4)-ZxZpE*VS1(2,4))
VS0(3,11)=QCx*VS0(3,5)+WQx*VS1(3,5)-r1x2E*VR1(3,5)&
   +2D0*r1x2E*(VS0(3,2)-ZxZpE*VS1(3,2))
VS0(3,12)=QCx*VS0(3,6)+WQx*VS1(3,6)-r1x2E*VR1(3,6)&
   +r1x2E*(VS0(3,3)-ZxZpE*VS1(3,3))
VS0(3,13)=QCy*VS0(3,6)+WQy*VS1(3,6)&
   +r1x2E*(VS0(3,2)-ZxZpE*VS1(3,2))&
   +HfxZpE*VS1(1,6)
VS0(3,14)=QCy*VS0(3,7)+WQy*VS1(3,7)&
   +2D0*r1x2E*(VS0(3,3)-ZxZpE*VS1(3,3))&
   +HfxZpE*VS1(1,7)
VS0(3,15)=QCx*VS0(3,8)+WQx*VS1(3,8)-r1x2E*VR1(3,8)&
   +r1x2E*(VS0(3,4)-ZxZpE*VS1(3,4))
VS0(3,16)=QCx*VS0(3,9)+WQx*VS1(3,9)-r1x2E*VR1(3,9)
VS0(3,17)=QCy*VS0(3,9)+WQy*VS1(3,9)&
   +r1x2E*(VS0(3,4)-ZxZpE*VS1(3,4))&
   +HfxZpE*VS1(1,9)
VS0(3,18)=QCz*VS0(3,8)+WQz*VS1(3,8)&
   +r1x2E*(VS0(3,2)-ZxZpE*VS1(3,2))
VS0(3,19)=QCz*VS0(3,9)+WQz*VS1(3,9)&
   +r1x2E*(VS0(3,3)-ZxZpE*VS1(3,3))
VS0(3,20)=QCz*VS0(3,10)+WQz*VS1(3,10)&
   +2D0*r1x2E*(VS0(3,4)-ZxZpE*VS1(3,4))
VS0(4,11)=QCx*VS0(4,5)+WQx*VS1(4,5)-r1x2E*VR1(4,5)&
   +2D0*r1x2E*(VS0(4,2)-ZxZpE*VS1(4,2))
VS0(4,12)=QCx*VS0(4,6)+WQx*VS1(4,6)-r1x2E*VR1(4,6)&
   +r1x2E*(VS0(4,3)-ZxZpE*VS1(4,3))
VS0(4,13)=QCy*VS0(4,6)+WQy*VS1(4,6)&
   +r1x2E*(VS0(4,2)-ZxZpE*VS1(4,2))
VS0(4,14)=QCy*VS0(4,7)+WQy*VS1(4,7)&
   +2D0*r1x2E*(VS0(4,3)-ZxZpE*VS1(4,3))
VS0(4,15)=QCx*VS0(4,8)+WQx*VS1(4,8)-r1x2E*VR1(4,8)&
   +r1x2E*(VS0(4,4)-ZxZpE*VS1(4,4))
VS0(4,16)=QCx*VS0(4,9)+WQx*VS1(4,9)-r1x2E*VR1(4,9)
VS0(4,17)=QCy*VS0(4,9)+WQy*VS1(4,9)&
   +r1x2E*(VS0(4,4)-ZxZpE*VS1(4,4))
VS0(4,18)=QCz*VS0(4,8)+WQz*VS1(4,8)&
   +r1x2E*(VS0(4,2)-ZxZpE*VS1(4,2))&
   +HfxZpE*VS1(1,8)
VS0(4,19)=QCz*VS0(4,9)+WQz*VS1(4,9)&
   +r1x2E*(VS0(4,3)-ZxZpE*VS1(4,3))&
   +HfxZpE*VS1(1,9)
VS0(4,20)=QCz*VS0(4,10)+WQz*VS1(4,10)&
   +2D0*r1x2E*(VS0(4,4)-ZxZpE*VS1(4,4))&
   +HfxZpE*VS1(1,10)
CASE(2)
VS0(2,11)=QCx*VS0(2,5)+WQx*VS1(2,5)&
   +2D0*r1x2E*(VS0(2,2)-ZxZpE*VS1(2,2))&
   +HfxZpE*VS1(1,5)
VS0(2,12)=QCx*VS0(2,6)+WQx*VS1(2,6)&
   +r1x2E*(VS0(2,3)-ZxZpE*VS1(2,3))&
   +HfxZpE*VS1(1,6)
VS0(2,13)=QCy*VS0(2,6)+WQy*VS1(2,6)-r1x2E*VR1(2,6)&
   +r1x2E*(VS0(2,2)-ZxZpE*VS1(2,2))
VS0(2,14)=QCy*VS0(2,7)+WQy*VS1(2,7)-r1x2E*VR1(2,7)&
   +2D0*r1x2E*(VS0(2,3)-ZxZpE*VS1(2,3))
VS0(2,15)=QCx*VS0(2,8)+WQx*VS1(2,8)&
   +r1x2E*(VS0(2,4)-ZxZpE*VS1(2,4))&
   +HfxZpE*VS1(1,8)
VS0(2,16)=QCx*VS0(2,9)+WQx*VS1(2,9)&
   +HfxZpE*VS1(1,9)
VS0(2,17)=QCy*VS0(2,9)+WQy*VS1(2,9)-r1x2E*VR1(2,9)&
   +r1x2E*(VS0(2,4)-ZxZpE*VS1(2,4))
VS0(2,18)=QCz*VS0(2,8)+WQz*VS1(2,8)&
   +r1x2E*(VS0(2,2)-ZxZpE*VS1(2,2))
VS0(2,19)=QCz*VS0(2,9)+WQz*VS1(2,9)&
   +r1x2E*(VS0(2,3)-ZxZpE*VS1(2,3))
VS0(2,20)=QCz*VS0(2,10)+WQz*VS1(2,10)&
   +2D0*r1x2E*(VS0(2,4)-ZxZpE*VS1(2,4))
VS0(3,11)=QCx*VS0(3,5)+WQx*VS1(3,5)&
   +2D0*r1x2E*(VS0(3,2)-ZxZpE*VS1(3,2))
VS0(3,12)=QCx*VS0(3,6)+WQx*VS1(3,6)&
   +r1x2E*(VS0(3,3)-ZxZpE*VS1(3,3))
VS0(3,13)=QCy*VS0(3,6)+WQy*VS1(3,6)-r1x2E*VR1(3,6)&
   +r1x2E*(VS0(3,2)-ZxZpE*VS1(3,2))&
   +HfxZpE*VS1(1,6)
VS0(3,14)=QCy*VS0(3,7)+WQy*VS1(3,7)-r1x2E*VR1(3,7)&
   +2D0*r1x2E*(VS0(3,3)-ZxZpE*VS1(3,3))&
   +HfxZpE*VS1(1,7)
VS0(3,15)=QCx*VS0(3,8)+WQx*VS1(3,8)&
   +r1x2E*(VS0(3,4)-ZxZpE*VS1(3,4))
VS0(3,16)=QCx*VS0(3,9)+WQx*VS1(3,9)
VS0(3,17)=QCy*VS0(3,9)+WQy*VS1(3,9)-r1x2E*VR1(3,9)&
   +r1x2E*(VS0(3,4)-ZxZpE*VS1(3,4))&
   +HfxZpE*VS1(1,9)
VS0(3,18)=QCz*VS0(3,8)+WQz*VS1(3,8)&
   +r1x2E*(VS0(3,2)-ZxZpE*VS1(3,2))
VS0(3,19)=QCz*VS0(3,9)+WQz*VS1(3,9)&
   +r1x2E*(VS0(3,3)-ZxZpE*VS1(3,3))
VS0(3,20)=QCz*VS0(3,10)+WQz*VS1(3,10)&
   +2D0*r1x2E*(VS0(3,4)-ZxZpE*VS1(3,4))
VS0(4,11)=QCx*VS0(4,5)+WQx*VS1(4,5)&
   +2D0*r1x2E*(VS0(4,2)-ZxZpE*VS1(4,2))
VS0(4,12)=QCx*VS0(4,6)+WQx*VS1(4,6)&
   +r1x2E*(VS0(4,3)-ZxZpE*VS1(4,3))
VS0(4,13)=QCy*VS0(4,6)+WQy*VS1(4,6)-r1x2E*VR1(4,6)&
   +r1x2E*(VS0(4,2)-ZxZpE*VS1(4,2))
VS0(4,14)=QCy*VS0(4,7)+WQy*VS1(4,7)-r1x2E*VR1(4,7)&
   +2D0*r1x2E*(VS0(4,3)-ZxZpE*VS1(4,3))
VS0(4,15)=QCx*VS0(4,8)+WQx*VS1(4,8)&
   +r1x2E*(VS0(4,4)-ZxZpE*VS1(4,4))
VS0(4,16)=QCx*VS0(4,9)+WQx*VS1(4,9)
VS0(4,17)=QCy*VS0(4,9)+WQy*VS1(4,9)-r1x2E*VR1(4,9)&
   +r1x2E*(VS0(4,4)-ZxZpE*VS1(4,4))
VS0(4,18)=QCz*VS0(4,8)+WQz*VS1(4,8)&
   +r1x2E*(VS0(4,2)-ZxZpE*VS1(4,2))&
   +HfxZpE*VS1(1,8)
VS0(4,19)=QCz*VS0(4,9)+WQz*VS1(4,9)&
   +r1x2E*(VS0(4,3)-ZxZpE*VS1(4,3))&
   +HfxZpE*VS1(1,9)
VS0(4,20)=QCz*VS0(4,10)+WQz*VS1(4,10)&
   +2D0*r1x2E*(VS0(4,4)-ZxZpE*VS1(4,4))&
   +HfxZpE*VS1(1,10)
CASE(3)
VS0(2,11)=QCx*VS0(2,5)+WQx*VS1(2,5)&
   +2D0*r1x2E*(VS0(2,2)-ZxZpE*VS1(2,2))&
   +HfxZpE*VS1(1,5)
VS0(2,12)=QCx*VS0(2,6)+WQx*VS1(2,6)&
   +r1x2E*(VS0(2,3)-ZxZpE*VS1(2,3))&
   +HfxZpE*VS1(1,6)
VS0(2,13)=QCy*VS0(2,6)+WQy*VS1(2,6)&
   +r1x2E*(VS0(2,2)-ZxZpE*VS1(2,2))
VS0(2,14)=QCy*VS0(2,7)+WQy*VS1(2,7)&
   +2D0*r1x2E*(VS0(2,3)-ZxZpE*VS1(2,3))
VS0(2,15)=QCx*VS0(2,8)+WQx*VS1(2,8)&
   +r1x2E*(VS0(2,4)-ZxZpE*VS1(2,4))&
   +HfxZpE*VS1(1,8)
VS0(2,16)=QCx*VS0(2,9)+WQx*VS1(2,9)&
   +HfxZpE*VS1(1,9)
VS0(2,17)=QCy*VS0(2,9)+WQy*VS1(2,9)&
   +r1x2E*(VS0(2,4)-ZxZpE*VS1(2,4))
VS0(2,18)=QCz*VS0(2,8)+WQz*VS1(2,8)-r1x2E*VR1(2,8)&
   +r1x2E*(VS0(2,2)-ZxZpE*VS1(2,2))
VS0(2,19)=QCz*VS0(2,9)+WQz*VS1(2,9)-r1x2E*VR1(2,9)&
   +r1x2E*(VS0(2,3)-ZxZpE*VS1(2,3))
VS0(2,20)=QCz*VS0(2,10)+WQz*VS1(2,10)-r1x2E*VR1(2,10)&
   +2D0*r1x2E*(VS0(2,4)-ZxZpE*VS1(2,4))
VS0(3,11)=QCx*VS0(3,5)+WQx*VS1(3,5)&
   +2D0*r1x2E*(VS0(3,2)-ZxZpE*VS1(3,2))
VS0(3,12)=QCx*VS0(3,6)+WQx*VS1(3,6)&
   +r1x2E*(VS0(3,3)-ZxZpE*VS1(3,3))
VS0(3,13)=QCy*VS0(3,6)+WQy*VS1(3,6)&
   +r1x2E*(VS0(3,2)-ZxZpE*VS1(3,2))&
   +HfxZpE*VS1(1,6)
VS0(3,14)=QCy*VS0(3,7)+WQy*VS1(3,7)&
   +2D0*r1x2E*(VS0(3,3)-ZxZpE*VS1(3,3))&
   +HfxZpE*VS1(1,7)
VS0(3,15)=QCx*VS0(3,8)+WQx*VS1(3,8)&
   +r1x2E*(VS0(3,4)-ZxZpE*VS1(3,4))
VS0(3,16)=QCx*VS0(3,9)+WQx*VS1(3,9)
VS0(3,17)=QCy*VS0(3,9)+WQy*VS1(3,9)&
   +r1x2E*(VS0(3,4)-ZxZpE*VS1(3,4))&
   +HfxZpE*VS1(1,9)
VS0(3,18)=QCz*VS0(3,8)+WQz*VS1(3,8)-r1x2E*VR1(3,8)&
   +r1x2E*(VS0(3,2)-ZxZpE*VS1(3,2))
VS0(3,19)=QCz*VS0(3,9)+WQz*VS1(3,9)-r1x2E*VR1(3,9)&
   +r1x2E*(VS0(3,3)-ZxZpE*VS1(3,3))
VS0(3,20)=QCz*VS0(3,10)+WQz*VS1(3,10)-r1x2E*VR1(3,10)&
   +2D0*r1x2E*(VS0(3,4)-ZxZpE*VS1(3,4))
VS0(4,11)=QCx*VS0(4,5)+WQx*VS1(4,5)&
   +2D0*r1x2E*(VS0(4,2)-ZxZpE*VS1(4,2))
VS0(4,12)=QCx*VS0(4,6)+WQx*VS1(4,6)&
   +r1x2E*(VS0(4,3)-ZxZpE*VS1(4,3))
VS0(4,13)=QCy*VS0(4,6)+WQy*VS1(4,6)&
   +r1x2E*(VS0(4,2)-ZxZpE*VS1(4,2))
VS0(4,14)=QCy*VS0(4,7)+WQy*VS1(4,7)&
   +2D0*r1x2E*(VS0(4,3)-ZxZpE*VS1(4,3))
VS0(4,15)=QCx*VS0(4,8)+WQx*VS1(4,8)&
   +r1x2E*(VS0(4,4)-ZxZpE*VS1(4,4))
VS0(4,16)=QCx*VS0(4,9)+WQx*VS1(4,9)
VS0(4,17)=QCy*VS0(4,9)+WQy*VS1(4,9)&
   +r1x2E*(VS0(4,4)-ZxZpE*VS1(4,4))
VS0(4,18)=QCz*VS0(4,8)+WQz*VS1(4,8)-r1x2E*VR1(4,8)&
   +r1x2E*(VS0(4,2)-ZxZpE*VS1(4,2))&
   +HfxZpE*VS1(1,8)
VS0(4,19)=QCz*VS0(4,9)+WQz*VS1(4,9)-r1x2E*VR1(4,9)&
   +r1x2E*(VS0(4,3)-ZxZpE*VS1(4,3))&
   +HfxZpE*VS1(1,9)
VS0(4,20)=QCz*VS0(4,10)+WQz*VS1(4,10)-r1x2E*VR1(4,10)&
   +2D0*r1x2E*(VS0(4,4)-ZxZpE*VS1(4,4))&
   +HfxZpE*VS1(1,10)
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRp0f0'
STOP
END SELECT
END SUBROUTINE MVRRp0f0
