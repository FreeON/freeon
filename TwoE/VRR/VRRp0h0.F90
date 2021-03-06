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
   SUBROUTINE VRRp0h0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2E*VRR0(2,12)
      V(2)=r1x2E*ZxZpE*VRR1(2,12)
      V(3)=r1x2E*VRR0(2,13)
      V(4)=HfxZpE*VRR1(1,23)
      V(5)=r1x2E*ZxZpE*VRR1(2,13)
      V(6)=r1x2E*VRR0(2,15)
      V(7)=r1x2E*ZxZpE*VRR1(2,15)
      V(8)=r1x2E*VRR0(2,16)
      V(9)=2.D0*V(8)
      V(10)=r1x2E*ZxZpE*VRR1(2,16)
      V(11)=-2.D0*V(10)
      V(12)=r1x2E*VRR0(2,17)
      V(13)=HfxZpE*VRR1(1,28)
      V(14)=r1x2E*ZxZpE*VRR1(2,17)
      V(15)=r1x2E*VRR0(2,18)
      V(16)=HfxZpE*VRR1(1,30)
      V(17)=r1x2E*ZxZpE*VRR1(2,18)
      V(18)=r1x2E*VRR0(2,19)
      V(19)=HfxZpE*VRR1(1,31)
      V(20)=r1x2E*ZxZpE*VRR1(2,19)
      V(21)=r1x2E*VRR0(3,12)
      V(22)=r1x2E*ZxZpE*VRR1(3,12)
      V(23)=r1x2E*VRR0(3,13)
      V(24)=r1x2E*ZxZpE*VRR1(3,13)
      V(25)=r1x2E*VRR0(3,15)
      V(26)=r1x2E*ZxZpE*VRR1(3,15)
      V(27)=r1x2E*VRR0(3,16)
      V(28)=2.D0*V(27)
      V(29)=r1x2E*ZxZpE*VRR1(3,16)
      V(30)=-2.D0*V(29)
      V(31)=r1x2E*VRR0(3,17)
      V(32)=r1x2E*ZxZpE*VRR1(3,17)
      V(33)=r1x2E*VRR0(3,18)
      V(34)=r1x2E*ZxZpE*VRR1(3,18)
      V(35)=r1x2E*VRR0(3,19)
      V(36)=r1x2E*ZxZpE*VRR1(3,19)
      V(37)=HfxZpE*VRR1(1,32)
      V(38)=r1x2E*VRR0(4,12)
      V(39)=r1x2E*ZxZpE*VRR1(4,12)
      V(40)=r1x2E*VRR0(4,13)
      V(41)=r1x2E*ZxZpE*VRR1(4,13)
      V(42)=r1x2E*VRR0(4,15)
      V(43)=r1x2E*ZxZpE*VRR1(4,15)
      V(44)=r1x2E*VRR0(4,16)
      V(45)=2.D0*V(44)
      V(46)=r1x2E*ZxZpE*VRR1(4,16)
      V(47)=-2.D0*V(46)
      V(48)=r1x2E*VRR0(4,17)
      V(49)=r1x2E*ZxZpE*VRR1(4,17)
      V(50)=r1x2E*VRR0(4,18)
      V(51)=r1x2E*ZxZpE*VRR1(4,18)
      V(52)=r1x2E*VRR0(4,19)
      V(53)=r1x2E*ZxZpE*VRR1(4,19)
      VRR0(2,36)=4.D0*r1x2E*VRR0(2,11)+QCx*VRR0(2,21)+HfxZpE*VRR1(1,21)-4.D0*r1x2E*ZxZpE*VRR1(2,11)+WQx*VRR1(2,21)
      VRR0(2,37)=3.D0*V(1)-3.D0*V(2)+QCx*VRR0(2,22)+HfxZpE*VRR1(1,22)+WQx*VRR1(2,22)
      VRR0(2,38)=2.D0*V(3)+V(4)-2.D0*V(5)+QCx*VRR0(2,23)+WQx*VRR1(2,23)
      VRR0(2,39)=2.D0*V(1)-2.D0*V(2)+QCy*VRR0(2,23)+WQy*VRR1(2,23)
      VRR0(2,40)=3.D0*V(3)-3.D0*V(5)+QCy*VRR0(2,24)+WQy*VRR1(2,24)
      VRR0(2,41)=4.D0*r1x2E*VRR0(2,14)+QCy*VRR0(2,25)-4.D0*r1x2E*ZxZpE*VRR1(2,14)+WQy*VRR1(2,25)
      VRR0(2,42)=3.D0*V(6)-3.D0*V(7)+QCx*VRR0(2,26)+HfxZpE*VRR1(1,26)+WQx*VRR1(2,26)
      VRR0(2,43)=V(9)+V(11)+QCx*VRR0(2,27)+HfxZpE*VRR1(1,27)+WQx*VRR1(2,27)
      VRR0(2,44)=V(12)+V(13)-V(14)+QCx*VRR0(2,28)+WQx*VRR1(2,28)
      VRR0(2,45)=V(9)+V(11)+QCy*VRR0(2,28)+WQy*VRR1(2,28)
      VRR0(2,46)=3.D0*V(12)-3.D0*V(14)+QCy*VRR0(2,29)+WQy*VRR1(2,29)
      VRR0(2,47)=2.D0*V(15)+V(16)-2.D0*V(17)+QCx*VRR0(2,30)+WQx*VRR1(2,30)
      VRR0(2,48)=V(18)+V(19)-V(20)+QCx*VRR0(2,31)+WQx*VRR1(2,31)
      VRR0(2,49)=V(15)-V(17)+QCy*VRR0(2,31)+WQy*VRR1(2,31)
      VRR0(2,50)=2.D0*V(18)-2.D0*V(20)+QCy*VRR0(2,32)+WQy*VRR1(2,32)
      VRR0(2,51)=2.D0*V(6)-2.D0*V(7)+QCz*VRR0(2,30)+WQz*VRR1(2,30)
      VRR0(2,52)=V(9)+V(11)+QCz*VRR0(2,31)+WQz*VRR1(2,31)
      VRR0(2,53)=2.D0*V(12)-2.D0*V(14)+QCz*VRR0(2,32)+WQz*VRR1(2,32)
      VRR0(2,54)=3.D0*V(15)-3.D0*V(17)+QCz*VRR0(2,33)+WQz*VRR1(2,33)
      VRR0(2,55)=3.D0*V(18)-3.D0*V(20)+QCz*VRR0(2,34)+WQz*VRR1(2,34)
      VRR0(2,56)=4.D0*r1x2E*VRR0(2,20)+QCz*VRR0(2,35)-4.D0*r1x2E*ZxZpE*VRR1(2,20)+WQz*VRR1(2,35)
      VRR0(3,36)=4.D0*r1x2E*VRR0(3,11)+QCx*VRR0(3,21)-4.D0*r1x2E*ZxZpE*VRR1(3,11)+WQx*VRR1(3,21)
      VRR0(3,37)=3.D0*V(21)-3.D0*V(22)+QCx*VRR0(3,22)+WQx*VRR1(3,22)
      VRR0(3,38)=2.D0*V(23)-2.D0*V(24)+QCx*VRR0(3,23)+WQx*VRR1(3,23)
      VRR0(3,39)=V(4)+2.D0*V(21)-2.D0*V(22)+QCy*VRR0(3,23)+WQy*VRR1(3,23)
      VRR0(3,40)=3.D0*V(23)-3.D0*V(24)+QCy*VRR0(3,24)+HfxZpE*VRR1(1,24)+WQy*VRR1(3,24)
      VRR0(3,41)=4.D0*r1x2E*VRR0(3,14)+QCy*VRR0(3,25)+HfxZpE*VRR1(1,25)-4.D0*r1x2E*ZxZpE*VRR1(3,14)+WQy*VRR1(3,25)
      VRR0(3,42)=3.D0*V(25)-3.D0*V(26)+QCx*VRR0(3,26)+WQx*VRR1(3,26)
      VRR0(3,43)=V(28)+V(30)+QCx*VRR0(3,27)+WQx*VRR1(3,27)
      VRR0(3,44)=V(31)-V(32)+QCx*VRR0(3,28)+WQx*VRR1(3,28)
      VRR0(3,45)=V(13)+V(28)+V(30)+QCy*VRR0(3,28)+WQy*VRR1(3,28)
      VRR0(3,46)=3.D0*V(31)-3.D0*V(32)+QCy*VRR0(3,29)+HfxZpE*VRR1(1,29)+WQy*VRR1(3,29)
      VRR0(3,47)=2.D0*V(33)-2.D0*V(34)+QCx*VRR0(3,30)+WQx*VRR1(3,30)
      VRR0(3,48)=V(35)-V(36)+QCx*VRR0(3,31)+WQx*VRR1(3,31)
      VRR0(3,49)=V(19)+V(33)-V(34)+QCy*VRR0(3,31)+WQy*VRR1(3,31)
      VRR0(3,50)=2.D0*V(35)-2.D0*V(36)+V(37)+QCy*VRR0(3,32)+WQy*VRR1(3,32)
      VRR0(3,51)=2.D0*V(25)-2.D0*V(26)+QCz*VRR0(3,30)+WQz*VRR1(3,30)
      VRR0(3,52)=V(28)+V(30)+QCz*VRR0(3,31)+WQz*VRR1(3,31)
      VRR0(3,53)=2.D0*V(31)-2.D0*V(32)+QCz*VRR0(3,32)+WQz*VRR1(3,32)
      VRR0(3,54)=3.D0*V(33)-3.D0*V(34)+QCz*VRR0(3,33)+WQz*VRR1(3,33)
      VRR0(3,55)=3.D0*V(35)-3.D0*V(36)+QCz*VRR0(3,34)+WQz*VRR1(3,34)
      VRR0(3,56)=4.D0*r1x2E*VRR0(3,20)+QCz*VRR0(3,35)-4.D0*r1x2E*ZxZpE*VRR1(3,20)+WQz*VRR1(3,35)
      VRR0(4,36)=4.D0*r1x2E*VRR0(4,11)+QCx*VRR0(4,21)-4.D0*r1x2E*ZxZpE*VRR1(4,11)+WQx*VRR1(4,21)
      VRR0(4,37)=3.D0*V(38)-3.D0*V(39)+QCx*VRR0(4,22)+WQx*VRR1(4,22)
      VRR0(4,38)=2.D0*V(40)-2.D0*V(41)+QCx*VRR0(4,23)+WQx*VRR1(4,23)
      VRR0(4,39)=2.D0*V(38)-2.D0*V(39)+QCy*VRR0(4,23)+WQy*VRR1(4,23)
      VRR0(4,40)=3.D0*V(40)-3.D0*V(41)+QCy*VRR0(4,24)+WQy*VRR1(4,24)
      VRR0(4,41)=4.D0*r1x2E*VRR0(4,14)+QCy*VRR0(4,25)-4.D0*r1x2E*ZxZpE*VRR1(4,14)+WQy*VRR1(4,25)
      VRR0(4,42)=3.D0*V(42)-3.D0*V(43)+QCx*VRR0(4,26)+WQx*VRR1(4,26)
      VRR0(4,43)=V(45)+V(47)+QCx*VRR0(4,27)+WQx*VRR1(4,27)
      VRR0(4,44)=V(48)-V(49)+QCx*VRR0(4,28)+WQx*VRR1(4,28)
      VRR0(4,45)=V(45)+V(47)+QCy*VRR0(4,28)+WQy*VRR1(4,28)
      VRR0(4,46)=3.D0*V(48)-3.D0*V(49)+QCy*VRR0(4,29)+WQy*VRR1(4,29)
      VRR0(4,47)=2.D0*V(50)-2.D0*V(51)+QCx*VRR0(4,30)+WQx*VRR1(4,30)
      VRR0(4,48)=V(52)-V(53)+QCx*VRR0(4,31)+WQx*VRR1(4,31)
      VRR0(4,49)=V(50)-V(51)+QCy*VRR0(4,31)+WQy*VRR1(4,31)
      VRR0(4,50)=2.D0*V(52)-2.D0*V(53)+QCy*VRR0(4,32)+WQy*VRR1(4,32)
      VRR0(4,51)=V(16)+2.D0*V(42)-2.D0*V(43)+QCz*VRR0(4,30)+WQz*VRR1(4,30)
      VRR0(4,52)=V(19)+V(45)+V(47)+QCz*VRR0(4,31)+WQz*VRR1(4,31)
      VRR0(4,53)=V(37)+2.D0*V(48)-2.D0*V(49)+QCz*VRR0(4,32)+WQz*VRR1(4,32)
      VRR0(4,54)=3.D0*V(50)-3.D0*V(51)+QCz*VRR0(4,33)+HfxZpE*VRR1(1,33)+WQz*VRR1(4,33)
      VRR0(4,55)=3.D0*V(52)-3.D0*V(53)+QCz*VRR0(4,34)+HfxZpE*VRR1(1,34)+WQz*VRR1(4,34)
      VRR0(4,56)=4.D0*r1x2E*VRR0(4,20)+QCz*VRR0(4,35)+HfxZpE*VRR1(1,35)-4.D0*r1x2E*ZxZpE*VRR1(4,20)+WQz*VRR1(4,35)
END SUBROUTINE VRRp0h0
SUBROUTINE MVRRp0h0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(2,36)=QCx*VS0(2,21)+WQx*VS1(2,21)-r1x2E*VR1(2,21)&
   +4D0*r1x2E*(VS0(2,11)-ZxZpE*VS1(2,11))&
   +HfxZpE*VS1(1,21)
VS0(2,37)=QCx*VS0(2,22)+WQx*VS1(2,22)-r1x2E*VR1(2,22)&
   +3D0*r1x2E*(VS0(2,12)-ZxZpE*VS1(2,12))&
   +HfxZpE*VS1(1,22)
VS0(2,38)=QCx*VS0(2,23)+WQx*VS1(2,23)-r1x2E*VR1(2,23)&
   +2D0*r1x2E*(VS0(2,13)-ZxZpE*VS1(2,13))&
   +HfxZpE*VS1(1,23)
VS0(2,39)=QCy*VS0(2,23)+WQy*VS1(2,23)&
   +2D0*r1x2E*(VS0(2,12)-ZxZpE*VS1(2,12))
VS0(2,40)=QCy*VS0(2,24)+WQy*VS1(2,24)&
   +3D0*r1x2E*(VS0(2,13)-ZxZpE*VS1(2,13))
VS0(2,41)=QCy*VS0(2,25)+WQy*VS1(2,25)&
   +4D0*r1x2E*(VS0(2,14)-ZxZpE*VS1(2,14))
VS0(2,42)=QCx*VS0(2,26)+WQx*VS1(2,26)-r1x2E*VR1(2,26)&
   +3D0*r1x2E*(VS0(2,15)-ZxZpE*VS1(2,15))&
   +HfxZpE*VS1(1,26)
VS0(2,43)=QCx*VS0(2,27)+WQx*VS1(2,27)-r1x2E*VR1(2,27)&
   +2D0*r1x2E*(VS0(2,16)-ZxZpE*VS1(2,16))&
   +HfxZpE*VS1(1,27)
VS0(2,44)=QCx*VS0(2,28)+WQx*VS1(2,28)-r1x2E*VR1(2,28)&
   +r1x2E*(VS0(2,17)-ZxZpE*VS1(2,17))&
   +HfxZpE*VS1(1,28)
VS0(2,45)=QCy*VS0(2,28)+WQy*VS1(2,28)&
   +2D0*r1x2E*(VS0(2,16)-ZxZpE*VS1(2,16))
VS0(2,46)=QCy*VS0(2,29)+WQy*VS1(2,29)&
   +3D0*r1x2E*(VS0(2,17)-ZxZpE*VS1(2,17))
VS0(2,47)=QCx*VS0(2,30)+WQx*VS1(2,30)-r1x2E*VR1(2,30)&
   +2D0*r1x2E*(VS0(2,18)-ZxZpE*VS1(2,18))&
   +HfxZpE*VS1(1,30)
VS0(2,48)=QCx*VS0(2,31)+WQx*VS1(2,31)-r1x2E*VR1(2,31)&
   +r1x2E*(VS0(2,19)-ZxZpE*VS1(2,19))&
   +HfxZpE*VS1(1,31)
VS0(2,49)=QCy*VS0(2,31)+WQy*VS1(2,31)&
   +r1x2E*(VS0(2,18)-ZxZpE*VS1(2,18))
VS0(2,50)=QCy*VS0(2,32)+WQy*VS1(2,32)&
   +2D0*r1x2E*(VS0(2,19)-ZxZpE*VS1(2,19))
VS0(2,51)=QCz*VS0(2,30)+WQz*VS1(2,30)&
   +2D0*r1x2E*(VS0(2,15)-ZxZpE*VS1(2,15))
VS0(2,52)=QCz*VS0(2,31)+WQz*VS1(2,31)&
   +2D0*r1x2E*(VS0(2,16)-ZxZpE*VS1(2,16))
VS0(2,53)=QCz*VS0(2,32)+WQz*VS1(2,32)&
   +2D0*r1x2E*(VS0(2,17)-ZxZpE*VS1(2,17))
VS0(2,54)=QCz*VS0(2,33)+WQz*VS1(2,33)&
   +3D0*r1x2E*(VS0(2,18)-ZxZpE*VS1(2,18))
VS0(2,55)=QCz*VS0(2,34)+WQz*VS1(2,34)&
   +3D0*r1x2E*(VS0(2,19)-ZxZpE*VS1(2,19))
VS0(2,56)=QCz*VS0(2,35)+WQz*VS1(2,35)&
   +4D0*r1x2E*(VS0(2,20)-ZxZpE*VS1(2,20))
VS0(3,36)=QCx*VS0(3,21)+WQx*VS1(3,21)-r1x2E*VR1(3,21)&
   +4D0*r1x2E*(VS0(3,11)-ZxZpE*VS1(3,11))
VS0(3,37)=QCx*VS0(3,22)+WQx*VS1(3,22)-r1x2E*VR1(3,22)&
   +3D0*r1x2E*(VS0(3,12)-ZxZpE*VS1(3,12))
VS0(3,38)=QCx*VS0(3,23)+WQx*VS1(3,23)-r1x2E*VR1(3,23)&
   +2D0*r1x2E*(VS0(3,13)-ZxZpE*VS1(3,13))
VS0(3,39)=QCy*VS0(3,23)+WQy*VS1(3,23)&
   +2D0*r1x2E*(VS0(3,12)-ZxZpE*VS1(3,12))&
   +HfxZpE*VS1(1,23)
VS0(3,40)=QCy*VS0(3,24)+WQy*VS1(3,24)&
   +3D0*r1x2E*(VS0(3,13)-ZxZpE*VS1(3,13))&
   +HfxZpE*VS1(1,24)
VS0(3,41)=QCy*VS0(3,25)+WQy*VS1(3,25)&
   +4D0*r1x2E*(VS0(3,14)-ZxZpE*VS1(3,14))&
   +HfxZpE*VS1(1,25)
VS0(3,42)=QCx*VS0(3,26)+WQx*VS1(3,26)-r1x2E*VR1(3,26)&
   +3D0*r1x2E*(VS0(3,15)-ZxZpE*VS1(3,15))
VS0(3,43)=QCx*VS0(3,27)+WQx*VS1(3,27)-r1x2E*VR1(3,27)&
   +2D0*r1x2E*(VS0(3,16)-ZxZpE*VS1(3,16))
VS0(3,44)=QCx*VS0(3,28)+WQx*VS1(3,28)-r1x2E*VR1(3,28)&
   +r1x2E*(VS0(3,17)-ZxZpE*VS1(3,17))
VS0(3,45)=QCy*VS0(3,28)+WQy*VS1(3,28)&
   +2D0*r1x2E*(VS0(3,16)-ZxZpE*VS1(3,16))&
   +HfxZpE*VS1(1,28)
VS0(3,46)=QCy*VS0(3,29)+WQy*VS1(3,29)&
   +3D0*r1x2E*(VS0(3,17)-ZxZpE*VS1(3,17))&
   +HfxZpE*VS1(1,29)
VS0(3,47)=QCx*VS0(3,30)+WQx*VS1(3,30)-r1x2E*VR1(3,30)&
   +2D0*r1x2E*(VS0(3,18)-ZxZpE*VS1(3,18))
VS0(3,48)=QCx*VS0(3,31)+WQx*VS1(3,31)-r1x2E*VR1(3,31)&
   +r1x2E*(VS0(3,19)-ZxZpE*VS1(3,19))
VS0(3,49)=QCy*VS0(3,31)+WQy*VS1(3,31)&
   +r1x2E*(VS0(3,18)-ZxZpE*VS1(3,18))&
   +HfxZpE*VS1(1,31)
VS0(3,50)=QCy*VS0(3,32)+WQy*VS1(3,32)&
   +2D0*r1x2E*(VS0(3,19)-ZxZpE*VS1(3,19))&
   +HfxZpE*VS1(1,32)
VS0(3,51)=QCz*VS0(3,30)+WQz*VS1(3,30)&
   +2D0*r1x2E*(VS0(3,15)-ZxZpE*VS1(3,15))
VS0(3,52)=QCz*VS0(3,31)+WQz*VS1(3,31)&
   +2D0*r1x2E*(VS0(3,16)-ZxZpE*VS1(3,16))
VS0(3,53)=QCz*VS0(3,32)+WQz*VS1(3,32)&
   +2D0*r1x2E*(VS0(3,17)-ZxZpE*VS1(3,17))
VS0(3,54)=QCz*VS0(3,33)+WQz*VS1(3,33)&
   +3D0*r1x2E*(VS0(3,18)-ZxZpE*VS1(3,18))
VS0(3,55)=QCz*VS0(3,34)+WQz*VS1(3,34)&
   +3D0*r1x2E*(VS0(3,19)-ZxZpE*VS1(3,19))
VS0(3,56)=QCz*VS0(3,35)+WQz*VS1(3,35)&
   +4D0*r1x2E*(VS0(3,20)-ZxZpE*VS1(3,20))
VS0(4,36)=QCx*VS0(4,21)+WQx*VS1(4,21)-r1x2E*VR1(4,21)&
   +4D0*r1x2E*(VS0(4,11)-ZxZpE*VS1(4,11))
VS0(4,37)=QCx*VS0(4,22)+WQx*VS1(4,22)-r1x2E*VR1(4,22)&
   +3D0*r1x2E*(VS0(4,12)-ZxZpE*VS1(4,12))
VS0(4,38)=QCx*VS0(4,23)+WQx*VS1(4,23)-r1x2E*VR1(4,23)&
   +2D0*r1x2E*(VS0(4,13)-ZxZpE*VS1(4,13))
VS0(4,39)=QCy*VS0(4,23)+WQy*VS1(4,23)&
   +2D0*r1x2E*(VS0(4,12)-ZxZpE*VS1(4,12))
VS0(4,40)=QCy*VS0(4,24)+WQy*VS1(4,24)&
   +3D0*r1x2E*(VS0(4,13)-ZxZpE*VS1(4,13))
VS0(4,41)=QCy*VS0(4,25)+WQy*VS1(4,25)&
   +4D0*r1x2E*(VS0(4,14)-ZxZpE*VS1(4,14))
VS0(4,42)=QCx*VS0(4,26)+WQx*VS1(4,26)-r1x2E*VR1(4,26)&
   +3D0*r1x2E*(VS0(4,15)-ZxZpE*VS1(4,15))
VS0(4,43)=QCx*VS0(4,27)+WQx*VS1(4,27)-r1x2E*VR1(4,27)&
   +2D0*r1x2E*(VS0(4,16)-ZxZpE*VS1(4,16))
VS0(4,44)=QCx*VS0(4,28)+WQx*VS1(4,28)-r1x2E*VR1(4,28)&
   +r1x2E*(VS0(4,17)-ZxZpE*VS1(4,17))
VS0(4,45)=QCy*VS0(4,28)+WQy*VS1(4,28)&
   +2D0*r1x2E*(VS0(4,16)-ZxZpE*VS1(4,16))
VS0(4,46)=QCy*VS0(4,29)+WQy*VS1(4,29)&
   +3D0*r1x2E*(VS0(4,17)-ZxZpE*VS1(4,17))
VS0(4,47)=QCx*VS0(4,30)+WQx*VS1(4,30)-r1x2E*VR1(4,30)&
   +2D0*r1x2E*(VS0(4,18)-ZxZpE*VS1(4,18))
VS0(4,48)=QCx*VS0(4,31)+WQx*VS1(4,31)-r1x2E*VR1(4,31)&
   +r1x2E*(VS0(4,19)-ZxZpE*VS1(4,19))
VS0(4,49)=QCy*VS0(4,31)+WQy*VS1(4,31)&
   +r1x2E*(VS0(4,18)-ZxZpE*VS1(4,18))
VS0(4,50)=QCy*VS0(4,32)+WQy*VS1(4,32)&
   +2D0*r1x2E*(VS0(4,19)-ZxZpE*VS1(4,19))
VS0(4,51)=QCz*VS0(4,30)+WQz*VS1(4,30)&
   +2D0*r1x2E*(VS0(4,15)-ZxZpE*VS1(4,15))&
   +HfxZpE*VS1(1,30)
VS0(4,52)=QCz*VS0(4,31)+WQz*VS1(4,31)&
   +2D0*r1x2E*(VS0(4,16)-ZxZpE*VS1(4,16))&
   +HfxZpE*VS1(1,31)
VS0(4,53)=QCz*VS0(4,32)+WQz*VS1(4,32)&
   +2D0*r1x2E*(VS0(4,17)-ZxZpE*VS1(4,17))&
   +HfxZpE*VS1(1,32)
VS0(4,54)=QCz*VS0(4,33)+WQz*VS1(4,33)&
   +3D0*r1x2E*(VS0(4,18)-ZxZpE*VS1(4,18))&
   +HfxZpE*VS1(1,33)
VS0(4,55)=QCz*VS0(4,34)+WQz*VS1(4,34)&
   +3D0*r1x2E*(VS0(4,19)-ZxZpE*VS1(4,19))&
   +HfxZpE*VS1(1,34)
VS0(4,56)=QCz*VS0(4,35)+WQz*VS1(4,35)&
   +4D0*r1x2E*(VS0(4,20)-ZxZpE*VS1(4,20))&
   +HfxZpE*VS1(1,35)
CASE(2)
VS0(2,36)=QCx*VS0(2,21)+WQx*VS1(2,21)&
   +4D0*r1x2E*(VS0(2,11)-ZxZpE*VS1(2,11))&
   +HfxZpE*VS1(1,21)
VS0(2,37)=QCx*VS0(2,22)+WQx*VS1(2,22)&
   +3D0*r1x2E*(VS0(2,12)-ZxZpE*VS1(2,12))&
   +HfxZpE*VS1(1,22)
VS0(2,38)=QCx*VS0(2,23)+WQx*VS1(2,23)&
   +2D0*r1x2E*(VS0(2,13)-ZxZpE*VS1(2,13))&
   +HfxZpE*VS1(1,23)
VS0(2,39)=QCy*VS0(2,23)+WQy*VS1(2,23)-r1x2E*VR1(2,23)&
   +2D0*r1x2E*(VS0(2,12)-ZxZpE*VS1(2,12))
VS0(2,40)=QCy*VS0(2,24)+WQy*VS1(2,24)-r1x2E*VR1(2,24)&
   +3D0*r1x2E*(VS0(2,13)-ZxZpE*VS1(2,13))
VS0(2,41)=QCy*VS0(2,25)+WQy*VS1(2,25)-r1x2E*VR1(2,25)&
   +4D0*r1x2E*(VS0(2,14)-ZxZpE*VS1(2,14))
VS0(2,42)=QCx*VS0(2,26)+WQx*VS1(2,26)&
   +3D0*r1x2E*(VS0(2,15)-ZxZpE*VS1(2,15))&
   +HfxZpE*VS1(1,26)
VS0(2,43)=QCx*VS0(2,27)+WQx*VS1(2,27)&
   +2D0*r1x2E*(VS0(2,16)-ZxZpE*VS1(2,16))&
   +HfxZpE*VS1(1,27)
VS0(2,44)=QCx*VS0(2,28)+WQx*VS1(2,28)&
   +r1x2E*(VS0(2,17)-ZxZpE*VS1(2,17))&
   +HfxZpE*VS1(1,28)
VS0(2,45)=QCy*VS0(2,28)+WQy*VS1(2,28)-r1x2E*VR1(2,28)&
   +2D0*r1x2E*(VS0(2,16)-ZxZpE*VS1(2,16))
VS0(2,46)=QCy*VS0(2,29)+WQy*VS1(2,29)-r1x2E*VR1(2,29)&
   +3D0*r1x2E*(VS0(2,17)-ZxZpE*VS1(2,17))
VS0(2,47)=QCx*VS0(2,30)+WQx*VS1(2,30)&
   +2D0*r1x2E*(VS0(2,18)-ZxZpE*VS1(2,18))&
   +HfxZpE*VS1(1,30)
VS0(2,48)=QCx*VS0(2,31)+WQx*VS1(2,31)&
   +r1x2E*(VS0(2,19)-ZxZpE*VS1(2,19))&
   +HfxZpE*VS1(1,31)
VS0(2,49)=QCy*VS0(2,31)+WQy*VS1(2,31)-r1x2E*VR1(2,31)&
   +r1x2E*(VS0(2,18)-ZxZpE*VS1(2,18))
VS0(2,50)=QCy*VS0(2,32)+WQy*VS1(2,32)-r1x2E*VR1(2,32)&
   +2D0*r1x2E*(VS0(2,19)-ZxZpE*VS1(2,19))
VS0(2,51)=QCz*VS0(2,30)+WQz*VS1(2,30)&
   +2D0*r1x2E*(VS0(2,15)-ZxZpE*VS1(2,15))
VS0(2,52)=QCz*VS0(2,31)+WQz*VS1(2,31)&
   +2D0*r1x2E*(VS0(2,16)-ZxZpE*VS1(2,16))
VS0(2,53)=QCz*VS0(2,32)+WQz*VS1(2,32)&
   +2D0*r1x2E*(VS0(2,17)-ZxZpE*VS1(2,17))
VS0(2,54)=QCz*VS0(2,33)+WQz*VS1(2,33)&
   +3D0*r1x2E*(VS0(2,18)-ZxZpE*VS1(2,18))
VS0(2,55)=QCz*VS0(2,34)+WQz*VS1(2,34)&
   +3D0*r1x2E*(VS0(2,19)-ZxZpE*VS1(2,19))
VS0(2,56)=QCz*VS0(2,35)+WQz*VS1(2,35)&
   +4D0*r1x2E*(VS0(2,20)-ZxZpE*VS1(2,20))
VS0(3,36)=QCx*VS0(3,21)+WQx*VS1(3,21)&
   +4D0*r1x2E*(VS0(3,11)-ZxZpE*VS1(3,11))
VS0(3,37)=QCx*VS0(3,22)+WQx*VS1(3,22)&
   +3D0*r1x2E*(VS0(3,12)-ZxZpE*VS1(3,12))
VS0(3,38)=QCx*VS0(3,23)+WQx*VS1(3,23)&
   +2D0*r1x2E*(VS0(3,13)-ZxZpE*VS1(3,13))
VS0(3,39)=QCy*VS0(3,23)+WQy*VS1(3,23)-r1x2E*VR1(3,23)&
   +2D0*r1x2E*(VS0(3,12)-ZxZpE*VS1(3,12))&
   +HfxZpE*VS1(1,23)
VS0(3,40)=QCy*VS0(3,24)+WQy*VS1(3,24)-r1x2E*VR1(3,24)&
   +3D0*r1x2E*(VS0(3,13)-ZxZpE*VS1(3,13))&
   +HfxZpE*VS1(1,24)
VS0(3,41)=QCy*VS0(3,25)+WQy*VS1(3,25)-r1x2E*VR1(3,25)&
   +4D0*r1x2E*(VS0(3,14)-ZxZpE*VS1(3,14))&
   +HfxZpE*VS1(1,25)
VS0(3,42)=QCx*VS0(3,26)+WQx*VS1(3,26)&
   +3D0*r1x2E*(VS0(3,15)-ZxZpE*VS1(3,15))
VS0(3,43)=QCx*VS0(3,27)+WQx*VS1(3,27)&
   +2D0*r1x2E*(VS0(3,16)-ZxZpE*VS1(3,16))
VS0(3,44)=QCx*VS0(3,28)+WQx*VS1(3,28)&
   +r1x2E*(VS0(3,17)-ZxZpE*VS1(3,17))
VS0(3,45)=QCy*VS0(3,28)+WQy*VS1(3,28)-r1x2E*VR1(3,28)&
   +2D0*r1x2E*(VS0(3,16)-ZxZpE*VS1(3,16))&
   +HfxZpE*VS1(1,28)
VS0(3,46)=QCy*VS0(3,29)+WQy*VS1(3,29)-r1x2E*VR1(3,29)&
   +3D0*r1x2E*(VS0(3,17)-ZxZpE*VS1(3,17))&
   +HfxZpE*VS1(1,29)
VS0(3,47)=QCx*VS0(3,30)+WQx*VS1(3,30)&
   +2D0*r1x2E*(VS0(3,18)-ZxZpE*VS1(3,18))
VS0(3,48)=QCx*VS0(3,31)+WQx*VS1(3,31)&
   +r1x2E*(VS0(3,19)-ZxZpE*VS1(3,19))
VS0(3,49)=QCy*VS0(3,31)+WQy*VS1(3,31)-r1x2E*VR1(3,31)&
   +r1x2E*(VS0(3,18)-ZxZpE*VS1(3,18))&
   +HfxZpE*VS1(1,31)
VS0(3,50)=QCy*VS0(3,32)+WQy*VS1(3,32)-r1x2E*VR1(3,32)&
   +2D0*r1x2E*(VS0(3,19)-ZxZpE*VS1(3,19))&
   +HfxZpE*VS1(1,32)
VS0(3,51)=QCz*VS0(3,30)+WQz*VS1(3,30)&
   +2D0*r1x2E*(VS0(3,15)-ZxZpE*VS1(3,15))
VS0(3,52)=QCz*VS0(3,31)+WQz*VS1(3,31)&
   +2D0*r1x2E*(VS0(3,16)-ZxZpE*VS1(3,16))
VS0(3,53)=QCz*VS0(3,32)+WQz*VS1(3,32)&
   +2D0*r1x2E*(VS0(3,17)-ZxZpE*VS1(3,17))
VS0(3,54)=QCz*VS0(3,33)+WQz*VS1(3,33)&
   +3D0*r1x2E*(VS0(3,18)-ZxZpE*VS1(3,18))
VS0(3,55)=QCz*VS0(3,34)+WQz*VS1(3,34)&
   +3D0*r1x2E*(VS0(3,19)-ZxZpE*VS1(3,19))
VS0(3,56)=QCz*VS0(3,35)+WQz*VS1(3,35)&
   +4D0*r1x2E*(VS0(3,20)-ZxZpE*VS1(3,20))
VS0(4,36)=QCx*VS0(4,21)+WQx*VS1(4,21)&
   +4D0*r1x2E*(VS0(4,11)-ZxZpE*VS1(4,11))
VS0(4,37)=QCx*VS0(4,22)+WQx*VS1(4,22)&
   +3D0*r1x2E*(VS0(4,12)-ZxZpE*VS1(4,12))
VS0(4,38)=QCx*VS0(4,23)+WQx*VS1(4,23)&
   +2D0*r1x2E*(VS0(4,13)-ZxZpE*VS1(4,13))
VS0(4,39)=QCy*VS0(4,23)+WQy*VS1(4,23)-r1x2E*VR1(4,23)&
   +2D0*r1x2E*(VS0(4,12)-ZxZpE*VS1(4,12))
VS0(4,40)=QCy*VS0(4,24)+WQy*VS1(4,24)-r1x2E*VR1(4,24)&
   +3D0*r1x2E*(VS0(4,13)-ZxZpE*VS1(4,13))
VS0(4,41)=QCy*VS0(4,25)+WQy*VS1(4,25)-r1x2E*VR1(4,25)&
   +4D0*r1x2E*(VS0(4,14)-ZxZpE*VS1(4,14))
VS0(4,42)=QCx*VS0(4,26)+WQx*VS1(4,26)&
   +3D0*r1x2E*(VS0(4,15)-ZxZpE*VS1(4,15))
VS0(4,43)=QCx*VS0(4,27)+WQx*VS1(4,27)&
   +2D0*r1x2E*(VS0(4,16)-ZxZpE*VS1(4,16))
VS0(4,44)=QCx*VS0(4,28)+WQx*VS1(4,28)&
   +r1x2E*(VS0(4,17)-ZxZpE*VS1(4,17))
VS0(4,45)=QCy*VS0(4,28)+WQy*VS1(4,28)-r1x2E*VR1(4,28)&
   +2D0*r1x2E*(VS0(4,16)-ZxZpE*VS1(4,16))
VS0(4,46)=QCy*VS0(4,29)+WQy*VS1(4,29)-r1x2E*VR1(4,29)&
   +3D0*r1x2E*(VS0(4,17)-ZxZpE*VS1(4,17))
VS0(4,47)=QCx*VS0(4,30)+WQx*VS1(4,30)&
   +2D0*r1x2E*(VS0(4,18)-ZxZpE*VS1(4,18))
VS0(4,48)=QCx*VS0(4,31)+WQx*VS1(4,31)&
   +r1x2E*(VS0(4,19)-ZxZpE*VS1(4,19))
VS0(4,49)=QCy*VS0(4,31)+WQy*VS1(4,31)-r1x2E*VR1(4,31)&
   +r1x2E*(VS0(4,18)-ZxZpE*VS1(4,18))
VS0(4,50)=QCy*VS0(4,32)+WQy*VS1(4,32)-r1x2E*VR1(4,32)&
   +2D0*r1x2E*(VS0(4,19)-ZxZpE*VS1(4,19))
VS0(4,51)=QCz*VS0(4,30)+WQz*VS1(4,30)&
   +2D0*r1x2E*(VS0(4,15)-ZxZpE*VS1(4,15))&
   +HfxZpE*VS1(1,30)
VS0(4,52)=QCz*VS0(4,31)+WQz*VS1(4,31)&
   +2D0*r1x2E*(VS0(4,16)-ZxZpE*VS1(4,16))&
   +HfxZpE*VS1(1,31)
VS0(4,53)=QCz*VS0(4,32)+WQz*VS1(4,32)&
   +2D0*r1x2E*(VS0(4,17)-ZxZpE*VS1(4,17))&
   +HfxZpE*VS1(1,32)
VS0(4,54)=QCz*VS0(4,33)+WQz*VS1(4,33)&
   +3D0*r1x2E*(VS0(4,18)-ZxZpE*VS1(4,18))&
   +HfxZpE*VS1(1,33)
VS0(4,55)=QCz*VS0(4,34)+WQz*VS1(4,34)&
   +3D0*r1x2E*(VS0(4,19)-ZxZpE*VS1(4,19))&
   +HfxZpE*VS1(1,34)
VS0(4,56)=QCz*VS0(4,35)+WQz*VS1(4,35)&
   +4D0*r1x2E*(VS0(4,20)-ZxZpE*VS1(4,20))&
   +HfxZpE*VS1(1,35)
CASE(3)
VS0(2,36)=QCx*VS0(2,21)+WQx*VS1(2,21)&
   +4D0*r1x2E*(VS0(2,11)-ZxZpE*VS1(2,11))&
   +HfxZpE*VS1(1,21)
VS0(2,37)=QCx*VS0(2,22)+WQx*VS1(2,22)&
   +3D0*r1x2E*(VS0(2,12)-ZxZpE*VS1(2,12))&
   +HfxZpE*VS1(1,22)
VS0(2,38)=QCx*VS0(2,23)+WQx*VS1(2,23)&
   +2D0*r1x2E*(VS0(2,13)-ZxZpE*VS1(2,13))&
   +HfxZpE*VS1(1,23)
VS0(2,39)=QCy*VS0(2,23)+WQy*VS1(2,23)&
   +2D0*r1x2E*(VS0(2,12)-ZxZpE*VS1(2,12))
VS0(2,40)=QCy*VS0(2,24)+WQy*VS1(2,24)&
   +3D0*r1x2E*(VS0(2,13)-ZxZpE*VS1(2,13))
VS0(2,41)=QCy*VS0(2,25)+WQy*VS1(2,25)&
   +4D0*r1x2E*(VS0(2,14)-ZxZpE*VS1(2,14))
VS0(2,42)=QCx*VS0(2,26)+WQx*VS1(2,26)&
   +3D0*r1x2E*(VS0(2,15)-ZxZpE*VS1(2,15))&
   +HfxZpE*VS1(1,26)
VS0(2,43)=QCx*VS0(2,27)+WQx*VS1(2,27)&
   +2D0*r1x2E*(VS0(2,16)-ZxZpE*VS1(2,16))&
   +HfxZpE*VS1(1,27)
VS0(2,44)=QCx*VS0(2,28)+WQx*VS1(2,28)&
   +r1x2E*(VS0(2,17)-ZxZpE*VS1(2,17))&
   +HfxZpE*VS1(1,28)
VS0(2,45)=QCy*VS0(2,28)+WQy*VS1(2,28)&
   +2D0*r1x2E*(VS0(2,16)-ZxZpE*VS1(2,16))
VS0(2,46)=QCy*VS0(2,29)+WQy*VS1(2,29)&
   +3D0*r1x2E*(VS0(2,17)-ZxZpE*VS1(2,17))
VS0(2,47)=QCx*VS0(2,30)+WQx*VS1(2,30)&
   +2D0*r1x2E*(VS0(2,18)-ZxZpE*VS1(2,18))&
   +HfxZpE*VS1(1,30)
VS0(2,48)=QCx*VS0(2,31)+WQx*VS1(2,31)&
   +r1x2E*(VS0(2,19)-ZxZpE*VS1(2,19))&
   +HfxZpE*VS1(1,31)
VS0(2,49)=QCy*VS0(2,31)+WQy*VS1(2,31)&
   +r1x2E*(VS0(2,18)-ZxZpE*VS1(2,18))
VS0(2,50)=QCy*VS0(2,32)+WQy*VS1(2,32)&
   +2D0*r1x2E*(VS0(2,19)-ZxZpE*VS1(2,19))
VS0(2,51)=QCz*VS0(2,30)+WQz*VS1(2,30)-r1x2E*VR1(2,30)&
   +2D0*r1x2E*(VS0(2,15)-ZxZpE*VS1(2,15))
VS0(2,52)=QCz*VS0(2,31)+WQz*VS1(2,31)-r1x2E*VR1(2,31)&
   +2D0*r1x2E*(VS0(2,16)-ZxZpE*VS1(2,16))
VS0(2,53)=QCz*VS0(2,32)+WQz*VS1(2,32)-r1x2E*VR1(2,32)&
   +2D0*r1x2E*(VS0(2,17)-ZxZpE*VS1(2,17))
VS0(2,54)=QCz*VS0(2,33)+WQz*VS1(2,33)-r1x2E*VR1(2,33)&
   +3D0*r1x2E*(VS0(2,18)-ZxZpE*VS1(2,18))
VS0(2,55)=QCz*VS0(2,34)+WQz*VS1(2,34)-r1x2E*VR1(2,34)&
   +3D0*r1x2E*(VS0(2,19)-ZxZpE*VS1(2,19))
VS0(2,56)=QCz*VS0(2,35)+WQz*VS1(2,35)-r1x2E*VR1(2,35)&
   +4D0*r1x2E*(VS0(2,20)-ZxZpE*VS1(2,20))
VS0(3,36)=QCx*VS0(3,21)+WQx*VS1(3,21)&
   +4D0*r1x2E*(VS0(3,11)-ZxZpE*VS1(3,11))
VS0(3,37)=QCx*VS0(3,22)+WQx*VS1(3,22)&
   +3D0*r1x2E*(VS0(3,12)-ZxZpE*VS1(3,12))
VS0(3,38)=QCx*VS0(3,23)+WQx*VS1(3,23)&
   +2D0*r1x2E*(VS0(3,13)-ZxZpE*VS1(3,13))
VS0(3,39)=QCy*VS0(3,23)+WQy*VS1(3,23)&
   +2D0*r1x2E*(VS0(3,12)-ZxZpE*VS1(3,12))&
   +HfxZpE*VS1(1,23)
VS0(3,40)=QCy*VS0(3,24)+WQy*VS1(3,24)&
   +3D0*r1x2E*(VS0(3,13)-ZxZpE*VS1(3,13))&
   +HfxZpE*VS1(1,24)
VS0(3,41)=QCy*VS0(3,25)+WQy*VS1(3,25)&
   +4D0*r1x2E*(VS0(3,14)-ZxZpE*VS1(3,14))&
   +HfxZpE*VS1(1,25)
VS0(3,42)=QCx*VS0(3,26)+WQx*VS1(3,26)&
   +3D0*r1x2E*(VS0(3,15)-ZxZpE*VS1(3,15))
VS0(3,43)=QCx*VS0(3,27)+WQx*VS1(3,27)&
   +2D0*r1x2E*(VS0(3,16)-ZxZpE*VS1(3,16))
VS0(3,44)=QCx*VS0(3,28)+WQx*VS1(3,28)&
   +r1x2E*(VS0(3,17)-ZxZpE*VS1(3,17))
VS0(3,45)=QCy*VS0(3,28)+WQy*VS1(3,28)&
   +2D0*r1x2E*(VS0(3,16)-ZxZpE*VS1(3,16))&
   +HfxZpE*VS1(1,28)
VS0(3,46)=QCy*VS0(3,29)+WQy*VS1(3,29)&
   +3D0*r1x2E*(VS0(3,17)-ZxZpE*VS1(3,17))&
   +HfxZpE*VS1(1,29)
VS0(3,47)=QCx*VS0(3,30)+WQx*VS1(3,30)&
   +2D0*r1x2E*(VS0(3,18)-ZxZpE*VS1(3,18))
VS0(3,48)=QCx*VS0(3,31)+WQx*VS1(3,31)&
   +r1x2E*(VS0(3,19)-ZxZpE*VS1(3,19))
VS0(3,49)=QCy*VS0(3,31)+WQy*VS1(3,31)&
   +r1x2E*(VS0(3,18)-ZxZpE*VS1(3,18))&
   +HfxZpE*VS1(1,31)
VS0(3,50)=QCy*VS0(3,32)+WQy*VS1(3,32)&
   +2D0*r1x2E*(VS0(3,19)-ZxZpE*VS1(3,19))&
   +HfxZpE*VS1(1,32)
VS0(3,51)=QCz*VS0(3,30)+WQz*VS1(3,30)-r1x2E*VR1(3,30)&
   +2D0*r1x2E*(VS0(3,15)-ZxZpE*VS1(3,15))
VS0(3,52)=QCz*VS0(3,31)+WQz*VS1(3,31)-r1x2E*VR1(3,31)&
   +2D0*r1x2E*(VS0(3,16)-ZxZpE*VS1(3,16))
VS0(3,53)=QCz*VS0(3,32)+WQz*VS1(3,32)-r1x2E*VR1(3,32)&
   +2D0*r1x2E*(VS0(3,17)-ZxZpE*VS1(3,17))
VS0(3,54)=QCz*VS0(3,33)+WQz*VS1(3,33)-r1x2E*VR1(3,33)&
   +3D0*r1x2E*(VS0(3,18)-ZxZpE*VS1(3,18))
VS0(3,55)=QCz*VS0(3,34)+WQz*VS1(3,34)-r1x2E*VR1(3,34)&
   +3D0*r1x2E*(VS0(3,19)-ZxZpE*VS1(3,19))
VS0(3,56)=QCz*VS0(3,35)+WQz*VS1(3,35)-r1x2E*VR1(3,35)&
   +4D0*r1x2E*(VS0(3,20)-ZxZpE*VS1(3,20))
VS0(4,36)=QCx*VS0(4,21)+WQx*VS1(4,21)&
   +4D0*r1x2E*(VS0(4,11)-ZxZpE*VS1(4,11))
VS0(4,37)=QCx*VS0(4,22)+WQx*VS1(4,22)&
   +3D0*r1x2E*(VS0(4,12)-ZxZpE*VS1(4,12))
VS0(4,38)=QCx*VS0(4,23)+WQx*VS1(4,23)&
   +2D0*r1x2E*(VS0(4,13)-ZxZpE*VS1(4,13))
VS0(4,39)=QCy*VS0(4,23)+WQy*VS1(4,23)&
   +2D0*r1x2E*(VS0(4,12)-ZxZpE*VS1(4,12))
VS0(4,40)=QCy*VS0(4,24)+WQy*VS1(4,24)&
   +3D0*r1x2E*(VS0(4,13)-ZxZpE*VS1(4,13))
VS0(4,41)=QCy*VS0(4,25)+WQy*VS1(4,25)&
   +4D0*r1x2E*(VS0(4,14)-ZxZpE*VS1(4,14))
VS0(4,42)=QCx*VS0(4,26)+WQx*VS1(4,26)&
   +3D0*r1x2E*(VS0(4,15)-ZxZpE*VS1(4,15))
VS0(4,43)=QCx*VS0(4,27)+WQx*VS1(4,27)&
   +2D0*r1x2E*(VS0(4,16)-ZxZpE*VS1(4,16))
VS0(4,44)=QCx*VS0(4,28)+WQx*VS1(4,28)&
   +r1x2E*(VS0(4,17)-ZxZpE*VS1(4,17))
VS0(4,45)=QCy*VS0(4,28)+WQy*VS1(4,28)&
   +2D0*r1x2E*(VS0(4,16)-ZxZpE*VS1(4,16))
VS0(4,46)=QCy*VS0(4,29)+WQy*VS1(4,29)&
   +3D0*r1x2E*(VS0(4,17)-ZxZpE*VS1(4,17))
VS0(4,47)=QCx*VS0(4,30)+WQx*VS1(4,30)&
   +2D0*r1x2E*(VS0(4,18)-ZxZpE*VS1(4,18))
VS0(4,48)=QCx*VS0(4,31)+WQx*VS1(4,31)&
   +r1x2E*(VS0(4,19)-ZxZpE*VS1(4,19))
VS0(4,49)=QCy*VS0(4,31)+WQy*VS1(4,31)&
   +r1x2E*(VS0(4,18)-ZxZpE*VS1(4,18))
VS0(4,50)=QCy*VS0(4,32)+WQy*VS1(4,32)&
   +2D0*r1x2E*(VS0(4,19)-ZxZpE*VS1(4,19))
VS0(4,51)=QCz*VS0(4,30)+WQz*VS1(4,30)-r1x2E*VR1(4,30)&
   +2D0*r1x2E*(VS0(4,15)-ZxZpE*VS1(4,15))&
   +HfxZpE*VS1(1,30)
VS0(4,52)=QCz*VS0(4,31)+WQz*VS1(4,31)-r1x2E*VR1(4,31)&
   +2D0*r1x2E*(VS0(4,16)-ZxZpE*VS1(4,16))&
   +HfxZpE*VS1(1,31)
VS0(4,53)=QCz*VS0(4,32)+WQz*VS1(4,32)-r1x2E*VR1(4,32)&
   +2D0*r1x2E*(VS0(4,17)-ZxZpE*VS1(4,17))&
   +HfxZpE*VS1(1,32)
VS0(4,54)=QCz*VS0(4,33)+WQz*VS1(4,33)-r1x2E*VR1(4,33)&
   +3D0*r1x2E*(VS0(4,18)-ZxZpE*VS1(4,18))&
   +HfxZpE*VS1(1,33)
VS0(4,55)=QCz*VS0(4,34)+WQz*VS1(4,34)-r1x2E*VR1(4,34)&
   +3D0*r1x2E*(VS0(4,19)-ZxZpE*VS1(4,19))&
   +HfxZpE*VS1(1,34)
VS0(4,56)=QCz*VS0(4,35)+WQz*VS1(4,35)-r1x2E*VR1(4,35)&
   +4D0*r1x2E*(VS0(4,20)-ZxZpE*VS1(4,20))&
   +HfxZpE*VS1(1,35)
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRp0h0'
STOP
END SELECT
END SUBROUTINE MVRRp0h0
