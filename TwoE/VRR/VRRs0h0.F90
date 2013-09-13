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
   SUBROUTINE VRRs0h0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2E*VRR0(1,12)
      V(2)=r1x2E*ZxZpE*VRR1(1,12)
      V(3)=r1x2E*VRR0(1,13)
      V(4)=r1x2E*ZxZpE*VRR1(1,13)
      V(5)=r1x2E*VRR0(1,15)
      V(6)=r1x2E*ZxZpE*VRR1(1,15)
      V(7)=r1x2E*VRR0(1,16)
      V(8)=2.D0*V(7)
      V(9)=r1x2E*ZxZpE*VRR1(1,16)
      V(10)=-2.D0*V(9)
      V(11)=r1x2E*VRR0(1,17)
      V(12)=r1x2E*ZxZpE*VRR1(1,17)
      V(13)=r1x2E*VRR0(1,18)
      V(14)=r1x2E*ZxZpE*VRR1(1,18)
      V(15)=r1x2E*VRR0(1,19)
      V(16)=r1x2E*ZxZpE*VRR1(1,19)
      VRR0(1,36)=4.D0*r1x2E*VRR0(1,11)+QCx*VRR0(1,21)-4.D0*r1x2E*ZxZpE*VRR1(1,11)+WQx*VRR1(1,21)
      VRR0(1,37)=3.D0*V(1)-3.D0*V(2)+QCx*VRR0(1,22)+WQx*VRR1(1,22)
      VRR0(1,38)=2.D0*V(3)-2.D0*V(4)+QCx*VRR0(1,23)+WQx*VRR1(1,23)
      VRR0(1,39)=2.D0*V(1)-2.D0*V(2)+QCy*VRR0(1,23)+WQy*VRR1(1,23)
      VRR0(1,40)=3.D0*V(3)-3.D0*V(4)+QCy*VRR0(1,24)+WQy*VRR1(1,24)
      VRR0(1,41)=4.D0*r1x2E*VRR0(1,14)+QCy*VRR0(1,25)-4.D0*r1x2E*ZxZpE*VRR1(1,14)+WQy*VRR1(1,25)
      VRR0(1,42)=3.D0*V(5)-3.D0*V(6)+QCx*VRR0(1,26)+WQx*VRR1(1,26)
      VRR0(1,43)=V(8)+V(10)+QCx*VRR0(1,27)+WQx*VRR1(1,27)
      VRR0(1,44)=V(11)-V(12)+QCx*VRR0(1,28)+WQx*VRR1(1,28)
      VRR0(1,45)=V(8)+V(10)+QCy*VRR0(1,28)+WQy*VRR1(1,28)
      VRR0(1,46)=3.D0*V(11)-3.D0*V(12)+QCy*VRR0(1,29)+WQy*VRR1(1,29)
      VRR0(1,47)=2.D0*V(13)-2.D0*V(14)+QCx*VRR0(1,30)+WQx*VRR1(1,30)
      VRR0(1,48)=V(15)-V(16)+QCx*VRR0(1,31)+WQx*VRR1(1,31)
      VRR0(1,49)=V(13)-V(14)+QCy*VRR0(1,31)+WQy*VRR1(1,31)
      VRR0(1,50)=2.D0*V(15)-2.D0*V(16)+QCy*VRR0(1,32)+WQy*VRR1(1,32)
      VRR0(1,51)=2.D0*V(5)-2.D0*V(6)+QCz*VRR0(1,30)+WQz*VRR1(1,30)
      VRR0(1,52)=V(8)+V(10)+QCz*VRR0(1,31)+WQz*VRR1(1,31)
      VRR0(1,53)=2.D0*V(11)-2.D0*V(12)+QCz*VRR0(1,32)+WQz*VRR1(1,32)
      VRR0(1,54)=3.D0*V(13)-3.D0*V(14)+QCz*VRR0(1,33)+WQz*VRR1(1,33)
      VRR0(1,55)=3.D0*V(15)-3.D0*V(16)+QCz*VRR0(1,34)+WQz*VRR1(1,34)
      VRR0(1,56)=4.D0*r1x2E*VRR0(1,20)+QCz*VRR0(1,35)-4.D0*r1x2E*ZxZpE*VRR1(1,20)+WQz*VRR1(1,35)
END SUBROUTINE VRRs0h0
SUBROUTINE MVRRs0h0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(1,36)=QCx*VS0(1,21)+WQx*VS1(1,21)-r1x2E*VR1(1,21)&
   +4D0*r1x2E*(VS0(1,11)-ZxZpE*VS1(1,11))
VS0(1,37)=QCx*VS0(1,22)+WQx*VS1(1,22)-r1x2E*VR1(1,22)&
   +3D0*r1x2E*(VS0(1,12)-ZxZpE*VS1(1,12))
VS0(1,38)=QCx*VS0(1,23)+WQx*VS1(1,23)-r1x2E*VR1(1,23)&
   +2D0*r1x2E*(VS0(1,13)-ZxZpE*VS1(1,13))
VS0(1,39)=QCy*VS0(1,23)+WQy*VS1(1,23)&
   +2D0*r1x2E*(VS0(1,12)-ZxZpE*VS1(1,12))
VS0(1,40)=QCy*VS0(1,24)+WQy*VS1(1,24)&
   +3D0*r1x2E*(VS0(1,13)-ZxZpE*VS1(1,13))
VS0(1,41)=QCy*VS0(1,25)+WQy*VS1(1,25)&
   +4D0*r1x2E*(VS0(1,14)-ZxZpE*VS1(1,14))
VS0(1,42)=QCx*VS0(1,26)+WQx*VS1(1,26)-r1x2E*VR1(1,26)&
   +3D0*r1x2E*(VS0(1,15)-ZxZpE*VS1(1,15))
VS0(1,43)=QCx*VS0(1,27)+WQx*VS1(1,27)-r1x2E*VR1(1,27)&
   +2D0*r1x2E*(VS0(1,16)-ZxZpE*VS1(1,16))
VS0(1,44)=QCx*VS0(1,28)+WQx*VS1(1,28)-r1x2E*VR1(1,28)&
   +r1x2E*(VS0(1,17)-ZxZpE*VS1(1,17))
VS0(1,45)=QCy*VS0(1,28)+WQy*VS1(1,28)&
   +2D0*r1x2E*(VS0(1,16)-ZxZpE*VS1(1,16))
VS0(1,46)=QCy*VS0(1,29)+WQy*VS1(1,29)&
   +3D0*r1x2E*(VS0(1,17)-ZxZpE*VS1(1,17))
VS0(1,47)=QCx*VS0(1,30)+WQx*VS1(1,30)-r1x2E*VR1(1,30)&
   +2D0*r1x2E*(VS0(1,18)-ZxZpE*VS1(1,18))
VS0(1,48)=QCx*VS0(1,31)+WQx*VS1(1,31)-r1x2E*VR1(1,31)&
   +r1x2E*(VS0(1,19)-ZxZpE*VS1(1,19))
VS0(1,49)=QCy*VS0(1,31)+WQy*VS1(1,31)&
   +r1x2E*(VS0(1,18)-ZxZpE*VS1(1,18))
VS0(1,50)=QCy*VS0(1,32)+WQy*VS1(1,32)&
   +2D0*r1x2E*(VS0(1,19)-ZxZpE*VS1(1,19))
VS0(1,51)=QCz*VS0(1,30)+WQz*VS1(1,30)&
   +2D0*r1x2E*(VS0(1,15)-ZxZpE*VS1(1,15))
VS0(1,52)=QCz*VS0(1,31)+WQz*VS1(1,31)&
   +2D0*r1x2E*(VS0(1,16)-ZxZpE*VS1(1,16))
VS0(1,53)=QCz*VS0(1,32)+WQz*VS1(1,32)&
   +2D0*r1x2E*(VS0(1,17)-ZxZpE*VS1(1,17))
VS0(1,54)=QCz*VS0(1,33)+WQz*VS1(1,33)&
   +3D0*r1x2E*(VS0(1,18)-ZxZpE*VS1(1,18))
VS0(1,55)=QCz*VS0(1,34)+WQz*VS1(1,34)&
   +3D0*r1x2E*(VS0(1,19)-ZxZpE*VS1(1,19))
VS0(1,56)=QCz*VS0(1,35)+WQz*VS1(1,35)&
   +4D0*r1x2E*(VS0(1,20)-ZxZpE*VS1(1,20))
CASE(2)
VS0(1,36)=QCx*VS0(1,21)+WQx*VS1(1,21)&
   +4D0*r1x2E*(VS0(1,11)-ZxZpE*VS1(1,11))
VS0(1,37)=QCx*VS0(1,22)+WQx*VS1(1,22)&
   +3D0*r1x2E*(VS0(1,12)-ZxZpE*VS1(1,12))
VS0(1,38)=QCx*VS0(1,23)+WQx*VS1(1,23)&
   +2D0*r1x2E*(VS0(1,13)-ZxZpE*VS1(1,13))
VS0(1,39)=QCy*VS0(1,23)+WQy*VS1(1,23)-r1x2E*VR1(1,23)&
   +2D0*r1x2E*(VS0(1,12)-ZxZpE*VS1(1,12))
VS0(1,40)=QCy*VS0(1,24)+WQy*VS1(1,24)-r1x2E*VR1(1,24)&
   +3D0*r1x2E*(VS0(1,13)-ZxZpE*VS1(1,13))
VS0(1,41)=QCy*VS0(1,25)+WQy*VS1(1,25)-r1x2E*VR1(1,25)&
   +4D0*r1x2E*(VS0(1,14)-ZxZpE*VS1(1,14))
VS0(1,42)=QCx*VS0(1,26)+WQx*VS1(1,26)&
   +3D0*r1x2E*(VS0(1,15)-ZxZpE*VS1(1,15))
VS0(1,43)=QCx*VS0(1,27)+WQx*VS1(1,27)&
   +2D0*r1x2E*(VS0(1,16)-ZxZpE*VS1(1,16))
VS0(1,44)=QCx*VS0(1,28)+WQx*VS1(1,28)&
   +r1x2E*(VS0(1,17)-ZxZpE*VS1(1,17))
VS0(1,45)=QCy*VS0(1,28)+WQy*VS1(1,28)-r1x2E*VR1(1,28)&
   +2D0*r1x2E*(VS0(1,16)-ZxZpE*VS1(1,16))
VS0(1,46)=QCy*VS0(1,29)+WQy*VS1(1,29)-r1x2E*VR1(1,29)&
   +3D0*r1x2E*(VS0(1,17)-ZxZpE*VS1(1,17))
VS0(1,47)=QCx*VS0(1,30)+WQx*VS1(1,30)&
   +2D0*r1x2E*(VS0(1,18)-ZxZpE*VS1(1,18))
VS0(1,48)=QCx*VS0(1,31)+WQx*VS1(1,31)&
   +r1x2E*(VS0(1,19)-ZxZpE*VS1(1,19))
VS0(1,49)=QCy*VS0(1,31)+WQy*VS1(1,31)-r1x2E*VR1(1,31)&
   +r1x2E*(VS0(1,18)-ZxZpE*VS1(1,18))
VS0(1,50)=QCy*VS0(1,32)+WQy*VS1(1,32)-r1x2E*VR1(1,32)&
   +2D0*r1x2E*(VS0(1,19)-ZxZpE*VS1(1,19))
VS0(1,51)=QCz*VS0(1,30)+WQz*VS1(1,30)&
   +2D0*r1x2E*(VS0(1,15)-ZxZpE*VS1(1,15))
VS0(1,52)=QCz*VS0(1,31)+WQz*VS1(1,31)&
   +2D0*r1x2E*(VS0(1,16)-ZxZpE*VS1(1,16))
VS0(1,53)=QCz*VS0(1,32)+WQz*VS1(1,32)&
   +2D0*r1x2E*(VS0(1,17)-ZxZpE*VS1(1,17))
VS0(1,54)=QCz*VS0(1,33)+WQz*VS1(1,33)&
   +3D0*r1x2E*(VS0(1,18)-ZxZpE*VS1(1,18))
VS0(1,55)=QCz*VS0(1,34)+WQz*VS1(1,34)&
   +3D0*r1x2E*(VS0(1,19)-ZxZpE*VS1(1,19))
VS0(1,56)=QCz*VS0(1,35)+WQz*VS1(1,35)&
   +4D0*r1x2E*(VS0(1,20)-ZxZpE*VS1(1,20))
CASE(3)
VS0(1,36)=QCx*VS0(1,21)+WQx*VS1(1,21)&
   +4D0*r1x2E*(VS0(1,11)-ZxZpE*VS1(1,11))
VS0(1,37)=QCx*VS0(1,22)+WQx*VS1(1,22)&
   +3D0*r1x2E*(VS0(1,12)-ZxZpE*VS1(1,12))
VS0(1,38)=QCx*VS0(1,23)+WQx*VS1(1,23)&
   +2D0*r1x2E*(VS0(1,13)-ZxZpE*VS1(1,13))
VS0(1,39)=QCy*VS0(1,23)+WQy*VS1(1,23)&
   +2D0*r1x2E*(VS0(1,12)-ZxZpE*VS1(1,12))
VS0(1,40)=QCy*VS0(1,24)+WQy*VS1(1,24)&
   +3D0*r1x2E*(VS0(1,13)-ZxZpE*VS1(1,13))
VS0(1,41)=QCy*VS0(1,25)+WQy*VS1(1,25)&
   +4D0*r1x2E*(VS0(1,14)-ZxZpE*VS1(1,14))
VS0(1,42)=QCx*VS0(1,26)+WQx*VS1(1,26)&
   +3D0*r1x2E*(VS0(1,15)-ZxZpE*VS1(1,15))
VS0(1,43)=QCx*VS0(1,27)+WQx*VS1(1,27)&
   +2D0*r1x2E*(VS0(1,16)-ZxZpE*VS1(1,16))
VS0(1,44)=QCx*VS0(1,28)+WQx*VS1(1,28)&
   +r1x2E*(VS0(1,17)-ZxZpE*VS1(1,17))
VS0(1,45)=QCy*VS0(1,28)+WQy*VS1(1,28)&
   +2D0*r1x2E*(VS0(1,16)-ZxZpE*VS1(1,16))
VS0(1,46)=QCy*VS0(1,29)+WQy*VS1(1,29)&
   +3D0*r1x2E*(VS0(1,17)-ZxZpE*VS1(1,17))
VS0(1,47)=QCx*VS0(1,30)+WQx*VS1(1,30)&
   +2D0*r1x2E*(VS0(1,18)-ZxZpE*VS1(1,18))
VS0(1,48)=QCx*VS0(1,31)+WQx*VS1(1,31)&
   +r1x2E*(VS0(1,19)-ZxZpE*VS1(1,19))
VS0(1,49)=QCy*VS0(1,31)+WQy*VS1(1,31)&
   +r1x2E*(VS0(1,18)-ZxZpE*VS1(1,18))
VS0(1,50)=QCy*VS0(1,32)+WQy*VS1(1,32)&
   +2D0*r1x2E*(VS0(1,19)-ZxZpE*VS1(1,19))
VS0(1,51)=QCz*VS0(1,30)+WQz*VS1(1,30)-r1x2E*VR1(1,30)&
   +2D0*r1x2E*(VS0(1,15)-ZxZpE*VS1(1,15))
VS0(1,52)=QCz*VS0(1,31)+WQz*VS1(1,31)-r1x2E*VR1(1,31)&
   +2D0*r1x2E*(VS0(1,16)-ZxZpE*VS1(1,16))
VS0(1,53)=QCz*VS0(1,32)+WQz*VS1(1,32)-r1x2E*VR1(1,32)&
   +2D0*r1x2E*(VS0(1,17)-ZxZpE*VS1(1,17))
VS0(1,54)=QCz*VS0(1,33)+WQz*VS1(1,33)-r1x2E*VR1(1,33)&
   +3D0*r1x2E*(VS0(1,18)-ZxZpE*VS1(1,18))
VS0(1,55)=QCz*VS0(1,34)+WQz*VS1(1,34)-r1x2E*VR1(1,34)&
   +3D0*r1x2E*(VS0(1,19)-ZxZpE*VS1(1,19))
VS0(1,56)=QCz*VS0(1,35)+WQz*VS1(1,35)-r1x2E*VR1(1,35)&
   +4D0*r1x2E*(VS0(1,20)-ZxZpE*VS1(1,20))
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRs0h0'
STOP
END SELECT
END SUBROUTINE MVRRs0h0
