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
   SUBROUTINE VRRp0d0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2E*VRR0(2,1)
      V(2)=r1x2E*ZxZpE*VRR1(2,1)
      V(3)=-V(2)
      V(4)=HfxZpE*VRR1(1,3)
      V(5)=HfxZpE*VRR1(1,4)
      V(6)=r1x2E*VRR0(3,1)
      V(7)=r1x2E*ZxZpE*VRR1(3,1)
      V(8)=-V(7)
      V(9)=r1x2E*VRR0(4,1)
      V(10)=r1x2E*ZxZpE*VRR1(4,1)
      V(11)=-V(10)
      VRR0(2,5)=V(1)+V(3)+QCx*VRR0(2,2)+HfxZpE*VRR1(1,2)+WQx*VRR1(2,2)
      VRR0(2,6)=V(4)+QCx*VRR0(2,3)+WQx*VRR1(2,3)
      VRR0(2,7)=V(1)+V(3)+QCy*VRR0(2,3)+WQy*VRR1(2,3)
      VRR0(2,8)=V(5)+QCx*VRR0(2,4)+WQx*VRR1(2,4)
      VRR0(2,9)=QCy*VRR0(2,4)+WQy*VRR1(2,4)
      VRR0(2,10)=V(1)+V(3)+QCz*VRR0(2,4)+WQz*VRR1(2,4)
      VRR0(3,5)=V(6)+V(8)+QCx*VRR0(3,2)+WQx*VRR1(3,2)
      VRR0(3,6)=QCx*VRR0(3,3)+WQx*VRR1(3,3)
      VRR0(3,7)=V(4)+V(6)+V(8)+QCy*VRR0(3,3)+WQy*VRR1(3,3)
      VRR0(3,8)=QCx*VRR0(3,4)+WQx*VRR1(3,4)
      VRR0(3,9)=V(5)+QCy*VRR0(3,4)+WQy*VRR1(3,4)
      VRR0(3,10)=V(6)+V(8)+QCz*VRR0(3,4)+WQz*VRR1(3,4)
      VRR0(4,5)=V(9)+V(11)+QCx*VRR0(4,2)+WQx*VRR1(4,2)
      VRR0(4,6)=QCx*VRR0(4,3)+WQx*VRR1(4,3)
      VRR0(4,7)=V(9)+V(11)+QCy*VRR0(4,3)+WQy*VRR1(4,3)
      VRR0(4,8)=QCx*VRR0(4,4)+WQx*VRR1(4,4)
      VRR0(4,9)=QCy*VRR0(4,4)+WQy*VRR1(4,4)
      VRR0(4,10)=V(5)+V(9)+V(11)+QCz*VRR0(4,4)+WQz*VRR1(4,4)
END SUBROUTINE VRRp0d0
SUBROUTINE MVRRp0d0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(2,5)=QCx*VS0(2,2)+WQx*VS1(2,2)-r1x2E*VR1(2,2)&
   +r1x2E*(VS0(2,1)-ZxZpE*VS1(2,1))&
   +HfxZpE*VS1(1,2)
VS0(2,6)=QCx*VS0(2,3)+WQx*VS1(2,3)-r1x2E*VR1(2,3)&
   +HfxZpE*VS1(1,3)
VS0(2,7)=QCy*VS0(2,3)+WQy*VS1(2,3)&
   +r1x2E*(VS0(2,1)-ZxZpE*VS1(2,1))
VS0(2,8)=QCx*VS0(2,4)+WQx*VS1(2,4)-r1x2E*VR1(2,4)&
   +HfxZpE*VS1(1,4)
VS0(2,9)=QCy*VS0(2,4)+WQy*VS1(2,4)
VS0(2,10)=QCz*VS0(2,4)+WQz*VS1(2,4)&
   +r1x2E*(VS0(2,1)-ZxZpE*VS1(2,1))
VS0(3,5)=QCx*VS0(3,2)+WQx*VS1(3,2)-r1x2E*VR1(3,2)&
   +r1x2E*(VS0(3,1)-ZxZpE*VS1(3,1))
VS0(3,6)=QCx*VS0(3,3)+WQx*VS1(3,3)-r1x2E*VR1(3,3)
VS0(3,7)=QCy*VS0(3,3)+WQy*VS1(3,3)&
   +r1x2E*(VS0(3,1)-ZxZpE*VS1(3,1))&
   +HfxZpE*VS1(1,3)
VS0(3,8)=QCx*VS0(3,4)+WQx*VS1(3,4)-r1x2E*VR1(3,4)
VS0(3,9)=QCy*VS0(3,4)+WQy*VS1(3,4)&
   +HfxZpE*VS1(1,4)
VS0(3,10)=QCz*VS0(3,4)+WQz*VS1(3,4)&
   +r1x2E*(VS0(3,1)-ZxZpE*VS1(3,1))
VS0(4,5)=QCx*VS0(4,2)+WQx*VS1(4,2)-r1x2E*VR1(4,2)&
   +r1x2E*(VS0(4,1)-ZxZpE*VS1(4,1))
VS0(4,6)=QCx*VS0(4,3)+WQx*VS1(4,3)-r1x2E*VR1(4,3)
VS0(4,7)=QCy*VS0(4,3)+WQy*VS1(4,3)&
   +r1x2E*(VS0(4,1)-ZxZpE*VS1(4,1))
VS0(4,8)=QCx*VS0(4,4)+WQx*VS1(4,4)-r1x2E*VR1(4,4)
VS0(4,9)=QCy*VS0(4,4)+WQy*VS1(4,4)
VS0(4,10)=QCz*VS0(4,4)+WQz*VS1(4,4)&
   +r1x2E*(VS0(4,1)-ZxZpE*VS1(4,1))&
   +HfxZpE*VS1(1,4)
CASE(2)
VS0(2,5)=QCx*VS0(2,2)+WQx*VS1(2,2)&
   +r1x2E*(VS0(2,1)-ZxZpE*VS1(2,1))&
   +HfxZpE*VS1(1,2)
VS0(2,6)=QCx*VS0(2,3)+WQx*VS1(2,3)&
   +HfxZpE*VS1(1,3)
VS0(2,7)=QCy*VS0(2,3)+WQy*VS1(2,3)-r1x2E*VR1(2,3)&
   +r1x2E*(VS0(2,1)-ZxZpE*VS1(2,1))
VS0(2,8)=QCx*VS0(2,4)+WQx*VS1(2,4)&
   +HfxZpE*VS1(1,4)
VS0(2,9)=QCy*VS0(2,4)+WQy*VS1(2,4)-r1x2E*VR1(2,4)
VS0(2,10)=QCz*VS0(2,4)+WQz*VS1(2,4)&
   +r1x2E*(VS0(2,1)-ZxZpE*VS1(2,1))
VS0(3,5)=QCx*VS0(3,2)+WQx*VS1(3,2)&
   +r1x2E*(VS0(3,1)-ZxZpE*VS1(3,1))
VS0(3,6)=QCx*VS0(3,3)+WQx*VS1(3,3)
VS0(3,7)=QCy*VS0(3,3)+WQy*VS1(3,3)-r1x2E*VR1(3,3)&
   +r1x2E*(VS0(3,1)-ZxZpE*VS1(3,1))&
   +HfxZpE*VS1(1,3)
VS0(3,8)=QCx*VS0(3,4)+WQx*VS1(3,4)
VS0(3,9)=QCy*VS0(3,4)+WQy*VS1(3,4)-r1x2E*VR1(3,4)&
   +HfxZpE*VS1(1,4)
VS0(3,10)=QCz*VS0(3,4)+WQz*VS1(3,4)&
   +r1x2E*(VS0(3,1)-ZxZpE*VS1(3,1))
VS0(4,5)=QCx*VS0(4,2)+WQx*VS1(4,2)&
   +r1x2E*(VS0(4,1)-ZxZpE*VS1(4,1))
VS0(4,6)=QCx*VS0(4,3)+WQx*VS1(4,3)
VS0(4,7)=QCy*VS0(4,3)+WQy*VS1(4,3)-r1x2E*VR1(4,3)&
   +r1x2E*(VS0(4,1)-ZxZpE*VS1(4,1))
VS0(4,8)=QCx*VS0(4,4)+WQx*VS1(4,4)
VS0(4,9)=QCy*VS0(4,4)+WQy*VS1(4,4)-r1x2E*VR1(4,4)
VS0(4,10)=QCz*VS0(4,4)+WQz*VS1(4,4)&
   +r1x2E*(VS0(4,1)-ZxZpE*VS1(4,1))&
   +HfxZpE*VS1(1,4)
CASE(3)
VS0(2,5)=QCx*VS0(2,2)+WQx*VS1(2,2)&
   +r1x2E*(VS0(2,1)-ZxZpE*VS1(2,1))&
   +HfxZpE*VS1(1,2)
VS0(2,6)=QCx*VS0(2,3)+WQx*VS1(2,3)&
   +HfxZpE*VS1(1,3)
VS0(2,7)=QCy*VS0(2,3)+WQy*VS1(2,3)&
   +r1x2E*(VS0(2,1)-ZxZpE*VS1(2,1))
VS0(2,8)=QCx*VS0(2,4)+WQx*VS1(2,4)&
   +HfxZpE*VS1(1,4)
VS0(2,9)=QCy*VS0(2,4)+WQy*VS1(2,4)
VS0(2,10)=QCz*VS0(2,4)+WQz*VS1(2,4)-r1x2E*VR1(2,4)&
   +r1x2E*(VS0(2,1)-ZxZpE*VS1(2,1))
VS0(3,5)=QCx*VS0(3,2)+WQx*VS1(3,2)&
   +r1x2E*(VS0(3,1)-ZxZpE*VS1(3,1))
VS0(3,6)=QCx*VS0(3,3)+WQx*VS1(3,3)
VS0(3,7)=QCy*VS0(3,3)+WQy*VS1(3,3)&
   +r1x2E*(VS0(3,1)-ZxZpE*VS1(3,1))&
   +HfxZpE*VS1(1,3)
VS0(3,8)=QCx*VS0(3,4)+WQx*VS1(3,4)
VS0(3,9)=QCy*VS0(3,4)+WQy*VS1(3,4)&
   +HfxZpE*VS1(1,4)
VS0(3,10)=QCz*VS0(3,4)+WQz*VS1(3,4)-r1x2E*VR1(3,4)&
   +r1x2E*(VS0(3,1)-ZxZpE*VS1(3,1))
VS0(4,5)=QCx*VS0(4,2)+WQx*VS1(4,2)&
   +r1x2E*(VS0(4,1)-ZxZpE*VS1(4,1))
VS0(4,6)=QCx*VS0(4,3)+WQx*VS1(4,3)
VS0(4,7)=QCy*VS0(4,3)+WQy*VS1(4,3)&
   +r1x2E*(VS0(4,1)-ZxZpE*VS1(4,1))
VS0(4,8)=QCx*VS0(4,4)+WQx*VS1(4,4)
VS0(4,9)=QCy*VS0(4,4)+WQy*VS1(4,4)
VS0(4,10)=QCz*VS0(4,4)+WQz*VS1(4,4)-r1x2E*VR1(4,4)&
   +r1x2E*(VS0(4,1)-ZxZpE*VS1(4,1))&
   +HfxZpE*VS1(1,4)
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRp0d0'
STOP
END SELECT
END SUBROUTINE MVRRp0d0
