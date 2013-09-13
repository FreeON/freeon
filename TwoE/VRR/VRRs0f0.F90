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
   SUBROUTINE VRRs0f0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2E*VRR0(1,2)
      V(2)=r1x2E*ZxZpE*VRR1(1,2)
      V(3)=r1x2E*VRR0(1,3)
      V(4)=r1x2E*ZxZpE*VRR1(1,3)
      V(5)=-V(4)
      V(6)=-V(2)
      V(7)=r1x2E*VRR0(1,4)
      V(8)=r1x2E*ZxZpE*VRR1(1,4)
      V(9)=-V(8)
      VRR0(1,11)=2.D0*V(1)-2.D0*V(2)+QCx*VRR0(1,5)+WQx*VRR1(1,5)
      VRR0(1,12)=V(3)+V(5)+QCx*VRR0(1,6)+WQx*VRR1(1,6)
      VRR0(1,13)=V(1)+V(6)+QCy*VRR0(1,6)+WQy*VRR1(1,6)
      VRR0(1,14)=2.D0*V(3)-2.D0*V(4)+QCy*VRR0(1,7)+WQy*VRR1(1,7)
      VRR0(1,15)=V(7)+V(9)+QCx*VRR0(1,8)+WQx*VRR1(1,8)
      VRR0(1,16)=QCx*VRR0(1,9)+WQx*VRR1(1,9)
      VRR0(1,17)=V(7)+V(9)+QCy*VRR0(1,9)+WQy*VRR1(1,9)
      VRR0(1,18)=V(1)+V(6)+QCz*VRR0(1,8)+WQz*VRR1(1,8)
      VRR0(1,19)=V(3)+V(5)+QCz*VRR0(1,9)+WQz*VRR1(1,9)
      VRR0(1,20)=2.D0*V(7)-2.D0*V(8)+QCz*VRR0(1,10)+WQz*VRR1(1,10)
END SUBROUTINE VRRs0f0
SUBROUTINE MVRRs0f0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(1,11)=QCx*VS0(1,5)+WQx*VS1(1,5)-r1x2E*VR1(1,5)&
   +2D0*r1x2E*(VS0(1,2)-ZxZpE*VS1(1,2))
VS0(1,12)=QCx*VS0(1,6)+WQx*VS1(1,6)-r1x2E*VR1(1,6)&
   +r1x2E*(VS0(1,3)-ZxZpE*VS1(1,3))
VS0(1,13)=QCy*VS0(1,6)+WQy*VS1(1,6)&
   +r1x2E*(VS0(1,2)-ZxZpE*VS1(1,2))
VS0(1,14)=QCy*VS0(1,7)+WQy*VS1(1,7)&
   +2D0*r1x2E*(VS0(1,3)-ZxZpE*VS1(1,3))
VS0(1,15)=QCx*VS0(1,8)+WQx*VS1(1,8)-r1x2E*VR1(1,8)&
   +r1x2E*(VS0(1,4)-ZxZpE*VS1(1,4))
VS0(1,16)=QCx*VS0(1,9)+WQx*VS1(1,9)-r1x2E*VR1(1,9)
VS0(1,17)=QCy*VS0(1,9)+WQy*VS1(1,9)&
   +r1x2E*(VS0(1,4)-ZxZpE*VS1(1,4))
VS0(1,18)=QCz*VS0(1,8)+WQz*VS1(1,8)&
   +r1x2E*(VS0(1,2)-ZxZpE*VS1(1,2))
VS0(1,19)=QCz*VS0(1,9)+WQz*VS1(1,9)&
   +r1x2E*(VS0(1,3)-ZxZpE*VS1(1,3))
VS0(1,20)=QCz*VS0(1,10)+WQz*VS1(1,10)&
   +2D0*r1x2E*(VS0(1,4)-ZxZpE*VS1(1,4))
CASE(2)
VS0(1,11)=QCx*VS0(1,5)+WQx*VS1(1,5)&
   +2D0*r1x2E*(VS0(1,2)-ZxZpE*VS1(1,2))
VS0(1,12)=QCx*VS0(1,6)+WQx*VS1(1,6)&
   +r1x2E*(VS0(1,3)-ZxZpE*VS1(1,3))
VS0(1,13)=QCy*VS0(1,6)+WQy*VS1(1,6)-r1x2E*VR1(1,6)&
   +r1x2E*(VS0(1,2)-ZxZpE*VS1(1,2))
VS0(1,14)=QCy*VS0(1,7)+WQy*VS1(1,7)-r1x2E*VR1(1,7)&
   +2D0*r1x2E*(VS0(1,3)-ZxZpE*VS1(1,3))
VS0(1,15)=QCx*VS0(1,8)+WQx*VS1(1,8)&
   +r1x2E*(VS0(1,4)-ZxZpE*VS1(1,4))
VS0(1,16)=QCx*VS0(1,9)+WQx*VS1(1,9)
VS0(1,17)=QCy*VS0(1,9)+WQy*VS1(1,9)-r1x2E*VR1(1,9)&
   +r1x2E*(VS0(1,4)-ZxZpE*VS1(1,4))
VS0(1,18)=QCz*VS0(1,8)+WQz*VS1(1,8)&
   +r1x2E*(VS0(1,2)-ZxZpE*VS1(1,2))
VS0(1,19)=QCz*VS0(1,9)+WQz*VS1(1,9)&
   +r1x2E*(VS0(1,3)-ZxZpE*VS1(1,3))
VS0(1,20)=QCz*VS0(1,10)+WQz*VS1(1,10)&
   +2D0*r1x2E*(VS0(1,4)-ZxZpE*VS1(1,4))
CASE(3)
VS0(1,11)=QCx*VS0(1,5)+WQx*VS1(1,5)&
   +2D0*r1x2E*(VS0(1,2)-ZxZpE*VS1(1,2))
VS0(1,12)=QCx*VS0(1,6)+WQx*VS1(1,6)&
   +r1x2E*(VS0(1,3)-ZxZpE*VS1(1,3))
VS0(1,13)=QCy*VS0(1,6)+WQy*VS1(1,6)&
   +r1x2E*(VS0(1,2)-ZxZpE*VS1(1,2))
VS0(1,14)=QCy*VS0(1,7)+WQy*VS1(1,7)&
   +2D0*r1x2E*(VS0(1,3)-ZxZpE*VS1(1,3))
VS0(1,15)=QCx*VS0(1,8)+WQx*VS1(1,8)&
   +r1x2E*(VS0(1,4)-ZxZpE*VS1(1,4))
VS0(1,16)=QCx*VS0(1,9)+WQx*VS1(1,9)
VS0(1,17)=QCy*VS0(1,9)+WQy*VS1(1,9)&
   +r1x2E*(VS0(1,4)-ZxZpE*VS1(1,4))
VS0(1,18)=QCz*VS0(1,8)+WQz*VS1(1,8)-r1x2E*VR1(1,8)&
   +r1x2E*(VS0(1,2)-ZxZpE*VS1(1,2))
VS0(1,19)=QCz*VS0(1,9)+WQz*VS1(1,9)-r1x2E*VR1(1,9)&
   +r1x2E*(VS0(1,3)-ZxZpE*VS1(1,3))
VS0(1,20)=QCz*VS0(1,10)+WQz*VS1(1,10)-r1x2E*VR1(1,10)&
   +2D0*r1x2E*(VS0(1,4)-ZxZpE*VS1(1,4))
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRs0f0'
STOP
END SELECT
END SUBROUTINE MVRRs0f0
