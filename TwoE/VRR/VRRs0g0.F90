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
   SUBROUTINE VRRs0g0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2E*VRR0(1,6)
      V(2)=2.D0*V(1)
      V(3)=r1x2E*ZxZpE*VRR1(1,6)
      V(4)=-2.D0*V(3)
      V(5)=r1x2E*VRR0(1,7)
      V(6)=r1x2E*ZxZpE*VRR1(1,7)
      V(7)=r1x2E*VRR0(1,8)
      V(8)=2.D0*V(7)
      V(9)=r1x2E*ZxZpE*VRR1(1,8)
      V(10)=-2.D0*V(9)
      V(11)=r1x2E*VRR0(1,9)
      V(12)=r1x2E*ZxZpE*VRR1(1,9)
      V(13)=2.D0*V(11)
      V(14)=-2.D0*V(12)
      V(15)=r1x2E*VRR0(1,10)
      V(16)=r1x2E*ZxZpE*VRR1(1,10)
      V(17)=-V(16)
      VRR0(1,21)=3.D0*r1x2E*VRR0(1,5)+QCx*VRR0(1,11)-3.D0*r1x2E*ZxZpE*VRR1(1,5)+WQx*VRR1(1,11)
      VRR0(1,22)=V(2)+V(4)+QCx*VRR0(1,12)+WQx*VRR1(1,12)
      VRR0(1,23)=V(5)-V(6)+QCx*VRR0(1,13)+WQx*VRR1(1,13)
      VRR0(1,24)=V(2)+V(4)+QCy*VRR0(1,13)+WQy*VRR1(1,13)
      VRR0(1,25)=3.D0*V(5)-3.D0*V(6)+QCy*VRR0(1,14)+WQy*VRR1(1,14)
      VRR0(1,26)=V(8)+V(10)+QCx*VRR0(1,15)+WQx*VRR1(1,15)
      VRR0(1,27)=V(11)-V(12)+QCx*VRR0(1,16)+WQx*VRR1(1,16)
      VRR0(1,28)=V(7)-V(9)+QCy*VRR0(1,16)+WQy*VRR1(1,16)
      VRR0(1,29)=V(13)+V(14)+QCy*VRR0(1,17)+WQy*VRR1(1,17)
      VRR0(1,30)=V(15)+V(17)+QCx*VRR0(1,18)+WQx*VRR1(1,18)
      VRR0(1,31)=V(1)-V(3)+QCz*VRR0(1,16)+WQz*VRR1(1,16)
      VRR0(1,32)=V(15)+V(17)+QCy*VRR0(1,19)+WQy*VRR1(1,19)
      VRR0(1,33)=V(8)+V(10)+QCz*VRR0(1,18)+WQz*VRR1(1,18)
      VRR0(1,34)=V(13)+V(14)+QCz*VRR0(1,19)+WQz*VRR1(1,19)
      VRR0(1,35)=3.D0*V(15)-3.D0*V(16)+QCz*VRR0(1,20)+WQz*VRR1(1,20)
END SUBROUTINE VRRs0g0
SUBROUTINE MVRRs0g0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(1,21)=QCx*VS0(1,11)+WQx*VS1(1,11)-r1x2E*VR1(1,11)&
   +3D0*r1x2E*(VS0(1,5)-ZxZpE*VS1(1,5))
VS0(1,22)=QCx*VS0(1,12)+WQx*VS1(1,12)-r1x2E*VR1(1,12)&
   +2D0*r1x2E*(VS0(1,6)-ZxZpE*VS1(1,6))
VS0(1,23)=QCx*VS0(1,13)+WQx*VS1(1,13)-r1x2E*VR1(1,13)&
   +r1x2E*(VS0(1,7)-ZxZpE*VS1(1,7))
VS0(1,24)=QCy*VS0(1,13)+WQy*VS1(1,13)&
   +2D0*r1x2E*(VS0(1,6)-ZxZpE*VS1(1,6))
VS0(1,25)=QCy*VS0(1,14)+WQy*VS1(1,14)&
   +3D0*r1x2E*(VS0(1,7)-ZxZpE*VS1(1,7))
VS0(1,26)=QCx*VS0(1,15)+WQx*VS1(1,15)-r1x2E*VR1(1,15)&
   +2D0*r1x2E*(VS0(1,8)-ZxZpE*VS1(1,8))
VS0(1,27)=QCx*VS0(1,16)+WQx*VS1(1,16)-r1x2E*VR1(1,16)&
   +r1x2E*(VS0(1,9)-ZxZpE*VS1(1,9))
VS0(1,28)=QCy*VS0(1,16)+WQy*VS1(1,16)&
   +r1x2E*(VS0(1,8)-ZxZpE*VS1(1,8))
VS0(1,29)=QCy*VS0(1,17)+WQy*VS1(1,17)&
   +2D0*r1x2E*(VS0(1,9)-ZxZpE*VS1(1,9))
VS0(1,30)=QCx*VS0(1,18)+WQx*VS1(1,18)-r1x2E*VR1(1,18)&
   +r1x2E*(VS0(1,10)-ZxZpE*VS1(1,10))
VS0(1,31)=QCz*VS0(1,16)+WQz*VS1(1,16)&
   +r1x2E*(VS0(1,6)-ZxZpE*VS1(1,6))
VS0(1,32)=QCy*VS0(1,19)+WQy*VS1(1,19)&
   +r1x2E*(VS0(1,10)-ZxZpE*VS1(1,10))
VS0(1,33)=QCz*VS0(1,18)+WQz*VS1(1,18)&
   +2D0*r1x2E*(VS0(1,8)-ZxZpE*VS1(1,8))
VS0(1,34)=QCz*VS0(1,19)+WQz*VS1(1,19)&
   +2D0*r1x2E*(VS0(1,9)-ZxZpE*VS1(1,9))
VS0(1,35)=QCz*VS0(1,20)+WQz*VS1(1,20)&
   +3D0*r1x2E*(VS0(1,10)-ZxZpE*VS1(1,10))
CASE(2)
VS0(1,21)=QCx*VS0(1,11)+WQx*VS1(1,11)&
   +3D0*r1x2E*(VS0(1,5)-ZxZpE*VS1(1,5))
VS0(1,22)=QCx*VS0(1,12)+WQx*VS1(1,12)&
   +2D0*r1x2E*(VS0(1,6)-ZxZpE*VS1(1,6))
VS0(1,23)=QCx*VS0(1,13)+WQx*VS1(1,13)&
   +r1x2E*(VS0(1,7)-ZxZpE*VS1(1,7))
VS0(1,24)=QCy*VS0(1,13)+WQy*VS1(1,13)-r1x2E*VR1(1,13)&
   +2D0*r1x2E*(VS0(1,6)-ZxZpE*VS1(1,6))
VS0(1,25)=QCy*VS0(1,14)+WQy*VS1(1,14)-r1x2E*VR1(1,14)&
   +3D0*r1x2E*(VS0(1,7)-ZxZpE*VS1(1,7))
VS0(1,26)=QCx*VS0(1,15)+WQx*VS1(1,15)&
   +2D0*r1x2E*(VS0(1,8)-ZxZpE*VS1(1,8))
VS0(1,27)=QCx*VS0(1,16)+WQx*VS1(1,16)&
   +r1x2E*(VS0(1,9)-ZxZpE*VS1(1,9))
VS0(1,28)=QCy*VS0(1,16)+WQy*VS1(1,16)-r1x2E*VR1(1,16)&
   +r1x2E*(VS0(1,8)-ZxZpE*VS1(1,8))
VS0(1,29)=QCy*VS0(1,17)+WQy*VS1(1,17)-r1x2E*VR1(1,17)&
   +2D0*r1x2E*(VS0(1,9)-ZxZpE*VS1(1,9))
VS0(1,30)=QCx*VS0(1,18)+WQx*VS1(1,18)&
   +r1x2E*(VS0(1,10)-ZxZpE*VS1(1,10))
VS0(1,31)=QCz*VS0(1,16)+WQz*VS1(1,16)&
   +r1x2E*(VS0(1,6)-ZxZpE*VS1(1,6))
VS0(1,32)=QCy*VS0(1,19)+WQy*VS1(1,19)-r1x2E*VR1(1,19)&
   +r1x2E*(VS0(1,10)-ZxZpE*VS1(1,10))
VS0(1,33)=QCz*VS0(1,18)+WQz*VS1(1,18)&
   +2D0*r1x2E*(VS0(1,8)-ZxZpE*VS1(1,8))
VS0(1,34)=QCz*VS0(1,19)+WQz*VS1(1,19)&
   +2D0*r1x2E*(VS0(1,9)-ZxZpE*VS1(1,9))
VS0(1,35)=QCz*VS0(1,20)+WQz*VS1(1,20)&
   +3D0*r1x2E*(VS0(1,10)-ZxZpE*VS1(1,10))
CASE(3)
VS0(1,21)=QCx*VS0(1,11)+WQx*VS1(1,11)&
   +3D0*r1x2E*(VS0(1,5)-ZxZpE*VS1(1,5))
VS0(1,22)=QCx*VS0(1,12)+WQx*VS1(1,12)&
   +2D0*r1x2E*(VS0(1,6)-ZxZpE*VS1(1,6))
VS0(1,23)=QCx*VS0(1,13)+WQx*VS1(1,13)&
   +r1x2E*(VS0(1,7)-ZxZpE*VS1(1,7))
VS0(1,24)=QCy*VS0(1,13)+WQy*VS1(1,13)&
   +2D0*r1x2E*(VS0(1,6)-ZxZpE*VS1(1,6))
VS0(1,25)=QCy*VS0(1,14)+WQy*VS1(1,14)&
   +3D0*r1x2E*(VS0(1,7)-ZxZpE*VS1(1,7))
VS0(1,26)=QCx*VS0(1,15)+WQx*VS1(1,15)&
   +2D0*r1x2E*(VS0(1,8)-ZxZpE*VS1(1,8))
VS0(1,27)=QCx*VS0(1,16)+WQx*VS1(1,16)&
   +r1x2E*(VS0(1,9)-ZxZpE*VS1(1,9))
VS0(1,28)=QCy*VS0(1,16)+WQy*VS1(1,16)&
   +r1x2E*(VS0(1,8)-ZxZpE*VS1(1,8))
VS0(1,29)=QCy*VS0(1,17)+WQy*VS1(1,17)&
   +2D0*r1x2E*(VS0(1,9)-ZxZpE*VS1(1,9))
VS0(1,30)=QCx*VS0(1,18)+WQx*VS1(1,18)&
   +r1x2E*(VS0(1,10)-ZxZpE*VS1(1,10))
VS0(1,31)=QCz*VS0(1,16)+WQz*VS1(1,16)-r1x2E*VR1(1,16)&
   +r1x2E*(VS0(1,6)-ZxZpE*VS1(1,6))
VS0(1,32)=QCy*VS0(1,19)+WQy*VS1(1,19)&
   +r1x2E*(VS0(1,10)-ZxZpE*VS1(1,10))
VS0(1,33)=QCz*VS0(1,18)+WQz*VS1(1,18)-r1x2E*VR1(1,18)&
   +2D0*r1x2E*(VS0(1,8)-ZxZpE*VS1(1,8))
VS0(1,34)=QCz*VS0(1,19)+WQz*VS1(1,19)-r1x2E*VR1(1,19)&
   +2D0*r1x2E*(VS0(1,9)-ZxZpE*VS1(1,9))
VS0(1,35)=QCz*VS0(1,20)+WQz*VS1(1,20)-r1x2E*VR1(1,20)&
   +3D0*r1x2E*(VS0(1,10)-ZxZpE*VS1(1,10))
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRs0g0'
STOP
END SELECT
END SUBROUTINE MVRRs0g0
