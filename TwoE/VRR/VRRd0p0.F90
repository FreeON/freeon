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
   SUBROUTINE VRRd0p0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2Z*VRR0(1,2)
      V(2)=ExZpE*r1x2Z*VRR1(1,2)
      V(3)=-V(2)
      V(4)=HfxZpE*VRR1(2,1)
      V(5)=r1x2Z*VRR0(1,3)
      V(6)=ExZpE*r1x2Z*VRR1(1,3)
      V(7)=-V(6)
      V(8)=r1x2Z*VRR0(1,4)
      V(9)=ExZpE*r1x2Z*VRR1(1,4)
      V(10)=-V(9)
      V(11)=HfxZpE*VRR1(3,1)
      V(12)=HfxZpE*VRR1(4,1)
      VRR0(5,2)=V(1)+V(3)+V(4)+PAx*VRR0(2,2)+WPx*VRR1(2,2)
      VRR0(5,3)=V(5)+V(7)+PAx*VRR0(2,3)+WPx*VRR1(2,3)
      VRR0(5,4)=V(8)+V(10)+PAx*VRR0(2,4)+WPx*VRR1(2,4)
      VRR0(6,2)=V(11)+QCx*VRR0(6,1)+WQx*VRR1(6,1)
      VRR0(6,3)=V(4)+QCy*VRR0(6,1)+WQy*VRR1(6,1)
      VRR0(6,4)=QCz*VRR0(6,1)+WQz*VRR1(6,1)
      VRR0(7,2)=V(1)+V(3)+PAy*VRR0(3,2)+WPy*VRR1(3,2)
      VRR0(7,3)=V(5)+V(7)+V(11)+PAy*VRR0(3,3)+WPy*VRR1(3,3)
      VRR0(7,4)=V(8)+V(10)+PAy*VRR0(3,4)+WPy*VRR1(3,4)
      VRR0(8,2)=V(12)+QCx*VRR0(8,1)+WQx*VRR1(8,1)
      VRR0(8,3)=QCy*VRR0(8,1)+WQy*VRR1(8,1)
      VRR0(8,4)=V(4)+QCz*VRR0(8,1)+WQz*VRR1(8,1)
      VRR0(9,2)=QCx*VRR0(9,1)+WQx*VRR1(9,1)
      VRR0(9,3)=V(12)+QCy*VRR0(9,1)+WQy*VRR1(9,1)
      VRR0(9,4)=V(11)+QCz*VRR0(9,1)+WQz*VRR1(9,1)
      VRR0(10,2)=V(1)+V(3)+PAz*VRR0(4,2)+WPz*VRR1(4,2)
      VRR0(10,3)=V(5)+V(7)+PAz*VRR0(4,3)+WPz*VRR1(4,3)
      VRR0(10,4)=V(8)+V(10)+V(12)+PAz*VRR0(4,4)+WPz*VRR1(4,4)
END SUBROUTINE VRRd0p0
SUBROUTINE MVRRd0p0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(5,2)=PAx*VS0(2,2)+WPx*VS1(2,2)+r1x2Z*VR1(2,2)&
   +r1x2Z*(VS0(1,2)-ExZpE*VS1(1,2))&
   +HfxZpE*VS1(2,1)
VS0(5,3)=PAx*VS0(2,3)+WPx*VS1(2,3)+r1x2Z*VR1(2,3)&
   +r1x2Z*(VS0(1,3)-ExZpE*VS1(1,3))
VS0(5,4)=PAx*VS0(2,4)+WPx*VS1(2,4)+r1x2Z*VR1(2,4)&
   +r1x2Z*(VS0(1,4)-ExZpE*VS1(1,4))
VS0(6,2)=PAx*VS0(3,2)+WPx*VS1(3,2)+r1x2Z*VR1(3,2)&
   +HfxZpE*VS1(3,1)
VS0(6,3)=PAx*VS0(3,3)+WPx*VS1(3,3)+r1x2Z*VR1(3,3)
VS0(6,4)=PAx*VS0(3,4)+WPx*VS1(3,4)+r1x2Z*VR1(3,4)
VS0(7,2)=PAy*VS0(3,2)+WPy*VS1(3,2)&
   +r1x2Z*(VS0(1,2)-ExZpE*VS1(1,2))
VS0(7,3)=PAy*VS0(3,3)+WPy*VS1(3,3)&
   +r1x2Z*(VS0(1,3)-ExZpE*VS1(1,3))&
   +HfxZpE*VS1(3,1)
VS0(7,4)=PAy*VS0(3,4)+WPy*VS1(3,4)&
   +r1x2Z*(VS0(1,4)-ExZpE*VS1(1,4))
VS0(8,2)=PAx*VS0(4,2)+WPx*VS1(4,2)+r1x2Z*VR1(4,2)&
   +HfxZpE*VS1(4,1)
VS0(8,3)=PAx*VS0(4,3)+WPx*VS1(4,3)+r1x2Z*VR1(4,3)
VS0(8,4)=PAx*VS0(4,4)+WPx*VS1(4,4)+r1x2Z*VR1(4,4)
VS0(9,2)=PAy*VS0(4,2)+WPy*VS1(4,2)
VS0(9,3)=PAy*VS0(4,3)+WPy*VS1(4,3)&
   +HfxZpE*VS1(4,1)
VS0(9,4)=PAy*VS0(4,4)+WPy*VS1(4,4)
VS0(10,2)=PAz*VS0(4,2)+WPz*VS1(4,2)&
   +r1x2Z*(VS0(1,2)-ExZpE*VS1(1,2))
VS0(10,3)=PAz*VS0(4,3)+WPz*VS1(4,3)&
   +r1x2Z*(VS0(1,3)-ExZpE*VS1(1,3))
VS0(10,4)=PAz*VS0(4,4)+WPz*VS1(4,4)&
   +r1x2Z*(VS0(1,4)-ExZpE*VS1(1,4))&
   +HfxZpE*VS1(4,1)
CASE(2)
VS0(5,2)=PAx*VS0(2,2)+WPx*VS1(2,2)&
   +r1x2Z*(VS0(1,2)-ExZpE*VS1(1,2))&
   +HfxZpE*VS1(2,1)
VS0(5,3)=PAx*VS0(2,3)+WPx*VS1(2,3)&
   +r1x2Z*(VS0(1,3)-ExZpE*VS1(1,3))
VS0(5,4)=PAx*VS0(2,4)+WPx*VS1(2,4)&
   +r1x2Z*(VS0(1,4)-ExZpE*VS1(1,4))
VS0(6,2)=PAx*VS0(3,2)+WPx*VS1(3,2)&
   +HfxZpE*VS1(3,1)
VS0(6,3)=PAx*VS0(3,3)+WPx*VS1(3,3)
VS0(6,4)=PAx*VS0(3,4)+WPx*VS1(3,4)
VS0(7,2)=PAy*VS0(3,2)+WPy*VS1(3,2)+r1x2Z*VR1(3,2)&
   +r1x2Z*(VS0(1,2)-ExZpE*VS1(1,2))
VS0(7,3)=PAy*VS0(3,3)+WPy*VS1(3,3)+r1x2Z*VR1(3,3)&
   +r1x2Z*(VS0(1,3)-ExZpE*VS1(1,3))&
   +HfxZpE*VS1(3,1)
VS0(7,4)=PAy*VS0(3,4)+WPy*VS1(3,4)+r1x2Z*VR1(3,4)&
   +r1x2Z*(VS0(1,4)-ExZpE*VS1(1,4))
VS0(8,2)=PAx*VS0(4,2)+WPx*VS1(4,2)&
   +HfxZpE*VS1(4,1)
VS0(8,3)=PAx*VS0(4,3)+WPx*VS1(4,3)
VS0(8,4)=PAx*VS0(4,4)+WPx*VS1(4,4)
VS0(9,2)=PAy*VS0(4,2)+WPy*VS1(4,2)+r1x2Z*VR1(4,2)
VS0(9,3)=PAy*VS0(4,3)+WPy*VS1(4,3)+r1x2Z*VR1(4,3)&
   +HfxZpE*VS1(4,1)
VS0(9,4)=PAy*VS0(4,4)+WPy*VS1(4,4)+r1x2Z*VR1(4,4)
VS0(10,2)=PAz*VS0(4,2)+WPz*VS1(4,2)&
   +r1x2Z*(VS0(1,2)-ExZpE*VS1(1,2))
VS0(10,3)=PAz*VS0(4,3)+WPz*VS1(4,3)&
   +r1x2Z*(VS0(1,3)-ExZpE*VS1(1,3))
VS0(10,4)=PAz*VS0(4,4)+WPz*VS1(4,4)&
   +r1x2Z*(VS0(1,4)-ExZpE*VS1(1,4))&
   +HfxZpE*VS1(4,1)
CASE(3)
VS0(5,2)=PAx*VS0(2,2)+WPx*VS1(2,2)&
   +r1x2Z*(VS0(1,2)-ExZpE*VS1(1,2))&
   +HfxZpE*VS1(2,1)
VS0(5,3)=PAx*VS0(2,3)+WPx*VS1(2,3)&
   +r1x2Z*(VS0(1,3)-ExZpE*VS1(1,3))
VS0(5,4)=PAx*VS0(2,4)+WPx*VS1(2,4)&
   +r1x2Z*(VS0(1,4)-ExZpE*VS1(1,4))
VS0(6,2)=PAx*VS0(3,2)+WPx*VS1(3,2)&
   +HfxZpE*VS1(3,1)
VS0(6,3)=PAx*VS0(3,3)+WPx*VS1(3,3)
VS0(6,4)=PAx*VS0(3,4)+WPx*VS1(3,4)
VS0(7,2)=PAy*VS0(3,2)+WPy*VS1(3,2)&
   +r1x2Z*(VS0(1,2)-ExZpE*VS1(1,2))
VS0(7,3)=PAy*VS0(3,3)+WPy*VS1(3,3)&
   +r1x2Z*(VS0(1,3)-ExZpE*VS1(1,3))&
   +HfxZpE*VS1(3,1)
VS0(7,4)=PAy*VS0(3,4)+WPy*VS1(3,4)&
   +r1x2Z*(VS0(1,4)-ExZpE*VS1(1,4))
VS0(8,2)=PAx*VS0(4,2)+WPx*VS1(4,2)&
   +HfxZpE*VS1(4,1)
VS0(8,3)=PAx*VS0(4,3)+WPx*VS1(4,3)
VS0(8,4)=PAx*VS0(4,4)+WPx*VS1(4,4)
VS0(9,2)=PAy*VS0(4,2)+WPy*VS1(4,2)
VS0(9,3)=PAy*VS0(4,3)+WPy*VS1(4,3)&
   +HfxZpE*VS1(4,1)
VS0(9,4)=PAy*VS0(4,4)+WPy*VS1(4,4)
VS0(10,2)=PAz*VS0(4,2)+WPz*VS1(4,2)+r1x2Z*VR1(4,2)&
   +r1x2Z*(VS0(1,2)-ExZpE*VS1(1,2))
VS0(10,3)=PAz*VS0(4,3)+WPz*VS1(4,3)+r1x2Z*VR1(4,3)&
   +r1x2Z*(VS0(1,3)-ExZpE*VS1(1,3))
VS0(10,4)=PAz*VS0(4,4)+WPz*VS1(4,4)+r1x2Z*VR1(4,4)&
   +r1x2Z*(VS0(1,4)-ExZpE*VS1(1,4))&
   +HfxZpE*VS1(4,1)
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRd0p0'
STOP
END SELECT
END SUBROUTINE MVRRd0p0
