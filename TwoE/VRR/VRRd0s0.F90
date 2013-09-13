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
   SUBROUTINE VRRd0s0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2Z*VRR0(1,1)
      V(2)=ExZpE*r1x2Z*VRR1(1,1)
      V(3)=-V(2)
      VRR0(5,1)=V(1)+V(3)+PAx*VRR0(2,1)+WPx*VRR1(2,1)
      VRR0(6,1)=PAx*VRR0(3,1)+WPx*VRR1(3,1)
      VRR0(7,1)=V(1)+V(3)+PAy*VRR0(3,1)+WPy*VRR1(3,1)
      VRR0(8,1)=PAx*VRR0(4,1)+WPx*VRR1(4,1)
      VRR0(9,1)=PAy*VRR0(4,1)+WPy*VRR1(4,1)
      VRR0(10,1)=V(1)+V(3)+PAz*VRR0(4,1)+WPz*VRR1(4,1)
END SUBROUTINE VRRd0s0
SUBROUTINE MVRRd0s0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(5,1)=PAx*VS0(2,1)+WPx*VS1(2,1)+r1x2Z*VR1(2,1)&
   +r1x2Z*(VS0(1,1)-ExZpE*VS1(1,1))
VS0(6,1)=PAx*VS0(3,1)+WPx*VS1(3,1)+r1x2Z*VR1(3,1)
VS0(7,1)=PAy*VS0(3,1)+WPy*VS1(3,1)&
   +r1x2Z*(VS0(1,1)-ExZpE*VS1(1,1))
VS0(8,1)=PAx*VS0(4,1)+WPx*VS1(4,1)+r1x2Z*VR1(4,1)
VS0(9,1)=PAy*VS0(4,1)+WPy*VS1(4,1)
VS0(10,1)=PAz*VS0(4,1)+WPz*VS1(4,1)&
   +r1x2Z*(VS0(1,1)-ExZpE*VS1(1,1))
CASE(2)
VS0(5,1)=PAx*VS0(2,1)+WPx*VS1(2,1)&
   +r1x2Z*(VS0(1,1)-ExZpE*VS1(1,1))
VS0(6,1)=PAx*VS0(3,1)+WPx*VS1(3,1)
VS0(7,1)=PAy*VS0(3,1)+WPy*VS1(3,1)+r1x2Z*VR1(3,1)&
   +r1x2Z*(VS0(1,1)-ExZpE*VS1(1,1))
VS0(8,1)=PAx*VS0(4,1)+WPx*VS1(4,1)
VS0(9,1)=PAy*VS0(4,1)+WPy*VS1(4,1)+r1x2Z*VR1(4,1)
VS0(10,1)=PAz*VS0(4,1)+WPz*VS1(4,1)&
   +r1x2Z*(VS0(1,1)-ExZpE*VS1(1,1))
CASE(3)
VS0(5,1)=PAx*VS0(2,1)+WPx*VS1(2,1)&
   +r1x2Z*(VS0(1,1)-ExZpE*VS1(1,1))
VS0(6,1)=PAx*VS0(3,1)+WPx*VS1(3,1)
VS0(7,1)=PAy*VS0(3,1)+WPy*VS1(3,1)&
   +r1x2Z*(VS0(1,1)-ExZpE*VS1(1,1))
VS0(8,1)=PAx*VS0(4,1)+WPx*VS1(4,1)
VS0(9,1)=PAy*VS0(4,1)+WPy*VS1(4,1)
VS0(10,1)=PAz*VS0(4,1)+WPz*VS1(4,1)+r1x2Z*VR1(4,1)&
   +r1x2Z*(VS0(1,1)-ExZpE*VS1(1,1))
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRd0s0'
STOP
END SELECT
END SUBROUTINE MVRRd0s0
