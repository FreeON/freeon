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
   SUBROUTINE VRRh0s0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2Z*VRR0(12,1)
      V(2)=ExZpE*r1x2Z*VRR1(12,1)
      V(3)=r1x2Z*VRR0(13,1)
      V(4)=ExZpE*r1x2Z*VRR1(13,1)
      V(5)=r1x2Z*VRR0(15,1)
      V(6)=ExZpE*r1x2Z*VRR1(15,1)
      V(7)=r1x2Z*VRR0(16,1)
      V(8)=2.D0*V(7)
      V(9)=ExZpE*r1x2Z*VRR1(16,1)
      V(10)=-2.D0*V(9)
      V(11)=r1x2Z*VRR0(17,1)
      V(12)=ExZpE*r1x2Z*VRR1(17,1)
      V(13)=r1x2Z*VRR0(18,1)
      V(14)=ExZpE*r1x2Z*VRR1(18,1)
      V(15)=r1x2Z*VRR0(19,1)
      V(16)=ExZpE*r1x2Z*VRR1(19,1)
      VRR0(36,1)=4.D0*r1x2Z*VRR0(11,1)+PAx*VRR0(21,1)-4.D0*ExZpE*r1x2Z*VRR1(11,1)+WPx*VRR1(21,1)
      VRR0(37,1)=3.D0*V(1)-3.D0*V(2)+PAx*VRR0(22,1)+WPx*VRR1(22,1)
      VRR0(38,1)=2.D0*V(3)-2.D0*V(4)+PAx*VRR0(23,1)+WPx*VRR1(23,1)
      VRR0(39,1)=2.D0*V(1)-2.D0*V(2)+PAy*VRR0(23,1)+WPy*VRR1(23,1)
      VRR0(40,1)=3.D0*V(3)-3.D0*V(4)+PAy*VRR0(24,1)+WPy*VRR1(24,1)
      VRR0(41,1)=4.D0*r1x2Z*VRR0(14,1)+PAy*VRR0(25,1)-4.D0*ExZpE*r1x2Z*VRR1(14,1)+WPy*VRR1(25,1)
      VRR0(42,1)=3.D0*V(5)-3.D0*V(6)+PAx*VRR0(26,1)+WPx*VRR1(26,1)
      VRR0(43,1)=V(8)+V(10)+PAx*VRR0(27,1)+WPx*VRR1(27,1)
      VRR0(44,1)=V(11)-V(12)+PAx*VRR0(28,1)+WPx*VRR1(28,1)
      VRR0(45,1)=V(8)+V(10)+PAy*VRR0(28,1)+WPy*VRR1(28,1)
      VRR0(46,1)=3.D0*V(11)-3.D0*V(12)+PAy*VRR0(29,1)+WPy*VRR1(29,1)
      VRR0(47,1)=2.D0*V(13)-2.D0*V(14)+PAx*VRR0(30,1)+WPx*VRR1(30,1)
      VRR0(48,1)=V(15)-V(16)+PAx*VRR0(31,1)+WPx*VRR1(31,1)
      VRR0(49,1)=V(13)-V(14)+PAy*VRR0(31,1)+WPy*VRR1(31,1)
      VRR0(50,1)=2.D0*V(15)-2.D0*V(16)+PAy*VRR0(32,1)+WPy*VRR1(32,1)
      VRR0(51,1)=2.D0*V(5)-2.D0*V(6)+PAz*VRR0(30,1)+WPz*VRR1(30,1)
      VRR0(52,1)=V(8)+V(10)+PAz*VRR0(31,1)+WPz*VRR1(31,1)
      VRR0(53,1)=2.D0*V(11)-2.D0*V(12)+PAz*VRR0(32,1)+WPz*VRR1(32,1)
      VRR0(54,1)=3.D0*V(13)-3.D0*V(14)+PAz*VRR0(33,1)+WPz*VRR1(33,1)
      VRR0(55,1)=3.D0*V(15)-3.D0*V(16)+PAz*VRR0(34,1)+WPz*VRR1(34,1)
      VRR0(56,1)=4.D0*r1x2Z*VRR0(20,1)+PAz*VRR0(35,1)-4.D0*ExZpE*r1x2Z*VRR1(20,1)+WPz*VRR1(35,1)
END SUBROUTINE VRRh0s0
SUBROUTINE MVRRh0s0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(36,1)=PAx*VS0(21,1)+WPx*VS1(21,1)+r1x2Z*VR1(21,1)&
   +4D0*r1x2Z*(VS0(11,1)-ExZpE*VS1(11,1))
VS0(37,1)=PAx*VS0(22,1)+WPx*VS1(22,1)+r1x2Z*VR1(22,1)&
   +3D0*r1x2Z*(VS0(12,1)-ExZpE*VS1(12,1))
VS0(38,1)=PAx*VS0(23,1)+WPx*VS1(23,1)+r1x2Z*VR1(23,1)&
   +2D0*r1x2Z*(VS0(13,1)-ExZpE*VS1(13,1))
VS0(39,1)=PAy*VS0(23,1)+WPy*VS1(23,1)&
   +2D0*r1x2Z*(VS0(12,1)-ExZpE*VS1(12,1))
VS0(40,1)=PAy*VS0(24,1)+WPy*VS1(24,1)&
   +3D0*r1x2Z*(VS0(13,1)-ExZpE*VS1(13,1))
VS0(41,1)=PAy*VS0(25,1)+WPy*VS1(25,1)&
   +4D0*r1x2Z*(VS0(14,1)-ExZpE*VS1(14,1))
VS0(42,1)=PAx*VS0(26,1)+WPx*VS1(26,1)+r1x2Z*VR1(26,1)&
   +3D0*r1x2Z*(VS0(15,1)-ExZpE*VS1(15,1))
VS0(43,1)=PAx*VS0(27,1)+WPx*VS1(27,1)+r1x2Z*VR1(27,1)&
   +2D0*r1x2Z*(VS0(16,1)-ExZpE*VS1(16,1))
VS0(44,1)=PAx*VS0(28,1)+WPx*VS1(28,1)+r1x2Z*VR1(28,1)&
   +r1x2Z*(VS0(17,1)-ExZpE*VS1(17,1))
VS0(45,1)=PAy*VS0(28,1)+WPy*VS1(28,1)&
   +2D0*r1x2Z*(VS0(16,1)-ExZpE*VS1(16,1))
VS0(46,1)=PAy*VS0(29,1)+WPy*VS1(29,1)&
   +3D0*r1x2Z*(VS0(17,1)-ExZpE*VS1(17,1))
VS0(47,1)=PAx*VS0(30,1)+WPx*VS1(30,1)+r1x2Z*VR1(30,1)&
   +2D0*r1x2Z*(VS0(18,1)-ExZpE*VS1(18,1))
VS0(48,1)=PAx*VS0(31,1)+WPx*VS1(31,1)+r1x2Z*VR1(31,1)&
   +r1x2Z*(VS0(19,1)-ExZpE*VS1(19,1))
VS0(49,1)=PAy*VS0(31,1)+WPy*VS1(31,1)&
   +r1x2Z*(VS0(18,1)-ExZpE*VS1(18,1))
VS0(50,1)=PAy*VS0(32,1)+WPy*VS1(32,1)&
   +2D0*r1x2Z*(VS0(19,1)-ExZpE*VS1(19,1))
VS0(51,1)=PAz*VS0(30,1)+WPz*VS1(30,1)&
   +2D0*r1x2Z*(VS0(15,1)-ExZpE*VS1(15,1))
VS0(52,1)=PAz*VS0(31,1)+WPz*VS1(31,1)&
   +2D0*r1x2Z*(VS0(16,1)-ExZpE*VS1(16,1))
VS0(53,1)=PAz*VS0(32,1)+WPz*VS1(32,1)&
   +2D0*r1x2Z*(VS0(17,1)-ExZpE*VS1(17,1))
VS0(54,1)=PAz*VS0(33,1)+WPz*VS1(33,1)&
   +3D0*r1x2Z*(VS0(18,1)-ExZpE*VS1(18,1))
VS0(55,1)=PAz*VS0(34,1)+WPz*VS1(34,1)&
   +3D0*r1x2Z*(VS0(19,1)-ExZpE*VS1(19,1))
VS0(56,1)=PAz*VS0(35,1)+WPz*VS1(35,1)&
   +4D0*r1x2Z*(VS0(20,1)-ExZpE*VS1(20,1))
CASE(2)
VS0(36,1)=PAx*VS0(21,1)+WPx*VS1(21,1)&
   +4D0*r1x2Z*(VS0(11,1)-ExZpE*VS1(11,1))
VS0(37,1)=PAx*VS0(22,1)+WPx*VS1(22,1)&
   +3D0*r1x2Z*(VS0(12,1)-ExZpE*VS1(12,1))
VS0(38,1)=PAx*VS0(23,1)+WPx*VS1(23,1)&
   +2D0*r1x2Z*(VS0(13,1)-ExZpE*VS1(13,1))
VS0(39,1)=PAy*VS0(23,1)+WPy*VS1(23,1)+r1x2Z*VR1(23,1)&
   +2D0*r1x2Z*(VS0(12,1)-ExZpE*VS1(12,1))
VS0(40,1)=PAy*VS0(24,1)+WPy*VS1(24,1)+r1x2Z*VR1(24,1)&
   +3D0*r1x2Z*(VS0(13,1)-ExZpE*VS1(13,1))
VS0(41,1)=PAy*VS0(25,1)+WPy*VS1(25,1)+r1x2Z*VR1(25,1)&
   +4D0*r1x2Z*(VS0(14,1)-ExZpE*VS1(14,1))
VS0(42,1)=PAx*VS0(26,1)+WPx*VS1(26,1)&
   +3D0*r1x2Z*(VS0(15,1)-ExZpE*VS1(15,1))
VS0(43,1)=PAx*VS0(27,1)+WPx*VS1(27,1)&
   +2D0*r1x2Z*(VS0(16,1)-ExZpE*VS1(16,1))
VS0(44,1)=PAx*VS0(28,1)+WPx*VS1(28,1)&
   +r1x2Z*(VS0(17,1)-ExZpE*VS1(17,1))
VS0(45,1)=PAy*VS0(28,1)+WPy*VS1(28,1)+r1x2Z*VR1(28,1)&
   +2D0*r1x2Z*(VS0(16,1)-ExZpE*VS1(16,1))
VS0(46,1)=PAy*VS0(29,1)+WPy*VS1(29,1)+r1x2Z*VR1(29,1)&
   +3D0*r1x2Z*(VS0(17,1)-ExZpE*VS1(17,1))
VS0(47,1)=PAx*VS0(30,1)+WPx*VS1(30,1)&
   +2D0*r1x2Z*(VS0(18,1)-ExZpE*VS1(18,1))
VS0(48,1)=PAx*VS0(31,1)+WPx*VS1(31,1)&
   +r1x2Z*(VS0(19,1)-ExZpE*VS1(19,1))
VS0(49,1)=PAy*VS0(31,1)+WPy*VS1(31,1)+r1x2Z*VR1(31,1)&
   +r1x2Z*(VS0(18,1)-ExZpE*VS1(18,1))
VS0(50,1)=PAy*VS0(32,1)+WPy*VS1(32,1)+r1x2Z*VR1(32,1)&
   +2D0*r1x2Z*(VS0(19,1)-ExZpE*VS1(19,1))
VS0(51,1)=PAz*VS0(30,1)+WPz*VS1(30,1)&
   +2D0*r1x2Z*(VS0(15,1)-ExZpE*VS1(15,1))
VS0(52,1)=PAz*VS0(31,1)+WPz*VS1(31,1)&
   +2D0*r1x2Z*(VS0(16,1)-ExZpE*VS1(16,1))
VS0(53,1)=PAz*VS0(32,1)+WPz*VS1(32,1)&
   +2D0*r1x2Z*(VS0(17,1)-ExZpE*VS1(17,1))
VS0(54,1)=PAz*VS0(33,1)+WPz*VS1(33,1)&
   +3D0*r1x2Z*(VS0(18,1)-ExZpE*VS1(18,1))
VS0(55,1)=PAz*VS0(34,1)+WPz*VS1(34,1)&
   +3D0*r1x2Z*(VS0(19,1)-ExZpE*VS1(19,1))
VS0(56,1)=PAz*VS0(35,1)+WPz*VS1(35,1)&
   +4D0*r1x2Z*(VS0(20,1)-ExZpE*VS1(20,1))
CASE(3)
VS0(36,1)=PAx*VS0(21,1)+WPx*VS1(21,1)&
   +4D0*r1x2Z*(VS0(11,1)-ExZpE*VS1(11,1))
VS0(37,1)=PAx*VS0(22,1)+WPx*VS1(22,1)&
   +3D0*r1x2Z*(VS0(12,1)-ExZpE*VS1(12,1))
VS0(38,1)=PAx*VS0(23,1)+WPx*VS1(23,1)&
   +2D0*r1x2Z*(VS0(13,1)-ExZpE*VS1(13,1))
VS0(39,1)=PAy*VS0(23,1)+WPy*VS1(23,1)&
   +2D0*r1x2Z*(VS0(12,1)-ExZpE*VS1(12,1))
VS0(40,1)=PAy*VS0(24,1)+WPy*VS1(24,1)&
   +3D0*r1x2Z*(VS0(13,1)-ExZpE*VS1(13,1))
VS0(41,1)=PAy*VS0(25,1)+WPy*VS1(25,1)&
   +4D0*r1x2Z*(VS0(14,1)-ExZpE*VS1(14,1))
VS0(42,1)=PAx*VS0(26,1)+WPx*VS1(26,1)&
   +3D0*r1x2Z*(VS0(15,1)-ExZpE*VS1(15,1))
VS0(43,1)=PAx*VS0(27,1)+WPx*VS1(27,1)&
   +2D0*r1x2Z*(VS0(16,1)-ExZpE*VS1(16,1))
VS0(44,1)=PAx*VS0(28,1)+WPx*VS1(28,1)&
   +r1x2Z*(VS0(17,1)-ExZpE*VS1(17,1))
VS0(45,1)=PAy*VS0(28,1)+WPy*VS1(28,1)&
   +2D0*r1x2Z*(VS0(16,1)-ExZpE*VS1(16,1))
VS0(46,1)=PAy*VS0(29,1)+WPy*VS1(29,1)&
   +3D0*r1x2Z*(VS0(17,1)-ExZpE*VS1(17,1))
VS0(47,1)=PAx*VS0(30,1)+WPx*VS1(30,1)&
   +2D0*r1x2Z*(VS0(18,1)-ExZpE*VS1(18,1))
VS0(48,1)=PAx*VS0(31,1)+WPx*VS1(31,1)&
   +r1x2Z*(VS0(19,1)-ExZpE*VS1(19,1))
VS0(49,1)=PAy*VS0(31,1)+WPy*VS1(31,1)&
   +r1x2Z*(VS0(18,1)-ExZpE*VS1(18,1))
VS0(50,1)=PAy*VS0(32,1)+WPy*VS1(32,1)&
   +2D0*r1x2Z*(VS0(19,1)-ExZpE*VS1(19,1))
VS0(51,1)=PAz*VS0(30,1)+WPz*VS1(30,1)+r1x2Z*VR1(30,1)&
   +2D0*r1x2Z*(VS0(15,1)-ExZpE*VS1(15,1))
VS0(52,1)=PAz*VS0(31,1)+WPz*VS1(31,1)+r1x2Z*VR1(31,1)&
   +2D0*r1x2Z*(VS0(16,1)-ExZpE*VS1(16,1))
VS0(53,1)=PAz*VS0(32,1)+WPz*VS1(32,1)+r1x2Z*VR1(32,1)&
   +2D0*r1x2Z*(VS0(17,1)-ExZpE*VS1(17,1))
VS0(54,1)=PAz*VS0(33,1)+WPz*VS1(33,1)+r1x2Z*VR1(33,1)&
   +3D0*r1x2Z*(VS0(18,1)-ExZpE*VS1(18,1))
VS0(55,1)=PAz*VS0(34,1)+WPz*VS1(34,1)+r1x2Z*VR1(34,1)&
   +3D0*r1x2Z*(VS0(19,1)-ExZpE*VS1(19,1))
VS0(56,1)=PAz*VS0(35,1)+WPz*VS1(35,1)+r1x2Z*VR1(35,1)&
   +4D0*r1x2Z*(VS0(20,1)-ExZpE*VS1(20,1))
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRh0s0'
STOP
END SELECT
END SUBROUTINE MVRRh0s0
