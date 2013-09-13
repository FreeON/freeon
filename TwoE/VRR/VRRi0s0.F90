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
   SUBROUTINE VRRi0s0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2Z*VRR0(23,1)
      V(2)=3.D0*V(1)
      V(3)=ExZpE*r1x2Z*VRR1(23,1)
      V(4)=-3.D0*V(3)
      V(5)=r1x2Z*VRR0(24,1)
      V(6)=ExZpE*r1x2Z*VRR1(24,1)
      V(7)=r1x2Z*VRR0(27,1)
      V(8)=ExZpE*r1x2Z*VRR1(27,1)
      V(9)=r1x2Z*VRR0(28,1)
      V(10)=2.D0*V(9)
      V(11)=ExZpE*r1x2Z*VRR1(28,1)
      V(12)=-2.D0*V(11)
      V(13)=2.D0*V(7)
      V(14)=-2.D0*V(8)
      V(15)=r1x2Z*VRR0(30,1)
      V(16)=3.D0*V(15)
      V(17)=ExZpE*r1x2Z*VRR1(30,1)
      V(18)=-3.D0*V(17)
      V(19)=r1x2Z*VRR0(31,1)
      V(20)=2.D0*V(19)
      V(21)=ExZpE*r1x2Z*VRR1(31,1)
      V(22)=-2.D0*V(21)
      V(23)=r1x2Z*VRR0(32,1)
      V(24)=ExZpE*r1x2Z*VRR1(32,1)
      V(25)=3.D0*V(23)
      V(26)=-3.D0*V(24)
      V(27)=r1x2Z*VRR0(33,1)
      V(28)=ExZpE*r1x2Z*VRR1(33,1)
      V(29)=r1x2Z*VRR0(34,1)
      V(30)=ExZpE*r1x2Z*VRR1(34,1)
      VRR0(57,1)=5.D0*r1x2Z*VRR0(21,1)+PAx*VRR0(36,1)-5.D0*ExZpE*r1x2Z*VRR1(21,1)+WPx*VRR1(36,1)
      VRR0(58,1)=4.D0*r1x2Z*VRR0(22,1)+PAx*VRR0(37,1)-4.D0*ExZpE*r1x2Z*VRR1(22,1)+WPx*VRR1(37,1)
      VRR0(59,1)=V(2)+V(4)+PAx*VRR0(38,1)+WPx*VRR1(38,1)
      VRR0(60,1)=2.D0*V(5)-2.D0*V(6)+PAx*VRR0(39,1)+WPx*VRR1(39,1)
      VRR0(61,1)=V(2)+V(4)+PAy*VRR0(39,1)+WPy*VRR1(39,1)
      VRR0(62,1)=4.D0*V(5)-4.D0*V(6)+PAy*VRR0(40,1)+WPy*VRR1(40,1)
      VRR0(63,1)=5.D0*r1x2Z*VRR0(25,1)+PAy*VRR0(41,1)-5.D0*ExZpE*r1x2Z*VRR1(25,1)+WPy*VRR1(41,1)
      VRR0(64,1)=4.D0*r1x2Z*VRR0(26,1)+PAx*VRR0(42,1)-4.D0*ExZpE*r1x2Z*VRR1(26,1)+WPx*VRR1(42,1)
      VRR0(65,1)=3.D0*V(7)-3.D0*V(8)+PAx*VRR0(43,1)+WPx*VRR1(43,1)
      VRR0(66,1)=V(10)+V(12)+PAx*VRR0(44,1)+WPx*VRR1(44,1)
      VRR0(67,1)=V(13)+V(14)+PAy*VRR0(44,1)+WPy*VRR1(44,1)
      VRR0(68,1)=3.D0*V(9)-3.D0*V(11)+PAy*VRR0(45,1)+WPy*VRR1(45,1)
      VRR0(69,1)=4.D0*r1x2Z*VRR0(29,1)+PAy*VRR0(46,1)-4.D0*ExZpE*r1x2Z*VRR1(29,1)+WPy*VRR1(46,1)
      VRR0(70,1)=V(16)+V(18)+PAx*VRR0(47,1)+WPx*VRR1(47,1)
      VRR0(71,1)=V(20)+V(22)+PAx*VRR0(48,1)+WPx*VRR1(48,1)
      VRR0(72,1)=V(23)-V(24)+PAx*VRR0(49,1)+WPx*VRR1(49,1)
      VRR0(73,1)=V(20)+V(22)+PAy*VRR0(49,1)+WPy*VRR1(49,1)
      VRR0(74,1)=V(25)+V(26)+PAy*VRR0(50,1)+WPy*VRR1(50,1)
      VRR0(75,1)=2.D0*V(27)-2.D0*V(28)+PAx*VRR0(51,1)+WPx*VRR1(51,1)
      VRR0(76,1)=V(13)+V(14)+PAz*VRR0(48,1)+WPz*VRR1(48,1)
      VRR0(77,1)=V(10)+V(12)+PAz*VRR0(49,1)+WPz*VRR1(49,1)
      VRR0(78,1)=2.D0*V(29)-2.D0*V(30)+PAy*VRR0(53,1)+WPy*VRR1(53,1)
      VRR0(79,1)=V(16)+V(18)+PAz*VRR0(51,1)+WPz*VRR1(51,1)
      VRR0(80,1)=3.D0*V(19)-3.D0*V(21)+PAz*VRR0(52,1)+WPz*VRR1(52,1)
      VRR0(81,1)=V(25)+V(26)+PAz*VRR0(53,1)+WPz*VRR1(53,1)
      VRR0(82,1)=4.D0*V(27)-4.D0*V(28)+PAz*VRR0(54,1)+WPz*VRR1(54,1)
      VRR0(83,1)=4.D0*V(29)-4.D0*V(30)+PAz*VRR0(55,1)+WPz*VRR1(55,1)
      VRR0(84,1)=5.D0*r1x2Z*VRR0(35,1)+PAz*VRR0(56,1)-5.D0*ExZpE*r1x2Z*VRR1(35,1)+WPz*VRR1(56,1)
END SUBROUTINE VRRi0s0
SUBROUTINE MVRRi0s0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(57,1)=PAx*VS0(36,1)+WPx*VS1(36,1)+r1x2Z*VR1(36,1)&
   +5D0*r1x2Z*(VS0(21,1)-ExZpE*VS1(21,1))
VS0(58,1)=PAx*VS0(37,1)+WPx*VS1(37,1)+r1x2Z*VR1(37,1)&
   +4D0*r1x2Z*(VS0(22,1)-ExZpE*VS1(22,1))
VS0(59,1)=PAx*VS0(38,1)+WPx*VS1(38,1)+r1x2Z*VR1(38,1)&
   +3D0*r1x2Z*(VS0(23,1)-ExZpE*VS1(23,1))
VS0(60,1)=PAx*VS0(39,1)+WPx*VS1(39,1)+r1x2Z*VR1(39,1)&
   +2D0*r1x2Z*(VS0(24,1)-ExZpE*VS1(24,1))
VS0(61,1)=PAy*VS0(39,1)+WPy*VS1(39,1)&
   +3D0*r1x2Z*(VS0(23,1)-ExZpE*VS1(23,1))
VS0(62,1)=PAy*VS0(40,1)+WPy*VS1(40,1)&
   +4D0*r1x2Z*(VS0(24,1)-ExZpE*VS1(24,1))
VS0(63,1)=PAy*VS0(41,1)+WPy*VS1(41,1)&
   +5D0*r1x2Z*(VS0(25,1)-ExZpE*VS1(25,1))
VS0(64,1)=PAx*VS0(42,1)+WPx*VS1(42,1)+r1x2Z*VR1(42,1)&
   +4D0*r1x2Z*(VS0(26,1)-ExZpE*VS1(26,1))
VS0(65,1)=PAx*VS0(43,1)+WPx*VS1(43,1)+r1x2Z*VR1(43,1)&
   +3D0*r1x2Z*(VS0(27,1)-ExZpE*VS1(27,1))
VS0(66,1)=PAx*VS0(44,1)+WPx*VS1(44,1)+r1x2Z*VR1(44,1)&
   +2D0*r1x2Z*(VS0(28,1)-ExZpE*VS1(28,1))
VS0(67,1)=PAy*VS0(44,1)+WPy*VS1(44,1)&
   +2D0*r1x2Z*(VS0(27,1)-ExZpE*VS1(27,1))
VS0(68,1)=PAy*VS0(45,1)+WPy*VS1(45,1)&
   +3D0*r1x2Z*(VS0(28,1)-ExZpE*VS1(28,1))
VS0(69,1)=PAy*VS0(46,1)+WPy*VS1(46,1)&
   +4D0*r1x2Z*(VS0(29,1)-ExZpE*VS1(29,1))
VS0(70,1)=PAx*VS0(47,1)+WPx*VS1(47,1)+r1x2Z*VR1(47,1)&
   +3D0*r1x2Z*(VS0(30,1)-ExZpE*VS1(30,1))
VS0(71,1)=PAx*VS0(48,1)+WPx*VS1(48,1)+r1x2Z*VR1(48,1)&
   +2D0*r1x2Z*(VS0(31,1)-ExZpE*VS1(31,1))
VS0(72,1)=PAx*VS0(49,1)+WPx*VS1(49,1)+r1x2Z*VR1(49,1)&
   +r1x2Z*(VS0(32,1)-ExZpE*VS1(32,1))
VS0(73,1)=PAy*VS0(49,1)+WPy*VS1(49,1)&
   +2D0*r1x2Z*(VS0(31,1)-ExZpE*VS1(31,1))
VS0(74,1)=PAy*VS0(50,1)+WPy*VS1(50,1)&
   +3D0*r1x2Z*(VS0(32,1)-ExZpE*VS1(32,1))
VS0(75,1)=PAx*VS0(51,1)+WPx*VS1(51,1)+r1x2Z*VR1(51,1)&
   +2D0*r1x2Z*(VS0(33,1)-ExZpE*VS1(33,1))
VS0(76,1)=PAz*VS0(48,1)+WPz*VS1(48,1)&
   +2D0*r1x2Z*(VS0(27,1)-ExZpE*VS1(27,1))
VS0(77,1)=PAz*VS0(49,1)+WPz*VS1(49,1)&
   +2D0*r1x2Z*(VS0(28,1)-ExZpE*VS1(28,1))
VS0(78,1)=PAy*VS0(53,1)+WPy*VS1(53,1)&
   +2D0*r1x2Z*(VS0(34,1)-ExZpE*VS1(34,1))
VS0(79,1)=PAz*VS0(51,1)+WPz*VS1(51,1)&
   +3D0*r1x2Z*(VS0(30,1)-ExZpE*VS1(30,1))
VS0(80,1)=PAz*VS0(52,1)+WPz*VS1(52,1)&
   +3D0*r1x2Z*(VS0(31,1)-ExZpE*VS1(31,1))
VS0(81,1)=PAz*VS0(53,1)+WPz*VS1(53,1)&
   +3D0*r1x2Z*(VS0(32,1)-ExZpE*VS1(32,1))
VS0(82,1)=PAz*VS0(54,1)+WPz*VS1(54,1)&
   +4D0*r1x2Z*(VS0(33,1)-ExZpE*VS1(33,1))
VS0(83,1)=PAz*VS0(55,1)+WPz*VS1(55,1)&
   +4D0*r1x2Z*(VS0(34,1)-ExZpE*VS1(34,1))
VS0(84,1)=PAz*VS0(56,1)+WPz*VS1(56,1)&
   +5D0*r1x2Z*(VS0(35,1)-ExZpE*VS1(35,1))
CASE(2)
VS0(57,1)=PAx*VS0(36,1)+WPx*VS1(36,1)&
   +5D0*r1x2Z*(VS0(21,1)-ExZpE*VS1(21,1))
VS0(58,1)=PAx*VS0(37,1)+WPx*VS1(37,1)&
   +4D0*r1x2Z*(VS0(22,1)-ExZpE*VS1(22,1))
VS0(59,1)=PAx*VS0(38,1)+WPx*VS1(38,1)&
   +3D0*r1x2Z*(VS0(23,1)-ExZpE*VS1(23,1))
VS0(60,1)=PAx*VS0(39,1)+WPx*VS1(39,1)&
   +2D0*r1x2Z*(VS0(24,1)-ExZpE*VS1(24,1))
VS0(61,1)=PAy*VS0(39,1)+WPy*VS1(39,1)+r1x2Z*VR1(39,1)&
   +3D0*r1x2Z*(VS0(23,1)-ExZpE*VS1(23,1))
VS0(62,1)=PAy*VS0(40,1)+WPy*VS1(40,1)+r1x2Z*VR1(40,1)&
   +4D0*r1x2Z*(VS0(24,1)-ExZpE*VS1(24,1))
VS0(63,1)=PAy*VS0(41,1)+WPy*VS1(41,1)+r1x2Z*VR1(41,1)&
   +5D0*r1x2Z*(VS0(25,1)-ExZpE*VS1(25,1))
VS0(64,1)=PAx*VS0(42,1)+WPx*VS1(42,1)&
   +4D0*r1x2Z*(VS0(26,1)-ExZpE*VS1(26,1))
VS0(65,1)=PAx*VS0(43,1)+WPx*VS1(43,1)&
   +3D0*r1x2Z*(VS0(27,1)-ExZpE*VS1(27,1))
VS0(66,1)=PAx*VS0(44,1)+WPx*VS1(44,1)&
   +2D0*r1x2Z*(VS0(28,1)-ExZpE*VS1(28,1))
VS0(67,1)=PAy*VS0(44,1)+WPy*VS1(44,1)+r1x2Z*VR1(44,1)&
   +2D0*r1x2Z*(VS0(27,1)-ExZpE*VS1(27,1))
VS0(68,1)=PAy*VS0(45,1)+WPy*VS1(45,1)+r1x2Z*VR1(45,1)&
   +3D0*r1x2Z*(VS0(28,1)-ExZpE*VS1(28,1))
VS0(69,1)=PAy*VS0(46,1)+WPy*VS1(46,1)+r1x2Z*VR1(46,1)&
   +4D0*r1x2Z*(VS0(29,1)-ExZpE*VS1(29,1))
VS0(70,1)=PAx*VS0(47,1)+WPx*VS1(47,1)&
   +3D0*r1x2Z*(VS0(30,1)-ExZpE*VS1(30,1))
VS0(71,1)=PAx*VS0(48,1)+WPx*VS1(48,1)&
   +2D0*r1x2Z*(VS0(31,1)-ExZpE*VS1(31,1))
VS0(72,1)=PAx*VS0(49,1)+WPx*VS1(49,1)&
   +r1x2Z*(VS0(32,1)-ExZpE*VS1(32,1))
VS0(73,1)=PAy*VS0(49,1)+WPy*VS1(49,1)+r1x2Z*VR1(49,1)&
   +2D0*r1x2Z*(VS0(31,1)-ExZpE*VS1(31,1))
VS0(74,1)=PAy*VS0(50,1)+WPy*VS1(50,1)+r1x2Z*VR1(50,1)&
   +3D0*r1x2Z*(VS0(32,1)-ExZpE*VS1(32,1))
VS0(75,1)=PAx*VS0(51,1)+WPx*VS1(51,1)&
   +2D0*r1x2Z*(VS0(33,1)-ExZpE*VS1(33,1))
VS0(76,1)=PAz*VS0(48,1)+WPz*VS1(48,1)&
   +2D0*r1x2Z*(VS0(27,1)-ExZpE*VS1(27,1))
VS0(77,1)=PAz*VS0(49,1)+WPz*VS1(49,1)&
   +2D0*r1x2Z*(VS0(28,1)-ExZpE*VS1(28,1))
VS0(78,1)=PAy*VS0(53,1)+WPy*VS1(53,1)+r1x2Z*VR1(53,1)&
   +2D0*r1x2Z*(VS0(34,1)-ExZpE*VS1(34,1))
VS0(79,1)=PAz*VS0(51,1)+WPz*VS1(51,1)&
   +3D0*r1x2Z*(VS0(30,1)-ExZpE*VS1(30,1))
VS0(80,1)=PAz*VS0(52,1)+WPz*VS1(52,1)&
   +3D0*r1x2Z*(VS0(31,1)-ExZpE*VS1(31,1))
VS0(81,1)=PAz*VS0(53,1)+WPz*VS1(53,1)&
   +3D0*r1x2Z*(VS0(32,1)-ExZpE*VS1(32,1))
VS0(82,1)=PAz*VS0(54,1)+WPz*VS1(54,1)&
   +4D0*r1x2Z*(VS0(33,1)-ExZpE*VS1(33,1))
VS0(83,1)=PAz*VS0(55,1)+WPz*VS1(55,1)&
   +4D0*r1x2Z*(VS0(34,1)-ExZpE*VS1(34,1))
VS0(84,1)=PAz*VS0(56,1)+WPz*VS1(56,1)&
   +5D0*r1x2Z*(VS0(35,1)-ExZpE*VS1(35,1))
CASE(3)
VS0(57,1)=PAx*VS0(36,1)+WPx*VS1(36,1)&
   +5D0*r1x2Z*(VS0(21,1)-ExZpE*VS1(21,1))
VS0(58,1)=PAx*VS0(37,1)+WPx*VS1(37,1)&
   +4D0*r1x2Z*(VS0(22,1)-ExZpE*VS1(22,1))
VS0(59,1)=PAx*VS0(38,1)+WPx*VS1(38,1)&
   +3D0*r1x2Z*(VS0(23,1)-ExZpE*VS1(23,1))
VS0(60,1)=PAx*VS0(39,1)+WPx*VS1(39,1)&
   +2D0*r1x2Z*(VS0(24,1)-ExZpE*VS1(24,1))
VS0(61,1)=PAy*VS0(39,1)+WPy*VS1(39,1)&
   +3D0*r1x2Z*(VS0(23,1)-ExZpE*VS1(23,1))
VS0(62,1)=PAy*VS0(40,1)+WPy*VS1(40,1)&
   +4D0*r1x2Z*(VS0(24,1)-ExZpE*VS1(24,1))
VS0(63,1)=PAy*VS0(41,1)+WPy*VS1(41,1)&
   +5D0*r1x2Z*(VS0(25,1)-ExZpE*VS1(25,1))
VS0(64,1)=PAx*VS0(42,1)+WPx*VS1(42,1)&
   +4D0*r1x2Z*(VS0(26,1)-ExZpE*VS1(26,1))
VS0(65,1)=PAx*VS0(43,1)+WPx*VS1(43,1)&
   +3D0*r1x2Z*(VS0(27,1)-ExZpE*VS1(27,1))
VS0(66,1)=PAx*VS0(44,1)+WPx*VS1(44,1)&
   +2D0*r1x2Z*(VS0(28,1)-ExZpE*VS1(28,1))
VS0(67,1)=PAy*VS0(44,1)+WPy*VS1(44,1)&
   +2D0*r1x2Z*(VS0(27,1)-ExZpE*VS1(27,1))
VS0(68,1)=PAy*VS0(45,1)+WPy*VS1(45,1)&
   +3D0*r1x2Z*(VS0(28,1)-ExZpE*VS1(28,1))
VS0(69,1)=PAy*VS0(46,1)+WPy*VS1(46,1)&
   +4D0*r1x2Z*(VS0(29,1)-ExZpE*VS1(29,1))
VS0(70,1)=PAx*VS0(47,1)+WPx*VS1(47,1)&
   +3D0*r1x2Z*(VS0(30,1)-ExZpE*VS1(30,1))
VS0(71,1)=PAx*VS0(48,1)+WPx*VS1(48,1)&
   +2D0*r1x2Z*(VS0(31,1)-ExZpE*VS1(31,1))
VS0(72,1)=PAx*VS0(49,1)+WPx*VS1(49,1)&
   +r1x2Z*(VS0(32,1)-ExZpE*VS1(32,1))
VS0(73,1)=PAy*VS0(49,1)+WPy*VS1(49,1)&
   +2D0*r1x2Z*(VS0(31,1)-ExZpE*VS1(31,1))
VS0(74,1)=PAy*VS0(50,1)+WPy*VS1(50,1)&
   +3D0*r1x2Z*(VS0(32,1)-ExZpE*VS1(32,1))
VS0(75,1)=PAx*VS0(51,1)+WPx*VS1(51,1)&
   +2D0*r1x2Z*(VS0(33,1)-ExZpE*VS1(33,1))
VS0(76,1)=PAz*VS0(48,1)+WPz*VS1(48,1)+r1x2Z*VR1(48,1)&
   +2D0*r1x2Z*(VS0(27,1)-ExZpE*VS1(27,1))
VS0(77,1)=PAz*VS0(49,1)+WPz*VS1(49,1)+r1x2Z*VR1(49,1)&
   +2D0*r1x2Z*(VS0(28,1)-ExZpE*VS1(28,1))
VS0(78,1)=PAy*VS0(53,1)+WPy*VS1(53,1)&
   +2D0*r1x2Z*(VS0(34,1)-ExZpE*VS1(34,1))
VS0(79,1)=PAz*VS0(51,1)+WPz*VS1(51,1)+r1x2Z*VR1(51,1)&
   +3D0*r1x2Z*(VS0(30,1)-ExZpE*VS1(30,1))
VS0(80,1)=PAz*VS0(52,1)+WPz*VS1(52,1)+r1x2Z*VR1(52,1)&
   +3D0*r1x2Z*(VS0(31,1)-ExZpE*VS1(31,1))
VS0(81,1)=PAz*VS0(53,1)+WPz*VS1(53,1)+r1x2Z*VR1(53,1)&
   +3D0*r1x2Z*(VS0(32,1)-ExZpE*VS1(32,1))
VS0(82,1)=PAz*VS0(54,1)+WPz*VS1(54,1)+r1x2Z*VR1(54,1)&
   +4D0*r1x2Z*(VS0(33,1)-ExZpE*VS1(33,1))
VS0(83,1)=PAz*VS0(55,1)+WPz*VS1(55,1)+r1x2Z*VR1(55,1)&
   +4D0*r1x2Z*(VS0(34,1)-ExZpE*VS1(34,1))
VS0(84,1)=PAz*VS0(56,1)+WPz*VS1(56,1)+r1x2Z*VR1(56,1)&
   +5D0*r1x2Z*(VS0(35,1)-ExZpE*VS1(35,1))
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRi0s0'
STOP
END SELECT
END SUBROUTINE MVRRi0s0
