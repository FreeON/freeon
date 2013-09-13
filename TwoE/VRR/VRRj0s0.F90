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
   SUBROUTINE VRRj0s0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2Z*VRR0(38,1)
      V(2)=ExZpE*r1x2Z*VRR1(38,1)
      V(3)=r1x2Z*VRR0(39,1)
      V(4)=ExZpE*r1x2Z*VRR1(39,1)
      V(5)=r1x2Z*VRR0(44,1)
      V(6)=3.D0*V(5)
      V(7)=ExZpE*r1x2Z*VRR1(44,1)
      V(8)=-3.D0*V(7)
      V(9)=r1x2Z*VRR0(45,1)
      V(10)=ExZpE*r1x2Z*VRR1(45,1)
      V(11)=r1x2Z*VRR0(47,1)
      V(12)=ExZpE*r1x2Z*VRR1(47,1)
      V(13)=r1x2Z*VRR0(48,1)
      V(14)=3.D0*V(13)
      V(15)=ExZpE*r1x2Z*VRR1(48,1)
      V(16)=-3.D0*V(15)
      V(17)=r1x2Z*VRR0(49,1)
      V(18)=ExZpE*r1x2Z*VRR1(49,1)
      V(19)=3.D0*V(17)
      V(20)=-3.D0*V(18)
      V(21)=r1x2Z*VRR0(50,1)
      V(22)=ExZpE*r1x2Z*VRR1(50,1)
      V(23)=r1x2Z*VRR0(51,1)
      V(24)=ExZpE*r1x2Z*VRR1(51,1)
      V(25)=r1x2Z*VRR0(52,1)
      V(26)=2.D0*V(25)
      V(27)=ExZpE*r1x2Z*VRR1(52,1)
      V(28)=-2.D0*V(27)
      V(29)=r1x2Z*VRR0(53,1)
      V(30)=ExZpE*r1x2Z*VRR1(53,1)
      VRR0(85,1)=6.D0*r1x2Z*VRR0(36,1)+PAx*VRR0(57,1)-6.D0*ExZpE*r1x2Z*VRR1(36,1)+WPx*VRR1(57,1)
      VRR0(86,1)=5.D0*r1x2Z*VRR0(37,1)+PAx*VRR0(58,1)-5.D0*ExZpE*r1x2Z*VRR1(37,1)+WPx*VRR1(58,1)
      VRR0(87,1)=4.D0*V(1)-4.D0*V(2)+PAx*VRR0(59,1)+WPx*VRR1(59,1)
      VRR0(88,1)=3.D0*V(3)-3.D0*V(4)+PAx*VRR0(60,1)+WPx*VRR1(60,1)
      VRR0(89,1)=3.D0*V(1)-3.D0*V(2)+PAy*VRR0(60,1)+WPy*VRR1(60,1)
      VRR0(90,1)=4.D0*V(3)-4.D0*V(4)+PAy*VRR0(61,1)+WPy*VRR1(61,1)
      VRR0(91,1)=5.D0*r1x2Z*VRR0(40,1)+PAy*VRR0(62,1)-5.D0*ExZpE*r1x2Z*VRR1(40,1)+WPy*VRR1(62,1)
      VRR0(92,1)=6.D0*r1x2Z*VRR0(41,1)+PAy*VRR0(63,1)-6.D0*ExZpE*r1x2Z*VRR1(41,1)+WPy*VRR1(63,1)
      VRR0(93,1)=5.D0*r1x2Z*VRR0(42,1)+PAx*VRR0(64,1)-5.D0*ExZpE*r1x2Z*VRR1(42,1)+WPx*VRR1(64,1)
      VRR0(94,1)=4.D0*r1x2Z*VRR0(43,1)+PAx*VRR0(65,1)-4.D0*ExZpE*r1x2Z*VRR1(43,1)+WPx*VRR1(65,1)
      VRR0(95,1)=V(6)+V(8)+PAx*VRR0(66,1)+WPx*VRR1(66,1)
      VRR0(96,1)=2.D0*V(9)-2.D0*V(10)+PAx*VRR0(67,1)+WPx*VRR1(67,1)
      VRR0(97,1)=V(6)+V(8)+PAy*VRR0(67,1)+WPy*VRR1(67,1)
      VRR0(98,1)=4.D0*V(9)-4.D0*V(10)+PAy*VRR0(68,1)+WPy*VRR1(68,1)
      VRR0(99,1)=5.D0*r1x2Z*VRR0(46,1)+PAy*VRR0(69,1)-5.D0*ExZpE*r1x2Z*VRR1(46,1)+WPy*VRR1(69,1)
      VRR0(100,1)=4.D0*V(11)-4.D0*V(12)+PAx*VRR0(70,1)+WPx*VRR1(70,1)
      VRR0(101,1)=V(14)+V(16)+PAx*VRR0(71,1)+WPx*VRR1(71,1)
      VRR0(102,1)=2.D0*V(17)-2.D0*V(18)+PAx*VRR0(72,1)+WPx*VRR1(72,1)
      VRR0(103,1)=2.D0*V(13)-2.D0*V(15)+PAy*VRR0(72,1)+WPy*VRR1(72,1)
      VRR0(104,1)=V(19)+V(20)+PAy*VRR0(73,1)+WPy*VRR1(73,1)
      VRR0(105,1)=4.D0*V(21)-4.D0*V(22)+PAy*VRR0(74,1)+WPy*VRR1(74,1)
      VRR0(106,1)=3.D0*V(23)-3.D0*V(24)+PAx*VRR0(75,1)+WPx*VRR1(75,1)
      VRR0(107,1)=V(26)+V(28)+PAx*VRR0(76,1)+WPx*VRR1(76,1)
      VRR0(108,1)=2.D0*V(5)-2.D0*V(7)+PAz*VRR0(72,1)+WPz*VRR1(72,1)
      VRR0(109,1)=V(26)+V(28)+PAy*VRR0(77,1)+WPy*VRR1(77,1)
      VRR0(110,1)=3.D0*V(29)-3.D0*V(30)+PAy*VRR0(78,1)+WPy*VRR1(78,1)
      VRR0(111,1)=3.D0*V(11)-3.D0*V(12)+PAz*VRR0(75,1)+WPz*VRR1(75,1)
      VRR0(112,1)=V(14)+V(16)+PAz*VRR0(76,1)+WPz*VRR1(76,1)
      VRR0(113,1)=V(19)+V(20)+PAz*VRR0(77,1)+WPz*VRR1(77,1)
      VRR0(114,1)=3.D0*V(21)-3.D0*V(22)+PAz*VRR0(78,1)+WPz*VRR1(78,1)
      VRR0(115,1)=4.D0*V(23)-4.D0*V(24)+PAz*VRR0(79,1)+WPz*VRR1(79,1)
      VRR0(116,1)=4.D0*V(25)-4.D0*V(27)+PAz*VRR0(80,1)+WPz*VRR1(80,1)
      VRR0(117,1)=4.D0*V(29)-4.D0*V(30)+PAz*VRR0(81,1)+WPz*VRR1(81,1)
      VRR0(118,1)=5.D0*r1x2Z*VRR0(54,1)+PAz*VRR0(82,1)-5.D0*ExZpE*r1x2Z*VRR1(54,1)+WPz*VRR1(82,1)
      VRR0(119,1)=5.D0*r1x2Z*VRR0(55,1)+PAz*VRR0(83,1)-5.D0*ExZpE*r1x2Z*VRR1(55,1)+WPz*VRR1(83,1)
      VRR0(120,1)=6.D0*r1x2Z*VRR0(56,1)+PAz*VRR0(84,1)-6.D0*ExZpE*r1x2Z*VRR1(56,1)+WPz*VRR1(84,1)
END SUBROUTINE VRRj0s0
