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
   SUBROUTINE VRRh0p0(LB,LK,VRR0,VRR1) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2Z*VRR0(12,2)
      V(2)=ExZpE*r1x2Z*VRR1(12,2)
      V(3)=r1x2Z*VRR0(12,3)
      V(4)=ExZpE*r1x2Z*VRR1(12,3)
      V(5)=r1x2Z*VRR0(12,4)
      V(6)=ExZpE*r1x2Z*VRR1(12,4)
      V(7)=r1x2Z*VRR0(13,2)
      V(8)=ExZpE*r1x2Z*VRR1(13,2)
      V(9)=HfxZpE*VRR1(23,1)
      V(10)=r1x2Z*VRR0(13,3)
      V(11)=ExZpE*r1x2Z*VRR1(13,3)
      V(12)=r1x2Z*VRR0(13,4)
      V(13)=ExZpE*r1x2Z*VRR1(13,4)
      V(14)=r1x2Z*VRR0(15,2)
      V(15)=ExZpE*r1x2Z*VRR1(15,2)
      V(16)=r1x2Z*VRR0(15,3)
      V(17)=ExZpE*r1x2Z*VRR1(15,3)
      V(18)=r1x2Z*VRR0(15,4)
      V(19)=ExZpE*r1x2Z*VRR1(15,4)
      V(20)=r1x2Z*VRR0(16,2)
      V(21)=2.D0*V(20)
      V(22)=ExZpE*r1x2Z*VRR1(16,2)
      V(23)=-2.D0*V(22)
      V(24)=r1x2Z*VRR0(16,3)
      V(25)=2.D0*V(24)
      V(26)=ExZpE*r1x2Z*VRR1(16,3)
      V(27)=-2.D0*V(26)
      V(28)=r1x2Z*VRR0(16,4)
      V(29)=2.D0*V(28)
      V(30)=ExZpE*r1x2Z*VRR1(16,4)
      V(31)=-2.D0*V(30)
      V(32)=r1x2Z*VRR0(17,2)
      V(33)=ExZpE*r1x2Z*VRR1(17,2)
      V(34)=HfxZpE*VRR1(28,1)
      V(35)=r1x2Z*VRR0(17,3)
      V(36)=ExZpE*r1x2Z*VRR1(17,3)
      V(37)=r1x2Z*VRR0(17,4)
      V(38)=ExZpE*r1x2Z*VRR1(17,4)
      V(39)=r1x2Z*VRR0(18,2)
      V(40)=ExZpE*r1x2Z*VRR1(18,2)
      V(41)=HfxZpE*VRR1(30,1)
      V(42)=r1x2Z*VRR0(18,3)
      V(43)=ExZpE*r1x2Z*VRR1(18,3)
      V(44)=r1x2Z*VRR0(18,4)
      V(45)=ExZpE*r1x2Z*VRR1(18,4)
      V(46)=r1x2Z*VRR0(19,2)
      V(47)=ExZpE*r1x2Z*VRR1(19,2)
      V(48)=HfxZpE*VRR1(31,1)
      V(49)=r1x2Z*VRR0(19,3)
      V(50)=ExZpE*r1x2Z*VRR1(19,3)
      V(51)=r1x2Z*VRR0(19,4)
      V(52)=ExZpE*r1x2Z*VRR1(19,4)
      V(53)=HfxZpE*VRR1(32,1)
      VRR0(36,2)=4.D0*r1x2Z*VRR0(11,2)+PAx*VRR0(21,2)-4.D0*ExZpE*r1x2Z*VRR1(11,2)+HfxZpE*VRR1(21,1)+WPx*VRR1(21,2)
      VRR0(36,3)=4.D0*r1x2Z*VRR0(11,3)+PAx*VRR0(21,3)-4.D0*ExZpE*r1x2Z*VRR1(11,3)+WPx*VRR1(21,3)
      VRR0(36,4)=4.D0*r1x2Z*VRR0(11,4)+PAx*VRR0(21,4)-4.D0*ExZpE*r1x2Z*VRR1(11,4)+WPx*VRR1(21,4)
      VRR0(37,2)=3.D0*V(1)-3.D0*V(2)+PAx*VRR0(22,2)+HfxZpE*VRR1(22,1)+WPx*VRR1(22,2)
      VRR0(37,3)=3.D0*V(3)-3.D0*V(4)+PAx*VRR0(22,3)+WPx*VRR1(22,3)
      VRR0(37,4)=3.D0*V(5)-3.D0*V(6)+PAx*VRR0(22,4)+WPx*VRR1(22,4)
      VRR0(38,2)=2.D0*V(7)-2.D0*V(8)+V(9)+PAx*VRR0(23,2)+WPx*VRR1(23,2)
      VRR0(38,3)=2.D0*V(10)-2.D0*V(11)+PAx*VRR0(23,3)+WPx*VRR1(23,3)
      VRR0(38,4)=2.D0*V(12)-2.D0*V(13)+PAx*VRR0(23,4)+WPx*VRR1(23,4)
      VRR0(39,2)=2.D0*V(1)-2.D0*V(2)+PAy*VRR0(23,2)+WPy*VRR1(23,2)
      VRR0(39,3)=2.D0*V(3)-2.D0*V(4)+V(9)+PAy*VRR0(23,3)+WPy*VRR1(23,3)
      VRR0(39,4)=2.D0*V(5)-2.D0*V(6)+PAy*VRR0(23,4)+WPy*VRR1(23,4)
      VRR0(40,2)=3.D0*V(7)-3.D0*V(8)+PAy*VRR0(24,2)+WPy*VRR1(24,2)
      VRR0(40,3)=3.D0*V(10)-3.D0*V(11)+PAy*VRR0(24,3)+HfxZpE*VRR1(24,1)+WPy*VRR1(24,3)
      VRR0(40,4)=3.D0*V(12)-3.D0*V(13)+PAy*VRR0(24,4)+WPy*VRR1(24,4)
      VRR0(41,2)=4.D0*r1x2Z*VRR0(14,2)+PAy*VRR0(25,2)-4.D0*ExZpE*r1x2Z*VRR1(14,2)+WPy*VRR1(25,2)
      VRR0(41,3)=4.D0*r1x2Z*VRR0(14,3)+PAy*VRR0(25,3)-4.D0*ExZpE*r1x2Z*VRR1(14,3)+HfxZpE*VRR1(25,1)+WPy*VRR1(25,3)
      VRR0(41,4)=4.D0*r1x2Z*VRR0(14,4)+PAy*VRR0(25,4)-4.D0*ExZpE*r1x2Z*VRR1(14,4)+WPy*VRR1(25,4)
      VRR0(42,2)=3.D0*V(14)-3.D0*V(15)+PAx*VRR0(26,2)+HfxZpE*VRR1(26,1)+WPx*VRR1(26,2)
      VRR0(42,3)=3.D0*V(16)-3.D0*V(17)+PAx*VRR0(26,3)+WPx*VRR1(26,3)
      VRR0(42,4)=3.D0*V(18)-3.D0*V(19)+PAx*VRR0(26,4)+WPx*VRR1(26,4)
      VRR0(43,2)=V(21)+V(23)+PAx*VRR0(27,2)+HfxZpE*VRR1(27,1)+WPx*VRR1(27,2)
      VRR0(43,3)=V(25)+V(27)+PAx*VRR0(27,3)+WPx*VRR1(27,3)
      VRR0(43,4)=V(29)+V(31)+PAx*VRR0(27,4)+WPx*VRR1(27,4)
      VRR0(44,2)=V(32)-V(33)+V(34)+PAx*VRR0(28,2)+WPx*VRR1(28,2)
      VRR0(44,3)=V(35)-V(36)+PAx*VRR0(28,3)+WPx*VRR1(28,3)
      VRR0(44,4)=V(37)-V(38)+PAx*VRR0(28,4)+WPx*VRR1(28,4)
      VRR0(45,2)=V(21)+V(23)+PAy*VRR0(28,2)+WPy*VRR1(28,2)
      VRR0(45,3)=V(25)+V(27)+V(34)+PAy*VRR0(28,3)+WPy*VRR1(28,3)
      VRR0(45,4)=V(29)+V(31)+PAy*VRR0(28,4)+WPy*VRR1(28,4)
      VRR0(46,2)=3.D0*V(32)-3.D0*V(33)+PAy*VRR0(29,2)+WPy*VRR1(29,2)
      VRR0(46,3)=3.D0*V(35)-3.D0*V(36)+PAy*VRR0(29,3)+HfxZpE*VRR1(29,1)+WPy*VRR1(29,3)
      VRR0(46,4)=3.D0*V(37)-3.D0*V(38)+PAy*VRR0(29,4)+WPy*VRR1(29,4)
      VRR0(47,2)=2.D0*V(39)-2.D0*V(40)+V(41)+PAx*VRR0(30,2)+WPx*VRR1(30,2)
      VRR0(47,3)=2.D0*V(42)-2.D0*V(43)+PAx*VRR0(30,3)+WPx*VRR1(30,3)
      VRR0(47,4)=2.D0*V(44)-2.D0*V(45)+PAx*VRR0(30,4)+WPx*VRR1(30,4)
      VRR0(48,2)=V(46)-V(47)+V(48)+PAx*VRR0(31,2)+WPx*VRR1(31,2)
      VRR0(48,3)=V(49)-V(50)+PAx*VRR0(31,3)+WPx*VRR1(31,3)
      VRR0(48,4)=V(51)-V(52)+PAx*VRR0(31,4)+WPx*VRR1(31,4)
      VRR0(49,2)=V(39)-V(40)+PAy*VRR0(31,2)+WPy*VRR1(31,2)
      VRR0(49,3)=V(42)-V(43)+V(48)+PAy*VRR0(31,3)+WPy*VRR1(31,3)
      VRR0(49,4)=V(44)-V(45)+PAy*VRR0(31,4)+WPy*VRR1(31,4)
      VRR0(50,2)=2.D0*V(46)-2.D0*V(47)+PAy*VRR0(32,2)+WPy*VRR1(32,2)
      VRR0(50,3)=2.D0*V(49)-2.D0*V(50)+V(53)+PAy*VRR0(32,3)+WPy*VRR1(32,3)
      VRR0(50,4)=2.D0*V(51)-2.D0*V(52)+PAy*VRR0(32,4)+WPy*VRR1(32,4)
      VRR0(51,2)=2.D0*V(14)-2.D0*V(15)+PAz*VRR0(30,2)+WPz*VRR1(30,2)
      VRR0(51,3)=2.D0*V(16)-2.D0*V(17)+PAz*VRR0(30,3)+WPz*VRR1(30,3)
      VRR0(51,4)=2.D0*V(18)-2.D0*V(19)+V(41)+PAz*VRR0(30,4)+WPz*VRR1(30,4)
      VRR0(52,2)=V(21)+V(23)+PAz*VRR0(31,2)+WPz*VRR1(31,2)
      VRR0(52,3)=V(25)+V(27)+PAz*VRR0(31,3)+WPz*VRR1(31,3)
      VRR0(52,4)=V(29)+V(31)+V(48)+PAz*VRR0(31,4)+WPz*VRR1(31,4)
      VRR0(53,2)=2.D0*V(32)-2.D0*V(33)+PAz*VRR0(32,2)+WPz*VRR1(32,2)
      VRR0(53,3)=2.D0*V(35)-2.D0*V(36)+PAz*VRR0(32,3)+WPz*VRR1(32,3)
      VRR0(53,4)=2.D0*V(37)-2.D0*V(38)+V(53)+PAz*VRR0(32,4)+WPz*VRR1(32,4)
      VRR0(54,2)=3.D0*V(39)-3.D0*V(40)+PAz*VRR0(33,2)+WPz*VRR1(33,2)
      VRR0(54,3)=3.D0*V(42)-3.D0*V(43)+PAz*VRR0(33,3)+WPz*VRR1(33,3)
      VRR0(54,4)=3.D0*V(44)-3.D0*V(45)+PAz*VRR0(33,4)+HfxZpE*VRR1(33,1)+WPz*VRR1(33,4)
      VRR0(55,2)=3.D0*V(46)-3.D0*V(47)+PAz*VRR0(34,2)+WPz*VRR1(34,2)
      VRR0(55,3)=3.D0*V(49)-3.D0*V(50)+PAz*VRR0(34,3)+WPz*VRR1(34,3)
      VRR0(55,4)=3.D0*V(51)-3.D0*V(52)+PAz*VRR0(34,4)+HfxZpE*VRR1(34,1)+WPz*VRR1(34,4)
      VRR0(56,2)=4.D0*r1x2Z*VRR0(20,2)+PAz*VRR0(35,2)-4.D0*ExZpE*r1x2Z*VRR1(20,2)+WPz*VRR1(35,2)
      VRR0(56,3)=4.D0*r1x2Z*VRR0(20,3)+PAz*VRR0(35,3)-4.D0*ExZpE*r1x2Z*VRR1(20,3)+WPz*VRR1(35,3)
      VRR0(56,4)=4.D0*r1x2Z*VRR0(20,4)+PAz*VRR0(35,4)-4.D0*ExZpE*r1x2Z*VRR1(20,4)+HfxZpE*VRR1(35,1)+WPz*VRR1(35,4)
END SUBROUTINE VRRh0p0
SUBROUTINE MVRRh0p0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(36,2)=PAx*VS0(21,2)+WPx*VS1(21,2)+r1x2Z*VR1(21,2)&
   +4D0*r1x2Z*(VS0(11,2)-ExZpE*VS1(11,2))&
   +HfxZpE*VS1(21,1)
VS0(36,3)=PAx*VS0(21,3)+WPx*VS1(21,3)+r1x2Z*VR1(21,3)&
   +4D0*r1x2Z*(VS0(11,3)-ExZpE*VS1(11,3))
VS0(36,4)=PAx*VS0(21,4)+WPx*VS1(21,4)+r1x2Z*VR1(21,4)&
   +4D0*r1x2Z*(VS0(11,4)-ExZpE*VS1(11,4))
VS0(37,2)=PAx*VS0(22,2)+WPx*VS1(22,2)+r1x2Z*VR1(22,2)&
   +3D0*r1x2Z*(VS0(12,2)-ExZpE*VS1(12,2))&
   +HfxZpE*VS1(22,1)
VS0(37,3)=PAx*VS0(22,3)+WPx*VS1(22,3)+r1x2Z*VR1(22,3)&
   +3D0*r1x2Z*(VS0(12,3)-ExZpE*VS1(12,3))
VS0(37,4)=PAx*VS0(22,4)+WPx*VS1(22,4)+r1x2Z*VR1(22,4)&
   +3D0*r1x2Z*(VS0(12,4)-ExZpE*VS1(12,4))
VS0(38,2)=PAx*VS0(23,2)+WPx*VS1(23,2)+r1x2Z*VR1(23,2)&
   +2D0*r1x2Z*(VS0(13,2)-ExZpE*VS1(13,2))&
   +HfxZpE*VS1(23,1)
VS0(38,3)=PAx*VS0(23,3)+WPx*VS1(23,3)+r1x2Z*VR1(23,3)&
   +2D0*r1x2Z*(VS0(13,3)-ExZpE*VS1(13,3))
VS0(38,4)=PAx*VS0(23,4)+WPx*VS1(23,4)+r1x2Z*VR1(23,4)&
   +2D0*r1x2Z*(VS0(13,4)-ExZpE*VS1(13,4))
VS0(39,2)=PAy*VS0(23,2)+WPy*VS1(23,2)&
   +2D0*r1x2Z*(VS0(12,2)-ExZpE*VS1(12,2))
VS0(39,3)=PAy*VS0(23,3)+WPy*VS1(23,3)&
   +2D0*r1x2Z*(VS0(12,3)-ExZpE*VS1(12,3))&
   +HfxZpE*VS1(23,1)
VS0(39,4)=PAy*VS0(23,4)+WPy*VS1(23,4)&
   +2D0*r1x2Z*(VS0(12,4)-ExZpE*VS1(12,4))
VS0(40,2)=PAy*VS0(24,2)+WPy*VS1(24,2)&
   +3D0*r1x2Z*(VS0(13,2)-ExZpE*VS1(13,2))
VS0(40,3)=PAy*VS0(24,3)+WPy*VS1(24,3)&
   +3D0*r1x2Z*(VS0(13,3)-ExZpE*VS1(13,3))&
   +HfxZpE*VS1(24,1)
VS0(40,4)=PAy*VS0(24,4)+WPy*VS1(24,4)&
   +3D0*r1x2Z*(VS0(13,4)-ExZpE*VS1(13,4))
VS0(41,2)=PAy*VS0(25,2)+WPy*VS1(25,2)&
   +4D0*r1x2Z*(VS0(14,2)-ExZpE*VS1(14,2))
VS0(41,3)=PAy*VS0(25,3)+WPy*VS1(25,3)&
   +4D0*r1x2Z*(VS0(14,3)-ExZpE*VS1(14,3))&
   +HfxZpE*VS1(25,1)
VS0(41,4)=PAy*VS0(25,4)+WPy*VS1(25,4)&
   +4D0*r1x2Z*(VS0(14,4)-ExZpE*VS1(14,4))
VS0(42,2)=PAx*VS0(26,2)+WPx*VS1(26,2)+r1x2Z*VR1(26,2)&
   +3D0*r1x2Z*(VS0(15,2)-ExZpE*VS1(15,2))&
   +HfxZpE*VS1(26,1)
VS0(42,3)=PAx*VS0(26,3)+WPx*VS1(26,3)+r1x2Z*VR1(26,3)&
   +3D0*r1x2Z*(VS0(15,3)-ExZpE*VS1(15,3))
VS0(42,4)=PAx*VS0(26,4)+WPx*VS1(26,4)+r1x2Z*VR1(26,4)&
   +3D0*r1x2Z*(VS0(15,4)-ExZpE*VS1(15,4))
VS0(43,2)=PAx*VS0(27,2)+WPx*VS1(27,2)+r1x2Z*VR1(27,2)&
   +2D0*r1x2Z*(VS0(16,2)-ExZpE*VS1(16,2))&
   +HfxZpE*VS1(27,1)
VS0(43,3)=PAx*VS0(27,3)+WPx*VS1(27,3)+r1x2Z*VR1(27,3)&
   +2D0*r1x2Z*(VS0(16,3)-ExZpE*VS1(16,3))
VS0(43,4)=PAx*VS0(27,4)+WPx*VS1(27,4)+r1x2Z*VR1(27,4)&
   +2D0*r1x2Z*(VS0(16,4)-ExZpE*VS1(16,4))
VS0(44,2)=PAx*VS0(28,2)+WPx*VS1(28,2)+r1x2Z*VR1(28,2)&
   +r1x2Z*(VS0(17,2)-ExZpE*VS1(17,2))&
   +HfxZpE*VS1(28,1)
VS0(44,3)=PAx*VS0(28,3)+WPx*VS1(28,3)+r1x2Z*VR1(28,3)&
   +r1x2Z*(VS0(17,3)-ExZpE*VS1(17,3))
VS0(44,4)=PAx*VS0(28,4)+WPx*VS1(28,4)+r1x2Z*VR1(28,4)&
   +r1x2Z*(VS0(17,4)-ExZpE*VS1(17,4))
VS0(45,2)=PAy*VS0(28,2)+WPy*VS1(28,2)&
   +2D0*r1x2Z*(VS0(16,2)-ExZpE*VS1(16,2))
VS0(45,3)=PAy*VS0(28,3)+WPy*VS1(28,3)&
   +2D0*r1x2Z*(VS0(16,3)-ExZpE*VS1(16,3))&
   +HfxZpE*VS1(28,1)
VS0(45,4)=PAy*VS0(28,4)+WPy*VS1(28,4)&
   +2D0*r1x2Z*(VS0(16,4)-ExZpE*VS1(16,4))
VS0(46,2)=PAy*VS0(29,2)+WPy*VS1(29,2)&
   +3D0*r1x2Z*(VS0(17,2)-ExZpE*VS1(17,2))
VS0(46,3)=PAy*VS0(29,3)+WPy*VS1(29,3)&
   +3D0*r1x2Z*(VS0(17,3)-ExZpE*VS1(17,3))&
   +HfxZpE*VS1(29,1)
VS0(46,4)=PAy*VS0(29,4)+WPy*VS1(29,4)&
   +3D0*r1x2Z*(VS0(17,4)-ExZpE*VS1(17,4))
VS0(47,2)=PAx*VS0(30,2)+WPx*VS1(30,2)+r1x2Z*VR1(30,2)&
   +2D0*r1x2Z*(VS0(18,2)-ExZpE*VS1(18,2))&
   +HfxZpE*VS1(30,1)
VS0(47,3)=PAx*VS0(30,3)+WPx*VS1(30,3)+r1x2Z*VR1(30,3)&
   +2D0*r1x2Z*(VS0(18,3)-ExZpE*VS1(18,3))
VS0(47,4)=PAx*VS0(30,4)+WPx*VS1(30,4)+r1x2Z*VR1(30,4)&
   +2D0*r1x2Z*(VS0(18,4)-ExZpE*VS1(18,4))
VS0(48,2)=PAx*VS0(31,2)+WPx*VS1(31,2)+r1x2Z*VR1(31,2)&
   +r1x2Z*(VS0(19,2)-ExZpE*VS1(19,2))&
   +HfxZpE*VS1(31,1)
VS0(48,3)=PAx*VS0(31,3)+WPx*VS1(31,3)+r1x2Z*VR1(31,3)&
   +r1x2Z*(VS0(19,3)-ExZpE*VS1(19,3))
VS0(48,4)=PAx*VS0(31,4)+WPx*VS1(31,4)+r1x2Z*VR1(31,4)&
   +r1x2Z*(VS0(19,4)-ExZpE*VS1(19,4))
VS0(49,2)=PAy*VS0(31,2)+WPy*VS1(31,2)&
   +r1x2Z*(VS0(18,2)-ExZpE*VS1(18,2))
VS0(49,3)=PAy*VS0(31,3)+WPy*VS1(31,3)&
   +r1x2Z*(VS0(18,3)-ExZpE*VS1(18,3))&
   +HfxZpE*VS1(31,1)
VS0(49,4)=PAy*VS0(31,4)+WPy*VS1(31,4)&
   +r1x2Z*(VS0(18,4)-ExZpE*VS1(18,4))
VS0(50,2)=PAy*VS0(32,2)+WPy*VS1(32,2)&
   +2D0*r1x2Z*(VS0(19,2)-ExZpE*VS1(19,2))
VS0(50,3)=PAy*VS0(32,3)+WPy*VS1(32,3)&
   +2D0*r1x2Z*(VS0(19,3)-ExZpE*VS1(19,3))&
   +HfxZpE*VS1(32,1)
VS0(50,4)=PAy*VS0(32,4)+WPy*VS1(32,4)&
   +2D0*r1x2Z*(VS0(19,4)-ExZpE*VS1(19,4))
VS0(51,2)=PAz*VS0(30,2)+WPz*VS1(30,2)&
   +2D0*r1x2Z*(VS0(15,2)-ExZpE*VS1(15,2))
VS0(51,3)=PAz*VS0(30,3)+WPz*VS1(30,3)&
   +2D0*r1x2Z*(VS0(15,3)-ExZpE*VS1(15,3))
VS0(51,4)=PAz*VS0(30,4)+WPz*VS1(30,4)&
   +2D0*r1x2Z*(VS0(15,4)-ExZpE*VS1(15,4))&
   +HfxZpE*VS1(30,1)
VS0(52,2)=PAz*VS0(31,2)+WPz*VS1(31,2)&
   +2D0*r1x2Z*(VS0(16,2)-ExZpE*VS1(16,2))
VS0(52,3)=PAz*VS0(31,3)+WPz*VS1(31,3)&
   +2D0*r1x2Z*(VS0(16,3)-ExZpE*VS1(16,3))
VS0(52,4)=PAz*VS0(31,4)+WPz*VS1(31,4)&
   +2D0*r1x2Z*(VS0(16,4)-ExZpE*VS1(16,4))&
   +HfxZpE*VS1(31,1)
VS0(53,2)=PAz*VS0(32,2)+WPz*VS1(32,2)&
   +2D0*r1x2Z*(VS0(17,2)-ExZpE*VS1(17,2))
VS0(53,3)=PAz*VS0(32,3)+WPz*VS1(32,3)&
   +2D0*r1x2Z*(VS0(17,3)-ExZpE*VS1(17,3))
VS0(53,4)=PAz*VS0(32,4)+WPz*VS1(32,4)&
   +2D0*r1x2Z*(VS0(17,4)-ExZpE*VS1(17,4))&
   +HfxZpE*VS1(32,1)
VS0(54,2)=PAz*VS0(33,2)+WPz*VS1(33,2)&
   +3D0*r1x2Z*(VS0(18,2)-ExZpE*VS1(18,2))
VS0(54,3)=PAz*VS0(33,3)+WPz*VS1(33,3)&
   +3D0*r1x2Z*(VS0(18,3)-ExZpE*VS1(18,3))
VS0(54,4)=PAz*VS0(33,4)+WPz*VS1(33,4)&
   +3D0*r1x2Z*(VS0(18,4)-ExZpE*VS1(18,4))&
   +HfxZpE*VS1(33,1)
VS0(55,2)=PAz*VS0(34,2)+WPz*VS1(34,2)&
   +3D0*r1x2Z*(VS0(19,2)-ExZpE*VS1(19,2))
VS0(55,3)=PAz*VS0(34,3)+WPz*VS1(34,3)&
   +3D0*r1x2Z*(VS0(19,3)-ExZpE*VS1(19,3))
VS0(55,4)=PAz*VS0(34,4)+WPz*VS1(34,4)&
   +3D0*r1x2Z*(VS0(19,4)-ExZpE*VS1(19,4))&
   +HfxZpE*VS1(34,1)
VS0(56,2)=PAz*VS0(35,2)+WPz*VS1(35,2)&
   +4D0*r1x2Z*(VS0(20,2)-ExZpE*VS1(20,2))
VS0(56,3)=PAz*VS0(35,3)+WPz*VS1(35,3)&
   +4D0*r1x2Z*(VS0(20,3)-ExZpE*VS1(20,3))
VS0(56,4)=PAz*VS0(35,4)+WPz*VS1(35,4)&
   +4D0*r1x2Z*(VS0(20,4)-ExZpE*VS1(20,4))&
   +HfxZpE*VS1(35,1)
CASE(2)
VS0(36,2)=PAx*VS0(21,2)+WPx*VS1(21,2)&
   +4D0*r1x2Z*(VS0(11,2)-ExZpE*VS1(11,2))&
   +HfxZpE*VS1(21,1)
VS0(36,3)=PAx*VS0(21,3)+WPx*VS1(21,3)&
   +4D0*r1x2Z*(VS0(11,3)-ExZpE*VS1(11,3))
VS0(36,4)=PAx*VS0(21,4)+WPx*VS1(21,4)&
   +4D0*r1x2Z*(VS0(11,4)-ExZpE*VS1(11,4))
VS0(37,2)=PAx*VS0(22,2)+WPx*VS1(22,2)&
   +3D0*r1x2Z*(VS0(12,2)-ExZpE*VS1(12,2))&
   +HfxZpE*VS1(22,1)
VS0(37,3)=PAx*VS0(22,3)+WPx*VS1(22,3)&
   +3D0*r1x2Z*(VS0(12,3)-ExZpE*VS1(12,3))
VS0(37,4)=PAx*VS0(22,4)+WPx*VS1(22,4)&
   +3D0*r1x2Z*(VS0(12,4)-ExZpE*VS1(12,4))
VS0(38,2)=PAx*VS0(23,2)+WPx*VS1(23,2)&
   +2D0*r1x2Z*(VS0(13,2)-ExZpE*VS1(13,2))&
   +HfxZpE*VS1(23,1)
VS0(38,3)=PAx*VS0(23,3)+WPx*VS1(23,3)&
   +2D0*r1x2Z*(VS0(13,3)-ExZpE*VS1(13,3))
VS0(38,4)=PAx*VS0(23,4)+WPx*VS1(23,4)&
   +2D0*r1x2Z*(VS0(13,4)-ExZpE*VS1(13,4))
VS0(39,2)=PAy*VS0(23,2)+WPy*VS1(23,2)+r1x2Z*VR1(23,2)&
   +2D0*r1x2Z*(VS0(12,2)-ExZpE*VS1(12,2))
VS0(39,3)=PAy*VS0(23,3)+WPy*VS1(23,3)+r1x2Z*VR1(23,3)&
   +2D0*r1x2Z*(VS0(12,3)-ExZpE*VS1(12,3))&
   +HfxZpE*VS1(23,1)
VS0(39,4)=PAy*VS0(23,4)+WPy*VS1(23,4)+r1x2Z*VR1(23,4)&
   +2D0*r1x2Z*(VS0(12,4)-ExZpE*VS1(12,4))
VS0(40,2)=PAy*VS0(24,2)+WPy*VS1(24,2)+r1x2Z*VR1(24,2)&
   +3D0*r1x2Z*(VS0(13,2)-ExZpE*VS1(13,2))
VS0(40,3)=PAy*VS0(24,3)+WPy*VS1(24,3)+r1x2Z*VR1(24,3)&
   +3D0*r1x2Z*(VS0(13,3)-ExZpE*VS1(13,3))&
   +HfxZpE*VS1(24,1)
VS0(40,4)=PAy*VS0(24,4)+WPy*VS1(24,4)+r1x2Z*VR1(24,4)&
   +3D0*r1x2Z*(VS0(13,4)-ExZpE*VS1(13,4))
VS0(41,2)=PAy*VS0(25,2)+WPy*VS1(25,2)+r1x2Z*VR1(25,2)&
   +4D0*r1x2Z*(VS0(14,2)-ExZpE*VS1(14,2))
VS0(41,3)=PAy*VS0(25,3)+WPy*VS1(25,3)+r1x2Z*VR1(25,3)&
   +4D0*r1x2Z*(VS0(14,3)-ExZpE*VS1(14,3))&
   +HfxZpE*VS1(25,1)
VS0(41,4)=PAy*VS0(25,4)+WPy*VS1(25,4)+r1x2Z*VR1(25,4)&
   +4D0*r1x2Z*(VS0(14,4)-ExZpE*VS1(14,4))
VS0(42,2)=PAx*VS0(26,2)+WPx*VS1(26,2)&
   +3D0*r1x2Z*(VS0(15,2)-ExZpE*VS1(15,2))&
   +HfxZpE*VS1(26,1)
VS0(42,3)=PAx*VS0(26,3)+WPx*VS1(26,3)&
   +3D0*r1x2Z*(VS0(15,3)-ExZpE*VS1(15,3))
VS0(42,4)=PAx*VS0(26,4)+WPx*VS1(26,4)&
   +3D0*r1x2Z*(VS0(15,4)-ExZpE*VS1(15,4))
VS0(43,2)=PAx*VS0(27,2)+WPx*VS1(27,2)&
   +2D0*r1x2Z*(VS0(16,2)-ExZpE*VS1(16,2))&
   +HfxZpE*VS1(27,1)
VS0(43,3)=PAx*VS0(27,3)+WPx*VS1(27,3)&
   +2D0*r1x2Z*(VS0(16,3)-ExZpE*VS1(16,3))
VS0(43,4)=PAx*VS0(27,4)+WPx*VS1(27,4)&
   +2D0*r1x2Z*(VS0(16,4)-ExZpE*VS1(16,4))
VS0(44,2)=PAx*VS0(28,2)+WPx*VS1(28,2)&
   +r1x2Z*(VS0(17,2)-ExZpE*VS1(17,2))&
   +HfxZpE*VS1(28,1)
VS0(44,3)=PAx*VS0(28,3)+WPx*VS1(28,3)&
   +r1x2Z*(VS0(17,3)-ExZpE*VS1(17,3))
VS0(44,4)=PAx*VS0(28,4)+WPx*VS1(28,4)&
   +r1x2Z*(VS0(17,4)-ExZpE*VS1(17,4))
VS0(45,2)=PAy*VS0(28,2)+WPy*VS1(28,2)+r1x2Z*VR1(28,2)&
   +2D0*r1x2Z*(VS0(16,2)-ExZpE*VS1(16,2))
VS0(45,3)=PAy*VS0(28,3)+WPy*VS1(28,3)+r1x2Z*VR1(28,3)&
   +2D0*r1x2Z*(VS0(16,3)-ExZpE*VS1(16,3))&
   +HfxZpE*VS1(28,1)
VS0(45,4)=PAy*VS0(28,4)+WPy*VS1(28,4)+r1x2Z*VR1(28,4)&
   +2D0*r1x2Z*(VS0(16,4)-ExZpE*VS1(16,4))
VS0(46,2)=PAy*VS0(29,2)+WPy*VS1(29,2)+r1x2Z*VR1(29,2)&
   +3D0*r1x2Z*(VS0(17,2)-ExZpE*VS1(17,2))
VS0(46,3)=PAy*VS0(29,3)+WPy*VS1(29,3)+r1x2Z*VR1(29,3)&
   +3D0*r1x2Z*(VS0(17,3)-ExZpE*VS1(17,3))&
   +HfxZpE*VS1(29,1)
VS0(46,4)=PAy*VS0(29,4)+WPy*VS1(29,4)+r1x2Z*VR1(29,4)&
   +3D0*r1x2Z*(VS0(17,4)-ExZpE*VS1(17,4))
VS0(47,2)=PAx*VS0(30,2)+WPx*VS1(30,2)&
   +2D0*r1x2Z*(VS0(18,2)-ExZpE*VS1(18,2))&
   +HfxZpE*VS1(30,1)
VS0(47,3)=PAx*VS0(30,3)+WPx*VS1(30,3)&
   +2D0*r1x2Z*(VS0(18,3)-ExZpE*VS1(18,3))
VS0(47,4)=PAx*VS0(30,4)+WPx*VS1(30,4)&
   +2D0*r1x2Z*(VS0(18,4)-ExZpE*VS1(18,4))
VS0(48,2)=PAx*VS0(31,2)+WPx*VS1(31,2)&
   +r1x2Z*(VS0(19,2)-ExZpE*VS1(19,2))&
   +HfxZpE*VS1(31,1)
VS0(48,3)=PAx*VS0(31,3)+WPx*VS1(31,3)&
   +r1x2Z*(VS0(19,3)-ExZpE*VS1(19,3))
VS0(48,4)=PAx*VS0(31,4)+WPx*VS1(31,4)&
   +r1x2Z*(VS0(19,4)-ExZpE*VS1(19,4))
VS0(49,2)=PAy*VS0(31,2)+WPy*VS1(31,2)+r1x2Z*VR1(31,2)&
   +r1x2Z*(VS0(18,2)-ExZpE*VS1(18,2))
VS0(49,3)=PAy*VS0(31,3)+WPy*VS1(31,3)+r1x2Z*VR1(31,3)&
   +r1x2Z*(VS0(18,3)-ExZpE*VS1(18,3))&
   +HfxZpE*VS1(31,1)
VS0(49,4)=PAy*VS0(31,4)+WPy*VS1(31,4)+r1x2Z*VR1(31,4)&
   +r1x2Z*(VS0(18,4)-ExZpE*VS1(18,4))
VS0(50,2)=PAy*VS0(32,2)+WPy*VS1(32,2)+r1x2Z*VR1(32,2)&
   +2D0*r1x2Z*(VS0(19,2)-ExZpE*VS1(19,2))
VS0(50,3)=PAy*VS0(32,3)+WPy*VS1(32,3)+r1x2Z*VR1(32,3)&
   +2D0*r1x2Z*(VS0(19,3)-ExZpE*VS1(19,3))&
   +HfxZpE*VS1(32,1)
VS0(50,4)=PAy*VS0(32,4)+WPy*VS1(32,4)+r1x2Z*VR1(32,4)&
   +2D0*r1x2Z*(VS0(19,4)-ExZpE*VS1(19,4))
VS0(51,2)=PAz*VS0(30,2)+WPz*VS1(30,2)&
   +2D0*r1x2Z*(VS0(15,2)-ExZpE*VS1(15,2))
VS0(51,3)=PAz*VS0(30,3)+WPz*VS1(30,3)&
   +2D0*r1x2Z*(VS0(15,3)-ExZpE*VS1(15,3))
VS0(51,4)=PAz*VS0(30,4)+WPz*VS1(30,4)&
   +2D0*r1x2Z*(VS0(15,4)-ExZpE*VS1(15,4))&
   +HfxZpE*VS1(30,1)
VS0(52,2)=PAz*VS0(31,2)+WPz*VS1(31,2)&
   +2D0*r1x2Z*(VS0(16,2)-ExZpE*VS1(16,2))
VS0(52,3)=PAz*VS0(31,3)+WPz*VS1(31,3)&
   +2D0*r1x2Z*(VS0(16,3)-ExZpE*VS1(16,3))
VS0(52,4)=PAz*VS0(31,4)+WPz*VS1(31,4)&
   +2D0*r1x2Z*(VS0(16,4)-ExZpE*VS1(16,4))&
   +HfxZpE*VS1(31,1)
VS0(53,2)=PAz*VS0(32,2)+WPz*VS1(32,2)&
   +2D0*r1x2Z*(VS0(17,2)-ExZpE*VS1(17,2))
VS0(53,3)=PAz*VS0(32,3)+WPz*VS1(32,3)&
   +2D0*r1x2Z*(VS0(17,3)-ExZpE*VS1(17,3))
VS0(53,4)=PAz*VS0(32,4)+WPz*VS1(32,4)&
   +2D0*r1x2Z*(VS0(17,4)-ExZpE*VS1(17,4))&
   +HfxZpE*VS1(32,1)
VS0(54,2)=PAz*VS0(33,2)+WPz*VS1(33,2)&
   +3D0*r1x2Z*(VS0(18,2)-ExZpE*VS1(18,2))
VS0(54,3)=PAz*VS0(33,3)+WPz*VS1(33,3)&
   +3D0*r1x2Z*(VS0(18,3)-ExZpE*VS1(18,3))
VS0(54,4)=PAz*VS0(33,4)+WPz*VS1(33,4)&
   +3D0*r1x2Z*(VS0(18,4)-ExZpE*VS1(18,4))&
   +HfxZpE*VS1(33,1)
VS0(55,2)=PAz*VS0(34,2)+WPz*VS1(34,2)&
   +3D0*r1x2Z*(VS0(19,2)-ExZpE*VS1(19,2))
VS0(55,3)=PAz*VS0(34,3)+WPz*VS1(34,3)&
   +3D0*r1x2Z*(VS0(19,3)-ExZpE*VS1(19,3))
VS0(55,4)=PAz*VS0(34,4)+WPz*VS1(34,4)&
   +3D0*r1x2Z*(VS0(19,4)-ExZpE*VS1(19,4))&
   +HfxZpE*VS1(34,1)
VS0(56,2)=PAz*VS0(35,2)+WPz*VS1(35,2)&
   +4D0*r1x2Z*(VS0(20,2)-ExZpE*VS1(20,2))
VS0(56,3)=PAz*VS0(35,3)+WPz*VS1(35,3)&
   +4D0*r1x2Z*(VS0(20,3)-ExZpE*VS1(20,3))
VS0(56,4)=PAz*VS0(35,4)+WPz*VS1(35,4)&
   +4D0*r1x2Z*(VS0(20,4)-ExZpE*VS1(20,4))&
   +HfxZpE*VS1(35,1)
CASE(3)
VS0(36,2)=PAx*VS0(21,2)+WPx*VS1(21,2)&
   +4D0*r1x2Z*(VS0(11,2)-ExZpE*VS1(11,2))&
   +HfxZpE*VS1(21,1)
VS0(36,3)=PAx*VS0(21,3)+WPx*VS1(21,3)&
   +4D0*r1x2Z*(VS0(11,3)-ExZpE*VS1(11,3))
VS0(36,4)=PAx*VS0(21,4)+WPx*VS1(21,4)&
   +4D0*r1x2Z*(VS0(11,4)-ExZpE*VS1(11,4))
VS0(37,2)=PAx*VS0(22,2)+WPx*VS1(22,2)&
   +3D0*r1x2Z*(VS0(12,2)-ExZpE*VS1(12,2))&
   +HfxZpE*VS1(22,1)
VS0(37,3)=PAx*VS0(22,3)+WPx*VS1(22,3)&
   +3D0*r1x2Z*(VS0(12,3)-ExZpE*VS1(12,3))
VS0(37,4)=PAx*VS0(22,4)+WPx*VS1(22,4)&
   +3D0*r1x2Z*(VS0(12,4)-ExZpE*VS1(12,4))
VS0(38,2)=PAx*VS0(23,2)+WPx*VS1(23,2)&
   +2D0*r1x2Z*(VS0(13,2)-ExZpE*VS1(13,2))&
   +HfxZpE*VS1(23,1)
VS0(38,3)=PAx*VS0(23,3)+WPx*VS1(23,3)&
   +2D0*r1x2Z*(VS0(13,3)-ExZpE*VS1(13,3))
VS0(38,4)=PAx*VS0(23,4)+WPx*VS1(23,4)&
   +2D0*r1x2Z*(VS0(13,4)-ExZpE*VS1(13,4))
VS0(39,2)=PAy*VS0(23,2)+WPy*VS1(23,2)&
   +2D0*r1x2Z*(VS0(12,2)-ExZpE*VS1(12,2))
VS0(39,3)=PAy*VS0(23,3)+WPy*VS1(23,3)&
   +2D0*r1x2Z*(VS0(12,3)-ExZpE*VS1(12,3))&
   +HfxZpE*VS1(23,1)
VS0(39,4)=PAy*VS0(23,4)+WPy*VS1(23,4)&
   +2D0*r1x2Z*(VS0(12,4)-ExZpE*VS1(12,4))
VS0(40,2)=PAy*VS0(24,2)+WPy*VS1(24,2)&
   +3D0*r1x2Z*(VS0(13,2)-ExZpE*VS1(13,2))
VS0(40,3)=PAy*VS0(24,3)+WPy*VS1(24,3)&
   +3D0*r1x2Z*(VS0(13,3)-ExZpE*VS1(13,3))&
   +HfxZpE*VS1(24,1)
VS0(40,4)=PAy*VS0(24,4)+WPy*VS1(24,4)&
   +3D0*r1x2Z*(VS0(13,4)-ExZpE*VS1(13,4))
VS0(41,2)=PAy*VS0(25,2)+WPy*VS1(25,2)&
   +4D0*r1x2Z*(VS0(14,2)-ExZpE*VS1(14,2))
VS0(41,3)=PAy*VS0(25,3)+WPy*VS1(25,3)&
   +4D0*r1x2Z*(VS0(14,3)-ExZpE*VS1(14,3))&
   +HfxZpE*VS1(25,1)
VS0(41,4)=PAy*VS0(25,4)+WPy*VS1(25,4)&
   +4D0*r1x2Z*(VS0(14,4)-ExZpE*VS1(14,4))
VS0(42,2)=PAx*VS0(26,2)+WPx*VS1(26,2)&
   +3D0*r1x2Z*(VS0(15,2)-ExZpE*VS1(15,2))&
   +HfxZpE*VS1(26,1)
VS0(42,3)=PAx*VS0(26,3)+WPx*VS1(26,3)&
   +3D0*r1x2Z*(VS0(15,3)-ExZpE*VS1(15,3))
VS0(42,4)=PAx*VS0(26,4)+WPx*VS1(26,4)&
   +3D0*r1x2Z*(VS0(15,4)-ExZpE*VS1(15,4))
VS0(43,2)=PAx*VS0(27,2)+WPx*VS1(27,2)&
   +2D0*r1x2Z*(VS0(16,2)-ExZpE*VS1(16,2))&
   +HfxZpE*VS1(27,1)
VS0(43,3)=PAx*VS0(27,3)+WPx*VS1(27,3)&
   +2D0*r1x2Z*(VS0(16,3)-ExZpE*VS1(16,3))
VS0(43,4)=PAx*VS0(27,4)+WPx*VS1(27,4)&
   +2D0*r1x2Z*(VS0(16,4)-ExZpE*VS1(16,4))
VS0(44,2)=PAx*VS0(28,2)+WPx*VS1(28,2)&
   +r1x2Z*(VS0(17,2)-ExZpE*VS1(17,2))&
   +HfxZpE*VS1(28,1)
VS0(44,3)=PAx*VS0(28,3)+WPx*VS1(28,3)&
   +r1x2Z*(VS0(17,3)-ExZpE*VS1(17,3))
VS0(44,4)=PAx*VS0(28,4)+WPx*VS1(28,4)&
   +r1x2Z*(VS0(17,4)-ExZpE*VS1(17,4))
VS0(45,2)=PAy*VS0(28,2)+WPy*VS1(28,2)&
   +2D0*r1x2Z*(VS0(16,2)-ExZpE*VS1(16,2))
VS0(45,3)=PAy*VS0(28,3)+WPy*VS1(28,3)&
   +2D0*r1x2Z*(VS0(16,3)-ExZpE*VS1(16,3))&
   +HfxZpE*VS1(28,1)
VS0(45,4)=PAy*VS0(28,4)+WPy*VS1(28,4)&
   +2D0*r1x2Z*(VS0(16,4)-ExZpE*VS1(16,4))
VS0(46,2)=PAy*VS0(29,2)+WPy*VS1(29,2)&
   +3D0*r1x2Z*(VS0(17,2)-ExZpE*VS1(17,2))
VS0(46,3)=PAy*VS0(29,3)+WPy*VS1(29,3)&
   +3D0*r1x2Z*(VS0(17,3)-ExZpE*VS1(17,3))&
   +HfxZpE*VS1(29,1)
VS0(46,4)=PAy*VS0(29,4)+WPy*VS1(29,4)&
   +3D0*r1x2Z*(VS0(17,4)-ExZpE*VS1(17,4))
VS0(47,2)=PAx*VS0(30,2)+WPx*VS1(30,2)&
   +2D0*r1x2Z*(VS0(18,2)-ExZpE*VS1(18,2))&
   +HfxZpE*VS1(30,1)
VS0(47,3)=PAx*VS0(30,3)+WPx*VS1(30,3)&
   +2D0*r1x2Z*(VS0(18,3)-ExZpE*VS1(18,3))
VS0(47,4)=PAx*VS0(30,4)+WPx*VS1(30,4)&
   +2D0*r1x2Z*(VS0(18,4)-ExZpE*VS1(18,4))
VS0(48,2)=PAx*VS0(31,2)+WPx*VS1(31,2)&
   +r1x2Z*(VS0(19,2)-ExZpE*VS1(19,2))&
   +HfxZpE*VS1(31,1)
VS0(48,3)=PAx*VS0(31,3)+WPx*VS1(31,3)&
   +r1x2Z*(VS0(19,3)-ExZpE*VS1(19,3))
VS0(48,4)=PAx*VS0(31,4)+WPx*VS1(31,4)&
   +r1x2Z*(VS0(19,4)-ExZpE*VS1(19,4))
VS0(49,2)=PAy*VS0(31,2)+WPy*VS1(31,2)&
   +r1x2Z*(VS0(18,2)-ExZpE*VS1(18,2))
VS0(49,3)=PAy*VS0(31,3)+WPy*VS1(31,3)&
   +r1x2Z*(VS0(18,3)-ExZpE*VS1(18,3))&
   +HfxZpE*VS1(31,1)
VS0(49,4)=PAy*VS0(31,4)+WPy*VS1(31,4)&
   +r1x2Z*(VS0(18,4)-ExZpE*VS1(18,4))
VS0(50,2)=PAy*VS0(32,2)+WPy*VS1(32,2)&
   +2D0*r1x2Z*(VS0(19,2)-ExZpE*VS1(19,2))
VS0(50,3)=PAy*VS0(32,3)+WPy*VS1(32,3)&
   +2D0*r1x2Z*(VS0(19,3)-ExZpE*VS1(19,3))&
   +HfxZpE*VS1(32,1)
VS0(50,4)=PAy*VS0(32,4)+WPy*VS1(32,4)&
   +2D0*r1x2Z*(VS0(19,4)-ExZpE*VS1(19,4))
VS0(51,2)=PAz*VS0(30,2)+WPz*VS1(30,2)+r1x2Z*VR1(30,2)&
   +2D0*r1x2Z*(VS0(15,2)-ExZpE*VS1(15,2))
VS0(51,3)=PAz*VS0(30,3)+WPz*VS1(30,3)+r1x2Z*VR1(30,3)&
   +2D0*r1x2Z*(VS0(15,3)-ExZpE*VS1(15,3))
VS0(51,4)=PAz*VS0(30,4)+WPz*VS1(30,4)+r1x2Z*VR1(30,4)&
   +2D0*r1x2Z*(VS0(15,4)-ExZpE*VS1(15,4))&
   +HfxZpE*VS1(30,1)
VS0(52,2)=PAz*VS0(31,2)+WPz*VS1(31,2)+r1x2Z*VR1(31,2)&
   +2D0*r1x2Z*(VS0(16,2)-ExZpE*VS1(16,2))
VS0(52,3)=PAz*VS0(31,3)+WPz*VS1(31,3)+r1x2Z*VR1(31,3)&
   +2D0*r1x2Z*(VS0(16,3)-ExZpE*VS1(16,3))
VS0(52,4)=PAz*VS0(31,4)+WPz*VS1(31,4)+r1x2Z*VR1(31,4)&
   +2D0*r1x2Z*(VS0(16,4)-ExZpE*VS1(16,4))&
   +HfxZpE*VS1(31,1)
VS0(53,2)=PAz*VS0(32,2)+WPz*VS1(32,2)+r1x2Z*VR1(32,2)&
   +2D0*r1x2Z*(VS0(17,2)-ExZpE*VS1(17,2))
VS0(53,3)=PAz*VS0(32,3)+WPz*VS1(32,3)+r1x2Z*VR1(32,3)&
   +2D0*r1x2Z*(VS0(17,3)-ExZpE*VS1(17,3))
VS0(53,4)=PAz*VS0(32,4)+WPz*VS1(32,4)+r1x2Z*VR1(32,4)&
   +2D0*r1x2Z*(VS0(17,4)-ExZpE*VS1(17,4))&
   +HfxZpE*VS1(32,1)
VS0(54,2)=PAz*VS0(33,2)+WPz*VS1(33,2)+r1x2Z*VR1(33,2)&
   +3D0*r1x2Z*(VS0(18,2)-ExZpE*VS1(18,2))
VS0(54,3)=PAz*VS0(33,3)+WPz*VS1(33,3)+r1x2Z*VR1(33,3)&
   +3D0*r1x2Z*(VS0(18,3)-ExZpE*VS1(18,3))
VS0(54,4)=PAz*VS0(33,4)+WPz*VS1(33,4)+r1x2Z*VR1(33,4)&
   +3D0*r1x2Z*(VS0(18,4)-ExZpE*VS1(18,4))&
   +HfxZpE*VS1(33,1)
VS0(55,2)=PAz*VS0(34,2)+WPz*VS1(34,2)+r1x2Z*VR1(34,2)&
   +3D0*r1x2Z*(VS0(19,2)-ExZpE*VS1(19,2))
VS0(55,3)=PAz*VS0(34,3)+WPz*VS1(34,3)+r1x2Z*VR1(34,3)&
   +3D0*r1x2Z*(VS0(19,3)-ExZpE*VS1(19,3))
VS0(55,4)=PAz*VS0(34,4)+WPz*VS1(34,4)+r1x2Z*VR1(34,4)&
   +3D0*r1x2Z*(VS0(19,4)-ExZpE*VS1(19,4))&
   +HfxZpE*VS1(34,1)
VS0(56,2)=PAz*VS0(35,2)+WPz*VS1(35,2)+r1x2Z*VR1(35,2)&
   +4D0*r1x2Z*(VS0(20,2)-ExZpE*VS1(20,2))
VS0(56,3)=PAz*VS0(35,3)+WPz*VS1(35,3)+r1x2Z*VR1(35,3)&
   +4D0*r1x2Z*(VS0(20,3)-ExZpE*VS1(20,3))
VS0(56,4)=PAz*VS0(35,4)+WPz*VS1(35,4)+r1x2Z*VR1(35,4)&
   +4D0*r1x2Z*(VS0(20,4)-ExZpE*VS1(20,4))&
   +HfxZpE*VS1(35,1)
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRh0p0'
STOP
END SELECT
END SUBROUTINE MVRRh0p0