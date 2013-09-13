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
   SUBROUTINE VRRs0i0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2E*VRR0(1,23)
      V(2)=3.D0*V(1)
      V(3)=r1x2E*ZxZpE*VRR1(1,23)
      V(4)=-3.D0*V(3)
      V(5)=r1x2E*VRR0(1,24)
      V(6)=r1x2E*ZxZpE*VRR1(1,24)
      V(7)=r1x2E*VRR0(1,27)
      V(8)=r1x2E*ZxZpE*VRR1(1,27)
      V(9)=r1x2E*VRR0(1,28)
      V(10)=2.D0*V(9)
      V(11)=r1x2E*ZxZpE*VRR1(1,28)
      V(12)=-2.D0*V(11)
      V(13)=2.D0*V(7)
      V(14)=-2.D0*V(8)
      V(15)=r1x2E*VRR0(1,30)
      V(16)=3.D0*V(15)
      V(17)=r1x2E*ZxZpE*VRR1(1,30)
      V(18)=-3.D0*V(17)
      V(19)=r1x2E*VRR0(1,31)
      V(20)=2.D0*V(19)
      V(21)=r1x2E*ZxZpE*VRR1(1,31)
      V(22)=-2.D0*V(21)
      V(23)=r1x2E*VRR0(1,32)
      V(24)=r1x2E*ZxZpE*VRR1(1,32)
      V(25)=3.D0*V(23)
      V(26)=-3.D0*V(24)
      V(27)=r1x2E*VRR0(1,33)
      V(28)=r1x2E*ZxZpE*VRR1(1,33)
      V(29)=r1x2E*VRR0(1,34)
      V(30)=r1x2E*ZxZpE*VRR1(1,34)
      VRR0(1,57)=5.D0*r1x2E*VRR0(1,21)+QCx*VRR0(1,36)-5.D0*r1x2E*ZxZpE*VRR1(1,21)+WQx*VRR1(1,36)
      VRR0(1,58)=4.D0*r1x2E*VRR0(1,22)+QCx*VRR0(1,37)-4.D0*r1x2E*ZxZpE*VRR1(1,22)+WQx*VRR1(1,37)
      VRR0(1,59)=V(2)+V(4)+QCx*VRR0(1,38)+WQx*VRR1(1,38)
      VRR0(1,60)=2.D0*V(5)-2.D0*V(6)+QCx*VRR0(1,39)+WQx*VRR1(1,39)
      VRR0(1,61)=V(2)+V(4)+QCy*VRR0(1,39)+WQy*VRR1(1,39)
      VRR0(1,62)=4.D0*V(5)-4.D0*V(6)+QCy*VRR0(1,40)+WQy*VRR1(1,40)
      VRR0(1,63)=5.D0*r1x2E*VRR0(1,25)+QCy*VRR0(1,41)-5.D0*r1x2E*ZxZpE*VRR1(1,25)+WQy*VRR1(1,41)
      VRR0(1,64)=4.D0*r1x2E*VRR0(1,26)+QCx*VRR0(1,42)-4.D0*r1x2E*ZxZpE*VRR1(1,26)+WQx*VRR1(1,42)
      VRR0(1,65)=3.D0*V(7)-3.D0*V(8)+QCx*VRR0(1,43)+WQx*VRR1(1,43)
      VRR0(1,66)=V(10)+V(12)+QCx*VRR0(1,44)+WQx*VRR1(1,44)
      VRR0(1,67)=V(13)+V(14)+QCy*VRR0(1,44)+WQy*VRR1(1,44)
      VRR0(1,68)=3.D0*V(9)-3.D0*V(11)+QCy*VRR0(1,45)+WQy*VRR1(1,45)
      VRR0(1,69)=4.D0*r1x2E*VRR0(1,29)+QCy*VRR0(1,46)-4.D0*r1x2E*ZxZpE*VRR1(1,29)+WQy*VRR1(1,46)
      VRR0(1,70)=V(16)+V(18)+QCx*VRR0(1,47)+WQx*VRR1(1,47)
      VRR0(1,71)=V(20)+V(22)+QCx*VRR0(1,48)+WQx*VRR1(1,48)
      VRR0(1,72)=V(23)-V(24)+QCx*VRR0(1,49)+WQx*VRR1(1,49)
      VRR0(1,73)=V(20)+V(22)+QCy*VRR0(1,49)+WQy*VRR1(1,49)
      VRR0(1,74)=V(25)+V(26)+QCy*VRR0(1,50)+WQy*VRR1(1,50)
      VRR0(1,75)=2.D0*V(27)-2.D0*V(28)+QCx*VRR0(1,51)+WQx*VRR1(1,51)
      VRR0(1,76)=V(13)+V(14)+QCz*VRR0(1,48)+WQz*VRR1(1,48)
      VRR0(1,77)=V(10)+V(12)+QCz*VRR0(1,49)+WQz*VRR1(1,49)
      VRR0(1,78)=2.D0*V(29)-2.D0*V(30)+QCy*VRR0(1,53)+WQy*VRR1(1,53)
      VRR0(1,79)=V(16)+V(18)+QCz*VRR0(1,51)+WQz*VRR1(1,51)
      VRR0(1,80)=3.D0*V(19)-3.D0*V(21)+QCz*VRR0(1,52)+WQz*VRR1(1,52)
      VRR0(1,81)=V(25)+V(26)+QCz*VRR0(1,53)+WQz*VRR1(1,53)
      VRR0(1,82)=4.D0*V(27)-4.D0*V(28)+QCz*VRR0(1,54)+WQz*VRR1(1,54)
      VRR0(1,83)=4.D0*V(29)-4.D0*V(30)+QCz*VRR0(1,55)+WQz*VRR1(1,55)
      VRR0(1,84)=5.D0*r1x2E*VRR0(1,35)+QCz*VRR0(1,56)-5.D0*r1x2E*ZxZpE*VRR1(1,35)+WQz*VRR1(1,56)
END SUBROUTINE VRRs0i0
SUBROUTINE MVRRs0i0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(1,57)=QCx*VS0(1,36)+WQx*VS1(1,36)-r1x2E*VR1(1,36)&
   +5D0*r1x2E*(VS0(1,21)-ZxZpE*VS1(1,21))
VS0(1,58)=QCx*VS0(1,37)+WQx*VS1(1,37)-r1x2E*VR1(1,37)&
   +4D0*r1x2E*(VS0(1,22)-ZxZpE*VS1(1,22))
VS0(1,59)=QCx*VS0(1,38)+WQx*VS1(1,38)-r1x2E*VR1(1,38)&
   +3D0*r1x2E*(VS0(1,23)-ZxZpE*VS1(1,23))
VS0(1,60)=QCx*VS0(1,39)+WQx*VS1(1,39)-r1x2E*VR1(1,39)&
   +2D0*r1x2E*(VS0(1,24)-ZxZpE*VS1(1,24))
VS0(1,61)=QCy*VS0(1,39)+WQy*VS1(1,39)&
   +3D0*r1x2E*(VS0(1,23)-ZxZpE*VS1(1,23))
VS0(1,62)=QCy*VS0(1,40)+WQy*VS1(1,40)&
   +4D0*r1x2E*(VS0(1,24)-ZxZpE*VS1(1,24))
VS0(1,63)=QCy*VS0(1,41)+WQy*VS1(1,41)&
   +5D0*r1x2E*(VS0(1,25)-ZxZpE*VS1(1,25))
VS0(1,64)=QCx*VS0(1,42)+WQx*VS1(1,42)-r1x2E*VR1(1,42)&
   +4D0*r1x2E*(VS0(1,26)-ZxZpE*VS1(1,26))
VS0(1,65)=QCx*VS0(1,43)+WQx*VS1(1,43)-r1x2E*VR1(1,43)&
   +3D0*r1x2E*(VS0(1,27)-ZxZpE*VS1(1,27))
VS0(1,66)=QCx*VS0(1,44)+WQx*VS1(1,44)-r1x2E*VR1(1,44)&
   +2D0*r1x2E*(VS0(1,28)-ZxZpE*VS1(1,28))
VS0(1,67)=QCy*VS0(1,44)+WQy*VS1(1,44)&
   +2D0*r1x2E*(VS0(1,27)-ZxZpE*VS1(1,27))
VS0(1,68)=QCy*VS0(1,45)+WQy*VS1(1,45)&
   +3D0*r1x2E*(VS0(1,28)-ZxZpE*VS1(1,28))
VS0(1,69)=QCy*VS0(1,46)+WQy*VS1(1,46)&
   +4D0*r1x2E*(VS0(1,29)-ZxZpE*VS1(1,29))
VS0(1,70)=QCx*VS0(1,47)+WQx*VS1(1,47)-r1x2E*VR1(1,47)&
   +3D0*r1x2E*(VS0(1,30)-ZxZpE*VS1(1,30))
VS0(1,71)=QCx*VS0(1,48)+WQx*VS1(1,48)-r1x2E*VR1(1,48)&
   +2D0*r1x2E*(VS0(1,31)-ZxZpE*VS1(1,31))
VS0(1,72)=QCx*VS0(1,49)+WQx*VS1(1,49)-r1x2E*VR1(1,49)&
   +r1x2E*(VS0(1,32)-ZxZpE*VS1(1,32))
VS0(1,73)=QCy*VS0(1,49)+WQy*VS1(1,49)&
   +2D0*r1x2E*(VS0(1,31)-ZxZpE*VS1(1,31))
VS0(1,74)=QCy*VS0(1,50)+WQy*VS1(1,50)&
   +3D0*r1x2E*(VS0(1,32)-ZxZpE*VS1(1,32))
VS0(1,75)=QCx*VS0(1,51)+WQx*VS1(1,51)-r1x2E*VR1(1,51)&
   +2D0*r1x2E*(VS0(1,33)-ZxZpE*VS1(1,33))
VS0(1,76)=QCz*VS0(1,48)+WQz*VS1(1,48)&
   +2D0*r1x2E*(VS0(1,27)-ZxZpE*VS1(1,27))
VS0(1,77)=QCz*VS0(1,49)+WQz*VS1(1,49)&
   +2D0*r1x2E*(VS0(1,28)-ZxZpE*VS1(1,28))
VS0(1,78)=QCy*VS0(1,53)+WQy*VS1(1,53)&
   +2D0*r1x2E*(VS0(1,34)-ZxZpE*VS1(1,34))
VS0(1,79)=QCz*VS0(1,51)+WQz*VS1(1,51)&
   +3D0*r1x2E*(VS0(1,30)-ZxZpE*VS1(1,30))
VS0(1,80)=QCz*VS0(1,52)+WQz*VS1(1,52)&
   +3D0*r1x2E*(VS0(1,31)-ZxZpE*VS1(1,31))
VS0(1,81)=QCz*VS0(1,53)+WQz*VS1(1,53)&
   +3D0*r1x2E*(VS0(1,32)-ZxZpE*VS1(1,32))
VS0(1,82)=QCz*VS0(1,54)+WQz*VS1(1,54)&
   +4D0*r1x2E*(VS0(1,33)-ZxZpE*VS1(1,33))
VS0(1,83)=QCz*VS0(1,55)+WQz*VS1(1,55)&
   +4D0*r1x2E*(VS0(1,34)-ZxZpE*VS1(1,34))
VS0(1,84)=QCz*VS0(1,56)+WQz*VS1(1,56)&
   +5D0*r1x2E*(VS0(1,35)-ZxZpE*VS1(1,35))
CASE(2)
VS0(1,57)=QCx*VS0(1,36)+WQx*VS1(1,36)&
   +5D0*r1x2E*(VS0(1,21)-ZxZpE*VS1(1,21))
VS0(1,58)=QCx*VS0(1,37)+WQx*VS1(1,37)&
   +4D0*r1x2E*(VS0(1,22)-ZxZpE*VS1(1,22))
VS0(1,59)=QCx*VS0(1,38)+WQx*VS1(1,38)&
   +3D0*r1x2E*(VS0(1,23)-ZxZpE*VS1(1,23))
VS0(1,60)=QCx*VS0(1,39)+WQx*VS1(1,39)&
   +2D0*r1x2E*(VS0(1,24)-ZxZpE*VS1(1,24))
VS0(1,61)=QCy*VS0(1,39)+WQy*VS1(1,39)-r1x2E*VR1(1,39)&
   +3D0*r1x2E*(VS0(1,23)-ZxZpE*VS1(1,23))
VS0(1,62)=QCy*VS0(1,40)+WQy*VS1(1,40)-r1x2E*VR1(1,40)&
   +4D0*r1x2E*(VS0(1,24)-ZxZpE*VS1(1,24))
VS0(1,63)=QCy*VS0(1,41)+WQy*VS1(1,41)-r1x2E*VR1(1,41)&
   +5D0*r1x2E*(VS0(1,25)-ZxZpE*VS1(1,25))
VS0(1,64)=QCx*VS0(1,42)+WQx*VS1(1,42)&
   +4D0*r1x2E*(VS0(1,26)-ZxZpE*VS1(1,26))
VS0(1,65)=QCx*VS0(1,43)+WQx*VS1(1,43)&
   +3D0*r1x2E*(VS0(1,27)-ZxZpE*VS1(1,27))
VS0(1,66)=QCx*VS0(1,44)+WQx*VS1(1,44)&
   +2D0*r1x2E*(VS0(1,28)-ZxZpE*VS1(1,28))
VS0(1,67)=QCy*VS0(1,44)+WQy*VS1(1,44)-r1x2E*VR1(1,44)&
   +2D0*r1x2E*(VS0(1,27)-ZxZpE*VS1(1,27))
VS0(1,68)=QCy*VS0(1,45)+WQy*VS1(1,45)-r1x2E*VR1(1,45)&
   +3D0*r1x2E*(VS0(1,28)-ZxZpE*VS1(1,28))
VS0(1,69)=QCy*VS0(1,46)+WQy*VS1(1,46)-r1x2E*VR1(1,46)&
   +4D0*r1x2E*(VS0(1,29)-ZxZpE*VS1(1,29))
VS0(1,70)=QCx*VS0(1,47)+WQx*VS1(1,47)&
   +3D0*r1x2E*(VS0(1,30)-ZxZpE*VS1(1,30))
VS0(1,71)=QCx*VS0(1,48)+WQx*VS1(1,48)&
   +2D0*r1x2E*(VS0(1,31)-ZxZpE*VS1(1,31))
VS0(1,72)=QCx*VS0(1,49)+WQx*VS1(1,49)&
   +r1x2E*(VS0(1,32)-ZxZpE*VS1(1,32))
VS0(1,73)=QCy*VS0(1,49)+WQy*VS1(1,49)-r1x2E*VR1(1,49)&
   +2D0*r1x2E*(VS0(1,31)-ZxZpE*VS1(1,31))
VS0(1,74)=QCy*VS0(1,50)+WQy*VS1(1,50)-r1x2E*VR1(1,50)&
   +3D0*r1x2E*(VS0(1,32)-ZxZpE*VS1(1,32))
VS0(1,75)=QCx*VS0(1,51)+WQx*VS1(1,51)&
   +2D0*r1x2E*(VS0(1,33)-ZxZpE*VS1(1,33))
VS0(1,76)=QCz*VS0(1,48)+WQz*VS1(1,48)&
   +2D0*r1x2E*(VS0(1,27)-ZxZpE*VS1(1,27))
VS0(1,77)=QCz*VS0(1,49)+WQz*VS1(1,49)&
   +2D0*r1x2E*(VS0(1,28)-ZxZpE*VS1(1,28))
VS0(1,78)=QCy*VS0(1,53)+WQy*VS1(1,53)-r1x2E*VR1(1,53)&
   +2D0*r1x2E*(VS0(1,34)-ZxZpE*VS1(1,34))
VS0(1,79)=QCz*VS0(1,51)+WQz*VS1(1,51)&
   +3D0*r1x2E*(VS0(1,30)-ZxZpE*VS1(1,30))
VS0(1,80)=QCz*VS0(1,52)+WQz*VS1(1,52)&
   +3D0*r1x2E*(VS0(1,31)-ZxZpE*VS1(1,31))
VS0(1,81)=QCz*VS0(1,53)+WQz*VS1(1,53)&
   +3D0*r1x2E*(VS0(1,32)-ZxZpE*VS1(1,32))
VS0(1,82)=QCz*VS0(1,54)+WQz*VS1(1,54)&
   +4D0*r1x2E*(VS0(1,33)-ZxZpE*VS1(1,33))
VS0(1,83)=QCz*VS0(1,55)+WQz*VS1(1,55)&
   +4D0*r1x2E*(VS0(1,34)-ZxZpE*VS1(1,34))
VS0(1,84)=QCz*VS0(1,56)+WQz*VS1(1,56)&
   +5D0*r1x2E*(VS0(1,35)-ZxZpE*VS1(1,35))
CASE(3)
VS0(1,57)=QCx*VS0(1,36)+WQx*VS1(1,36)&
   +5D0*r1x2E*(VS0(1,21)-ZxZpE*VS1(1,21))
VS0(1,58)=QCx*VS0(1,37)+WQx*VS1(1,37)&
   +4D0*r1x2E*(VS0(1,22)-ZxZpE*VS1(1,22))
VS0(1,59)=QCx*VS0(1,38)+WQx*VS1(1,38)&
   +3D0*r1x2E*(VS0(1,23)-ZxZpE*VS1(1,23))
VS0(1,60)=QCx*VS0(1,39)+WQx*VS1(1,39)&
   +2D0*r1x2E*(VS0(1,24)-ZxZpE*VS1(1,24))
VS0(1,61)=QCy*VS0(1,39)+WQy*VS1(1,39)&
   +3D0*r1x2E*(VS0(1,23)-ZxZpE*VS1(1,23))
VS0(1,62)=QCy*VS0(1,40)+WQy*VS1(1,40)&
   +4D0*r1x2E*(VS0(1,24)-ZxZpE*VS1(1,24))
VS0(1,63)=QCy*VS0(1,41)+WQy*VS1(1,41)&
   +5D0*r1x2E*(VS0(1,25)-ZxZpE*VS1(1,25))
VS0(1,64)=QCx*VS0(1,42)+WQx*VS1(1,42)&
   +4D0*r1x2E*(VS0(1,26)-ZxZpE*VS1(1,26))
VS0(1,65)=QCx*VS0(1,43)+WQx*VS1(1,43)&
   +3D0*r1x2E*(VS0(1,27)-ZxZpE*VS1(1,27))
VS0(1,66)=QCx*VS0(1,44)+WQx*VS1(1,44)&
   +2D0*r1x2E*(VS0(1,28)-ZxZpE*VS1(1,28))
VS0(1,67)=QCy*VS0(1,44)+WQy*VS1(1,44)&
   +2D0*r1x2E*(VS0(1,27)-ZxZpE*VS1(1,27))
VS0(1,68)=QCy*VS0(1,45)+WQy*VS1(1,45)&
   +3D0*r1x2E*(VS0(1,28)-ZxZpE*VS1(1,28))
VS0(1,69)=QCy*VS0(1,46)+WQy*VS1(1,46)&
   +4D0*r1x2E*(VS0(1,29)-ZxZpE*VS1(1,29))
VS0(1,70)=QCx*VS0(1,47)+WQx*VS1(1,47)&
   +3D0*r1x2E*(VS0(1,30)-ZxZpE*VS1(1,30))
VS0(1,71)=QCx*VS0(1,48)+WQx*VS1(1,48)&
   +2D0*r1x2E*(VS0(1,31)-ZxZpE*VS1(1,31))
VS0(1,72)=QCx*VS0(1,49)+WQx*VS1(1,49)&
   +r1x2E*(VS0(1,32)-ZxZpE*VS1(1,32))
VS0(1,73)=QCy*VS0(1,49)+WQy*VS1(1,49)&
   +2D0*r1x2E*(VS0(1,31)-ZxZpE*VS1(1,31))
VS0(1,74)=QCy*VS0(1,50)+WQy*VS1(1,50)&
   +3D0*r1x2E*(VS0(1,32)-ZxZpE*VS1(1,32))
VS0(1,75)=QCx*VS0(1,51)+WQx*VS1(1,51)&
   +2D0*r1x2E*(VS0(1,33)-ZxZpE*VS1(1,33))
VS0(1,76)=QCz*VS0(1,48)+WQz*VS1(1,48)-r1x2E*VR1(1,48)&
   +2D0*r1x2E*(VS0(1,27)-ZxZpE*VS1(1,27))
VS0(1,77)=QCz*VS0(1,49)+WQz*VS1(1,49)-r1x2E*VR1(1,49)&
   +2D0*r1x2E*(VS0(1,28)-ZxZpE*VS1(1,28))
VS0(1,78)=QCy*VS0(1,53)+WQy*VS1(1,53)&
   +2D0*r1x2E*(VS0(1,34)-ZxZpE*VS1(1,34))
VS0(1,79)=QCz*VS0(1,51)+WQz*VS1(1,51)-r1x2E*VR1(1,51)&
   +3D0*r1x2E*(VS0(1,30)-ZxZpE*VS1(1,30))
VS0(1,80)=QCz*VS0(1,52)+WQz*VS1(1,52)-r1x2E*VR1(1,52)&
   +3D0*r1x2E*(VS0(1,31)-ZxZpE*VS1(1,31))
VS0(1,81)=QCz*VS0(1,53)+WQz*VS1(1,53)-r1x2E*VR1(1,53)&
   +3D0*r1x2E*(VS0(1,32)-ZxZpE*VS1(1,32))
VS0(1,82)=QCz*VS0(1,54)+WQz*VS1(1,54)-r1x2E*VR1(1,54)&
   +4D0*r1x2E*(VS0(1,33)-ZxZpE*VS1(1,33))
VS0(1,83)=QCz*VS0(1,55)+WQz*VS1(1,55)-r1x2E*VR1(1,55)&
   +4D0*r1x2E*(VS0(1,34)-ZxZpE*VS1(1,34))
VS0(1,84)=QCz*VS0(1,56)+WQz*VS1(1,56)-r1x2E*VR1(1,56)&
   +5D0*r1x2E*(VS0(1,35)-ZxZpE*VS1(1,35))
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRs0i0'
STOP
END SELECT
END SUBROUTINE MVRRs0i0
