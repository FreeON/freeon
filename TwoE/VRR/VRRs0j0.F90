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
   SUBROUTINE VRRs0j0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2E*VRR0(1,38)
      V(2)=r1x2E*ZxZpE*VRR1(1,38)
      V(3)=r1x2E*VRR0(1,39)
      V(4)=r1x2E*ZxZpE*VRR1(1,39)
      V(5)=r1x2E*VRR0(1,44)
      V(6)=3.D0*V(5)
      V(7)=r1x2E*ZxZpE*VRR1(1,44)
      V(8)=-3.D0*V(7)
      V(9)=r1x2E*VRR0(1,45)
      V(10)=r1x2E*ZxZpE*VRR1(1,45)
      V(11)=r1x2E*VRR0(1,47)
      V(12)=r1x2E*ZxZpE*VRR1(1,47)
      V(13)=r1x2E*VRR0(1,48)
      V(14)=3.D0*V(13)
      V(15)=r1x2E*ZxZpE*VRR1(1,48)
      V(16)=-3.D0*V(15)
      V(17)=r1x2E*VRR0(1,49)
      V(18)=r1x2E*ZxZpE*VRR1(1,49)
      V(19)=3.D0*V(17)
      V(20)=-3.D0*V(18)
      V(21)=r1x2E*VRR0(1,50)
      V(22)=r1x2E*ZxZpE*VRR1(1,50)
      V(23)=r1x2E*VRR0(1,51)
      V(24)=r1x2E*ZxZpE*VRR1(1,51)
      V(25)=r1x2E*VRR0(1,52)
      V(26)=2.D0*V(25)
      V(27)=r1x2E*ZxZpE*VRR1(1,52)
      V(28)=-2.D0*V(27)
      V(29)=r1x2E*VRR0(1,53)
      V(30)=r1x2E*ZxZpE*VRR1(1,53)
      VRR0(1,85)=6.D0*r1x2E*VRR0(1,36)+QCx*VRR0(1,57)-6.D0*r1x2E*ZxZpE*VRR1(1,36)+WQx*VRR1(1,57)
      VRR0(1,86)=5.D0*r1x2E*VRR0(1,37)+QCx*VRR0(1,58)-5.D0*r1x2E*ZxZpE*VRR1(1,37)+WQx*VRR1(1,58)
      VRR0(1,87)=4.D0*V(1)-4.D0*V(2)+QCx*VRR0(1,59)+WQx*VRR1(1,59)
      VRR0(1,88)=3.D0*V(3)-3.D0*V(4)+QCx*VRR0(1,60)+WQx*VRR1(1,60)
      VRR0(1,89)=3.D0*V(1)-3.D0*V(2)+QCy*VRR0(1,60)+WQy*VRR1(1,60)
      VRR0(1,90)=4.D0*V(3)-4.D0*V(4)+QCy*VRR0(1,61)+WQy*VRR1(1,61)
      VRR0(1,91)=5.D0*r1x2E*VRR0(1,40)+QCy*VRR0(1,62)-5.D0*r1x2E*ZxZpE*VRR1(1,40)+WQy*VRR1(1,62)
      VRR0(1,92)=6.D0*r1x2E*VRR0(1,41)+QCy*VRR0(1,63)-6.D0*r1x2E*ZxZpE*VRR1(1,41)+WQy*VRR1(1,63)
      VRR0(1,93)=5.D0*r1x2E*VRR0(1,42)+QCx*VRR0(1,64)-5.D0*r1x2E*ZxZpE*VRR1(1,42)+WQx*VRR1(1,64)
      VRR0(1,94)=4.D0*r1x2E*VRR0(1,43)+QCx*VRR0(1,65)-4.D0*r1x2E*ZxZpE*VRR1(1,43)+WQx*VRR1(1,65)
      VRR0(1,95)=V(6)+V(8)+QCx*VRR0(1,66)+WQx*VRR1(1,66)
      VRR0(1,96)=2.D0*V(9)-2.D0*V(10)+QCx*VRR0(1,67)+WQx*VRR1(1,67)
      VRR0(1,97)=V(6)+V(8)+QCy*VRR0(1,67)+WQy*VRR1(1,67)
      VRR0(1,98)=4.D0*V(9)-4.D0*V(10)+QCy*VRR0(1,68)+WQy*VRR1(1,68)
      VRR0(1,99)=5.D0*r1x2E*VRR0(1,46)+QCy*VRR0(1,69)-5.D0*r1x2E*ZxZpE*VRR1(1,46)+WQy*VRR1(1,69)
      VRR0(1,100)=4.D0*V(11)-4.D0*V(12)+QCx*VRR0(1,70)+WQx*VRR1(1,70)
      VRR0(1,101)=V(14)+V(16)+QCx*VRR0(1,71)+WQx*VRR1(1,71)
      VRR0(1,102)=2.D0*V(17)-2.D0*V(18)+QCx*VRR0(1,72)+WQx*VRR1(1,72)
      VRR0(1,103)=2.D0*V(13)-2.D0*V(15)+QCy*VRR0(1,72)+WQy*VRR1(1,72)
      VRR0(1,104)=V(19)+V(20)+QCy*VRR0(1,73)+WQy*VRR1(1,73)
      VRR0(1,105)=4.D0*V(21)-4.D0*V(22)+QCy*VRR0(1,74)+WQy*VRR1(1,74)
      VRR0(1,106)=3.D0*V(23)-3.D0*V(24)+QCx*VRR0(1,75)+WQx*VRR1(1,75)
      VRR0(1,107)=V(26)+V(28)+QCx*VRR0(1,76)+WQx*VRR1(1,76)
      VRR0(1,108)=2.D0*V(5)-2.D0*V(7)+QCz*VRR0(1,72)+WQz*VRR1(1,72)
      VRR0(1,109)=V(26)+V(28)+QCy*VRR0(1,77)+WQy*VRR1(1,77)
      VRR0(1,110)=3.D0*V(29)-3.D0*V(30)+QCy*VRR0(1,78)+WQy*VRR1(1,78)
      VRR0(1,111)=3.D0*V(11)-3.D0*V(12)+QCz*VRR0(1,75)+WQz*VRR1(1,75)
      VRR0(1,112)=V(14)+V(16)+QCz*VRR0(1,76)+WQz*VRR1(1,76)
      VRR0(1,113)=V(19)+V(20)+QCz*VRR0(1,77)+WQz*VRR1(1,77)
      VRR0(1,114)=3.D0*V(21)-3.D0*V(22)+QCz*VRR0(1,78)+WQz*VRR1(1,78)
      VRR0(1,115)=4.D0*V(23)-4.D0*V(24)+QCz*VRR0(1,79)+WQz*VRR1(1,79)
      VRR0(1,116)=4.D0*V(25)-4.D0*V(27)+QCz*VRR0(1,80)+WQz*VRR1(1,80)
      VRR0(1,117)=4.D0*V(29)-4.D0*V(30)+QCz*VRR0(1,81)+WQz*VRR1(1,81)
      VRR0(1,118)=5.D0*r1x2E*VRR0(1,54)+QCz*VRR0(1,82)-5.D0*r1x2E*ZxZpE*VRR1(1,54)+WQz*VRR1(1,82)
      VRR0(1,119)=5.D0*r1x2E*VRR0(1,55)+QCz*VRR0(1,83)-5.D0*r1x2E*ZxZpE*VRR1(1,55)+WQz*VRR1(1,83)
      VRR0(1,120)=6.D0*r1x2E*VRR0(1,56)+QCz*VRR0(1,84)-6.D0*r1x2E*ZxZpE*VRR1(1,56)+WQz*VRR1(1,84)
END SUBROUTINE VRRs0j0
