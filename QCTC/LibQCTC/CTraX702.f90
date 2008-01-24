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
      SUBROUTINE CTraX702(Q)
         USE PoleNodeType
         USE PoleGlobals
         USE ERIGlobals
         TYPE(PoleNode) :: Q
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)
        W(2)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)
        CTwo=W(1)+W(2)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)-Cpq(6)*Q%C(3)
        W(1)=-(Cpq(4)*Q%C(2))-Cpq(7)*Q%C(4)-Cpq(8)*Q%C(5)
        W(2)=-(Spq(4)*Q%S(2))-Spq(7)*Q%S(4)-Spq(8)*Q%S(5)
        CTwo=W(1)+W(2)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)+Cpq(3)*Q%C(2)
        W(2)=-(Cpq(5)*Q%C(2))-Cpq(7)*Q%C(3)+Cpq(6)*Q%C(4)
        W(3)=-(Cpq(8)*Q%C(4))+Cpq(7)*Q%C(5)-Cpq(9)*Q%C(5)
        W(4)=-(Spq(5)*Q%S(2))-Spq(8)*Q%S(4)+Spq(7)*Q%S(5)-Spq(9)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)-Q%C(3)*Spq(7)
        W(2)=-(Q%C(4)*Spq(8))-Q%C(5)*(Spq(7)+Spq(9))
        W(3)=(Cpq(3)+Cpq(5))*Q%S(2)+(Cpq(6)+Cpq(8))*Q%S(4)+(Cpq(7)+Cpq(9))*Q%S(5)
        STwo=W(1)+W(2)+W(3)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)+Cpq(6)*Q%C(1)+Cpq(10)*Q%C(3)
        W(1)=Cpq(7)*Q%C(2)+Cpq(11)*Q%C(4)+Cpq(12)*Q%C(5)
        W(2)=Spq(7)*Q%S(2)+Spq(11)*Q%S(4)+Spq(12)*Q%S(5)
        CTwo=W(1)+W(2)
        SPKetC(3)=SPKetC(3)+COne+CTwo*Two
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        W(1)=Cpq(4)*Q%C(0)+Cpq(7)*Q%C(1)+(-Cpq(6)+Cpq(8))*Q%C(2)+Cpq(11)*Q%C(3)
        W(2)=(-Cpq(10)+Cpq(12))*Q%C(4)+(-Cpq(11)+Cpq(13))*Q%C(5)
        W(3)=Spq(8)*Q%S(2)+Spq(12)*Q%S(4)+(-Spq(11)+Spq(13))*Q%S(5)
        CTwo=W(1)+W(2)+W(3)
        W(1)=Q%C(0)*Spq(4)+Q%C(1)*Spq(7)+Q%C(2)*Spq(8)+Q%C(3)*Spq(11)
        W(2)=Q%C(4)*Spq(12)+Q%C(5)*(Spq(11)+Spq(13))
        W(3)=-((Cpq(6)+Cpq(8))*Q%S(2))
        W(4)=-((Cpq(10)+Cpq(12))*Q%S(4))-(Cpq(11)+Cpq(13))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        W(1)=Cpq(5)*Q%C(0)+Cpq(8)*Q%C(1)+(-Cpq(7)+Cpq(9))*Q%C(2)+Cpq(12)*Q%C(3)
        W(2)=(-Cpq(11)+Cpq(13))*Q%C(4)+(Cpq(10)+Cpq(14))*Q%C(5)
        W(3)=(Spq(7)+Spq(9))*Q%S(2)+(Spq(11)+Spq(13))*Q%S(4)+Spq(14)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)
        W(1)=Q%C(0)*Spq(5)+Q%C(1)*Spq(8)+Q%C(2)*(-Spq(7)+Spq(9))+Q%C(3)*Spq(12)
        W(2)=Q%C(4)*(-Spq(11)+Spq(13))+Q%C(5)*Spq(14)
        W(3)=-((Cpq(7)+Cpq(9))*Q%S(2))
        W(4)=-((Cpq(11)+Cpq(13))*Q%S(4))+(Cpq(10)-Cpq(14))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
!       -------------------------------- Lp=3, Mp = 0-------------------------------------
        COne=-(Cpq(6)*Q%C(0))-Cpq(10)*Q%C(1)-Cpq(15)*Q%C(3)
        W(1)=-(Cpq(11)*Q%C(2))-Cpq(16)*Q%C(4)-Cpq(17)*Q%C(5)
        W(2)=-(Spq(11)*Q%S(2))-Spq(16)*Q%S(4)-Spq(17)*Q%S(5)
        CTwo=W(1)+W(2)
        SPKetC(6)=SPKetC(6)+COne+CTwo*Two
!       -------------------------------- Lp=3, Mp = 1-------------------------------------
        W(1)=-(Cpq(7)*Q%C(0))-Cpq(11)*Q%C(1)+Cpq(10)*Q%C(2)
        W(2)=-(Cpq(12)*Q%C(2))-Cpq(16)*Q%C(3)+Cpq(15)*Q%C(4)
        W(3)=-(Cpq(17)*Q%C(4))+Cpq(16)*Q%C(5)-Cpq(18)*Q%C(5)
        W(4)=-(Spq(12)*Q%S(2))-Spq(17)*Q%S(4)+Spq(16)*Q%S(5)-Spq(18)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(7))-Q%C(1)*Spq(11)-Q%C(2)*Spq(12)-Q%C(3)*Spq(16)
        W(2)=-(Q%C(4)*Spq(17))-Q%C(5)*(Spq(16)+Spq(18))
        W(3)=(Cpq(10)+Cpq(12))*Q%S(2)+(Cpq(15)+Cpq(17))*Q%S(4)+(Cpq(16)+Cpq(18))*Q%S(5)
        STwo=W(1)+W(2)+W(3)
        SPKetC(7)=SPKetC(7)+CTwo*Two
        SPKetS(7)=SPKetS(7)+STwo*Two
!       -------------------------------- Lp=3, Mp = 2-------------------------------------
        W(1)=-(Cpq(8)*Q%C(0))-Cpq(12)*Q%C(1)+Cpq(11)*Q%C(2)
        W(2)=-(Cpq(13)*Q%C(2))-Cpq(17)*Q%C(3)+Cpq(16)*Q%C(4)-Cpq(18)*Q%C(4)
        W(3)=-(Cpq(15)*Q%C(5))-Cpq(19)*Q%C(5)-Spq(11)*Q%S(2)
        W(4)=-(Spq(13)*Q%S(2))-Spq(16)*Q%S(4)-Spq(18)*Q%S(4)-Spq(19)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(8))-Q%C(1)*Spq(12)
        W(2)=Q%C(2)*(Spq(11)-Spq(13))-Q%C(3)*Spq(17)
        W(3)=Q%C(4)*(Spq(16)-Spq(18))-Q%C(5)*Spq(19)
        W(4)=(Cpq(11)+Cpq(13))*Q%S(2)+(Cpq(16)+Cpq(18))*Q%S(4)+(-Cpq(15)+Cpq(19))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(8)=SPKetC(8)+CTwo*Two
        SPKetS(8)=SPKetS(8)+STwo*Two
!       -------------------------------- Lp=3, Mp = 3-------------------------------------
        W(1)=-(Cpq(9)*Q%C(0))-Cpq(13)*Q%C(1)+Cpq(12)*Q%C(2)
        W(2)=-(Cpq(14)*Q%C(2))-Cpq(18)*Q%C(3)+Cpq(17)*Q%C(4)-Cpq(19)*Q%C(4)
        W(3)=-(Cpq(16)*Q%C(5))-Cpq(20)*Q%C(5)-Spq(12)*Q%S(2)-Spq(14)*Q%S(2)
        W(4)=-(Spq(17)*Q%S(4))-Spq(19)*Q%S(4)+Spq(16)*Q%S(5)-Spq(20)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(9))-Q%C(1)*Spq(13)
        W(2)=Q%C(2)*(Spq(12)-Spq(14))-Q%C(3)*Spq(18)
        W(3)=Q%C(4)*(Spq(17)-Spq(19))-Q%C(5)*(Spq(16)+Spq(20))
        W(4)=(Cpq(12)+Cpq(14))*Q%S(2)+(Cpq(17)+Cpq(19))*Q%S(4)+(-Cpq(16)+Cpq(20))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(9)=SPKetC(9)+CTwo*Two
        SPKetS(9)=SPKetS(9)+STwo*Two
!       -------------------------------- Lp=4, Mp = 0-------------------------------------
        COne=Cpq(10)*Q%C(0)+Cpq(15)*Q%C(1)+Cpq(21)*Q%C(3)
        W(1)=Cpq(16)*Q%C(2)+Cpq(22)*Q%C(4)+Cpq(23)*Q%C(5)
        W(2)=Spq(16)*Q%S(2)+Spq(22)*Q%S(4)+Spq(23)*Q%S(5)
        CTwo=W(1)+W(2)
        SPKetC(10)=SPKetC(10)+COne+CTwo*Two
!       -------------------------------- Lp=4, Mp = 1-------------------------------------
        W(1)=Cpq(11)*Q%C(0)+Cpq(16)*Q%C(1)+(-Cpq(15)+Cpq(17))*Q%C(2)+Cpq(22)*Q%C(3)
        W(2)=(-Cpq(21)+Cpq(23))*Q%C(4)+(-Cpq(22)+Cpq(24))*Q%C(5)
        W(3)=Spq(17)*Q%S(2)+Spq(23)*Q%S(4)+(-Spq(22)+Spq(24))*Q%S(5)
        CTwo=W(1)+W(2)+W(3)
        W(1)=Q%C(0)*Spq(11)+Q%C(1)*Spq(16)+Q%C(2)*Spq(17)+Q%C(3)*Spq(22)
        W(2)=Q%C(4)*Spq(23)+Q%C(5)*(Spq(22)+Spq(24))
        W(3)=-((Cpq(15)+Cpq(17))*Q%S(2))
        W(4)=-((Cpq(21)+Cpq(23))*Q%S(4))-(Cpq(22)+Cpq(24))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(11)=SPKetC(11)+CTwo*Two
        SPKetS(11)=SPKetS(11)+STwo*Two
!       -------------------------------- Lp=4, Mp = 2-------------------------------------
        W(1)=Cpq(12)*Q%C(0)+Cpq(17)*Q%C(1)+(-Cpq(16)+Cpq(18))*Q%C(2)+Cpq(23)*Q%C(3)
        W(2)=(-Cpq(22)+Cpq(24))*Q%C(4)+(Cpq(21)+Cpq(25))*Q%C(5)
        W(3)=(Spq(16)+Spq(18))*Q%S(2)+(Spq(22)+Spq(24))*Q%S(4)+Spq(25)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)
        W(1)=Q%C(0)*Spq(12)+Q%C(1)*Spq(17)+Q%C(2)*(-Spq(16)+Spq(18))+Q%C(3)*Spq(23)
        W(2)=Q%C(4)*(-Spq(22)+Spq(24))+Q%C(5)*Spq(25)
        W(3)=-((Cpq(16)+Cpq(18))*Q%S(2))
        W(4)=-((Cpq(22)+Cpq(24))*Q%S(4))+(Cpq(21)-Cpq(25))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(12)=SPKetC(12)+CTwo*Two
        SPKetS(12)=SPKetS(12)+STwo*Two
!       -------------------------------- Lp=4, Mp = 3-------------------------------------
        W(1)=Cpq(13)*Q%C(0)+Cpq(18)*Q%C(1)+(-Cpq(17)+Cpq(19))*Q%C(2)+Cpq(24)*Q%C(3)
        W(2)=(-Cpq(23)+Cpq(25))*Q%C(4)+(Cpq(22)+Cpq(26))*Q%C(5)
        W(3)=(Spq(17)+Spq(19))*Q%S(2)+(Spq(23)+Spq(25))*Q%S(4)+(-Spq(22)+Spq(26))*Q%S(5)
        CTwo=W(1)+W(2)+W(3)
        W(1)=Q%C(0)*Spq(13)+Q%C(1)*Spq(18)+Q%C(2)*(-Spq(17)+Spq(19))+Q%C(3)*Spq(24)
        W(2)=Q%C(4)*(-Spq(23)+Spq(25))+Q%C(5)*(Spq(22)+Spq(26))
        W(3)=-((Cpq(17)+Cpq(19))*Q%S(2))
        W(4)=-((Cpq(23)+Cpq(25))*Q%S(4))+(Cpq(22)-Cpq(26))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(13)=SPKetC(13)+CTwo*Two
        SPKetS(13)=SPKetS(13)+STwo*Two
!       -------------------------------- Lp=4, Mp = 4-------------------------------------
        W(1)=Cpq(14)*Q%C(0)+Cpq(19)*Q%C(1)+(-Cpq(18)+Cpq(20))*Q%C(2)+Cpq(25)*Q%C(3)
        W(2)=(-Cpq(24)+Cpq(26))*Q%C(4)+(Cpq(23)+Cpq(27))*Q%C(5)
        W(3)=(Spq(18)+Spq(20))*Q%S(2)+(Spq(24)+Spq(26))*Q%S(4)+(-Spq(23)+Spq(27))*Q%S(5)
        CTwo=W(1)+W(2)+W(3)
        W(1)=Q%C(0)*Spq(14)+Q%C(1)*Spq(19)+Q%C(2)*(-Spq(18)+Spq(20))+Q%C(3)*Spq(25)
        W(2)=Q%C(4)*(-Spq(24)+Spq(26))+Q%C(5)*(Spq(23)+Spq(27))
        W(3)=-((Cpq(18)+Cpq(20))*Q%S(2))
        W(4)=-((Cpq(24)+Cpq(26))*Q%S(4))+(Cpq(23)-Cpq(27))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(14)=SPKetC(14)+CTwo*Two
        SPKetS(14)=SPKetS(14)+STwo*Two
!       -------------------------------- Lp=5, Mp = 0-------------------------------------
        COne=-(Cpq(15)*Q%C(0))-Cpq(21)*Q%C(1)-Cpq(28)*Q%C(3)
        W(1)=-(Cpq(22)*Q%C(2))-Cpq(29)*Q%C(4)-Cpq(30)*Q%C(5)
        W(2)=-(Spq(22)*Q%S(2))-Spq(29)*Q%S(4)-Spq(30)*Q%S(5)
        CTwo=W(1)+W(2)
        SPKetC(15)=SPKetC(15)+COne+CTwo*Two
!       -------------------------------- Lp=5, Mp = 1-------------------------------------
        W(1)=-(Cpq(16)*Q%C(0))-Cpq(22)*Q%C(1)+Cpq(21)*Q%C(2)
        W(2)=-(Cpq(23)*Q%C(2))-Cpq(29)*Q%C(3)+Cpq(28)*Q%C(4)
        W(3)=-(Cpq(30)*Q%C(4))+Cpq(29)*Q%C(5)-Cpq(31)*Q%C(5)
        W(4)=-(Spq(23)*Q%S(2))-Spq(30)*Q%S(4)+Spq(29)*Q%S(5)-Spq(31)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(16))-Q%C(1)*Spq(22)-Q%C(2)*Spq(23)-Q%C(3)*Spq(29)
        W(2)=-(Q%C(4)*Spq(30))-Q%C(5)*(Spq(29)+Spq(31))
        W(3)=(Cpq(21)+Cpq(23))*Q%S(2)+(Cpq(28)+Cpq(30))*Q%S(4)+(Cpq(29)+Cpq(31))*Q%S(5)
        STwo=W(1)+W(2)+W(3)
        SPKetC(16)=SPKetC(16)+CTwo*Two
        SPKetS(16)=SPKetS(16)+STwo*Two
!       -------------------------------- Lp=5, Mp = 2-------------------------------------
        W(1)=-(Cpq(17)*Q%C(0))-Cpq(23)*Q%C(1)+Cpq(22)*Q%C(2)
        W(2)=-(Cpq(24)*Q%C(2))-Cpq(30)*Q%C(3)+Cpq(29)*Q%C(4)-Cpq(31)*Q%C(4)
        W(3)=-(Cpq(28)*Q%C(5))-Cpq(32)*Q%C(5)-Spq(22)*Q%S(2)
        W(4)=-(Spq(24)*Q%S(2))-Spq(29)*Q%S(4)-Spq(31)*Q%S(4)-Spq(32)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(17))-Q%C(1)*Spq(23)
        W(2)=Q%C(2)*(Spq(22)-Spq(24))-Q%C(3)*Spq(30)
        W(3)=Q%C(4)*(Spq(29)-Spq(31))-Q%C(5)*Spq(32)
        W(4)=(Cpq(22)+Cpq(24))*Q%S(2)+(Cpq(29)+Cpq(31))*Q%S(4)+(-Cpq(28)+Cpq(32))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(17)=SPKetC(17)+CTwo*Two
        SPKetS(17)=SPKetS(17)+STwo*Two
!       -------------------------------- Lp=5, Mp = 3-------------------------------------
        W(1)=-(Cpq(18)*Q%C(0))-Cpq(24)*Q%C(1)+Cpq(23)*Q%C(2)
        W(2)=-(Cpq(25)*Q%C(2))-Cpq(31)*Q%C(3)+Cpq(30)*Q%C(4)-Cpq(32)*Q%C(4)
        W(3)=-(Cpq(29)*Q%C(5))-Cpq(33)*Q%C(5)-Spq(23)*Q%S(2)-Spq(25)*Q%S(2)
        W(4)=-(Spq(30)*Q%S(4))-Spq(32)*Q%S(4)+Spq(29)*Q%S(5)-Spq(33)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(18))-Q%C(1)*Spq(24)
        W(2)=Q%C(2)*(Spq(23)-Spq(25))-Q%C(3)*Spq(31)
        W(3)=Q%C(4)*(Spq(30)-Spq(32))-Q%C(5)*(Spq(29)+Spq(33))
        W(4)=(Cpq(23)+Cpq(25))*Q%S(2)+(Cpq(30)+Cpq(32))*Q%S(4)+(-Cpq(29)+Cpq(33))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(18)=SPKetC(18)+CTwo*Two
        SPKetS(18)=SPKetS(18)+STwo*Two
!       -------------------------------- Lp=5, Mp = 4-------------------------------------
        W(1)=-(Cpq(19)*Q%C(0))-Cpq(25)*Q%C(1)+Cpq(24)*Q%C(2)
        W(2)=-(Cpq(26)*Q%C(2))-Cpq(32)*Q%C(3)+Cpq(31)*Q%C(4)-Cpq(33)*Q%C(4)
        W(3)=-(Cpq(30)*Q%C(5))-Cpq(34)*Q%C(5)-Spq(24)*Q%S(2)-Spq(26)*Q%S(2)
        W(4)=-(Spq(31)*Q%S(4))-Spq(33)*Q%S(4)+Spq(30)*Q%S(5)-Spq(34)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(19))-Q%C(1)*Spq(25)
        W(2)=Q%C(2)*(Spq(24)-Spq(26))-Q%C(3)*Spq(32)
        W(3)=Q%C(4)*(Spq(31)-Spq(33))-Q%C(5)*(Spq(30)+Spq(34))
        W(4)=(Cpq(24)+Cpq(26))*Q%S(2)+(Cpq(31)+Cpq(33))*Q%S(4)+(-Cpq(30)+Cpq(34))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(19)=SPKetC(19)+CTwo*Two
        SPKetS(19)=SPKetS(19)+STwo*Two
!       -------------------------------- Lp=5, Mp = 5-------------------------------------
        W(1)=-(Cpq(20)*Q%C(0))-Cpq(26)*Q%C(1)+Cpq(25)*Q%C(2)
        W(2)=-(Cpq(27)*Q%C(2))-Cpq(33)*Q%C(3)+Cpq(32)*Q%C(4)-Cpq(34)*Q%C(4)
        W(3)=-(Cpq(31)*Q%C(5))-Cpq(35)*Q%C(5)-Spq(25)*Q%S(2)-Spq(27)*Q%S(2)
        W(4)=-(Spq(32)*Q%S(4))-Spq(34)*Q%S(4)+Spq(31)*Q%S(5)-Spq(35)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(20))-Q%C(1)*Spq(26)
        W(2)=Q%C(2)*(Spq(25)-Spq(27))-Q%C(3)*Spq(33)
        W(3)=Q%C(4)*(Spq(32)-Spq(34))-Q%C(5)*(Spq(31)+Spq(35))
        W(4)=(Cpq(25)+Cpq(27))*Q%S(2)+(Cpq(32)+Cpq(34))*Q%S(4)+(-Cpq(31)+Cpq(35))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(20)=SPKetC(20)+CTwo*Two
        SPKetS(20)=SPKetS(20)+STwo*Two
!       -------------------------------- Lp=6, Mp = 0-------------------------------------
        COne=Cpq(21)*Q%C(0)+Cpq(28)*Q%C(1)+Cpq(36)*Q%C(3)
        W(1)=Cpq(29)*Q%C(2)+Cpq(37)*Q%C(4)+Cpq(38)*Q%C(5)
        W(2)=Spq(29)*Q%S(2)+Spq(37)*Q%S(4)+Spq(38)*Q%S(5)
        CTwo=W(1)+W(2)
        SPKetC(21)=SPKetC(21)+COne+CTwo*Two
!       -------------------------------- Lp=6, Mp = 1-------------------------------------
        W(1)=Cpq(22)*Q%C(0)+Cpq(29)*Q%C(1)+(-Cpq(28)+Cpq(30))*Q%C(2)+Cpq(37)*Q%C(3)
        W(2)=(-Cpq(36)+Cpq(38))*Q%C(4)+(-Cpq(37)+Cpq(39))*Q%C(5)
        W(3)=Spq(30)*Q%S(2)+Spq(38)*Q%S(4)+(-Spq(37)+Spq(39))*Q%S(5)
        CTwo=W(1)+W(2)+W(3)
        W(1)=Q%C(0)*Spq(22)+Q%C(1)*Spq(29)+Q%C(2)*Spq(30)+Q%C(3)*Spq(37)
        W(2)=Q%C(4)*Spq(38)+Q%C(5)*(Spq(37)+Spq(39))
        W(3)=-((Cpq(28)+Cpq(30))*Q%S(2))
        W(4)=-((Cpq(36)+Cpq(38))*Q%S(4))-(Cpq(37)+Cpq(39))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(22)=SPKetC(22)+CTwo*Two
        SPKetS(22)=SPKetS(22)+STwo*Two
!       -------------------------------- Lp=6, Mp = 2-------------------------------------
        W(1)=Cpq(23)*Q%C(0)+Cpq(30)*Q%C(1)+(-Cpq(29)+Cpq(31))*Q%C(2)+Cpq(38)*Q%C(3)
        W(2)=(-Cpq(37)+Cpq(39))*Q%C(4)+(Cpq(36)+Cpq(40))*Q%C(5)
        W(3)=(Spq(29)+Spq(31))*Q%S(2)+(Spq(37)+Spq(39))*Q%S(4)+Spq(40)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)
        W(1)=Q%C(0)*Spq(23)+Q%C(1)*Spq(30)+Q%C(2)*(-Spq(29)+Spq(31))+Q%C(3)*Spq(38)
        W(2)=Q%C(4)*(-Spq(37)+Spq(39))+Q%C(5)*Spq(40)
        W(3)=-((Cpq(29)+Cpq(31))*Q%S(2))
        W(4)=-((Cpq(37)+Cpq(39))*Q%S(4))+(Cpq(36)-Cpq(40))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(23)=SPKetC(23)+CTwo*Two
        SPKetS(23)=SPKetS(23)+STwo*Two
!       -------------------------------- Lp=6, Mp = 3-------------------------------------
        W(1)=Cpq(24)*Q%C(0)+Cpq(31)*Q%C(1)+(-Cpq(30)+Cpq(32))*Q%C(2)+Cpq(39)*Q%C(3)
        W(2)=(-Cpq(38)+Cpq(40))*Q%C(4)+(Cpq(37)+Cpq(41))*Q%C(5)
        W(3)=(Spq(30)+Spq(32))*Q%S(2)+(Spq(38)+Spq(40))*Q%S(4)+(-Spq(37)+Spq(41))*Q%S(5)
        CTwo=W(1)+W(2)+W(3)
        W(1)=Q%C(0)*Spq(24)+Q%C(1)*Spq(31)+Q%C(2)*(-Spq(30)+Spq(32))+Q%C(3)*Spq(39)
        W(2)=Q%C(4)*(-Spq(38)+Spq(40))+Q%C(5)*(Spq(37)+Spq(41))
        W(3)=-((Cpq(30)+Cpq(32))*Q%S(2))
        W(4)=-((Cpq(38)+Cpq(40))*Q%S(4))+(Cpq(37)-Cpq(41))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(24)=SPKetC(24)+CTwo*Two
        SPKetS(24)=SPKetS(24)+STwo*Two
!       -------------------------------- Lp=6, Mp = 4-------------------------------------
        W(1)=Cpq(25)*Q%C(0)+Cpq(32)*Q%C(1)+(-Cpq(31)+Cpq(33))*Q%C(2)+Cpq(40)*Q%C(3)
        W(2)=(-Cpq(39)+Cpq(41))*Q%C(4)+(Cpq(38)+Cpq(42))*Q%C(5)
        W(3)=(Spq(31)+Spq(33))*Q%S(2)+(Spq(39)+Spq(41))*Q%S(4)+(-Spq(38)+Spq(42))*Q%S(5)
        CTwo=W(1)+W(2)+W(3)
        W(1)=Q%C(0)*Spq(25)+Q%C(1)*Spq(32)+Q%C(2)*(-Spq(31)+Spq(33))+Q%C(3)*Spq(40)
        W(2)=Q%C(4)*(-Spq(39)+Spq(41))+Q%C(5)*(Spq(38)+Spq(42))
        W(3)=-((Cpq(31)+Cpq(33))*Q%S(2))
        W(4)=-((Cpq(39)+Cpq(41))*Q%S(4))+(Cpq(38)-Cpq(42))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(25)=SPKetC(25)+CTwo*Two
        SPKetS(25)=SPKetS(25)+STwo*Two
!       -------------------------------- Lp=6, Mp = 5-------------------------------------
        W(1)=Cpq(26)*Q%C(0)+Cpq(33)*Q%C(1)+(-Cpq(32)+Cpq(34))*Q%C(2)+Cpq(41)*Q%C(3)
        W(2)=(-Cpq(40)+Cpq(42))*Q%C(4)+(Cpq(39)+Cpq(43))*Q%C(5)
        W(3)=(Spq(32)+Spq(34))*Q%S(2)+(Spq(40)+Spq(42))*Q%S(4)+(-Spq(39)+Spq(43))*Q%S(5)
        CTwo=W(1)+W(2)+W(3)
        W(1)=Q%C(0)*Spq(26)+Q%C(1)*Spq(33)+Q%C(2)*(-Spq(32)+Spq(34))+Q%C(3)*Spq(41)
        W(2)=Q%C(4)*(-Spq(40)+Spq(42))+Q%C(5)*(Spq(39)+Spq(43))
        W(3)=-((Cpq(32)+Cpq(34))*Q%S(2))
        W(4)=-((Cpq(40)+Cpq(42))*Q%S(4))+(Cpq(39)-Cpq(43))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(26)=SPKetC(26)+CTwo*Two
        SPKetS(26)=SPKetS(26)+STwo*Two
!       -------------------------------- Lp=6, Mp = 6-------------------------------------
        W(1)=Cpq(27)*Q%C(0)+Cpq(34)*Q%C(1)+(-Cpq(33)+Cpq(35))*Q%C(2)+Cpq(42)*Q%C(3)
        W(2)=(-Cpq(41)+Cpq(43))*Q%C(4)+(Cpq(40)+Cpq(44))*Q%C(5)
        W(3)=(Spq(33)+Spq(35))*Q%S(2)+(Spq(41)+Spq(43))*Q%S(4)+(-Spq(40)+Spq(44))*Q%S(5)
        CTwo=W(1)+W(2)+W(3)
        W(1)=Q%C(0)*Spq(27)+Q%C(1)*Spq(34)+Q%C(2)*(-Spq(33)+Spq(35))+Q%C(3)*Spq(42)
        W(2)=Q%C(4)*(-Spq(41)+Spq(43))+Q%C(5)*(Spq(40)+Spq(44))
        W(3)=-((Cpq(33)+Cpq(35))*Q%S(2))
        W(4)=-((Cpq(41)+Cpq(43))*Q%S(4))+(Cpq(40)-Cpq(44))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(27)=SPKetC(27)+CTwo*Two
        SPKetS(27)=SPKetS(27)+STwo*Two
!       -------------------------------- Lp=7, Mp = 0-------------------------------------
        COne=-(Cpq(28)*Q%C(0))-Cpq(36)*Q%C(1)-Cpq(45)*Q%C(3)
        W(1)=-(Cpq(37)*Q%C(2))-Cpq(46)*Q%C(4)-Cpq(47)*Q%C(5)
        W(2)=-(Spq(37)*Q%S(2))-Spq(46)*Q%S(4)-Spq(47)*Q%S(5)
        CTwo=W(1)+W(2)
        SPKetC(28)=SPKetC(28)+COne+CTwo*Two
!       -------------------------------- Lp=7, Mp = 1-------------------------------------
        W(1)=-(Cpq(29)*Q%C(0))-Cpq(37)*Q%C(1)+Cpq(36)*Q%C(2)
        W(2)=-(Cpq(38)*Q%C(2))-Cpq(46)*Q%C(3)+Cpq(45)*Q%C(4)
        W(3)=-(Cpq(47)*Q%C(4))+Cpq(46)*Q%C(5)-Cpq(48)*Q%C(5)
        W(4)=-(Spq(38)*Q%S(2))-Spq(47)*Q%S(4)+Spq(46)*Q%S(5)-Spq(48)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(29))-Q%C(1)*Spq(37)-Q%C(2)*Spq(38)-Q%C(3)*Spq(46)
        W(2)=-(Q%C(4)*Spq(47))-Q%C(5)*(Spq(46)+Spq(48))
        W(3)=(Cpq(36)+Cpq(38))*Q%S(2)+(Cpq(45)+Cpq(47))*Q%S(4)+(Cpq(46)+Cpq(48))*Q%S(5)
        STwo=W(1)+W(2)+W(3)
        SPKetC(29)=SPKetC(29)+CTwo*Two
        SPKetS(29)=SPKetS(29)+STwo*Two
!       -------------------------------- Lp=7, Mp = 2-------------------------------------
        W(1)=-(Cpq(30)*Q%C(0))-Cpq(38)*Q%C(1)+Cpq(37)*Q%C(2)
        W(2)=-(Cpq(39)*Q%C(2))-Cpq(47)*Q%C(3)+Cpq(46)*Q%C(4)-Cpq(48)*Q%C(4)
        W(3)=-(Cpq(45)*Q%C(5))-Cpq(49)*Q%C(5)-Spq(37)*Q%S(2)
        W(4)=-(Spq(39)*Q%S(2))-Spq(46)*Q%S(4)-Spq(48)*Q%S(4)-Spq(49)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(30))-Q%C(1)*Spq(38)
        W(2)=Q%C(2)*(Spq(37)-Spq(39))-Q%C(3)*Spq(47)
        W(3)=Q%C(4)*(Spq(46)-Spq(48))-Q%C(5)*Spq(49)
        W(4)=(Cpq(37)+Cpq(39))*Q%S(2)+(Cpq(46)+Cpq(48))*Q%S(4)+(-Cpq(45)+Cpq(49))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(30)=SPKetC(30)+CTwo*Two
        SPKetS(30)=SPKetS(30)+STwo*Two
!       -------------------------------- Lp=7, Mp = 3-------------------------------------
        W(1)=-(Cpq(31)*Q%C(0))-Cpq(39)*Q%C(1)+Cpq(38)*Q%C(2)
        W(2)=-(Cpq(40)*Q%C(2))-Cpq(48)*Q%C(3)+Cpq(47)*Q%C(4)-Cpq(49)*Q%C(4)
        W(3)=-(Cpq(46)*Q%C(5))-Cpq(50)*Q%C(5)-Spq(38)*Q%S(2)-Spq(40)*Q%S(2)
        W(4)=-(Spq(47)*Q%S(4))-Spq(49)*Q%S(4)+Spq(46)*Q%S(5)-Spq(50)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(31))-Q%C(1)*Spq(39)
        W(2)=Q%C(2)*(Spq(38)-Spq(40))-Q%C(3)*Spq(48)
        W(3)=Q%C(4)*(Spq(47)-Spq(49))-Q%C(5)*(Spq(46)+Spq(50))
        W(4)=(Cpq(38)+Cpq(40))*Q%S(2)+(Cpq(47)+Cpq(49))*Q%S(4)+(-Cpq(46)+Cpq(50))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(31)=SPKetC(31)+CTwo*Two
        SPKetS(31)=SPKetS(31)+STwo*Two
!       -------------------------------- Lp=7, Mp = 4-------------------------------------
        W(1)=-(Cpq(32)*Q%C(0))-Cpq(40)*Q%C(1)+Cpq(39)*Q%C(2)
        W(2)=-(Cpq(41)*Q%C(2))-Cpq(49)*Q%C(3)+Cpq(48)*Q%C(4)-Cpq(50)*Q%C(4)
        W(3)=-(Cpq(47)*Q%C(5))-Cpq(51)*Q%C(5)-Spq(39)*Q%S(2)-Spq(41)*Q%S(2)
        W(4)=-(Spq(48)*Q%S(4))-Spq(50)*Q%S(4)+Spq(47)*Q%S(5)-Spq(51)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(32))-Q%C(1)*Spq(40)
        W(2)=Q%C(2)*(Spq(39)-Spq(41))-Q%C(3)*Spq(49)
        W(3)=Q%C(4)*(Spq(48)-Spq(50))-Q%C(5)*(Spq(47)+Spq(51))
        W(4)=(Cpq(39)+Cpq(41))*Q%S(2)+(Cpq(48)+Cpq(50))*Q%S(4)+(-Cpq(47)+Cpq(51))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(32)=SPKetC(32)+CTwo*Two
        SPKetS(32)=SPKetS(32)+STwo*Two
!       -------------------------------- Lp=7, Mp = 5-------------------------------------
        W(1)=-(Cpq(33)*Q%C(0))-Cpq(41)*Q%C(1)+Cpq(40)*Q%C(2)
        W(2)=-(Cpq(42)*Q%C(2))-Cpq(50)*Q%C(3)+Cpq(49)*Q%C(4)-Cpq(51)*Q%C(4)
        W(3)=-(Cpq(48)*Q%C(5))-Cpq(52)*Q%C(5)-Spq(40)*Q%S(2)-Spq(42)*Q%S(2)
        W(4)=-(Spq(49)*Q%S(4))-Spq(51)*Q%S(4)+Spq(48)*Q%S(5)-Spq(52)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(33))-Q%C(1)*Spq(41)
        W(2)=Q%C(2)*(Spq(40)-Spq(42))-Q%C(3)*Spq(50)
        W(3)=Q%C(4)*(Spq(49)-Spq(51))-Q%C(5)*(Spq(48)+Spq(52))
        W(4)=(Cpq(40)+Cpq(42))*Q%S(2)+(Cpq(49)+Cpq(51))*Q%S(4)+(-Cpq(48)+Cpq(52))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(33)=SPKetC(33)+CTwo*Two
        SPKetS(33)=SPKetS(33)+STwo*Two
!       -------------------------------- Lp=7, Mp = 6-------------------------------------
        W(1)=-(Cpq(34)*Q%C(0))-Cpq(42)*Q%C(1)+Cpq(41)*Q%C(2)
        W(2)=-(Cpq(43)*Q%C(2))-Cpq(51)*Q%C(3)+Cpq(50)*Q%C(4)-Cpq(52)*Q%C(4)
        W(3)=-(Cpq(49)*Q%C(5))-Cpq(53)*Q%C(5)-Spq(41)*Q%S(2)-Spq(43)*Q%S(2)
        W(4)=-(Spq(50)*Q%S(4))-Spq(52)*Q%S(4)+Spq(49)*Q%S(5)-Spq(53)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(34))-Q%C(1)*Spq(42)
        W(2)=Q%C(2)*(Spq(41)-Spq(43))-Q%C(3)*Spq(51)
        W(3)=Q%C(4)*(Spq(50)-Spq(52))-Q%C(5)*(Spq(49)+Spq(53))
        W(4)=(Cpq(41)+Cpq(43))*Q%S(2)+(Cpq(50)+Cpq(52))*Q%S(4)+(-Cpq(49)+Cpq(53))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(34)=SPKetC(34)+CTwo*Two
        SPKetS(34)=SPKetS(34)+STwo*Two
!       -------------------------------- Lp=7, Mp = 7-------------------------------------
        W(1)=-(Cpq(35)*Q%C(0))-Cpq(43)*Q%C(1)+Cpq(42)*Q%C(2)
        W(2)=-(Cpq(44)*Q%C(2))-Cpq(52)*Q%C(3)+Cpq(51)*Q%C(4)-Cpq(53)*Q%C(4)
        W(3)=-(Cpq(50)*Q%C(5))-Cpq(54)*Q%C(5)-Spq(42)*Q%S(2)-Spq(44)*Q%S(2)
        W(4)=-(Spq(51)*Q%S(4))-Spq(53)*Q%S(4)+Spq(50)*Q%S(5)-Spq(54)*Q%S(5)
        CTwo=W(1)+W(2)+W(3)+W(4)
        W(1)=-(Q%C(0)*Spq(35))-Q%C(1)*Spq(43)
        W(2)=Q%C(2)*(Spq(42)-Spq(44))-Q%C(3)*Spq(52)
        W(3)=Q%C(4)*(Spq(51)-Spq(53))-Q%C(5)*(Spq(50)+Spq(54))
        W(4)=(Cpq(42)+Cpq(44))*Q%S(2)+(Cpq(51)+Cpq(53))*Q%S(4)+(-Cpq(50)+Cpq(54))*Q%S(5)
        STwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(35)=SPKetC(35)+CTwo*Two
        SPKetS(35)=SPKetS(35)+STwo*Two
      END SUBROUTINE CTraX702
