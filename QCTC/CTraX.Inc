      SELECT CASE(LCode)
      CASE(0)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)
        SPKetC(0)=SPKetC(0)+COne
      CASE(1)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)
        CTwo=Cpq(2)*Q%C(2)+Spq(2)*Q%S(2)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
      CASE(2)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)
        W(2)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)
        CTwo=W(1)+W(2)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
      CASE(3)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)+Cpq(6)*Q%C(6)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)
        W(2)=Cpq(7)*Q%C(7)+Cpq(8)*Q%C(8)+Cpq(9)*Q%C(9)
        W(3)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)
        W(4)=Spq(7)*Q%S(7)+Spq(8)*Q%S(8)+Spq(9)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
      CASE(4)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)+Cpq(6)*Q%C(6)+Cpq(10)*Q%C(10)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)+Cpq(7)*Q%C(7)+Cpq(8)*Q%C(8)
        W(2)=Cpq(9)*Q%C(9)+Cpq(11)*Q%C(11)+Cpq(12)*Q%C(12)+Cpq(13)*Q%C(13)+Cpq(14)*Q%C(14)
        W(3)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)+Spq(7)*Q%S(7)+Spq(8)*Q%S(8)
        W(4)=Spq(9)*Q%S(9)+Spq(11)*Q%S(11)+Spq(12)*Q%S(12)+Spq(13)*Q%S(13)+Spq(14)*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
      CASE(5)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        W(1)=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)
        W(2)=Cpq(6)*Q%C(6)+Cpq(10)*Q%C(10)+Cpq(15)*Q%C(15)
        COne=W(1)+W(2)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)
        W(2)=Cpq(7)*Q%C(7)+Cpq(8)*Q%C(8)+Cpq(9)*Q%C(9)+Cpq(11)*Q%C(11)
        W(3)=Cpq(12)*Q%C(12)+Cpq(13)*Q%C(13)+Cpq(14)*Q%C(14)+Cpq(16)*Q%C(16)
        W(4)=Cpq(17)*Q%C(17)+Cpq(18)*Q%C(18)+Cpq(19)*Q%C(19)+Cpq(20)*Q%C(20)
        W(5)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)
        W(6)=Spq(7)*Q%S(7)+Spq(8)*Q%S(8)+Spq(9)*Q%S(9)+Spq(11)*Q%S(11)
        W(7)=Spq(12)*Q%S(12)+Spq(13)*Q%S(13)+Spq(14)*Q%S(14)+Spq(16)*Q%S(16)
        W(8)=Spq(17)*Q%S(17)+Spq(18)*Q%S(18)+Spq(19)*Q%S(19)+Spq(20)*Q%S(20)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
      CASE(6)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        W(1)=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)
        W(2)=Cpq(6)*Q%C(6)+Cpq(10)*Q%C(10)+Cpq(15)*Q%C(15)+Cpq(21)*Q%C(21)
        COne=W(1)+W(2)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)+Cpq(7)*Q%C(7)+Cpq(8)*Q%C(8)
        W(2)=Cpq(9)*Q%C(9)+Cpq(11)*Q%C(11)+Cpq(12)*Q%C(12)+Cpq(13)*Q%C(13)+Cpq(14)*Q%C(14)
        W(3)=Cpq(16)*Q%C(16)+Cpq(17)*Q%C(17)+Cpq(18)*Q%C(18)+Cpq(19)*Q%C(19)+Cpq(20)*Q%C(20)
        W(4)=Cpq(22)*Q%C(22)+Cpq(23)*Q%C(23)+Cpq(24)*Q%C(24)
        W(5)=Cpq(25)*Q%C(25)+Cpq(26)*Q%C(26)+Cpq(27)*Q%C(27)
        W(6)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)+Spq(7)*Q%S(7)+Spq(8)*Q%S(8)
        W(7)=Spq(9)*Q%S(9)+Spq(11)*Q%S(11)+Spq(12)*Q%S(12)+Spq(13)*Q%S(13)+Spq(14)*Q%S(14)
        W(8)=Spq(16)*Q%S(16)+Spq(17)*Q%S(17)+Spq(18)*Q%S(18)+Spq(19)*Q%S(19)+Spq(20)*Q%S(20)
        W(9)=Spq(22)*Q%S(22)+Spq(23)*Q%S(23)+Spq(24)*Q%S(24)
        W(10)=Spq(25)*Q%S(25)+Spq(26)*Q%S(26)+Spq(27)*Q%S(27)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
      CASE(7)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        W(1)=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)+Cpq(6)*Q%C(6)
        W(2)=Cpq(10)*Q%C(10)+Cpq(15)*Q%C(15)+Cpq(21)*Q%C(21)+Cpq(28)*Q%C(28)
        COne=W(1)+W(2)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)
        W(2)=Cpq(7)*Q%C(7)+Cpq(8)*Q%C(8)+Cpq(9)*Q%C(9)+Cpq(11)*Q%C(11)
        W(3)=Cpq(12)*Q%C(12)+Cpq(13)*Q%C(13)+Cpq(14)*Q%C(14)
        W(4)=Cpq(16)*Q%C(16)+Cpq(17)*Q%C(17)+Cpq(18)*Q%C(18)+Cpq(19)*Q%C(19)
        W(5)=Cpq(20)*Q%C(20)+Cpq(22)*Q%C(22)+Cpq(23)*Q%C(23)
        W(6)=Cpq(24)*Q%C(24)+Cpq(25)*Q%C(25)+Cpq(26)*Q%C(26)+Cpq(27)*Q%C(27)
        W(7)=Cpq(29)*Q%C(29)+Cpq(30)*Q%C(30)+Cpq(31)*Q%C(31)
        W(8)=Cpq(32)*Q%C(32)+Cpq(33)*Q%C(33)+Cpq(34)*Q%C(34)+Cpq(35)*Q%C(35)
        W(9)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)
        W(10)=Spq(7)*Q%S(7)+Spq(8)*Q%S(8)+Spq(9)*Q%S(9)+Spq(11)*Q%S(11)
        W(11)=Spq(12)*Q%S(12)+Spq(13)*Q%S(13)+Spq(14)*Q%S(14)
        W(12)=Spq(16)*Q%S(16)+Spq(17)*Q%S(17)+Spq(18)*Q%S(18)+Spq(19)*Q%S(19)
        W(13)=Spq(20)*Q%S(20)+Spq(22)*Q%S(22)+Spq(23)*Q%S(23)
        W(14)=Spq(24)*Q%S(24)+Spq(25)*Q%S(25)+Spq(26)*Q%S(26)+Spq(27)*Q%S(27)
        W(15)=Spq(29)*Q%S(29)+Spq(30)*Q%S(30)+Spq(31)*Q%S(31)
        W(16)=Spq(32)*Q%S(32)+Spq(33)*Q%S(33)+Spq(34)*Q%S(34)+Spq(35)*Q%S(35)
        W(17)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(18)=W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)+W(16)
        CTwo=W(17)+W(18)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
      CASE(100)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)
        SPKetC(0)=SPKetC(0)+COne
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))
        SPKetC(1)=SPKetC(1)+COne
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        CTwo=-(Cpq(2)*Q%C(0))
        STwo=-(Q%C(0)*Spq(2))
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
      CASE(101)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)
        CTwo=Cpq(2)*Q%C(2)+Spq(2)*Q%S(2)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)
        CTwo=-(Cpq(4)*Q%C(2))-Spq(4)*Q%S(2)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)
        W(2)=Cpq(3)*Q%C(2)-Cpq(5)*Q%C(2)-Spq(5)*Q%S(2)
        CTwo=W(1)+W(2)
        STwo=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)+(Cpq(3)+Cpq(5))*Q%S(2)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
      CASE(102)
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
      CASE(103)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)+Cpq(6)*Q%C(6)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)
        W(2)=Cpq(7)*Q%C(7)+Cpq(8)*Q%C(8)+Cpq(9)*Q%C(9)
        W(3)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)
        W(4)=Spq(7)*Q%S(7)+Spq(8)*Q%S(8)+Spq(9)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)-Cpq(6)*Q%C(3)-Cpq(10)*Q%C(6)
        W(1)=-(Cpq(4)*Q%C(2))-Cpq(7)*Q%C(4)-Cpq(8)*Q%C(5)
        W(2)=-(Cpq(11)*Q%C(7))-Cpq(12)*Q%C(8)-Cpq(13)*Q%C(9)
        W(3)=-(Spq(4)*Q%S(2))-Spq(7)*Q%S(4)-Spq(8)*Q%S(5)
        W(4)=-(Spq(11)*Q%S(7))-Spq(12)*Q%S(8)-Spq(13)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)+Cpq(3)*Q%C(2)
        W(2)=-(Cpq(5)*Q%C(2))-Cpq(7)*Q%C(3)+Cpq(6)*Q%C(4)
        W(3)=-(Cpq(8)*Q%C(4))+Cpq(7)*Q%C(5)-Cpq(9)*Q%C(5)
        W(4)=-(Cpq(11)*Q%C(6))+Cpq(10)*Q%C(7)-Cpq(12)*Q%C(7)
        W(5)=Cpq(11)*Q%C(8)-Cpq(13)*Q%C(8)+Cpq(12)*Q%C(9)
        W(6)=-(Cpq(14)*Q%C(9))-Spq(5)*Q%S(2)-Spq(8)*Q%S(4)
        W(7)=Spq(7)*Q%S(5)-Spq(9)*Q%S(5)-Spq(12)*Q%S(7)
        W(8)=Spq(11)*Q%S(8)-Spq(13)*Q%S(8)+Spq(12)*Q%S(9)-Spq(14)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(1)=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)-Q%C(3)*Spq(7)
        W(2)=-(Q%C(4)*Spq(8))-Q%C(5)*(Spq(7)+Spq(9))
        W(3)=-(Q%C(6)*Spq(11))-Q%C(7)*Spq(12)
        W(4)=-(Q%C(8)*(Spq(11)+Spq(13)))-Q%C(9)*(Spq(12)+Spq(14))
        W(5)=(Cpq(3)+Cpq(5))*Q%S(2)+(Cpq(6)+Cpq(8))*Q%S(4)
        W(6)=(Cpq(7)+Cpq(9))*Q%S(5)+(Cpq(10)+Cpq(12))*Q%S(7)
        W(7)=(Cpq(11)+Cpq(13))*Q%S(8)+(Cpq(12)+Cpq(14))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
      CASE(104)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)+Cpq(6)*Q%C(6)+Cpq(10)*Q%C(10)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)+Cpq(7)*Q%C(7)+Cpq(8)*Q%C(8)
        W(2)=Cpq(9)*Q%C(9)+Cpq(11)*Q%C(11)+Cpq(12)*Q%C(12)+Cpq(13)*Q%C(13)+Cpq(14)*Q%C(14)
        W(3)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)+Spq(7)*Q%S(7)+Spq(8)*Q%S(8)
        W(4)=Spq(9)*Q%S(9)+Spq(11)*Q%S(11)+Spq(12)*Q%S(12)+Spq(13)*Q%S(13)+Spq(14)*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        W(1)=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)
        W(2)=-(Cpq(6)*Q%C(3))-Cpq(10)*Q%C(6)-Cpq(15)*Q%C(10)
        COne=W(1)+W(2)
        W(1)=-(Cpq(4)*Q%C(2))-Cpq(7)*Q%C(4)
        W(2)=-(Cpq(8)*Q%C(5))-Cpq(11)*Q%C(7)-Cpq(12)*Q%C(8)
        W(3)=-(Cpq(13)*Q%C(9))-Cpq(16)*Q%C(11)
        W(4)=-(Cpq(17)*Q%C(12))-Cpq(18)*Q%C(13)-Cpq(19)*Q%C(14)
        W(5)=-(Spq(4)*Q%S(2))-Spq(7)*Q%S(4)
        W(6)=-(Spq(8)*Q%S(5))-Spq(11)*Q%S(7)-Spq(12)*Q%S(8)
        W(7)=-(Spq(13)*Q%S(9))-Spq(16)*Q%S(11)
        W(8)=-(Spq(17)*Q%S(12))-Spq(18)*Q%S(13)-Spq(19)*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)+(Cpq(3)-Cpq(5))*Q%C(2)
        W(2)=-(Cpq(7)*Q%C(3))+(Cpq(6)-Cpq(8))*Q%C(4)+(Cpq(7)-Cpq(9))*Q%C(5)
        W(3)=-(Cpq(11)*Q%C(6))+(Cpq(10)-Cpq(12))*Q%C(7)+(Cpq(11)-Cpq(13))*Q%C(8)
        W(4)=(Cpq(12)-Cpq(14))*Q%C(9)-Cpq(16)*Q%C(10)+(Cpq(15)-Cpq(17))*Q%C(11)
        W(5)=(Cpq(16)-Cpq(18))*Q%C(12)
        W(6)=(Cpq(17)-Cpq(19))*Q%C(13)+(Cpq(18)-Cpq(20))*Q%C(14)
        W(7)=-(Spq(5)*Q%S(2))-Spq(8)*Q%S(4)+(Spq(7)-Spq(9))*Q%S(5)
        W(8)=-(Spq(12)*Q%S(7))+(Spq(11)-Spq(13))*Q%S(8)+(Spq(12)-Spq(14))*Q%S(9)
        W(9)=-(Spq(17)*Q%S(11))+(Spq(16)-Spq(18))*Q%S(12)
        W(10)=(Spq(17)-Spq(19))*Q%S(13)+(Spq(18)-Spq(20))*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)
        W(1)=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)
        W(2)=-(Q%C(3)*Spq(7))-Q%C(4)*Spq(8)-Q%C(5)*(Spq(7)+Spq(9))
        W(3)=-(Q%C(6)*Spq(11))-Q%C(7)*Spq(12)-Q%C(8)*(Spq(11)+Spq(13))
        W(4)=-(Q%C(9)*(Spq(12)+Spq(14)))-Q%C(10)*Spq(16)-Q%C(11)*Spq(17)
        W(5)=-(Q%C(12)*(Spq(16)+Spq(18)))
        W(6)=-(Q%C(13)*(Spq(17)+Spq(19)))-Q%C(14)*(Spq(18)+Spq(20))
        W(7)=(Cpq(3)+Cpq(5))*Q%S(2)+(Cpq(6)+Cpq(8))*Q%S(4)+(Cpq(7)+Cpq(9))*Q%S(5)
        W(8)=(Cpq(10)+Cpq(12))*Q%S(7)+(Cpq(11)+Cpq(13))*Q%S(8)+(Cpq(12)+Cpq(14))*Q%S(9)
        W(9)=(Cpq(15)+Cpq(17))*Q%S(11)+(Cpq(16)+Cpq(18))*Q%S(12)
        W(10)=(Cpq(17)+Cpq(19))*Q%S(13)+(Cpq(18)+Cpq(20))*Q%S(14)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
      CASE(105)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        W(1)=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)
        W(2)=Cpq(6)*Q%C(6)+Cpq(10)*Q%C(10)+Cpq(15)*Q%C(15)
        COne=W(1)+W(2)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)
        W(2)=Cpq(7)*Q%C(7)+Cpq(8)*Q%C(8)+Cpq(9)*Q%C(9)+Cpq(11)*Q%C(11)
        W(3)=Cpq(12)*Q%C(12)+Cpq(13)*Q%C(13)+Cpq(14)*Q%C(14)+Cpq(16)*Q%C(16)
        W(4)=Cpq(17)*Q%C(17)+Cpq(18)*Q%C(18)+Cpq(19)*Q%C(19)+Cpq(20)*Q%C(20)
        W(5)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)
        W(6)=Spq(7)*Q%S(7)+Spq(8)*Q%S(8)+Spq(9)*Q%S(9)+Spq(11)*Q%S(11)
        W(7)=Spq(12)*Q%S(12)+Spq(13)*Q%S(13)+Spq(14)*Q%S(14)+Spq(16)*Q%S(16)
        W(8)=Spq(17)*Q%S(17)+Spq(18)*Q%S(18)+Spq(19)*Q%S(19)+Spq(20)*Q%S(20)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        W(1)=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)-Cpq(6)*Q%C(3)
        W(2)=-(Cpq(10)*Q%C(6))-Cpq(15)*Q%C(10)-Cpq(21)*Q%C(15)
        COne=W(1)+W(2)
        W(1)=-(Cpq(4)*Q%C(2))-Cpq(7)*Q%C(4)-Cpq(8)*Q%C(5)
        W(2)=-(Cpq(11)*Q%C(7))-Cpq(12)*Q%C(8)-Cpq(13)*Q%C(9)-Cpq(16)*Q%C(11)
        W(3)=-(Cpq(17)*Q%C(12))-Cpq(18)*Q%C(13)-Cpq(19)*Q%C(14)-Cpq(22)*Q%C(16)
        W(4)=-(Cpq(23)*Q%C(17))-Cpq(24)*Q%C(18)-Cpq(25)*Q%C(19)-Cpq(26)*Q%C(20)
        W(5)=-(Spq(4)*Q%S(2))-Spq(7)*Q%S(4)-Spq(8)*Q%S(5)
        W(6)=-(Spq(11)*Q%S(7))-Spq(12)*Q%S(8)-Spq(13)*Q%S(9)-Spq(16)*Q%S(11)
        W(7)=-(Spq(17)*Q%S(12))-Spq(18)*Q%S(13)-Spq(19)*Q%S(14)-Spq(22)*Q%S(16)
        W(8)=-(Spq(23)*Q%S(17))-Spq(24)*Q%S(18)-Spq(25)*Q%S(19)-Spq(26)*Q%S(20)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)
        W(2)=(Cpq(3)-Cpq(5))*Q%C(2)-Cpq(7)*Q%C(3)
        W(3)=(Cpq(6)-Cpq(8))*Q%C(4)+(Cpq(7)-Cpq(9))*Q%C(5)
        W(4)=-(Cpq(11)*Q%C(6))+(Cpq(10)-Cpq(12))*Q%C(7)+(Cpq(11)-Cpq(13))*Q%C(8)
        W(5)=(Cpq(12)-Cpq(14))*Q%C(9)-Cpq(16)*Q%C(10)
        W(6)=(Cpq(15)-Cpq(17))*Q%C(11)+(Cpq(16)-Cpq(18))*Q%C(12)
        W(7)=(Cpq(17)-Cpq(19))*Q%C(13)+(Cpq(18)-Cpq(20))*Q%C(14)
        W(8)=-(Cpq(22)*Q%C(15))+(Cpq(21)-Cpq(23))*Q%C(16)+(Cpq(22)-Cpq(24))*Q%C(17)
        W(9)=(Cpq(23)-Cpq(25))*Q%C(18)+(Cpq(24)-Cpq(26))*Q%C(19)
        W(10)=(Cpq(25)-Cpq(27))*Q%C(20)-Spq(5)*Q%S(2)
        W(11)=-(Spq(8)*Q%S(4))+(Spq(7)-Spq(9))*Q%S(5)
        W(12)=-(Spq(12)*Q%S(7))+(Spq(11)-Spq(13))*Q%S(8)+(Spq(12)-Spq(14))*Q%S(9)
        W(13)=-(Spq(17)*Q%S(11))+(Spq(16)-Spq(18))*Q%S(12)
        W(14)=(Spq(17)-Spq(19))*Q%S(13)+(Spq(18)-Spq(20))*Q%S(14)
        W(15)=-(Spq(23)*Q%S(16))+(Spq(22)-Spq(24))*Q%S(17)
        W(16)=(Spq(23)-Spq(25))*Q%S(18)
        W(17)=(Spq(24)-Spq(26))*Q%S(19)+(Spq(25)-Spq(27))*Q%S(20)
        W(18)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(19)=W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)+W(16)+W(17)
        CTwo=W(18)+W(19)
        W(1)=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)-Q%C(3)*Spq(7)
        W(2)=-(Q%C(4)*Spq(8))-Q%C(5)*(Spq(7)+Spq(9))
        W(3)=-(Q%C(6)*Spq(11))-Q%C(7)*Spq(12)-Q%C(8)*(Spq(11)+Spq(13))
        W(4)=-(Q%C(9)*(Spq(12)+Spq(14)))-Q%C(10)*Spq(16)
        W(5)=-(Q%C(11)*Spq(17))-Q%C(12)*(Spq(16)+Spq(18))
        W(6)=-(Q%C(13)*(Spq(17)+Spq(19)))-Q%C(14)*(Spq(18)+Spq(20))
        W(7)=-(Q%C(15)*Spq(22))-Q%C(16)*Spq(23)-Q%C(17)*(Spq(22)+Spq(24))
        W(8)=-(Q%C(18)*(Spq(23)+Spq(25)))-Q%C(19)*(Spq(24)+Spq(26))
        W(9)=-(Q%C(20)*(Spq(25)+Spq(27)))+(Cpq(3)+Cpq(5))*Q%S(2)
        W(10)=(Cpq(6)+Cpq(8))*Q%S(4)+(Cpq(7)+Cpq(9))*Q%S(5)
        W(11)=(Cpq(10)+Cpq(12))*Q%S(7)+(Cpq(11)+Cpq(13))*Q%S(8)+(Cpq(12)+Cpq(14))*Q%S(9)
        W(12)=(Cpq(15)+Cpq(17))*Q%S(11)+(Cpq(16)+Cpq(18))*Q%S(12)
        W(13)=(Cpq(17)+Cpq(19))*Q%S(13)+(Cpq(18)+Cpq(20))*Q%S(14)
        W(14)=(Cpq(21)+Cpq(23))*Q%S(16)+(Cpq(22)+Cpq(24))*Q%S(17)
        W(15)=(Cpq(23)+Cpq(25))*Q%S(18)+(Cpq(24)+Cpq(26))*Q%S(19)+(Cpq(25)+Cpq(27))*Q%S(20)
        W(16)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        W(17)=W(8)+W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)
        STwo=W(16)+W(17)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
      CASE(106)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        W(1)=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)
        W(2)=Cpq(6)*Q%C(6)+Cpq(10)*Q%C(10)+Cpq(15)*Q%C(15)+Cpq(21)*Q%C(21)
        COne=W(1)+W(2)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)+Cpq(7)*Q%C(7)+Cpq(8)*Q%C(8)
        W(2)=Cpq(9)*Q%C(9)+Cpq(11)*Q%C(11)+Cpq(12)*Q%C(12)+Cpq(13)*Q%C(13)+Cpq(14)*Q%C(14)
        W(3)=Cpq(16)*Q%C(16)+Cpq(17)*Q%C(17)+Cpq(18)*Q%C(18)+Cpq(19)*Q%C(19)+Cpq(20)*Q%C(20)
        W(4)=Cpq(22)*Q%C(22)+Cpq(23)*Q%C(23)+Cpq(24)*Q%C(24)
        W(5)=Cpq(25)*Q%C(25)+Cpq(26)*Q%C(26)+Cpq(27)*Q%C(27)
        W(6)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)+Spq(7)*Q%S(7)+Spq(8)*Q%S(8)
        W(7)=Spq(9)*Q%S(9)+Spq(11)*Q%S(11)+Spq(12)*Q%S(12)+Spq(13)*Q%S(13)+Spq(14)*Q%S(14)
        W(8)=Spq(16)*Q%S(16)+Spq(17)*Q%S(17)+Spq(18)*Q%S(18)+Spq(19)*Q%S(19)+Spq(20)*Q%S(20)
        W(9)=Spq(22)*Q%S(22)+Spq(23)*Q%S(23)+Spq(24)*Q%S(24)
        W(10)=Spq(25)*Q%S(25)+Spq(26)*Q%S(26)+Spq(27)*Q%S(27)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        W(1)=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)-Cpq(6)*Q%C(3)
        W(2)=-(Cpq(10)*Q%C(6))-Cpq(15)*Q%C(10)-Cpq(21)*Q%C(15)-Cpq(28)*Q%C(21)
        COne=W(1)+W(2)
        W(1)=-(Cpq(4)*Q%C(2))-Cpq(7)*Q%C(4)
        W(2)=-(Cpq(8)*Q%C(5))-Cpq(11)*Q%C(7)-Cpq(12)*Q%C(8)
        W(3)=-(Cpq(13)*Q%C(9))-Cpq(16)*Q%C(11)
        W(4)=-(Cpq(17)*Q%C(12))-Cpq(18)*Q%C(13)-Cpq(19)*Q%C(14)
        W(5)=-(Cpq(22)*Q%C(16))-Cpq(23)*Q%C(17)
        W(6)=-(Cpq(24)*Q%C(18))-Cpq(25)*Q%C(19)-Cpq(26)*Q%C(20)
        W(7)=-(Cpq(29)*Q%C(22))-Cpq(30)*Q%C(23)-Cpq(31)*Q%C(24)
        W(8)=-(Cpq(32)*Q%C(25))-Cpq(33)*Q%C(26)-Cpq(34)*Q%C(27)
        W(9)=-(Spq(4)*Q%S(2))-Spq(7)*Q%S(4)
        W(10)=-(Spq(8)*Q%S(5))-Spq(11)*Q%S(7)-Spq(12)*Q%S(8)
        W(11)=-(Spq(13)*Q%S(9))-Spq(16)*Q%S(11)
        W(12)=-(Spq(17)*Q%S(12))-Spq(18)*Q%S(13)-Spq(19)*Q%S(14)
        W(13)=-(Spq(22)*Q%S(16))-Spq(23)*Q%S(17)
        W(14)=-(Spq(24)*Q%S(18))-Spq(25)*Q%S(19)-Spq(26)*Q%S(20)
        W(15)=-(Spq(29)*Q%S(22))-Spq(30)*Q%S(23)-Spq(31)*Q%S(24)
        W(16)=-(Spq(32)*Q%S(25))-Spq(33)*Q%S(26)-Spq(34)*Q%S(27)
        W(17)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(18)=W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)+W(16)
        CTwo=W(17)+W(18)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)+(Cpq(3)-Cpq(5))*Q%C(2)
        W(2)=-(Cpq(7)*Q%C(3))+(Cpq(6)-Cpq(8))*Q%C(4)+(Cpq(7)-Cpq(9))*Q%C(5)
        W(3)=-(Cpq(11)*Q%C(6))+(Cpq(10)-Cpq(12))*Q%C(7)+(Cpq(11)-Cpq(13))*Q%C(8)
        W(4)=(Cpq(12)-Cpq(14))*Q%C(9)-Cpq(16)*Q%C(10)+(Cpq(15)-Cpq(17))*Q%C(11)
        W(5)=(Cpq(16)-Cpq(18))*Q%C(12)
        W(6)=(Cpq(17)-Cpq(19))*Q%C(13)+(Cpq(18)-Cpq(20))*Q%C(14)
        W(7)=-(Cpq(22)*Q%C(15))+(Cpq(21)-Cpq(23))*Q%C(16)+(Cpq(22)-Cpq(24))*Q%C(17)
        W(8)=(Cpq(23)-Cpq(25))*Q%C(18)
        W(9)=(Cpq(24)-Cpq(26))*Q%C(19)+(Cpq(25)-Cpq(27))*Q%C(20)
        W(10)=-(Cpq(29)*Q%C(21))+(Cpq(28)-Cpq(30))*Q%C(22)+(Cpq(29)-Cpq(31))*Q%C(23)
        W(11)=(Cpq(30)-Cpq(32))*Q%C(24)
        W(12)=(Cpq(31)-Cpq(33))*Q%C(25)+(Cpq(32)-Cpq(34))*Q%C(26)
        W(13)=(Cpq(33)-Cpq(35))*Q%C(27)-Spq(5)*Q%S(2)-Spq(8)*Q%S(4)
        W(14)=(Spq(7)-Spq(9))*Q%S(5)-Spq(12)*Q%S(7)+(Spq(11)-Spq(13))*Q%S(8)
        W(15)=(Spq(12)-Spq(14))*Q%S(9)-Spq(17)*Q%S(11)+(Spq(16)-Spq(18))*Q%S(12)
        W(16)=(Spq(17)-Spq(19))*Q%S(13)+(Spq(18)-Spq(20))*Q%S(14)-Spq(23)*Q%S(16)
        W(17)=(Spq(22)-Spq(24))*Q%S(17)
        W(18)=(Spq(23)-Spq(25))*Q%S(18)+(Spq(24)-Spq(26))*Q%S(19)
        W(19)=(Spq(25)-Spq(27))*Q%S(20)-Spq(30)*Q%S(22)+(Spq(29)-Spq(31))*Q%S(23)
        W(20)=(Spq(30)-Spq(32))*Q%S(24)+(Spq(31)-Spq(33))*Q%S(25)
        W(21)=(Spq(32)-Spq(34))*Q%S(26)+(Spq(33)-Spq(35))*Q%S(27)
        W(22)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)
        W(23)=W(11)+W(12)+W(13)+W(14)+W(15)+W(16)+W(17)+W(18)+W(19)+W(20)+W(21)
        CTwo=W(22)+W(23)
        W(1)=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)
        W(2)=-(Q%C(3)*Spq(7))-Q%C(4)*Spq(8)-Q%C(5)*(Spq(7)+Spq(9))
        W(3)=-(Q%C(6)*Spq(11))-Q%C(7)*Spq(12)-Q%C(8)*(Spq(11)+Spq(13))
        W(4)=-(Q%C(9)*(Spq(12)+Spq(14)))-Q%C(10)*Spq(16)-Q%C(11)*Spq(17)
        W(5)=-(Q%C(12)*(Spq(16)+Spq(18)))
        W(6)=-(Q%C(13)*(Spq(17)+Spq(19)))-Q%C(14)*(Spq(18)+Spq(20))
        W(7)=-(Q%C(15)*Spq(22))-Q%C(16)*Spq(23)-Q%C(17)*(Spq(22)+Spq(24))
        W(8)=-(Q%C(18)*(Spq(23)+Spq(25)))
        W(9)=-(Q%C(19)*(Spq(24)+Spq(26)))-Q%C(20)*(Spq(25)+Spq(27))
        W(10)=-(Q%C(21)*Spq(29))-Q%C(22)*Spq(30)-Q%C(23)*(Spq(29)+Spq(31))
        W(11)=-(Q%C(24)*(Spq(30)+Spq(32)))
        W(12)=-(Q%C(25)*(Spq(31)+Spq(33)))-Q%C(26)*(Spq(32)+Spq(34))
        W(13)=-(Q%C(27)*(Spq(33)+Spq(35)))+(Cpq(3)+Cpq(5))*Q%S(2)+(Cpq(6)+Cpq(8))*Q%S(4)
        W(14)=(Cpq(7)+Cpq(9))*Q%S(5)+(Cpq(10)+Cpq(12))*Q%S(7)+(Cpq(11)+Cpq(13))*Q%S(8)
        W(15)=(Cpq(12)+Cpq(14))*Q%S(9)+(Cpq(15)+Cpq(17))*Q%S(11)+(Cpq(16)+Cpq(18))*Q%S(12)
        W(16)=(Cpq(17)+Cpq(19))*Q%S(13)+(Cpq(18)+Cpq(20))*Q%S(14)+(Cpq(21)+Cpq(23))*Q%S(16)
        W(17)=(Cpq(22)+Cpq(24))*Q%S(17)+(Cpq(23)+Cpq(25))*Q%S(18)+(Cpq(24)+Cpq(26))*Q%S(19)
        W(18)=(Cpq(25)+Cpq(27))*Q%S(20)+(Cpq(28)+Cpq(30))*Q%S(22)+(Cpq(29)+Cpq(31))*Q%S(23)
        W(19)=(Cpq(30)+Cpq(32))*Q%S(24)+(Cpq(31)+Cpq(33))*Q%S(25)
        W(20)=(Cpq(32)+Cpq(34))*Q%S(26)+(Cpq(33)+Cpq(35))*Q%S(27)
        W(21)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)
        W(22)=W(11)+W(12)+W(13)+W(14)+W(15)+W(16)+W(17)+W(18)+W(19)+W(20)
        STwo=W(21)+W(22)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
      CASE(107)
         CALL CTraX107(Q)
      CASE(200)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)
        SPKetC(0)=SPKetC(0)+COne
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))
        SPKetC(1)=SPKetC(1)+COne
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        CTwo=-(Cpq(2)*Q%C(0))
        STwo=-(Q%C(0)*Spq(2))
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)
        SPKetC(3)=SPKetC(3)+COne
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        CTwo=Cpq(4)*Q%C(0)
        STwo=Q%C(0)*Spq(4)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        CTwo=Cpq(5)*Q%C(0)
        STwo=Q%C(0)*Spq(5)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
      CASE(201)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)
        CTwo=Cpq(2)*Q%C(2)+Spq(2)*Q%S(2)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)
        CTwo=-(Cpq(4)*Q%C(2))-Spq(4)*Q%S(2)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)
        W(2)=Cpq(3)*Q%C(2)-Cpq(5)*Q%C(2)-Spq(5)*Q%S(2)
        CTwo=W(1)+W(2)
        STwo=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)+(Cpq(3)+Cpq(5))*Q%S(2)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)+Cpq(6)*Q%C(1)
        CTwo=Cpq(7)*Q%C(2)+Spq(7)*Q%S(2)
        SPKetC(3)=SPKetC(3)+COne+CTwo*Two
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        CTwo=Cpq(4)*Q%C(0)+Cpq(7)*Q%C(1)+(-Cpq(6)+Cpq(8))*Q%C(2)+Spq(8)*Q%S(2)
        W(1)=Q%C(0)*Spq(4)+Q%C(1)*Spq(7)
        W(2)=Q%C(2)*Spq(8)-Cpq(6)*Q%S(2)-Cpq(8)*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        W(1)=Cpq(5)*Q%C(0)+Cpq(8)*Q%C(1)
        W(2)=(-Cpq(7)+Cpq(9))*Q%C(2)+(Spq(7)+Spq(9))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(5)+Q%C(1)*Spq(8)
        W(2)=Q%C(2)*(-Spq(7)+Spq(9))-(Cpq(7)+Cpq(9))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
      CASE(202)
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
      CASE(203)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)+Cpq(6)*Q%C(6)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)
        W(2)=Cpq(7)*Q%C(7)+Cpq(8)*Q%C(8)+Cpq(9)*Q%C(9)
        W(3)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)
        W(4)=Spq(7)*Q%S(7)+Spq(8)*Q%S(8)+Spq(9)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)-Cpq(6)*Q%C(3)-Cpq(10)*Q%C(6)
        W(1)=-(Cpq(4)*Q%C(2))-Cpq(7)*Q%C(4)-Cpq(8)*Q%C(5)
        W(2)=-(Cpq(11)*Q%C(7))-Cpq(12)*Q%C(8)-Cpq(13)*Q%C(9)
        W(3)=-(Spq(4)*Q%S(2))-Spq(7)*Q%S(4)-Spq(8)*Q%S(5)
        W(4)=-(Spq(11)*Q%S(7))-Spq(12)*Q%S(8)-Spq(13)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)+Cpq(3)*Q%C(2)
        W(2)=-(Cpq(5)*Q%C(2))-Cpq(7)*Q%C(3)+Cpq(6)*Q%C(4)
        W(3)=-(Cpq(8)*Q%C(4))+Cpq(7)*Q%C(5)-Cpq(9)*Q%C(5)
        W(4)=-(Cpq(11)*Q%C(6))+Cpq(10)*Q%C(7)-Cpq(12)*Q%C(7)
        W(5)=Cpq(11)*Q%C(8)-Cpq(13)*Q%C(8)+Cpq(12)*Q%C(9)
        W(6)=-(Cpq(14)*Q%C(9))-Spq(5)*Q%S(2)-Spq(8)*Q%S(4)
        W(7)=Spq(7)*Q%S(5)-Spq(9)*Q%S(5)-Spq(12)*Q%S(7)
        W(8)=Spq(11)*Q%S(8)-Spq(13)*Q%S(8)+Spq(12)*Q%S(9)-Spq(14)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(1)=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)-Q%C(3)*Spq(7)
        W(2)=-(Q%C(4)*Spq(8))-Q%C(5)*(Spq(7)+Spq(9))
        W(3)=-(Q%C(6)*Spq(11))-Q%C(7)*Spq(12)
        W(4)=-(Q%C(8)*(Spq(11)+Spq(13)))-Q%C(9)*(Spq(12)+Spq(14))
        W(5)=(Cpq(3)+Cpq(5))*Q%S(2)+(Cpq(6)+Cpq(8))*Q%S(4)
        W(6)=(Cpq(7)+Cpq(9))*Q%S(5)+(Cpq(10)+Cpq(12))*Q%S(7)
        W(7)=(Cpq(11)+Cpq(13))*Q%S(8)+(Cpq(12)+Cpq(14))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)+Cpq(6)*Q%C(1)+Cpq(10)*Q%C(3)+Cpq(15)*Q%C(6)
        W(1)=Cpq(7)*Q%C(2)+Cpq(11)*Q%C(4)+Cpq(12)*Q%C(5)
        W(2)=Cpq(16)*Q%C(7)+Cpq(17)*Q%C(8)+Cpq(18)*Q%C(9)
        W(3)=Spq(7)*Q%S(2)+Spq(11)*Q%S(4)+Spq(12)*Q%S(5)
        W(4)=Spq(16)*Q%S(7)+Spq(17)*Q%S(8)+Spq(18)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(3)=SPKetC(3)+COne+CTwo*Two
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        W(1)=Cpq(4)*Q%C(0)+Cpq(7)*Q%C(1)+(-Cpq(6)+Cpq(8))*Q%C(2)+Cpq(11)*Q%C(3)
        W(2)=(-Cpq(10)+Cpq(12))*Q%C(4)+(-Cpq(11)+Cpq(13))*Q%C(5)
        W(3)=Cpq(16)*Q%C(6)+(-Cpq(15)+Cpq(17))*Q%C(7)
        W(4)=(-Cpq(16)+Cpq(18))*Q%C(8)+(-Cpq(17)+Cpq(19))*Q%C(9)
        W(5)=Spq(8)*Q%S(2)+Spq(12)*Q%S(4)
        W(6)=(-Spq(11)+Spq(13))*Q%S(5)+Spq(17)*Q%S(7)
        W(7)=(-Spq(16)+Spq(18))*Q%S(8)+(-Spq(17)+Spq(19))*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        W(1)=Q%C(0)*Spq(4)+Q%C(1)*Spq(7)+Q%C(2)*Spq(8)+Q%C(3)*Spq(11)
        W(2)=Q%C(4)*Spq(12)+Q%C(5)*(Spq(11)+Spq(13))+Q%C(6)*Spq(16)+Q%C(7)*Spq(17)
        W(3)=Q%C(8)*(Spq(16)+Spq(18))+Q%C(9)*(Spq(17)+Spq(19))
        W(4)=-((Cpq(6)+Cpq(8))*Q%S(2))-(Cpq(10)+Cpq(12))*Q%S(4)
        W(5)=-((Cpq(11)+Cpq(13))*Q%S(5))-(Cpq(15)+Cpq(17))*Q%S(7)
        W(6)=-((Cpq(16)+Cpq(18))*Q%S(8))-(Cpq(17)+Cpq(19))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        W(1)=Cpq(5)*Q%C(0)+Cpq(8)*Q%C(1)+(-Cpq(7)+Cpq(9))*Q%C(2)+Cpq(12)*Q%C(3)
        W(2)=(-Cpq(11)+Cpq(13))*Q%C(4)+(Cpq(10)+Cpq(14))*Q%C(5)
        W(3)=Cpq(17)*Q%C(6)+(-Cpq(16)+Cpq(18))*Q%C(7)
        W(4)=(Cpq(15)+Cpq(19))*Q%C(8)+(Cpq(16)+Cpq(20))*Q%C(9)
        W(5)=(Spq(7)+Spq(9))*Q%S(2)+(Spq(11)+Spq(13))*Q%S(4)
        W(6)=Spq(14)*Q%S(5)+(Spq(16)+Spq(18))*Q%S(7)+Spq(19)*Q%S(8)+(Spq(16)+Spq(20))*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)
        W(1)=Q%C(0)*Spq(5)+Q%C(1)*Spq(8)+Q%C(2)*(-Spq(7)+Spq(9))+Q%C(3)*Spq(12)
        W(2)=Q%C(4)*(-Spq(11)+Spq(13))+Q%C(5)*Spq(14)
        W(3)=Q%C(6)*Spq(17)+Q%C(7)*(-Spq(16)+Spq(18))
        W(4)=Q%C(8)*Spq(19)+Q%C(9)*(-Spq(16)+Spq(20))
        W(5)=-((Cpq(7)+Cpq(9))*Q%S(2))-(Cpq(11)+Cpq(13))*Q%S(4)
        W(6)=(Cpq(10)-Cpq(14))*Q%S(5)-(Cpq(16)+Cpq(18))*Q%S(7)
        W(7)=(Cpq(15)-Cpq(19))*Q%S(8)+(Cpq(16)-Cpq(20))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
      CASE(204)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)+Cpq(6)*Q%C(6)+Cpq(10)*Q%C(10)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)+Cpq(7)*Q%C(7)+Cpq(8)*Q%C(8)
        W(2)=Cpq(9)*Q%C(9)+Cpq(11)*Q%C(11)+Cpq(12)*Q%C(12)+Cpq(13)*Q%C(13)+Cpq(14)*Q%C(14)
        W(3)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)+Spq(7)*Q%S(7)+Spq(8)*Q%S(8)
        W(4)=Spq(9)*Q%S(9)+Spq(11)*Q%S(11)+Spq(12)*Q%S(12)+Spq(13)*Q%S(13)+Spq(14)*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        W(1)=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)
        W(2)=-(Cpq(6)*Q%C(3))-Cpq(10)*Q%C(6)-Cpq(15)*Q%C(10)
        COne=W(1)+W(2)
        W(1)=-(Cpq(4)*Q%C(2))-Cpq(7)*Q%C(4)
        W(2)=-(Cpq(8)*Q%C(5))-Cpq(11)*Q%C(7)-Cpq(12)*Q%C(8)
        W(3)=-(Cpq(13)*Q%C(9))-Cpq(16)*Q%C(11)
        W(4)=-(Cpq(17)*Q%C(12))-Cpq(18)*Q%C(13)-Cpq(19)*Q%C(14)
        W(5)=-(Spq(4)*Q%S(2))-Spq(7)*Q%S(4)
        W(6)=-(Spq(8)*Q%S(5))-Spq(11)*Q%S(7)-Spq(12)*Q%S(8)
        W(7)=-(Spq(13)*Q%S(9))-Spq(16)*Q%S(11)
        W(8)=-(Spq(17)*Q%S(12))-Spq(18)*Q%S(13)-Spq(19)*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)+(Cpq(3)-Cpq(5))*Q%C(2)
        W(2)=-(Cpq(7)*Q%C(3))+(Cpq(6)-Cpq(8))*Q%C(4)+(Cpq(7)-Cpq(9))*Q%C(5)
        W(3)=-(Cpq(11)*Q%C(6))+(Cpq(10)-Cpq(12))*Q%C(7)+(Cpq(11)-Cpq(13))*Q%C(8)
        W(4)=(Cpq(12)-Cpq(14))*Q%C(9)-Cpq(16)*Q%C(10)+(Cpq(15)-Cpq(17))*Q%C(11)
        W(5)=(Cpq(16)-Cpq(18))*Q%C(12)
        W(6)=(Cpq(17)-Cpq(19))*Q%C(13)+(Cpq(18)-Cpq(20))*Q%C(14)
        W(7)=-(Spq(5)*Q%S(2))-Spq(8)*Q%S(4)+(Spq(7)-Spq(9))*Q%S(5)
        W(8)=-(Spq(12)*Q%S(7))+(Spq(11)-Spq(13))*Q%S(8)+(Spq(12)-Spq(14))*Q%S(9)
        W(9)=-(Spq(17)*Q%S(11))+(Spq(16)-Spq(18))*Q%S(12)
        W(10)=(Spq(17)-Spq(19))*Q%S(13)+(Spq(18)-Spq(20))*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)
        W(1)=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)
        W(2)=-(Q%C(3)*Spq(7))-Q%C(4)*Spq(8)-Q%C(5)*(Spq(7)+Spq(9))
        W(3)=-(Q%C(6)*Spq(11))-Q%C(7)*Spq(12)-Q%C(8)*(Spq(11)+Spq(13))
        W(4)=-(Q%C(9)*(Spq(12)+Spq(14)))-Q%C(10)*Spq(16)-Q%C(11)*Spq(17)
        W(5)=-(Q%C(12)*(Spq(16)+Spq(18)))
        W(6)=-(Q%C(13)*(Spq(17)+Spq(19)))-Q%C(14)*(Spq(18)+Spq(20))
        W(7)=(Cpq(3)+Cpq(5))*Q%S(2)+(Cpq(6)+Cpq(8))*Q%S(4)+(Cpq(7)+Cpq(9))*Q%S(5)
        W(8)=(Cpq(10)+Cpq(12))*Q%S(7)+(Cpq(11)+Cpq(13))*Q%S(8)+(Cpq(12)+Cpq(14))*Q%S(9)
        W(9)=(Cpq(15)+Cpq(17))*Q%S(11)+(Cpq(16)+Cpq(18))*Q%S(12)
        W(10)=(Cpq(17)+Cpq(19))*Q%S(13)+(Cpq(18)+Cpq(20))*Q%S(14)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)+Cpq(6)*Q%C(1)+Cpq(10)*Q%C(3)+Cpq(15)*Q%C(6)+Cpq(21)*Q%C(10)
        W(1)=Cpq(7)*Q%C(2)+Cpq(11)*Q%C(4)+Cpq(12)*Q%C(5)+Cpq(16)*Q%C(7)+Cpq(17)*Q%C(8)
        W(2)=Cpq(18)*Q%C(9)+Cpq(22)*Q%C(11)+Cpq(23)*Q%C(12)+Cpq(24)*Q%C(13)+Cpq(25)*Q%C(14)
        W(3)=Spq(7)*Q%S(2)+Spq(11)*Q%S(4)+Spq(12)*Q%S(5)+Spq(16)*Q%S(7)+Spq(17)*Q%S(8)
        W(4)=Spq(18)*Q%S(9)+Spq(22)*Q%S(11)+Spq(23)*Q%S(12)+Spq(24)*Q%S(13)+Spq(25)*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(3)=SPKetC(3)+COne+CTwo*Two
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        W(1)=Cpq(4)*Q%C(0)+Cpq(7)*Q%C(1)+(-Cpq(6)+Cpq(8))*Q%C(2)
        W(2)=Cpq(11)*Q%C(3)+(-Cpq(10)+Cpq(12))*Q%C(4)+(-Cpq(11)+Cpq(13))*Q%C(5)
        W(3)=Cpq(16)*Q%C(6)+(-Cpq(15)+Cpq(17))*Q%C(7)+(-Cpq(16)+Cpq(18))*Q%C(8)
        W(4)=(-Cpq(17)+Cpq(19))*Q%C(9)+Cpq(22)*Q%C(10)+(-Cpq(21)+Cpq(23))*Q%C(11)
        W(5)=(-Cpq(22)+Cpq(24))*Q%C(12)
        W(6)=(-Cpq(23)+Cpq(25))*Q%C(13)+(-Cpq(24)+Cpq(26))*Q%C(14)
        W(7)=Spq(8)*Q%S(2)+Spq(12)*Q%S(4)+(-Spq(11)+Spq(13))*Q%S(5)
        W(8)=Spq(17)*Q%S(7)+(-Spq(16)+Spq(18))*Q%S(8)+(-Spq(17)+Spq(19))*Q%S(9)
        W(9)=Spq(23)*Q%S(11)+(-Spq(22)+Spq(24))*Q%S(12)
        W(10)=(-Spq(23)+Spq(25))*Q%S(13)+(-Spq(24)+Spq(26))*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)
        W(1)=Q%C(0)*Spq(4)+Q%C(1)*Spq(7)+Q%C(2)*Spq(8)
        W(2)=Q%C(3)*Spq(11)+Q%C(4)*Spq(12)+Q%C(5)*(Spq(11)+Spq(13))
        W(3)=Q%C(6)*Spq(16)+Q%C(7)*Spq(17)+Q%C(8)*(Spq(16)+Spq(18))
        W(4)=Q%C(9)*(Spq(17)+Spq(19))+Q%C(10)*Spq(22)+Q%C(11)*Spq(23)
        W(5)=Q%C(12)*(Spq(22)+Spq(24))+Q%C(13)*(Spq(23)+Spq(25))+Q%C(14)*(Spq(24)+Spq(26))
        W(6)=-((Cpq(6)+Cpq(8))*Q%S(2))
        W(7)=-((Cpq(10)+Cpq(12))*Q%S(4))-(Cpq(11)+Cpq(13))*Q%S(5)
        W(8)=-((Cpq(15)+Cpq(17))*Q%S(7))
        W(9)=-((Cpq(16)+Cpq(18))*Q%S(8))-(Cpq(17)+Cpq(19))*Q%S(9)
        W(10)=-((Cpq(21)+Cpq(23))*Q%S(11))-(Cpq(22)+Cpq(24))*Q%S(12)
        W(11)=-((Cpq(23)+Cpq(25))*Q%S(13))-(Cpq(24)+Cpq(26))*Q%S(14)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)+W(11)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        W(1)=Cpq(5)*Q%C(0)+Cpq(8)*Q%C(1)+(-Cpq(7)+Cpq(9))*Q%C(2)
        W(2)=Cpq(12)*Q%C(3)+(-Cpq(11)+Cpq(13))*Q%C(4)+(Cpq(10)+Cpq(14))*Q%C(5)
        W(3)=Cpq(17)*Q%C(6)+(-Cpq(16)+Cpq(18))*Q%C(7)+(Cpq(15)+Cpq(19))*Q%C(8)
        W(4)=(Cpq(16)+Cpq(20))*Q%C(9)+Cpq(23)*Q%C(10)+(-Cpq(22)+Cpq(24))*Q%C(11)
        W(5)=(Cpq(21)+Cpq(25))*Q%C(12)+(Cpq(22)+Cpq(26))*Q%C(13)+(Cpq(23)+Cpq(27))*Q%C(14)
        W(6)=(Spq(7)+Spq(9))*Q%S(2)+(Spq(11)+Spq(13))*Q%S(4)+Spq(14)*Q%S(5)
        W(7)=(Spq(16)+Spq(18))*Q%S(7)+Spq(19)*Q%S(8)+(Spq(16)+Spq(20))*Q%S(9)
        W(8)=(Spq(22)+Spq(24))*Q%S(11)+Spq(25)*Q%S(12)
        W(9)=(Spq(22)+Spq(26))*Q%S(13)+(Spq(23)+Spq(27))*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)
        W(1)=Q%C(0)*Spq(5)+Q%C(1)*Spq(8)+Q%C(2)*(-Spq(7)+Spq(9))
        W(2)=Q%C(3)*Spq(12)+Q%C(4)*(-Spq(11)+Spq(13))+Q%C(5)*Spq(14)
        W(3)=Q%C(6)*Spq(17)+Q%C(7)*(-Spq(16)+Spq(18))+Q%C(8)*Spq(19)
        W(4)=Q%C(9)*(-Spq(16)+Spq(20))+Q%C(10)*Spq(23)+Q%C(11)*(-Spq(22)+Spq(24))
        W(5)=Q%C(12)*Spq(25)+Q%C(13)*(-Spq(22)+Spq(26))+Q%C(14)*(-Spq(23)+Spq(27))
        W(6)=-((Cpq(7)+Cpq(9))*Q%S(2))
        W(7)=-((Cpq(11)+Cpq(13))*Q%S(4))+(Cpq(10)-Cpq(14))*Q%S(5)
        W(8)=-((Cpq(16)+Cpq(18))*Q%S(7))
        W(9)=(Cpq(15)-Cpq(19))*Q%S(8)+(Cpq(16)-Cpq(20))*Q%S(9)
        W(10)=-((Cpq(22)+Cpq(24))*Q%S(11))+(Cpq(21)-Cpq(25))*Q%S(12)
        W(11)=(Cpq(22)-Cpq(26))*Q%S(13)+(Cpq(23)-Cpq(27))*Q%S(14)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)+W(11)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
      CASE(205)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        W(1)=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)
        W(2)=Cpq(6)*Q%C(6)+Cpq(10)*Q%C(10)+Cpq(15)*Q%C(15)
        COne=W(1)+W(2)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)
        W(2)=Cpq(7)*Q%C(7)+Cpq(8)*Q%C(8)+Cpq(9)*Q%C(9)+Cpq(11)*Q%C(11)
        W(3)=Cpq(12)*Q%C(12)+Cpq(13)*Q%C(13)+Cpq(14)*Q%C(14)+Cpq(16)*Q%C(16)
        W(4)=Cpq(17)*Q%C(17)+Cpq(18)*Q%C(18)+Cpq(19)*Q%C(19)+Cpq(20)*Q%C(20)
        W(5)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)
        W(6)=Spq(7)*Q%S(7)+Spq(8)*Q%S(8)+Spq(9)*Q%S(9)+Spq(11)*Q%S(11)
        W(7)=Spq(12)*Q%S(12)+Spq(13)*Q%S(13)+Spq(14)*Q%S(14)+Spq(16)*Q%S(16)
        W(8)=Spq(17)*Q%S(17)+Spq(18)*Q%S(18)+Spq(19)*Q%S(19)+Spq(20)*Q%S(20)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        W(1)=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)-Cpq(6)*Q%C(3)
        W(2)=-(Cpq(10)*Q%C(6))-Cpq(15)*Q%C(10)-Cpq(21)*Q%C(15)
        COne=W(1)+W(2)
        W(1)=-(Cpq(4)*Q%C(2))-Cpq(7)*Q%C(4)-Cpq(8)*Q%C(5)
        W(2)=-(Cpq(11)*Q%C(7))-Cpq(12)*Q%C(8)-Cpq(13)*Q%C(9)-Cpq(16)*Q%C(11)
        W(3)=-(Cpq(17)*Q%C(12))-Cpq(18)*Q%C(13)-Cpq(19)*Q%C(14)-Cpq(22)*Q%C(16)
        W(4)=-(Cpq(23)*Q%C(17))-Cpq(24)*Q%C(18)-Cpq(25)*Q%C(19)-Cpq(26)*Q%C(20)
        W(5)=-(Spq(4)*Q%S(2))-Spq(7)*Q%S(4)-Spq(8)*Q%S(5)
        W(6)=-(Spq(11)*Q%S(7))-Spq(12)*Q%S(8)-Spq(13)*Q%S(9)-Spq(16)*Q%S(11)
        W(7)=-(Spq(17)*Q%S(12))-Spq(18)*Q%S(13)-Spq(19)*Q%S(14)-Spq(22)*Q%S(16)
        W(8)=-(Spq(23)*Q%S(17))-Spq(24)*Q%S(18)-Spq(25)*Q%S(19)-Spq(26)*Q%S(20)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)
        W(2)=(Cpq(3)-Cpq(5))*Q%C(2)-Cpq(7)*Q%C(3)
        W(3)=(Cpq(6)-Cpq(8))*Q%C(4)+(Cpq(7)-Cpq(9))*Q%C(5)
        W(4)=-(Cpq(11)*Q%C(6))+(Cpq(10)-Cpq(12))*Q%C(7)+(Cpq(11)-Cpq(13))*Q%C(8)
        W(5)=(Cpq(12)-Cpq(14))*Q%C(9)-Cpq(16)*Q%C(10)
        W(6)=(Cpq(15)-Cpq(17))*Q%C(11)+(Cpq(16)-Cpq(18))*Q%C(12)
        W(7)=(Cpq(17)-Cpq(19))*Q%C(13)+(Cpq(18)-Cpq(20))*Q%C(14)
        W(8)=-(Cpq(22)*Q%C(15))+(Cpq(21)-Cpq(23))*Q%C(16)+(Cpq(22)-Cpq(24))*Q%C(17)
        W(9)=(Cpq(23)-Cpq(25))*Q%C(18)+(Cpq(24)-Cpq(26))*Q%C(19)
        W(10)=(Cpq(25)-Cpq(27))*Q%C(20)-Spq(5)*Q%S(2)
        W(11)=-(Spq(8)*Q%S(4))+(Spq(7)-Spq(9))*Q%S(5)
        W(12)=-(Spq(12)*Q%S(7))+(Spq(11)-Spq(13))*Q%S(8)+(Spq(12)-Spq(14))*Q%S(9)
        W(13)=-(Spq(17)*Q%S(11))+(Spq(16)-Spq(18))*Q%S(12)
        W(14)=(Spq(17)-Spq(19))*Q%S(13)+(Spq(18)-Spq(20))*Q%S(14)
        W(15)=-(Spq(23)*Q%S(16))+(Spq(22)-Spq(24))*Q%S(17)
        W(16)=(Spq(23)-Spq(25))*Q%S(18)
        W(17)=(Spq(24)-Spq(26))*Q%S(19)+(Spq(25)-Spq(27))*Q%S(20)
        W(18)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(19)=W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)+W(16)+W(17)
        CTwo=W(18)+W(19)
        W(1)=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)-Q%C(3)*Spq(7)
        W(2)=-(Q%C(4)*Spq(8))-Q%C(5)*(Spq(7)+Spq(9))
        W(3)=-(Q%C(6)*Spq(11))-Q%C(7)*Spq(12)-Q%C(8)*(Spq(11)+Spq(13))
        W(4)=-(Q%C(9)*(Spq(12)+Spq(14)))-Q%C(10)*Spq(16)
        W(5)=-(Q%C(11)*Spq(17))-Q%C(12)*(Spq(16)+Spq(18))
        W(6)=-(Q%C(13)*(Spq(17)+Spq(19)))-Q%C(14)*(Spq(18)+Spq(20))
        W(7)=-(Q%C(15)*Spq(22))-Q%C(16)*Spq(23)-Q%C(17)*(Spq(22)+Spq(24))
        W(8)=-(Q%C(18)*(Spq(23)+Spq(25)))-Q%C(19)*(Spq(24)+Spq(26))
        W(9)=-(Q%C(20)*(Spq(25)+Spq(27)))+(Cpq(3)+Cpq(5))*Q%S(2)
        W(10)=(Cpq(6)+Cpq(8))*Q%S(4)+(Cpq(7)+Cpq(9))*Q%S(5)
        W(11)=(Cpq(10)+Cpq(12))*Q%S(7)+(Cpq(11)+Cpq(13))*Q%S(8)+(Cpq(12)+Cpq(14))*Q%S(9)
        W(12)=(Cpq(15)+Cpq(17))*Q%S(11)+(Cpq(16)+Cpq(18))*Q%S(12)
        W(13)=(Cpq(17)+Cpq(19))*Q%S(13)+(Cpq(18)+Cpq(20))*Q%S(14)
        W(14)=(Cpq(21)+Cpq(23))*Q%S(16)+(Cpq(22)+Cpq(24))*Q%S(17)
        W(15)=(Cpq(23)+Cpq(25))*Q%S(18)+(Cpq(24)+Cpq(26))*Q%S(19)+(Cpq(25)+Cpq(27))*Q%S(20)
        W(16)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        W(17)=W(8)+W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)
        STwo=W(16)+W(17)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        W(1)=Cpq(3)*Q%C(0)+Cpq(6)*Q%C(1)+Cpq(10)*Q%C(3)
        W(2)=Cpq(15)*Q%C(6)+Cpq(21)*Q%C(10)+Cpq(28)*Q%C(15)
        COne=W(1)+W(2)
        W(1)=Cpq(7)*Q%C(2)+Cpq(11)*Q%C(4)+Cpq(12)*Q%C(5)
        W(2)=Cpq(16)*Q%C(7)+Cpq(17)*Q%C(8)+Cpq(18)*Q%C(9)+Cpq(22)*Q%C(11)
        W(3)=Cpq(23)*Q%C(12)+Cpq(24)*Q%C(13)+Cpq(25)*Q%C(14)+Cpq(29)*Q%C(16)
        W(4)=Cpq(30)*Q%C(17)+Cpq(31)*Q%C(18)+Cpq(32)*Q%C(19)+Cpq(33)*Q%C(20)
        W(5)=Spq(7)*Q%S(2)+Spq(11)*Q%S(4)+Spq(12)*Q%S(5)
        W(6)=Spq(16)*Q%S(7)+Spq(17)*Q%S(8)+Spq(18)*Q%S(9)+Spq(22)*Q%S(11)
        W(7)=Spq(23)*Q%S(12)+Spq(24)*Q%S(13)+Spq(25)*Q%S(14)+Spq(29)*Q%S(16)
        W(8)=Spq(30)*Q%S(17)+Spq(31)*Q%S(18)+Spq(32)*Q%S(19)+Spq(33)*Q%S(20)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        SPKetC(3)=SPKetC(3)+COne+CTwo*Two
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        W(1)=Cpq(4)*Q%C(0)+Cpq(7)*Q%C(1)+(-Cpq(6)+Cpq(8))*Q%C(2)+Cpq(11)*Q%C(3)
        W(2)=(-Cpq(10)+Cpq(12))*Q%C(4)+(-Cpq(11)+Cpq(13))*Q%C(5)
        W(3)=Cpq(16)*Q%C(6)+(-Cpq(15)+Cpq(17))*Q%C(7)+(-Cpq(16)+Cpq(18))*Q%C(8)
        W(4)=(-Cpq(17)+Cpq(19))*Q%C(9)+Cpq(22)*Q%C(10)
        W(5)=(-Cpq(21)+Cpq(23))*Q%C(11)+(-Cpq(22)+Cpq(24))*Q%C(12)
        W(6)=(-Cpq(23)+Cpq(25))*Q%C(13)+(-Cpq(24)+Cpq(26))*Q%C(14)
        W(7)=Cpq(29)*Q%C(15)+(-Cpq(28)+Cpq(30))*Q%C(16)+(-Cpq(29)+Cpq(31))*Q%C(17)
        W(8)=(-Cpq(30)+Cpq(32))*Q%C(18)+(-Cpq(31)+Cpq(33))*Q%C(19)
        W(9)=(-Cpq(32)+Cpq(34))*Q%C(20)+Spq(8)*Q%S(2)
        W(10)=Spq(12)*Q%S(4)+(-Spq(11)+Spq(13))*Q%S(5)
        W(11)=Spq(17)*Q%S(7)+(-Spq(16)+Spq(18))*Q%S(8)+(-Spq(17)+Spq(19))*Q%S(9)
        W(12)=Spq(23)*Q%S(11)+(-Spq(22)+Spq(24))*Q%S(12)
        W(13)=(-Spq(23)+Spq(25))*Q%S(13)+(-Spq(24)+Spq(26))*Q%S(14)
        W(14)=Spq(30)*Q%S(16)+(-Spq(29)+Spq(31))*Q%S(17)
        W(15)=(-Spq(30)+Spq(32))*Q%S(18)
        W(16)=(-Spq(31)+Spq(33))*Q%S(19)+(-Spq(32)+Spq(34))*Q%S(20)
        W(17)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(18)=W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)+W(16)
        CTwo=W(17)+W(18)
        W(1)=Q%C(0)*Spq(4)+Q%C(1)*Spq(7)+Q%C(2)*Spq(8)+Q%C(3)*Spq(11)
        W(2)=Q%C(4)*Spq(12)+Q%C(5)*(Spq(11)+Spq(13))
        W(3)=Q%C(6)*Spq(16)+Q%C(7)*Spq(17)+Q%C(8)*(Spq(16)+Spq(18))
        W(4)=Q%C(9)*(Spq(17)+Spq(19))+Q%C(10)*Spq(22)+Q%C(11)*Spq(23)+Q%C(12)*(Spq(22)+Spq(24))
        W(5)=Q%C(13)*(Spq(23)+Spq(25))+Q%C(14)*(Spq(24)+Spq(26))
        W(6)=Q%C(15)*Spq(29)+Q%C(16)*Spq(30)+Q%C(17)*(Spq(29)+Spq(31))
        W(7)=Q%C(18)*(Spq(30)+Spq(32))+Q%C(19)*(Spq(31)+Spq(33))
        W(8)=Q%C(20)*(Spq(32)+Spq(34))-(Cpq(6)+Cpq(8))*Q%S(2)
        W(9)=-((Cpq(10)+Cpq(12))*Q%S(4))-(Cpq(11)+Cpq(13))*Q%S(5)
        W(10)=-((Cpq(15)+Cpq(17))*Q%S(7))
        W(11)=-((Cpq(16)+Cpq(18))*Q%S(8))-(Cpq(17)+Cpq(19))*Q%S(9)
        W(12)=-((Cpq(21)+Cpq(23))*Q%S(11))-(Cpq(22)+Cpq(24))*Q%S(12)
        W(13)=-((Cpq(23)+Cpq(25))*Q%S(13))-(Cpq(24)+Cpq(26))*Q%S(14)
        W(14)=-((Cpq(28)+Cpq(30))*Q%S(16))-(Cpq(29)+Cpq(31))*Q%S(17)
        W(15)=-((Cpq(30)+Cpq(32))*Q%S(18))
        W(16)=-((Cpq(31)+Cpq(33))*Q%S(19))-(Cpq(32)+Cpq(34))*Q%S(20)
        W(17)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(18)=W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)+W(16)
        STwo=W(17)+W(18)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        W(1)=Cpq(5)*Q%C(0)+Cpq(8)*Q%C(1)+(-Cpq(7)+Cpq(9))*Q%C(2)+Cpq(12)*Q%C(3)
        W(2)=(-Cpq(11)+Cpq(13))*Q%C(4)+(Cpq(10)+Cpq(14))*Q%C(5)
        W(3)=Cpq(17)*Q%C(6)+(-Cpq(16)+Cpq(18))*Q%C(7)+(Cpq(15)+Cpq(19))*Q%C(8)
        W(4)=(Cpq(16)+Cpq(20))*Q%C(9)+Cpq(23)*Q%C(10)
        W(5)=(-Cpq(22)+Cpq(24))*Q%C(11)+(Cpq(21)+Cpq(25))*Q%C(12)
        W(6)=(Cpq(22)+Cpq(26))*Q%C(13)+(Cpq(23)+Cpq(27))*Q%C(14)
        W(7)=Cpq(30)*Q%C(15)+(-Cpq(29)+Cpq(31))*Q%C(16)+(Cpq(28)+Cpq(32))*Q%C(17)
        W(8)=(Cpq(29)+Cpq(33))*Q%C(18)+(Cpq(30)+Cpq(34))*Q%C(19)
        W(9)=(Cpq(31)+Cpq(35))*Q%C(20)+(Spq(7)+Spq(9))*Q%S(2)
        W(10)=(Spq(11)+Spq(13))*Q%S(4)+Spq(14)*Q%S(5)
        W(11)=(Spq(16)+Spq(18))*Q%S(7)+Spq(19)*Q%S(8)+(Spq(16)+Spq(20))*Q%S(9)
        W(12)=(Spq(22)+Spq(24))*Q%S(11)+Spq(25)*Q%S(12)
        W(13)=(Spq(22)+Spq(26))*Q%S(13)+(Spq(23)+Spq(27))*Q%S(14)
        W(14)=(Spq(29)+Spq(31))*Q%S(16)+Spq(32)*Q%S(17)
        W(15)=(Spq(29)+Spq(33))*Q%S(18)+(Spq(30)+Spq(34))*Q%S(19)+(Spq(31)+Spq(35))*Q%S(20)
        W(16)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        W(17)=W(8)+W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)
        CTwo=W(16)+W(17)
        W(1)=Q%C(0)*Spq(5)+Q%C(1)*Spq(8)+Q%C(2)*(-Spq(7)+Spq(9))+Q%C(3)*Spq(12)
        W(2)=Q%C(4)*(-Spq(11)+Spq(13))+Q%C(5)*Spq(14)
        W(3)=Q%C(6)*Spq(17)+Q%C(7)*(-Spq(16)+Spq(18))+Q%C(8)*Spq(19)
        W(4)=Q%C(9)*(-Spq(16)+Spq(20))+Q%C(10)*Spq(23)
        W(5)=Q%C(11)*(-Spq(22)+Spq(24))+Q%C(12)*Spq(25)
        W(6)=Q%C(13)*(-Spq(22)+Spq(26))+Q%C(14)*(-Spq(23)+Spq(27))
        W(7)=Q%C(15)*Spq(30)+Q%C(16)*(-Spq(29)+Spq(31))+Q%C(17)*Spq(32)
        W(8)=Q%C(18)*(-Spq(29)+Spq(33))+Q%C(19)*(-Spq(30)+Spq(34))
        W(9)=Q%C(20)*(-Spq(31)+Spq(35))-(Cpq(7)+Cpq(9))*Q%S(2)
        W(10)=-((Cpq(11)+Cpq(13))*Q%S(4))+(Cpq(10)-Cpq(14))*Q%S(5)
        W(11)=-((Cpq(16)+Cpq(18))*Q%S(7))
        W(12)=(Cpq(15)-Cpq(19))*Q%S(8)+(Cpq(16)-Cpq(20))*Q%S(9)
        W(13)=-((Cpq(22)+Cpq(24))*Q%S(11))+(Cpq(21)-Cpq(25))*Q%S(12)
        W(14)=(Cpq(22)-Cpq(26))*Q%S(13)+(Cpq(23)-Cpq(27))*Q%S(14)
        W(15)=-((Cpq(29)+Cpq(31))*Q%S(16))+(Cpq(28)-Cpq(32))*Q%S(17)
        W(16)=(Cpq(29)-Cpq(33))*Q%S(18)
        W(17)=(Cpq(30)-Cpq(34))*Q%S(19)+(Cpq(31)-Cpq(35))*Q%S(20)
        W(18)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(19)=W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)+W(16)+W(17)
        STwo=W(18)+W(19)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
      CASE(206)
         CALL CTraX206(Q)
      CASE(207)
         CALL CTraX207(Q)
      CASE(300)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)
        SPKetC(0)=SPKetC(0)+COne
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))
        SPKetC(1)=SPKetC(1)+COne
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        CTwo=-(Cpq(2)*Q%C(0))
        STwo=-(Q%C(0)*Spq(2))
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)
        SPKetC(3)=SPKetC(3)+COne
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        CTwo=Cpq(4)*Q%C(0)
        STwo=Q%C(0)*Spq(4)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        CTwo=Cpq(5)*Q%C(0)
        STwo=Q%C(0)*Spq(5)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
!       -------------------------------- Lp=3, Mp = 0-------------------------------------
        COne=-(Cpq(6)*Q%C(0))
        SPKetC(6)=SPKetC(6)+COne
!       -------------------------------- Lp=3, Mp = 1-------------------------------------
        CTwo=-(Cpq(7)*Q%C(0))
        STwo=-(Q%C(0)*Spq(7))
        SPKetC(7)=SPKetC(7)+CTwo*Two
        SPKetS(7)=SPKetS(7)+STwo*Two
!       -------------------------------- Lp=3, Mp = 2-------------------------------------
        CTwo=-(Cpq(8)*Q%C(0))
        STwo=-(Q%C(0)*Spq(8))
        SPKetC(8)=SPKetC(8)+CTwo*Two
        SPKetS(8)=SPKetS(8)+STwo*Two
!       -------------------------------- Lp=3, Mp = 3-------------------------------------
        CTwo=-(Cpq(9)*Q%C(0))
        STwo=-(Q%C(0)*Spq(9))
        SPKetC(9)=SPKetC(9)+CTwo*Two
        SPKetS(9)=SPKetS(9)+STwo*Two
      CASE(301)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)
        CTwo=Cpq(2)*Q%C(2)+Spq(2)*Q%S(2)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)
        CTwo=-(Cpq(4)*Q%C(2))-Spq(4)*Q%S(2)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)
        W(2)=Cpq(3)*Q%C(2)-Cpq(5)*Q%C(2)-Spq(5)*Q%S(2)
        CTwo=W(1)+W(2)
        STwo=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)+(Cpq(3)+Cpq(5))*Q%S(2)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)+Cpq(6)*Q%C(1)
        CTwo=Cpq(7)*Q%C(2)+Spq(7)*Q%S(2)
        SPKetC(3)=SPKetC(3)+COne+CTwo*Two
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        CTwo=Cpq(4)*Q%C(0)+Cpq(7)*Q%C(1)+(-Cpq(6)+Cpq(8))*Q%C(2)+Spq(8)*Q%S(2)
        W(1)=Q%C(0)*Spq(4)+Q%C(1)*Spq(7)
        W(2)=Q%C(2)*Spq(8)-Cpq(6)*Q%S(2)-Cpq(8)*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        W(1)=Cpq(5)*Q%C(0)+Cpq(8)*Q%C(1)
        W(2)=(-Cpq(7)+Cpq(9))*Q%C(2)+(Spq(7)+Spq(9))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(5)+Q%C(1)*Spq(8)
        W(2)=Q%C(2)*(-Spq(7)+Spq(9))-(Cpq(7)+Cpq(9))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
!       -------------------------------- Lp=3, Mp = 0-------------------------------------
        COne=-(Cpq(6)*Q%C(0))-Cpq(10)*Q%C(1)
        CTwo=-(Cpq(11)*Q%C(2))-Spq(11)*Q%S(2)
        SPKetC(6)=SPKetC(6)+COne+CTwo*Two
!       -------------------------------- Lp=3, Mp = 1-------------------------------------
        W(1)=-(Cpq(7)*Q%C(0))-Cpq(11)*Q%C(1)
        W(2)=Cpq(10)*Q%C(2)-Cpq(12)*Q%C(2)-Spq(12)*Q%S(2)
        CTwo=W(1)+W(2)
        STwo=-(Q%C(0)*Spq(7))-Q%C(1)*Spq(11)-Q%C(2)*Spq(12)+(Cpq(10)+Cpq(12))*Q%S(2)
        SPKetC(7)=SPKetC(7)+CTwo*Two
        SPKetS(7)=SPKetS(7)+STwo*Two
!       -------------------------------- Lp=3, Mp = 2-------------------------------------
        W(1)=-(Cpq(8)*Q%C(0))-Cpq(12)*Q%C(1)+Cpq(11)*Q%C(2)
        W(2)=-(Cpq(13)*Q%C(2))-Spq(11)*Q%S(2)-Spq(13)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(8))-Q%C(1)*Spq(12)
        W(2)=Q%C(2)*(Spq(11)-Spq(13))+(Cpq(11)+Cpq(13))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(8)=SPKetC(8)+CTwo*Two
        SPKetS(8)=SPKetS(8)+STwo*Two
!       -------------------------------- Lp=3, Mp = 3-------------------------------------
        W(1)=-(Cpq(9)*Q%C(0))-Cpq(13)*Q%C(1)+Cpq(12)*Q%C(2)
        W(2)=-(Cpq(14)*Q%C(2))-Spq(12)*Q%S(2)-Spq(14)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(9))-Q%C(1)*Spq(13)
        W(2)=Q%C(2)*(Spq(12)-Spq(14))+(Cpq(12)+Cpq(14))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(9)=SPKetC(9)+CTwo*Two
        SPKetS(9)=SPKetS(9)+STwo*Two
      CASE(302)
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
      CASE(303)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)+Cpq(6)*Q%C(6)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)
        W(2)=Cpq(7)*Q%C(7)+Cpq(8)*Q%C(8)+Cpq(9)*Q%C(9)
        W(3)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)
        W(4)=Spq(7)*Q%S(7)+Spq(8)*Q%S(8)+Spq(9)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)-Cpq(6)*Q%C(3)-Cpq(10)*Q%C(6)
        W(1)=-(Cpq(4)*Q%C(2))-Cpq(7)*Q%C(4)-Cpq(8)*Q%C(5)
        W(2)=-(Cpq(11)*Q%C(7))-Cpq(12)*Q%C(8)-Cpq(13)*Q%C(9)
        W(3)=-(Spq(4)*Q%S(2))-Spq(7)*Q%S(4)-Spq(8)*Q%S(5)
        W(4)=-(Spq(11)*Q%S(7))-Spq(12)*Q%S(8)-Spq(13)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)+Cpq(3)*Q%C(2)
        W(2)=-(Cpq(5)*Q%C(2))-Cpq(7)*Q%C(3)+Cpq(6)*Q%C(4)
        W(3)=-(Cpq(8)*Q%C(4))+Cpq(7)*Q%C(5)-Cpq(9)*Q%C(5)
        W(4)=-(Cpq(11)*Q%C(6))+Cpq(10)*Q%C(7)-Cpq(12)*Q%C(7)
        W(5)=Cpq(11)*Q%C(8)-Cpq(13)*Q%C(8)+Cpq(12)*Q%C(9)
        W(6)=-(Cpq(14)*Q%C(9))-Spq(5)*Q%S(2)-Spq(8)*Q%S(4)
        W(7)=Spq(7)*Q%S(5)-Spq(9)*Q%S(5)-Spq(12)*Q%S(7)
        W(8)=Spq(11)*Q%S(8)-Spq(13)*Q%S(8)+Spq(12)*Q%S(9)-Spq(14)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(1)=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)-Q%C(3)*Spq(7)
        W(2)=-(Q%C(4)*Spq(8))-Q%C(5)*(Spq(7)+Spq(9))
        W(3)=-(Q%C(6)*Spq(11))-Q%C(7)*Spq(12)
        W(4)=-(Q%C(8)*(Spq(11)+Spq(13)))-Q%C(9)*(Spq(12)+Spq(14))
        W(5)=(Cpq(3)+Cpq(5))*Q%S(2)+(Cpq(6)+Cpq(8))*Q%S(4)
        W(6)=(Cpq(7)+Cpq(9))*Q%S(5)+(Cpq(10)+Cpq(12))*Q%S(7)
        W(7)=(Cpq(11)+Cpq(13))*Q%S(8)+(Cpq(12)+Cpq(14))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)+Cpq(6)*Q%C(1)+Cpq(10)*Q%C(3)+Cpq(15)*Q%C(6)
        W(1)=Cpq(7)*Q%C(2)+Cpq(11)*Q%C(4)+Cpq(12)*Q%C(5)
        W(2)=Cpq(16)*Q%C(7)+Cpq(17)*Q%C(8)+Cpq(18)*Q%C(9)
        W(3)=Spq(7)*Q%S(2)+Spq(11)*Q%S(4)+Spq(12)*Q%S(5)
        W(4)=Spq(16)*Q%S(7)+Spq(17)*Q%S(8)+Spq(18)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(3)=SPKetC(3)+COne+CTwo*Two
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        W(1)=Cpq(4)*Q%C(0)+Cpq(7)*Q%C(1)+(-Cpq(6)+Cpq(8))*Q%C(2)+Cpq(11)*Q%C(3)
        W(2)=(-Cpq(10)+Cpq(12))*Q%C(4)+(-Cpq(11)+Cpq(13))*Q%C(5)
        W(3)=Cpq(16)*Q%C(6)+(-Cpq(15)+Cpq(17))*Q%C(7)
        W(4)=(-Cpq(16)+Cpq(18))*Q%C(8)+(-Cpq(17)+Cpq(19))*Q%C(9)
        W(5)=Spq(8)*Q%S(2)+Spq(12)*Q%S(4)
        W(6)=(-Spq(11)+Spq(13))*Q%S(5)+Spq(17)*Q%S(7)
        W(7)=(-Spq(16)+Spq(18))*Q%S(8)+(-Spq(17)+Spq(19))*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        W(1)=Q%C(0)*Spq(4)+Q%C(1)*Spq(7)+Q%C(2)*Spq(8)+Q%C(3)*Spq(11)
        W(2)=Q%C(4)*Spq(12)+Q%C(5)*(Spq(11)+Spq(13))+Q%C(6)*Spq(16)+Q%C(7)*Spq(17)
        W(3)=Q%C(8)*(Spq(16)+Spq(18))+Q%C(9)*(Spq(17)+Spq(19))
        W(4)=-((Cpq(6)+Cpq(8))*Q%S(2))-(Cpq(10)+Cpq(12))*Q%S(4)
        W(5)=-((Cpq(11)+Cpq(13))*Q%S(5))-(Cpq(15)+Cpq(17))*Q%S(7)
        W(6)=-((Cpq(16)+Cpq(18))*Q%S(8))-(Cpq(17)+Cpq(19))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        W(1)=Cpq(5)*Q%C(0)+Cpq(8)*Q%C(1)+(-Cpq(7)+Cpq(9))*Q%C(2)+Cpq(12)*Q%C(3)
        W(2)=(-Cpq(11)+Cpq(13))*Q%C(4)+(Cpq(10)+Cpq(14))*Q%C(5)
        W(3)=Cpq(17)*Q%C(6)+(-Cpq(16)+Cpq(18))*Q%C(7)
        W(4)=(Cpq(15)+Cpq(19))*Q%C(8)+(Cpq(16)+Cpq(20))*Q%C(9)
        W(5)=(Spq(7)+Spq(9))*Q%S(2)+(Spq(11)+Spq(13))*Q%S(4)
        W(6)=Spq(14)*Q%S(5)+(Spq(16)+Spq(18))*Q%S(7)+Spq(19)*Q%S(8)+(Spq(16)+Spq(20))*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)
        W(1)=Q%C(0)*Spq(5)+Q%C(1)*Spq(8)+Q%C(2)*(-Spq(7)+Spq(9))+Q%C(3)*Spq(12)
        W(2)=Q%C(4)*(-Spq(11)+Spq(13))+Q%C(5)*Spq(14)
        W(3)=Q%C(6)*Spq(17)+Q%C(7)*(-Spq(16)+Spq(18))
        W(4)=Q%C(8)*Spq(19)+Q%C(9)*(-Spq(16)+Spq(20))
        W(5)=-((Cpq(7)+Cpq(9))*Q%S(2))-(Cpq(11)+Cpq(13))*Q%S(4)
        W(6)=(Cpq(10)-Cpq(14))*Q%S(5)-(Cpq(16)+Cpq(18))*Q%S(7)
        W(7)=(Cpq(15)-Cpq(19))*Q%S(8)+(Cpq(16)-Cpq(20))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
!       -------------------------------- Lp=3, Mp = 0-------------------------------------
        COne=-(Cpq(6)*Q%C(0))-Cpq(10)*Q%C(1)-Cpq(15)*Q%C(3)-Cpq(21)*Q%C(6)
        W(1)=-(Cpq(11)*Q%C(2))-Cpq(16)*Q%C(4)-Cpq(17)*Q%C(5)
        W(2)=-(Cpq(22)*Q%C(7))-Cpq(23)*Q%C(8)-Cpq(24)*Q%C(9)
        W(3)=-(Spq(11)*Q%S(2))-Spq(16)*Q%S(4)-Spq(17)*Q%S(5)
        W(4)=-(Spq(22)*Q%S(7))-Spq(23)*Q%S(8)-Spq(24)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(6)=SPKetC(6)+COne+CTwo*Two
!       -------------------------------- Lp=3, Mp = 1-------------------------------------
        W(1)=-(Cpq(7)*Q%C(0))-Cpq(11)*Q%C(1)+Cpq(10)*Q%C(2)
        W(2)=-(Cpq(12)*Q%C(2))-Cpq(16)*Q%C(3)+Cpq(15)*Q%C(4)
        W(3)=-(Cpq(17)*Q%C(4))+Cpq(16)*Q%C(5)-Cpq(18)*Q%C(5)
        W(4)=-(Cpq(22)*Q%C(6))+Cpq(21)*Q%C(7)-Cpq(23)*Q%C(7)
        W(5)=Cpq(22)*Q%C(8)-Cpq(24)*Q%C(8)+Cpq(23)*Q%C(9)
        W(6)=-(Cpq(25)*Q%C(9))-Spq(12)*Q%S(2)-Spq(17)*Q%S(4)
        W(7)=Spq(16)*Q%S(5)-Spq(18)*Q%S(5)-Spq(23)*Q%S(7)
        W(8)=Spq(22)*Q%S(8)-Spq(24)*Q%S(8)+Spq(23)*Q%S(9)-Spq(25)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(1)=-(Q%C(0)*Spq(7))-Q%C(1)*Spq(11)-Q%C(2)*Spq(12)-Q%C(3)*Spq(16)
        W(2)=-(Q%C(4)*Spq(17))-Q%C(5)*(Spq(16)+Spq(18))
        W(3)=-(Q%C(6)*Spq(22))-Q%C(7)*Spq(23)
        W(4)=-(Q%C(8)*(Spq(22)+Spq(24)))-Q%C(9)*(Spq(23)+Spq(25))
        W(5)=(Cpq(10)+Cpq(12))*Q%S(2)+(Cpq(15)+Cpq(17))*Q%S(4)
        W(6)=(Cpq(16)+Cpq(18))*Q%S(5)+(Cpq(21)+Cpq(23))*Q%S(7)
        W(7)=(Cpq(22)+Cpq(24))*Q%S(8)+(Cpq(23)+Cpq(25))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        SPKetC(7)=SPKetC(7)+CTwo*Two
        SPKetS(7)=SPKetS(7)+STwo*Two
!       -------------------------------- Lp=3, Mp = 2-------------------------------------
        W(1)=-(Cpq(8)*Q%C(0))-Cpq(12)*Q%C(1)+Cpq(11)*Q%C(2)
        W(2)=-(Cpq(13)*Q%C(2))-Cpq(17)*Q%C(3)+Cpq(16)*Q%C(4)
        W(3)=-(Cpq(18)*Q%C(4))-Cpq(15)*Q%C(5)-Cpq(19)*Q%C(5)
        W(4)=-(Cpq(23)*Q%C(6))+Cpq(22)*Q%C(7)-Cpq(24)*Q%C(7)-Cpq(21)*Q%C(8)
        W(5)=-(Cpq(25)*Q%C(8))-Cpq(22)*Q%C(9)-Cpq(26)*Q%C(9)
        W(6)=-(Spq(11)*Q%S(2))-Spq(13)*Q%S(2)-Spq(16)*Q%S(4)
        W(7)=-(Spq(18)*Q%S(4))-Spq(19)*Q%S(5)-Spq(22)*Q%S(7)
        W(8)=-(Spq(24)*Q%S(7))-Spq(25)*Q%S(8)-Spq(22)*Q%S(9)-Spq(26)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(1)=-(Q%C(0)*Spq(8))-Q%C(1)*Spq(12)
        W(2)=Q%C(2)*(Spq(11)-Spq(13))-Q%C(3)*Spq(17)
        W(3)=Q%C(4)*(Spq(16)-Spq(18))-Q%C(5)*Spq(19)
        W(4)=-(Q%C(6)*Spq(23))+Q%C(7)*(Spq(22)-Spq(24))
        W(5)=-(Q%C(8)*Spq(25))+Q%C(9)*(Spq(22)-Spq(26))
        W(6)=(Cpq(11)+Cpq(13))*Q%S(2)+(Cpq(16)+Cpq(18))*Q%S(4)
        W(7)=(-Cpq(15)+Cpq(19))*Q%S(5)+(Cpq(22)+Cpq(24))*Q%S(7)
        W(8)=(-Cpq(21)+Cpq(25))*Q%S(8)+(-Cpq(22)+Cpq(26))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        SPKetC(8)=SPKetC(8)+CTwo*Two
        SPKetS(8)=SPKetS(8)+STwo*Two
!       -------------------------------- Lp=3, Mp = 3-------------------------------------
        W(1)=-(Cpq(9)*Q%C(0))-Cpq(13)*Q%C(1)+Cpq(12)*Q%C(2)
        W(2)=-(Cpq(14)*Q%C(2))-Cpq(18)*Q%C(3)+Cpq(17)*Q%C(4)
        W(3)=-(Cpq(19)*Q%C(4))-Cpq(16)*Q%C(5)-Cpq(20)*Q%C(5)
        W(4)=-(Cpq(24)*Q%C(6))+Cpq(23)*Q%C(7)-Cpq(25)*Q%C(7)-Cpq(22)*Q%C(8)
        W(5)=-(Cpq(26)*Q%C(8))+Cpq(21)*Q%C(9)-Cpq(27)*Q%C(9)
        W(6)=-(Spq(12)*Q%S(2))-Spq(14)*Q%S(2)-Spq(17)*Q%S(4)-Spq(19)*Q%S(4)
        W(7)=Spq(16)*Q%S(5)-Spq(20)*Q%S(5)-Spq(23)*Q%S(7)
        W(8)=-(Spq(25)*Q%S(7))+Spq(22)*Q%S(8)-Spq(26)*Q%S(8)-Spq(27)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(1)=-(Q%C(0)*Spq(9))-Q%C(1)*Spq(13)
        W(2)=Q%C(2)*(Spq(12)-Spq(14))-Q%C(3)*Spq(18)
        W(3)=Q%C(4)*(Spq(17)-Spq(19))-Q%C(5)*(Spq(16)+Spq(20))
        W(4)=-(Q%C(6)*Spq(24))+Q%C(7)*(Spq(23)-Spq(25))
        W(5)=-(Q%C(8)*(Spq(22)+Spq(26)))-Q%C(9)*Spq(27)
        W(6)=(Cpq(12)+Cpq(14))*Q%S(2)+(Cpq(17)+Cpq(19))*Q%S(4)
        W(7)=(-Cpq(16)+Cpq(20))*Q%S(5)+(Cpq(23)+Cpq(25))*Q%S(7)
        W(8)=(-Cpq(22)+Cpq(26))*Q%S(8)+(Cpq(21)+Cpq(27))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        SPKetC(9)=SPKetC(9)+CTwo*Two
        SPKetS(9)=SPKetS(9)+STwo*Two
      CASE(304)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)+Cpq(6)*Q%C(6)+Cpq(10)*Q%C(10)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)+Cpq(7)*Q%C(7)+Cpq(8)*Q%C(8)
        W(2)=Cpq(9)*Q%C(9)+Cpq(11)*Q%C(11)+Cpq(12)*Q%C(12)+Cpq(13)*Q%C(13)+Cpq(14)*Q%C(14)
        W(3)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)+Spq(7)*Q%S(7)+Spq(8)*Q%S(8)
        W(4)=Spq(9)*Q%S(9)+Spq(11)*Q%S(11)+Spq(12)*Q%S(12)+Spq(13)*Q%S(13)+Spq(14)*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        W(1)=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)
        W(2)=-(Cpq(6)*Q%C(3))-Cpq(10)*Q%C(6)-Cpq(15)*Q%C(10)
        COne=W(1)+W(2)
        W(1)=-(Cpq(4)*Q%C(2))-Cpq(7)*Q%C(4)
        W(2)=-(Cpq(8)*Q%C(5))-Cpq(11)*Q%C(7)-Cpq(12)*Q%C(8)
        W(3)=-(Cpq(13)*Q%C(9))-Cpq(16)*Q%C(11)
        W(4)=-(Cpq(17)*Q%C(12))-Cpq(18)*Q%C(13)-Cpq(19)*Q%C(14)
        W(5)=-(Spq(4)*Q%S(2))-Spq(7)*Q%S(4)
        W(6)=-(Spq(8)*Q%S(5))-Spq(11)*Q%S(7)-Spq(12)*Q%S(8)
        W(7)=-(Spq(13)*Q%S(9))-Spq(16)*Q%S(11)
        W(8)=-(Spq(17)*Q%S(12))-Spq(18)*Q%S(13)-Spq(19)*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)+(Cpq(3)-Cpq(5))*Q%C(2)
        W(2)=-(Cpq(7)*Q%C(3))+(Cpq(6)-Cpq(8))*Q%C(4)+(Cpq(7)-Cpq(9))*Q%C(5)
        W(3)=-(Cpq(11)*Q%C(6))+(Cpq(10)-Cpq(12))*Q%C(7)+(Cpq(11)-Cpq(13))*Q%C(8)
        W(4)=(Cpq(12)-Cpq(14))*Q%C(9)-Cpq(16)*Q%C(10)+(Cpq(15)-Cpq(17))*Q%C(11)
        W(5)=(Cpq(16)-Cpq(18))*Q%C(12)
        W(6)=(Cpq(17)-Cpq(19))*Q%C(13)+(Cpq(18)-Cpq(20))*Q%C(14)
        W(7)=-(Spq(5)*Q%S(2))-Spq(8)*Q%S(4)+(Spq(7)-Spq(9))*Q%S(5)
        W(8)=-(Spq(12)*Q%S(7))+(Spq(11)-Spq(13))*Q%S(8)+(Spq(12)-Spq(14))*Q%S(9)
        W(9)=-(Spq(17)*Q%S(11))+(Spq(16)-Spq(18))*Q%S(12)
        W(10)=(Spq(17)-Spq(19))*Q%S(13)+(Spq(18)-Spq(20))*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)
        W(1)=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)
        W(2)=-(Q%C(3)*Spq(7))-Q%C(4)*Spq(8)-Q%C(5)*(Spq(7)+Spq(9))
        W(3)=-(Q%C(6)*Spq(11))-Q%C(7)*Spq(12)-Q%C(8)*(Spq(11)+Spq(13))
        W(4)=-(Q%C(9)*(Spq(12)+Spq(14)))-Q%C(10)*Spq(16)-Q%C(11)*Spq(17)
        W(5)=-(Q%C(12)*(Spq(16)+Spq(18)))
        W(6)=-(Q%C(13)*(Spq(17)+Spq(19)))-Q%C(14)*(Spq(18)+Spq(20))
        W(7)=(Cpq(3)+Cpq(5))*Q%S(2)+(Cpq(6)+Cpq(8))*Q%S(4)+(Cpq(7)+Cpq(9))*Q%S(5)
        W(8)=(Cpq(10)+Cpq(12))*Q%S(7)+(Cpq(11)+Cpq(13))*Q%S(8)+(Cpq(12)+Cpq(14))*Q%S(9)
        W(9)=(Cpq(15)+Cpq(17))*Q%S(11)+(Cpq(16)+Cpq(18))*Q%S(12)
        W(10)=(Cpq(17)+Cpq(19))*Q%S(13)+(Cpq(18)+Cpq(20))*Q%S(14)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)+Cpq(6)*Q%C(1)+Cpq(10)*Q%C(3)+Cpq(15)*Q%C(6)+Cpq(21)*Q%C(10)
        W(1)=Cpq(7)*Q%C(2)+Cpq(11)*Q%C(4)+Cpq(12)*Q%C(5)+Cpq(16)*Q%C(7)+Cpq(17)*Q%C(8)
        W(2)=Cpq(18)*Q%C(9)+Cpq(22)*Q%C(11)+Cpq(23)*Q%C(12)+Cpq(24)*Q%C(13)+Cpq(25)*Q%C(14)
        W(3)=Spq(7)*Q%S(2)+Spq(11)*Q%S(4)+Spq(12)*Q%S(5)+Spq(16)*Q%S(7)+Spq(17)*Q%S(8)
        W(4)=Spq(18)*Q%S(9)+Spq(22)*Q%S(11)+Spq(23)*Q%S(12)+Spq(24)*Q%S(13)+Spq(25)*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(3)=SPKetC(3)+COne+CTwo*Two
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        W(1)=Cpq(4)*Q%C(0)+Cpq(7)*Q%C(1)+(-Cpq(6)+Cpq(8))*Q%C(2)
        W(2)=Cpq(11)*Q%C(3)+(-Cpq(10)+Cpq(12))*Q%C(4)+(-Cpq(11)+Cpq(13))*Q%C(5)
        W(3)=Cpq(16)*Q%C(6)+(-Cpq(15)+Cpq(17))*Q%C(7)+(-Cpq(16)+Cpq(18))*Q%C(8)
        W(4)=(-Cpq(17)+Cpq(19))*Q%C(9)+Cpq(22)*Q%C(10)+(-Cpq(21)+Cpq(23))*Q%C(11)
        W(5)=(-Cpq(22)+Cpq(24))*Q%C(12)
        W(6)=(-Cpq(23)+Cpq(25))*Q%C(13)+(-Cpq(24)+Cpq(26))*Q%C(14)
        W(7)=Spq(8)*Q%S(2)+Spq(12)*Q%S(4)+(-Spq(11)+Spq(13))*Q%S(5)
        W(8)=Spq(17)*Q%S(7)+(-Spq(16)+Spq(18))*Q%S(8)+(-Spq(17)+Spq(19))*Q%S(9)
        W(9)=Spq(23)*Q%S(11)+(-Spq(22)+Spq(24))*Q%S(12)
        W(10)=(-Spq(23)+Spq(25))*Q%S(13)+(-Spq(24)+Spq(26))*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)
        W(1)=Q%C(0)*Spq(4)+Q%C(1)*Spq(7)+Q%C(2)*Spq(8)
        W(2)=Q%C(3)*Spq(11)+Q%C(4)*Spq(12)+Q%C(5)*(Spq(11)+Spq(13))
        W(3)=Q%C(6)*Spq(16)+Q%C(7)*Spq(17)+Q%C(8)*(Spq(16)+Spq(18))
        W(4)=Q%C(9)*(Spq(17)+Spq(19))+Q%C(10)*Spq(22)+Q%C(11)*Spq(23)
        W(5)=Q%C(12)*(Spq(22)+Spq(24))+Q%C(13)*(Spq(23)+Spq(25))+Q%C(14)*(Spq(24)+Spq(26))
        W(6)=-((Cpq(6)+Cpq(8))*Q%S(2))
        W(7)=-((Cpq(10)+Cpq(12))*Q%S(4))-(Cpq(11)+Cpq(13))*Q%S(5)
        W(8)=-((Cpq(15)+Cpq(17))*Q%S(7))
        W(9)=-((Cpq(16)+Cpq(18))*Q%S(8))-(Cpq(17)+Cpq(19))*Q%S(9)
        W(10)=-((Cpq(21)+Cpq(23))*Q%S(11))-(Cpq(22)+Cpq(24))*Q%S(12)
        W(11)=-((Cpq(23)+Cpq(25))*Q%S(13))-(Cpq(24)+Cpq(26))*Q%S(14)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)+W(11)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        W(1)=Cpq(5)*Q%C(0)+Cpq(8)*Q%C(1)+(-Cpq(7)+Cpq(9))*Q%C(2)
        W(2)=Cpq(12)*Q%C(3)+(-Cpq(11)+Cpq(13))*Q%C(4)+(Cpq(10)+Cpq(14))*Q%C(5)
        W(3)=Cpq(17)*Q%C(6)+(-Cpq(16)+Cpq(18))*Q%C(7)+(Cpq(15)+Cpq(19))*Q%C(8)
        W(4)=(Cpq(16)+Cpq(20))*Q%C(9)+Cpq(23)*Q%C(10)+(-Cpq(22)+Cpq(24))*Q%C(11)
        W(5)=(Cpq(21)+Cpq(25))*Q%C(12)+(Cpq(22)+Cpq(26))*Q%C(13)+(Cpq(23)+Cpq(27))*Q%C(14)
        W(6)=(Spq(7)+Spq(9))*Q%S(2)+(Spq(11)+Spq(13))*Q%S(4)+Spq(14)*Q%S(5)
        W(7)=(Spq(16)+Spq(18))*Q%S(7)+Spq(19)*Q%S(8)+(Spq(16)+Spq(20))*Q%S(9)
        W(8)=(Spq(22)+Spq(24))*Q%S(11)+Spq(25)*Q%S(12)
        W(9)=(Spq(22)+Spq(26))*Q%S(13)+(Spq(23)+Spq(27))*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)
        W(1)=Q%C(0)*Spq(5)+Q%C(1)*Spq(8)+Q%C(2)*(-Spq(7)+Spq(9))
        W(2)=Q%C(3)*Spq(12)+Q%C(4)*(-Spq(11)+Spq(13))+Q%C(5)*Spq(14)
        W(3)=Q%C(6)*Spq(17)+Q%C(7)*(-Spq(16)+Spq(18))+Q%C(8)*Spq(19)
        W(4)=Q%C(9)*(-Spq(16)+Spq(20))+Q%C(10)*Spq(23)+Q%C(11)*(-Spq(22)+Spq(24))
        W(5)=Q%C(12)*Spq(25)+Q%C(13)*(-Spq(22)+Spq(26))+Q%C(14)*(-Spq(23)+Spq(27))
        W(6)=-((Cpq(7)+Cpq(9))*Q%S(2))
        W(7)=-((Cpq(11)+Cpq(13))*Q%S(4))+(Cpq(10)-Cpq(14))*Q%S(5)
        W(8)=-((Cpq(16)+Cpq(18))*Q%S(7))
        W(9)=(Cpq(15)-Cpq(19))*Q%S(8)+(Cpq(16)-Cpq(20))*Q%S(9)
        W(10)=-((Cpq(22)+Cpq(24))*Q%S(11))+(Cpq(21)-Cpq(25))*Q%S(12)
        W(11)=(Cpq(22)-Cpq(26))*Q%S(13)+(Cpq(23)-Cpq(27))*Q%S(14)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)+W(11)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
!       -------------------------------- Lp=3, Mp = 0-------------------------------------
        W(1)=-(Cpq(6)*Q%C(0))-Cpq(10)*Q%C(1)
        W(2)=-(Cpq(15)*Q%C(3))-Cpq(21)*Q%C(6)-Cpq(28)*Q%C(10)
        COne=W(1)+W(2)
        W(1)=-(Cpq(11)*Q%C(2))-Cpq(16)*Q%C(4)
        W(2)=-(Cpq(17)*Q%C(5))-Cpq(22)*Q%C(7)-Cpq(23)*Q%C(8)
        W(3)=-(Cpq(24)*Q%C(9))-Cpq(29)*Q%C(11)
        W(4)=-(Cpq(30)*Q%C(12))-Cpq(31)*Q%C(13)-Cpq(32)*Q%C(14)
        W(5)=-(Spq(11)*Q%S(2))-Spq(16)*Q%S(4)
        W(6)=-(Spq(17)*Q%S(5))-Spq(22)*Q%S(7)-Spq(23)*Q%S(8)
        W(7)=-(Spq(24)*Q%S(9))-Spq(29)*Q%S(11)
        W(8)=-(Spq(30)*Q%S(12))-Spq(31)*Q%S(13)-Spq(32)*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        SPKetC(6)=SPKetC(6)+COne+CTwo*Two
!       -------------------------------- Lp=3, Mp = 1-------------------------------------
        W(1)=-(Cpq(7)*Q%C(0))-Cpq(11)*Q%C(1)+(Cpq(10)-Cpq(12))*Q%C(2)
        W(2)=-(Cpq(16)*Q%C(3))+(Cpq(15)-Cpq(17))*Q%C(4)+(Cpq(16)-Cpq(18))*Q%C(5)
        W(3)=-(Cpq(22)*Q%C(6))+(Cpq(21)-Cpq(23))*Q%C(7)+(Cpq(22)-Cpq(24))*Q%C(8)
        W(4)=(Cpq(23)-Cpq(25))*Q%C(9)-Cpq(29)*Q%C(10)+(Cpq(28)-Cpq(30))*Q%C(11)
        W(5)=(Cpq(29)-Cpq(31))*Q%C(12)
        W(6)=(Cpq(30)-Cpq(32))*Q%C(13)+(Cpq(31)-Cpq(33))*Q%C(14)
        W(7)=-(Spq(12)*Q%S(2))-Spq(17)*Q%S(4)+(Spq(16)-Spq(18))*Q%S(5)
        W(8)=-(Spq(23)*Q%S(7))+(Spq(22)-Spq(24))*Q%S(8)+(Spq(23)-Spq(25))*Q%S(9)
        W(9)=-(Spq(30)*Q%S(11))+(Spq(29)-Spq(31))*Q%S(12)
        W(10)=(Spq(30)-Spq(32))*Q%S(13)+(Spq(31)-Spq(33))*Q%S(14)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)
        W(1)=-(Q%C(0)*Spq(7))-Q%C(1)*Spq(11)-Q%C(2)*Spq(12)
        W(2)=-(Q%C(3)*Spq(16))-Q%C(4)*Spq(17)-Q%C(5)*(Spq(16)+Spq(18))
        W(3)=-(Q%C(6)*Spq(22))-Q%C(7)*Spq(23)-Q%C(8)*(Spq(22)+Spq(24))
        W(4)=-(Q%C(9)*(Spq(23)+Spq(25)))-Q%C(10)*Spq(29)-Q%C(11)*Spq(30)
        W(5)=-(Q%C(12)*(Spq(29)+Spq(31)))
        W(6)=-(Q%C(13)*(Spq(30)+Spq(32)))-Q%C(14)*(Spq(31)+Spq(33))
        W(7)=(Cpq(10)+Cpq(12))*Q%S(2)+(Cpq(15)+Cpq(17))*Q%S(4)+(Cpq(16)+Cpq(18))*Q%S(5)
        W(8)=(Cpq(21)+Cpq(23))*Q%S(7)+(Cpq(22)+Cpq(24))*Q%S(8)+(Cpq(23)+Cpq(25))*Q%S(9)
        W(9)=(Cpq(28)+Cpq(30))*Q%S(11)+(Cpq(29)+Cpq(31))*Q%S(12)
        W(10)=(Cpq(30)+Cpq(32))*Q%S(13)+(Cpq(31)+Cpq(33))*Q%S(14)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)
        SPKetC(7)=SPKetC(7)+CTwo*Two
        SPKetS(7)=SPKetS(7)+STwo*Two
!       -------------------------------- Lp=3, Mp = 2-------------------------------------
        W(1)=-(Cpq(8)*Q%C(0))-Cpq(12)*Q%C(1)
        W(2)=Cpq(11)*Q%C(2)-Cpq(13)*Q%C(2)-Cpq(17)*Q%C(3)
        W(3)=Cpq(16)*Q%C(4)-Cpq(18)*Q%C(4)
        W(4)=-(Cpq(15)*Q%C(5))-Cpq(19)*Q%C(5)-Cpq(23)*Q%C(6)
        W(5)=Cpq(22)*Q%C(7)-Cpq(24)*Q%C(7)
        W(6)=-(Cpq(21)*Q%C(8))-Cpq(25)*Q%C(8)-Cpq(22)*Q%C(9)
        W(7)=-(Cpq(26)*Q%C(9))-Cpq(30)*Q%C(10)+Cpq(29)*Q%C(11)
        W(8)=-(Cpq(31)*Q%C(11))-Cpq(28)*Q%C(12)-Cpq(32)*Q%C(12)
        W(9)=-(Cpq(29)*Q%C(13))-Cpq(33)*Q%C(13)
        W(10)=-(Cpq(30)*Q%C(14))-Cpq(34)*Q%C(14)-Spq(11)*Q%S(2)
        W(11)=-(Spq(13)*Q%S(2))-Spq(16)*Q%S(4)
        W(12)=-(Spq(18)*Q%S(4))-Spq(19)*Q%S(5)-Spq(22)*Q%S(7)
        W(13)=-(Spq(24)*Q%S(7))-Spq(25)*Q%S(8)
        W(14)=-(Spq(22)*Q%S(9))-Spq(26)*Q%S(9)-Spq(29)*Q%S(11)
        W(15)=-(Spq(31)*Q%S(11))-Spq(32)*Q%S(12)-Spq(29)*Q%S(13)
        W(16)=-(Spq(33)*Q%S(13))-Spq(30)*Q%S(14)-Spq(34)*Q%S(14)
        W(17)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(18)=W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)+W(16)
        CTwo=W(17)+W(18)
        W(1)=-(Q%C(0)*Spq(8))-Q%C(1)*Spq(12)+Q%C(2)*(Spq(11)-Spq(13))
        W(2)=-(Q%C(3)*Spq(17))+Q%C(4)*(Spq(16)-Spq(18))-Q%C(5)*Spq(19)
        W(3)=-(Q%C(6)*Spq(23))+Q%C(7)*(Spq(22)-Spq(24))-Q%C(8)*Spq(25)
        W(4)=Q%C(9)*(Spq(22)-Spq(26))-Q%C(10)*Spq(30)+Q%C(11)*(Spq(29)-Spq(31))
        W(5)=-(Q%C(12)*Spq(32))+Q%C(13)*(Spq(29)-Spq(33))+Q%C(14)*(Spq(30)-Spq(34))
        W(6)=(Cpq(11)+Cpq(13))*Q%S(2)+(Cpq(16)+Cpq(18))*Q%S(4)+(-Cpq(15)+Cpq(19))*Q%S(5)
        W(7)=(Cpq(22)+Cpq(24))*Q%S(7)
        W(8)=(-Cpq(21)+Cpq(25))*Q%S(8)+(-Cpq(22)+Cpq(26))*Q%S(9)
        W(9)=(Cpq(29)+Cpq(31))*Q%S(11)+(-Cpq(28)+Cpq(32))*Q%S(12)
        W(10)=(-Cpq(29)+Cpq(33))*Q%S(13)+(-Cpq(30)+Cpq(34))*Q%S(14)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)+W(10)
        SPKetC(8)=SPKetC(8)+CTwo*Two
        SPKetS(8)=SPKetS(8)+STwo*Two
!       -------------------------------- Lp=3, Mp = 3-------------------------------------
        W(1)=-(Cpq(9)*Q%C(0))-Cpq(13)*Q%C(1)
        W(2)=Cpq(12)*Q%C(2)-Cpq(14)*Q%C(2)-Cpq(18)*Q%C(3)
        W(3)=Cpq(17)*Q%C(4)-Cpq(19)*Q%C(4)
        W(4)=-(Cpq(16)*Q%C(5))-Cpq(20)*Q%C(5)-Cpq(24)*Q%C(6)
        W(5)=Cpq(23)*Q%C(7)-Cpq(25)*Q%C(7)
        W(6)=-(Cpq(22)*Q%C(8))-Cpq(26)*Q%C(8)+Cpq(21)*Q%C(9)
        W(7)=-(Cpq(27)*Q%C(9))-Cpq(31)*Q%C(10)+Cpq(30)*Q%C(11)
        W(8)=-(Cpq(32)*Q%C(11))-Cpq(29)*Q%C(12)-Cpq(33)*Q%C(12)
        W(9)=Cpq(28)*Q%C(13)-Cpq(34)*Q%C(13)
        W(10)=Cpq(29)*Q%C(14)-Cpq(35)*Q%C(14)-Spq(12)*Q%S(2)
        W(11)=-(Spq(14)*Q%S(2))-Spq(17)*Q%S(4)-Spq(19)*Q%S(4)
        W(12)=Spq(16)*Q%S(5)-Spq(20)*Q%S(5)-Spq(23)*Q%S(7)
        W(13)=-(Spq(25)*Q%S(7))+Spq(22)*Q%S(8)
        W(14)=-(Spq(26)*Q%S(8))-Spq(27)*Q%S(9)-Spq(30)*Q%S(11)
        W(15)=-(Spq(32)*Q%S(11))+Spq(29)*Q%S(12)-Spq(33)*Q%S(12)
        W(16)=-(Spq(34)*Q%S(13))+Spq(29)*Q%S(14)-Spq(35)*Q%S(14)
        W(17)=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(18)=W(9)+W(10)+W(11)+W(12)+W(13)+W(14)+W(15)+W(16)
        CTwo=W(17)+W(18)
        W(1)=-(Q%C(0)*Spq(9))-Q%C(1)*Spq(13)+Q%C(2)*(Spq(12)-Spq(14))
        W(2)=-(Q%C(3)*Spq(18))+Q%C(4)*(Spq(17)-Spq(19))-Q%C(5)*(Spq(16)+Spq(20))
        W(3)=-(Q%C(6)*Spq(24))+Q%C(7)*(Spq(23)-Spq(25))-Q%C(8)*(Spq(22)+Spq(26))
        W(4)=-(Q%C(9)*Spq(27))-Q%C(10)*Spq(31)+Q%C(11)*(Spq(30)-Spq(32))
        W(5)=-(Q%C(12)*(Spq(29)+Spq(33)))-Q%C(13)*Spq(34)-Q%C(14)*(Spq(29)+Spq(35))
        W(6)=(Cpq(12)+Cpq(14))*Q%S(2)+(Cpq(17)+Cpq(19))*Q%S(4)+(-Cpq(16)+Cpq(20))*Q%S(5)
        W(7)=(Cpq(23)+Cpq(25))*Q%S(7)+(-Cpq(22)+Cpq(26))*Q%S(8)+(Cpq(21)+Cpq(27))*Q%S(9)
        W(8)=(Cpq(30)+Cpq(32))*Q%S(11)+(-Cpq(29)+Cpq(33))*Q%S(12)
        W(9)=(Cpq(28)+Cpq(34))*Q%S(13)+(Cpq(29)+Cpq(35))*Q%S(14)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)+W(9)
        SPKetC(9)=SPKetC(9)+CTwo*Two
        SPKetS(9)=SPKetS(9)+STwo*Two
      CASE(305)
         CALL CTraX305(Q)
      CASE(306)
         CALL CTraX306(Q)
      CASE(307)
         CALL CTraX307(Q)
      CASE(400)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)
        SPKetC(0)=SPKetC(0)+COne
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))
        SPKetC(1)=SPKetC(1)+COne
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        CTwo=-(Cpq(2)*Q%C(0))
        STwo=-(Q%C(0)*Spq(2))
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)
        SPKetC(3)=SPKetC(3)+COne
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        CTwo=Cpq(4)*Q%C(0)
        STwo=Q%C(0)*Spq(4)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        CTwo=Cpq(5)*Q%C(0)
        STwo=Q%C(0)*Spq(5)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
!       -------------------------------- Lp=3, Mp = 0-------------------------------------
        COne=-(Cpq(6)*Q%C(0))
        SPKetC(6)=SPKetC(6)+COne
!       -------------------------------- Lp=3, Mp = 1-------------------------------------
        CTwo=-(Cpq(7)*Q%C(0))
        STwo=-(Q%C(0)*Spq(7))
        SPKetC(7)=SPKetC(7)+CTwo*Two
        SPKetS(7)=SPKetS(7)+STwo*Two
!       -------------------------------- Lp=3, Mp = 2-------------------------------------
        CTwo=-(Cpq(8)*Q%C(0))
        STwo=-(Q%C(0)*Spq(8))
        SPKetC(8)=SPKetC(8)+CTwo*Two
        SPKetS(8)=SPKetS(8)+STwo*Two
!       -------------------------------- Lp=3, Mp = 3-------------------------------------
        CTwo=-(Cpq(9)*Q%C(0))
        STwo=-(Q%C(0)*Spq(9))
        SPKetC(9)=SPKetC(9)+CTwo*Two
        SPKetS(9)=SPKetS(9)+STwo*Two
!       -------------------------------- Lp=4, Mp = 0-------------------------------------
        COne=Cpq(10)*Q%C(0)
        SPKetC(10)=SPKetC(10)+COne
!       -------------------------------- Lp=4, Mp = 1-------------------------------------
        CTwo=Cpq(11)*Q%C(0)
        STwo=Q%C(0)*Spq(11)
        SPKetC(11)=SPKetC(11)+CTwo*Two
        SPKetS(11)=SPKetS(11)+STwo*Two
!       -------------------------------- Lp=4, Mp = 2-------------------------------------
        CTwo=Cpq(12)*Q%C(0)
        STwo=Q%C(0)*Spq(12)
        SPKetC(12)=SPKetC(12)+CTwo*Two
        SPKetS(12)=SPKetS(12)+STwo*Two
!       -------------------------------- Lp=4, Mp = 3-------------------------------------
        CTwo=Cpq(13)*Q%C(0)
        STwo=Q%C(0)*Spq(13)
        SPKetC(13)=SPKetC(13)+CTwo*Two
        SPKetS(13)=SPKetS(13)+STwo*Two
!       -------------------------------- Lp=4, Mp = 4-------------------------------------
        CTwo=Cpq(14)*Q%C(0)
        STwo=Q%C(0)*Spq(14)
        SPKetC(14)=SPKetC(14)+CTwo*Two
        SPKetS(14)=SPKetS(14)+STwo*Two
      CASE(401)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)
        CTwo=Cpq(2)*Q%C(2)+Spq(2)*Q%S(2)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)
        CTwo=-(Cpq(4)*Q%C(2))-Spq(4)*Q%S(2)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)
        W(2)=Cpq(3)*Q%C(2)-Cpq(5)*Q%C(2)-Spq(5)*Q%S(2)
        CTwo=W(1)+W(2)
        STwo=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)+(Cpq(3)+Cpq(5))*Q%S(2)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)+Cpq(6)*Q%C(1)
        CTwo=Cpq(7)*Q%C(2)+Spq(7)*Q%S(2)
        SPKetC(3)=SPKetC(3)+COne+CTwo*Two
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        CTwo=Cpq(4)*Q%C(0)+Cpq(7)*Q%C(1)+(-Cpq(6)+Cpq(8))*Q%C(2)+Spq(8)*Q%S(2)
        W(1)=Q%C(0)*Spq(4)+Q%C(1)*Spq(7)
        W(2)=Q%C(2)*Spq(8)-Cpq(6)*Q%S(2)-Cpq(8)*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        W(1)=Cpq(5)*Q%C(0)+Cpq(8)*Q%C(1)
        W(2)=(-Cpq(7)+Cpq(9))*Q%C(2)+(Spq(7)+Spq(9))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(5)+Q%C(1)*Spq(8)
        W(2)=Q%C(2)*(-Spq(7)+Spq(9))-(Cpq(7)+Cpq(9))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
!       -------------------------------- Lp=3, Mp = 0-------------------------------------
        COne=-(Cpq(6)*Q%C(0))-Cpq(10)*Q%C(1)
        CTwo=-(Cpq(11)*Q%C(2))-Spq(11)*Q%S(2)
        SPKetC(6)=SPKetC(6)+COne+CTwo*Two
!       -------------------------------- Lp=3, Mp = 1-------------------------------------
        W(1)=-(Cpq(7)*Q%C(0))-Cpq(11)*Q%C(1)
        W(2)=Cpq(10)*Q%C(2)-Cpq(12)*Q%C(2)-Spq(12)*Q%S(2)
        CTwo=W(1)+W(2)
        STwo=-(Q%C(0)*Spq(7))-Q%C(1)*Spq(11)-Q%C(2)*Spq(12)+(Cpq(10)+Cpq(12))*Q%S(2)
        SPKetC(7)=SPKetC(7)+CTwo*Two
        SPKetS(7)=SPKetS(7)+STwo*Two
!       -------------------------------- Lp=3, Mp = 2-------------------------------------
        W(1)=-(Cpq(8)*Q%C(0))-Cpq(12)*Q%C(1)+Cpq(11)*Q%C(2)
        W(2)=-(Cpq(13)*Q%C(2))-Spq(11)*Q%S(2)-Spq(13)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(8))-Q%C(1)*Spq(12)
        W(2)=Q%C(2)*(Spq(11)-Spq(13))+(Cpq(11)+Cpq(13))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(8)=SPKetC(8)+CTwo*Two
        SPKetS(8)=SPKetS(8)+STwo*Two
!       -------------------------------- Lp=3, Mp = 3-------------------------------------
        W(1)=-(Cpq(9)*Q%C(0))-Cpq(13)*Q%C(1)+Cpq(12)*Q%C(2)
        W(2)=-(Cpq(14)*Q%C(2))-Spq(12)*Q%S(2)-Spq(14)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(9))-Q%C(1)*Spq(13)
        W(2)=Q%C(2)*(Spq(12)-Spq(14))+(Cpq(12)+Cpq(14))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(9)=SPKetC(9)+CTwo*Two
        SPKetS(9)=SPKetS(9)+STwo*Two
!       -------------------------------- Lp=4, Mp = 0-------------------------------------
        COne=Cpq(10)*Q%C(0)+Cpq(15)*Q%C(1)
        CTwo=Cpq(16)*Q%C(2)+Spq(16)*Q%S(2)
        SPKetC(10)=SPKetC(10)+COne+CTwo*Two
!       -------------------------------- Lp=4, Mp = 1-------------------------------------
        CTwo=Cpq(11)*Q%C(0)+Cpq(16)*Q%C(1)+(-Cpq(15)+Cpq(17))*Q%C(2)+Spq(17)*Q%S(2)
        W(1)=Q%C(0)*Spq(11)+Q%C(1)*Spq(16)
        W(2)=Q%C(2)*Spq(17)-Cpq(15)*Q%S(2)-Cpq(17)*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(11)=SPKetC(11)+CTwo*Two
        SPKetS(11)=SPKetS(11)+STwo*Two
!       -------------------------------- Lp=4, Mp = 2-------------------------------------
        W(1)=Cpq(12)*Q%C(0)+Cpq(17)*Q%C(1)
        W(2)=(-Cpq(16)+Cpq(18))*Q%C(2)+(Spq(16)+Spq(18))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(12)+Q%C(1)*Spq(17)
        W(2)=Q%C(2)*(-Spq(16)+Spq(18))-(Cpq(16)+Cpq(18))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(12)=SPKetC(12)+CTwo*Two
        SPKetS(12)=SPKetS(12)+STwo*Two
!       -------------------------------- Lp=4, Mp = 3-------------------------------------
        W(1)=Cpq(13)*Q%C(0)+Cpq(18)*Q%C(1)
        W(2)=(-Cpq(17)+Cpq(19))*Q%C(2)+(Spq(17)+Spq(19))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(13)+Q%C(1)*Spq(18)
        W(2)=Q%C(2)*(-Spq(17)+Spq(19))-(Cpq(17)+Cpq(19))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(13)=SPKetC(13)+CTwo*Two
        SPKetS(13)=SPKetS(13)+STwo*Two
!       -------------------------------- Lp=4, Mp = 4-------------------------------------
        W(1)=Cpq(14)*Q%C(0)+Cpq(19)*Q%C(1)
        W(2)=(-Cpq(18)+Cpq(20))*Q%C(2)+(Spq(18)+Spq(20))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(14)+Q%C(1)*Spq(19)
        W(2)=Q%C(2)*(-Spq(18)+Spq(20))-(Cpq(18)+Cpq(20))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(14)=SPKetC(14)+CTwo*Two
        SPKetS(14)=SPKetS(14)+STwo*Two
      CASE(402)
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
      CASE(403)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)+Cpq(3)*Q%C(3)+Cpq(6)*Q%C(6)
        W(1)=Cpq(2)*Q%C(2)+Cpq(4)*Q%C(4)+Cpq(5)*Q%C(5)
        W(2)=Cpq(7)*Q%C(7)+Cpq(8)*Q%C(8)+Cpq(9)*Q%C(9)
        W(3)=Spq(2)*Q%S(2)+Spq(4)*Q%S(4)+Spq(5)*Q%S(5)
        W(4)=Spq(7)*Q%S(7)+Spq(8)*Q%S(8)+Spq(9)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)-Cpq(6)*Q%C(3)-Cpq(10)*Q%C(6)
        W(1)=-(Cpq(4)*Q%C(2))-Cpq(7)*Q%C(4)-Cpq(8)*Q%C(5)
        W(2)=-(Cpq(11)*Q%C(7))-Cpq(12)*Q%C(8)-Cpq(13)*Q%C(9)
        W(3)=-(Spq(4)*Q%S(2))-Spq(7)*Q%S(4)-Spq(8)*Q%S(5)
        W(4)=-(Spq(11)*Q%S(7))-Spq(12)*Q%S(8)-Spq(13)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)+Cpq(3)*Q%C(2)
        W(2)=-(Cpq(5)*Q%C(2))-Cpq(7)*Q%C(3)+Cpq(6)*Q%C(4)
        W(3)=-(Cpq(8)*Q%C(4))+Cpq(7)*Q%C(5)-Cpq(9)*Q%C(5)
        W(4)=-(Cpq(11)*Q%C(6))+Cpq(10)*Q%C(7)-Cpq(12)*Q%C(7)
        W(5)=Cpq(11)*Q%C(8)-Cpq(13)*Q%C(8)+Cpq(12)*Q%C(9)
        W(6)=-(Cpq(14)*Q%C(9))-Spq(5)*Q%S(2)-Spq(8)*Q%S(4)
        W(7)=Spq(7)*Q%S(5)-Spq(9)*Q%S(5)-Spq(12)*Q%S(7)
        W(8)=Spq(11)*Q%S(8)-Spq(13)*Q%S(8)+Spq(12)*Q%S(9)-Spq(14)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(1)=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)-Q%C(3)*Spq(7)
        W(2)=-(Q%C(4)*Spq(8))-Q%C(5)*(Spq(7)+Spq(9))
        W(3)=-(Q%C(6)*Spq(11))-Q%C(7)*Spq(12)
        W(4)=-(Q%C(8)*(Spq(11)+Spq(13)))-Q%C(9)*(Spq(12)+Spq(14))
        W(5)=(Cpq(3)+Cpq(5))*Q%S(2)+(Cpq(6)+Cpq(8))*Q%S(4)
        W(6)=(Cpq(7)+Cpq(9))*Q%S(5)+(Cpq(10)+Cpq(12))*Q%S(7)
        W(7)=(Cpq(11)+Cpq(13))*Q%S(8)+(Cpq(12)+Cpq(14))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)+Cpq(6)*Q%C(1)+Cpq(10)*Q%C(3)+Cpq(15)*Q%C(6)
        W(1)=Cpq(7)*Q%C(2)+Cpq(11)*Q%C(4)+Cpq(12)*Q%C(5)
        W(2)=Cpq(16)*Q%C(7)+Cpq(17)*Q%C(8)+Cpq(18)*Q%C(9)
        W(3)=Spq(7)*Q%S(2)+Spq(11)*Q%S(4)+Spq(12)*Q%S(5)
        W(4)=Spq(16)*Q%S(7)+Spq(17)*Q%S(8)+Spq(18)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(3)=SPKetC(3)+COne+CTwo*Two
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        W(1)=Cpq(4)*Q%C(0)+Cpq(7)*Q%C(1)+(-Cpq(6)+Cpq(8))*Q%C(2)+Cpq(11)*Q%C(3)
        W(2)=(-Cpq(10)+Cpq(12))*Q%C(4)+(-Cpq(11)+Cpq(13))*Q%C(5)
        W(3)=Cpq(16)*Q%C(6)+(-Cpq(15)+Cpq(17))*Q%C(7)
        W(4)=(-Cpq(16)+Cpq(18))*Q%C(8)+(-Cpq(17)+Cpq(19))*Q%C(9)
        W(5)=Spq(8)*Q%S(2)+Spq(12)*Q%S(4)
        W(6)=(-Spq(11)+Spq(13))*Q%S(5)+Spq(17)*Q%S(7)
        W(7)=(-Spq(16)+Spq(18))*Q%S(8)+(-Spq(17)+Spq(19))*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        W(1)=Q%C(0)*Spq(4)+Q%C(1)*Spq(7)+Q%C(2)*Spq(8)+Q%C(3)*Spq(11)
        W(2)=Q%C(4)*Spq(12)+Q%C(5)*(Spq(11)+Spq(13))+Q%C(6)*Spq(16)+Q%C(7)*Spq(17)
        W(3)=Q%C(8)*(Spq(16)+Spq(18))+Q%C(9)*(Spq(17)+Spq(19))
        W(4)=-((Cpq(6)+Cpq(8))*Q%S(2))-(Cpq(10)+Cpq(12))*Q%S(4)
        W(5)=-((Cpq(11)+Cpq(13))*Q%S(5))-(Cpq(15)+Cpq(17))*Q%S(7)
        W(6)=-((Cpq(16)+Cpq(18))*Q%S(8))-(Cpq(17)+Cpq(19))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        W(1)=Cpq(5)*Q%C(0)+Cpq(8)*Q%C(1)+(-Cpq(7)+Cpq(9))*Q%C(2)+Cpq(12)*Q%C(3)
        W(2)=(-Cpq(11)+Cpq(13))*Q%C(4)+(Cpq(10)+Cpq(14))*Q%C(5)
        W(3)=Cpq(17)*Q%C(6)+(-Cpq(16)+Cpq(18))*Q%C(7)
        W(4)=(Cpq(15)+Cpq(19))*Q%C(8)+(Cpq(16)+Cpq(20))*Q%C(9)
        W(5)=(Spq(7)+Spq(9))*Q%S(2)+(Spq(11)+Spq(13))*Q%S(4)
        W(6)=Spq(14)*Q%S(5)+(Spq(16)+Spq(18))*Q%S(7)+Spq(19)*Q%S(8)+(Spq(16)+Spq(20))*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)
        W(1)=Q%C(0)*Spq(5)+Q%C(1)*Spq(8)+Q%C(2)*(-Spq(7)+Spq(9))+Q%C(3)*Spq(12)
        W(2)=Q%C(4)*(-Spq(11)+Spq(13))+Q%C(5)*Spq(14)
        W(3)=Q%C(6)*Spq(17)+Q%C(7)*(-Spq(16)+Spq(18))
        W(4)=Q%C(8)*Spq(19)+Q%C(9)*(-Spq(16)+Spq(20))
        W(5)=-((Cpq(7)+Cpq(9))*Q%S(2))-(Cpq(11)+Cpq(13))*Q%S(4)
        W(6)=(Cpq(10)-Cpq(14))*Q%S(5)-(Cpq(16)+Cpq(18))*Q%S(7)
        W(7)=(Cpq(15)-Cpq(19))*Q%S(8)+(Cpq(16)-Cpq(20))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
!       -------------------------------- Lp=3, Mp = 0-------------------------------------
        COne=-(Cpq(6)*Q%C(0))-Cpq(10)*Q%C(1)-Cpq(15)*Q%C(3)-Cpq(21)*Q%C(6)
        W(1)=-(Cpq(11)*Q%C(2))-Cpq(16)*Q%C(4)-Cpq(17)*Q%C(5)
        W(2)=-(Cpq(22)*Q%C(7))-Cpq(23)*Q%C(8)-Cpq(24)*Q%C(9)
        W(3)=-(Spq(11)*Q%S(2))-Spq(16)*Q%S(4)-Spq(17)*Q%S(5)
        W(4)=-(Spq(22)*Q%S(7))-Spq(23)*Q%S(8)-Spq(24)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(6)=SPKetC(6)+COne+CTwo*Two
!       -------------------------------- Lp=3, Mp = 1-------------------------------------
        W(1)=-(Cpq(7)*Q%C(0))-Cpq(11)*Q%C(1)+Cpq(10)*Q%C(2)
        W(2)=-(Cpq(12)*Q%C(2))-Cpq(16)*Q%C(3)+Cpq(15)*Q%C(4)
        W(3)=-(Cpq(17)*Q%C(4))+Cpq(16)*Q%C(5)-Cpq(18)*Q%C(5)
        W(4)=-(Cpq(22)*Q%C(6))+Cpq(21)*Q%C(7)-Cpq(23)*Q%C(7)
        W(5)=Cpq(22)*Q%C(8)-Cpq(24)*Q%C(8)+Cpq(23)*Q%C(9)
        W(6)=-(Cpq(25)*Q%C(9))-Spq(12)*Q%S(2)-Spq(17)*Q%S(4)
        W(7)=Spq(16)*Q%S(5)-Spq(18)*Q%S(5)-Spq(23)*Q%S(7)
        W(8)=Spq(22)*Q%S(8)-Spq(24)*Q%S(8)+Spq(23)*Q%S(9)-Spq(25)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(1)=-(Q%C(0)*Spq(7))-Q%C(1)*Spq(11)-Q%C(2)*Spq(12)-Q%C(3)*Spq(16)
        W(2)=-(Q%C(4)*Spq(17))-Q%C(5)*(Spq(16)+Spq(18))
        W(3)=-(Q%C(6)*Spq(22))-Q%C(7)*Spq(23)
        W(4)=-(Q%C(8)*(Spq(22)+Spq(24)))-Q%C(9)*(Spq(23)+Spq(25))
        W(5)=(Cpq(10)+Cpq(12))*Q%S(2)+(Cpq(15)+Cpq(17))*Q%S(4)
        W(6)=(Cpq(16)+Cpq(18))*Q%S(5)+(Cpq(21)+Cpq(23))*Q%S(7)
        W(7)=(Cpq(22)+Cpq(24))*Q%S(8)+(Cpq(23)+Cpq(25))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        SPKetC(7)=SPKetC(7)+CTwo*Two
        SPKetS(7)=SPKetS(7)+STwo*Two
!       -------------------------------- Lp=3, Mp = 2-------------------------------------
        W(1)=-(Cpq(8)*Q%C(0))-Cpq(12)*Q%C(1)+Cpq(11)*Q%C(2)
        W(2)=-(Cpq(13)*Q%C(2))-Cpq(17)*Q%C(3)+Cpq(16)*Q%C(4)
        W(3)=-(Cpq(18)*Q%C(4))-Cpq(15)*Q%C(5)-Cpq(19)*Q%C(5)
        W(4)=-(Cpq(23)*Q%C(6))+Cpq(22)*Q%C(7)-Cpq(24)*Q%C(7)-Cpq(21)*Q%C(8)
        W(5)=-(Cpq(25)*Q%C(8))-Cpq(22)*Q%C(9)-Cpq(26)*Q%C(9)
        W(6)=-(Spq(11)*Q%S(2))-Spq(13)*Q%S(2)-Spq(16)*Q%S(4)
        W(7)=-(Spq(18)*Q%S(4))-Spq(19)*Q%S(5)-Spq(22)*Q%S(7)
        W(8)=-(Spq(24)*Q%S(7))-Spq(25)*Q%S(8)-Spq(22)*Q%S(9)-Spq(26)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(1)=-(Q%C(0)*Spq(8))-Q%C(1)*Spq(12)
        W(2)=Q%C(2)*(Spq(11)-Spq(13))-Q%C(3)*Spq(17)
        W(3)=Q%C(4)*(Spq(16)-Spq(18))-Q%C(5)*Spq(19)
        W(4)=-(Q%C(6)*Spq(23))+Q%C(7)*(Spq(22)-Spq(24))
        W(5)=-(Q%C(8)*Spq(25))+Q%C(9)*(Spq(22)-Spq(26))
        W(6)=(Cpq(11)+Cpq(13))*Q%S(2)+(Cpq(16)+Cpq(18))*Q%S(4)
        W(7)=(-Cpq(15)+Cpq(19))*Q%S(5)+(Cpq(22)+Cpq(24))*Q%S(7)
        W(8)=(-Cpq(21)+Cpq(25))*Q%S(8)+(-Cpq(22)+Cpq(26))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        SPKetC(8)=SPKetC(8)+CTwo*Two
        SPKetS(8)=SPKetS(8)+STwo*Two
!       -------------------------------- Lp=3, Mp = 3-------------------------------------
        W(1)=-(Cpq(9)*Q%C(0))-Cpq(13)*Q%C(1)+Cpq(12)*Q%C(2)
        W(2)=-(Cpq(14)*Q%C(2))-Cpq(18)*Q%C(3)+Cpq(17)*Q%C(4)
        W(3)=-(Cpq(19)*Q%C(4))-Cpq(16)*Q%C(5)-Cpq(20)*Q%C(5)
        W(4)=-(Cpq(24)*Q%C(6))+Cpq(23)*Q%C(7)-Cpq(25)*Q%C(7)-Cpq(22)*Q%C(8)
        W(5)=-(Cpq(26)*Q%C(8))+Cpq(21)*Q%C(9)-Cpq(27)*Q%C(9)
        W(6)=-(Spq(12)*Q%S(2))-Spq(14)*Q%S(2)-Spq(17)*Q%S(4)-Spq(19)*Q%S(4)
        W(7)=Spq(16)*Q%S(5)-Spq(20)*Q%S(5)-Spq(23)*Q%S(7)
        W(8)=-(Spq(25)*Q%S(7))+Spq(22)*Q%S(8)-Spq(26)*Q%S(8)-Spq(27)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        W(1)=-(Q%C(0)*Spq(9))-Q%C(1)*Spq(13)
        W(2)=Q%C(2)*(Spq(12)-Spq(14))-Q%C(3)*Spq(18)
        W(3)=Q%C(4)*(Spq(17)-Spq(19))-Q%C(5)*(Spq(16)+Spq(20))
        W(4)=-(Q%C(6)*Spq(24))+Q%C(7)*(Spq(23)-Spq(25))
        W(5)=-(Q%C(8)*(Spq(22)+Spq(26)))-Q%C(9)*Spq(27)
        W(6)=(Cpq(12)+Cpq(14))*Q%S(2)+(Cpq(17)+Cpq(19))*Q%S(4)
        W(7)=(-Cpq(16)+Cpq(20))*Q%S(5)+(Cpq(23)+Cpq(25))*Q%S(7)
        W(8)=(-Cpq(22)+Cpq(26))*Q%S(8)+(Cpq(21)+Cpq(27))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)+W(8)
        SPKetC(9)=SPKetC(9)+CTwo*Two
        SPKetS(9)=SPKetS(9)+STwo*Two
!       -------------------------------- Lp=4, Mp = 0-------------------------------------
        COne=Cpq(10)*Q%C(0)+Cpq(15)*Q%C(1)+Cpq(21)*Q%C(3)+Cpq(28)*Q%C(6)
        W(1)=Cpq(16)*Q%C(2)+Cpq(22)*Q%C(4)+Cpq(23)*Q%C(5)
        W(2)=Cpq(29)*Q%C(7)+Cpq(30)*Q%C(8)+Cpq(31)*Q%C(9)
        W(3)=Spq(16)*Q%S(2)+Spq(22)*Q%S(4)+Spq(23)*Q%S(5)
        W(4)=Spq(29)*Q%S(7)+Spq(30)*Q%S(8)+Spq(31)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)
        SPKetC(10)=SPKetC(10)+COne+CTwo*Two
!       -------------------------------- Lp=4, Mp = 1-------------------------------------
        W(1)=Cpq(11)*Q%C(0)+Cpq(16)*Q%C(1)+(-Cpq(15)+Cpq(17))*Q%C(2)+Cpq(22)*Q%C(3)
        W(2)=(-Cpq(21)+Cpq(23))*Q%C(4)+(-Cpq(22)+Cpq(24))*Q%C(5)
        W(3)=Cpq(29)*Q%C(6)+(-Cpq(28)+Cpq(30))*Q%C(7)
        W(4)=(-Cpq(29)+Cpq(31))*Q%C(8)+(-Cpq(30)+Cpq(32))*Q%C(9)
        W(5)=Spq(17)*Q%S(2)+Spq(23)*Q%S(4)
        W(6)=(-Spq(22)+Spq(24))*Q%S(5)+Spq(30)*Q%S(7)
        W(7)=(-Spq(29)+Spq(31))*Q%S(8)+(-Spq(30)+Spq(32))*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        W(1)=Q%C(0)*Spq(11)+Q%C(1)*Spq(16)+Q%C(2)*Spq(17)+Q%C(3)*Spq(22)
        W(2)=Q%C(4)*Spq(23)+Q%C(5)*(Spq(22)+Spq(24))+Q%C(6)*Spq(29)+Q%C(7)*Spq(30)
        W(3)=Q%C(8)*(Spq(29)+Spq(31))+Q%C(9)*(Spq(30)+Spq(32))
        W(4)=-((Cpq(15)+Cpq(17))*Q%S(2))-(Cpq(21)+Cpq(23))*Q%S(4)
        W(5)=-((Cpq(22)+Cpq(24))*Q%S(5))-(Cpq(28)+Cpq(30))*Q%S(7)
        W(6)=-((Cpq(29)+Cpq(31))*Q%S(8))-(Cpq(30)+Cpq(32))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)
        SPKetC(11)=SPKetC(11)+CTwo*Two
        SPKetS(11)=SPKetS(11)+STwo*Two
!       -------------------------------- Lp=4, Mp = 2-------------------------------------
        W(1)=Cpq(12)*Q%C(0)+Cpq(17)*Q%C(1)+(-Cpq(16)+Cpq(18))*Q%C(2)+Cpq(23)*Q%C(3)
        W(2)=(-Cpq(22)+Cpq(24))*Q%C(4)+(Cpq(21)+Cpq(25))*Q%C(5)
        W(3)=Cpq(30)*Q%C(6)+(-Cpq(29)+Cpq(31))*Q%C(7)
        W(4)=(Cpq(28)+Cpq(32))*Q%C(8)+(Cpq(29)+Cpq(33))*Q%C(9)
        W(5)=(Spq(16)+Spq(18))*Q%S(2)+(Spq(22)+Spq(24))*Q%S(4)
        W(6)=Spq(25)*Q%S(5)+(Spq(29)+Spq(31))*Q%S(7)+Spq(32)*Q%S(8)+(Spq(29)+Spq(33))*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)
        W(1)=Q%C(0)*Spq(12)+Q%C(1)*Spq(17)+Q%C(2)*(-Spq(16)+Spq(18))+Q%C(3)*Spq(23)
        W(2)=Q%C(4)*(-Spq(22)+Spq(24))+Q%C(5)*Spq(25)
        W(3)=Q%C(6)*Spq(30)+Q%C(7)*(-Spq(29)+Spq(31))
        W(4)=Q%C(8)*Spq(32)+Q%C(9)*(-Spq(29)+Spq(33))
        W(5)=-((Cpq(16)+Cpq(18))*Q%S(2))-(Cpq(22)+Cpq(24))*Q%S(4)
        W(6)=(Cpq(21)-Cpq(25))*Q%S(5)-(Cpq(29)+Cpq(31))*Q%S(7)
        W(7)=(Cpq(28)-Cpq(32))*Q%S(8)+(Cpq(29)-Cpq(33))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        SPKetC(12)=SPKetC(12)+CTwo*Two
        SPKetS(12)=SPKetS(12)+STwo*Two
!       -------------------------------- Lp=4, Mp = 3-------------------------------------
        W(1)=Cpq(13)*Q%C(0)+Cpq(18)*Q%C(1)+(-Cpq(17)+Cpq(19))*Q%C(2)+Cpq(24)*Q%C(3)
        W(2)=(-Cpq(23)+Cpq(25))*Q%C(4)+(Cpq(22)+Cpq(26))*Q%C(5)
        W(3)=Cpq(31)*Q%C(6)+(-Cpq(30)+Cpq(32))*Q%C(7)
        W(4)=(Cpq(29)+Cpq(33))*Q%C(8)+(-Cpq(28)+Cpq(34))*Q%C(9)
        W(5)=(Spq(17)+Spq(19))*Q%S(2)+(Spq(23)+Spq(25))*Q%S(4)
        W(6)=(-Spq(22)+Spq(26))*Q%S(5)+(Spq(30)+Spq(32))*Q%S(7)
        W(7)=(-Spq(29)+Spq(33))*Q%S(8)+Spq(34)*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        W(1)=Q%C(0)*Spq(13)+Q%C(1)*Spq(18)+Q%C(2)*(-Spq(17)+Spq(19))+Q%C(3)*Spq(24)
        W(2)=Q%C(4)*(-Spq(23)+Spq(25))+Q%C(5)*(Spq(22)+Spq(26))
        W(3)=Q%C(6)*Spq(31)+Q%C(7)*(-Spq(30)+Spq(32))
        W(4)=Q%C(8)*(Spq(29)+Spq(33))+Q%C(9)*Spq(34)
        W(5)=-((Cpq(17)+Cpq(19))*Q%S(2))-(Cpq(23)+Cpq(25))*Q%S(4)
        W(6)=(Cpq(22)-Cpq(26))*Q%S(5)-(Cpq(30)+Cpq(32))*Q%S(7)
        W(7)=(Cpq(29)-Cpq(33))*Q%S(8)-(Cpq(28)+Cpq(34))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        SPKetC(13)=SPKetC(13)+CTwo*Two
        SPKetS(13)=SPKetS(13)+STwo*Two
!       -------------------------------- Lp=4, Mp = 4-------------------------------------
        W(1)=Cpq(14)*Q%C(0)+Cpq(19)*Q%C(1)+(-Cpq(18)+Cpq(20))*Q%C(2)+Cpq(25)*Q%C(3)
        W(2)=(-Cpq(24)+Cpq(26))*Q%C(4)+(Cpq(23)+Cpq(27))*Q%C(5)
        W(3)=Cpq(32)*Q%C(6)+(-Cpq(31)+Cpq(33))*Q%C(7)
        W(4)=(Cpq(30)+Cpq(34))*Q%C(8)+(-Cpq(29)+Cpq(35))*Q%C(9)
        W(5)=(Spq(18)+Spq(20))*Q%S(2)+(Spq(24)+Spq(26))*Q%S(4)
        W(6)=(-Spq(23)+Spq(27))*Q%S(5)+(Spq(31)+Spq(33))*Q%S(7)
        W(7)=(-Spq(30)+Spq(34))*Q%S(8)+(Spq(29)+Spq(35))*Q%S(9)
        CTwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        W(1)=Q%C(0)*Spq(14)+Q%C(1)*Spq(19)+Q%C(2)*(-Spq(18)+Spq(20))+Q%C(3)*Spq(25)
        W(2)=Q%C(4)*(-Spq(24)+Spq(26))+Q%C(5)*(Spq(23)+Spq(27))
        W(3)=Q%C(6)*Spq(32)+Q%C(7)*(-Spq(31)+Spq(33))
        W(4)=Q%C(8)*(Spq(30)+Spq(34))+Q%C(9)*(-Spq(29)+Spq(35))
        W(5)=-((Cpq(18)+Cpq(20))*Q%S(2))-(Cpq(24)+Cpq(26))*Q%S(4)
        W(6)=(Cpq(23)-Cpq(27))*Q%S(5)-(Cpq(31)+Cpq(33))*Q%S(7)
        W(7)=(Cpq(30)-Cpq(34))*Q%S(8)-(Cpq(29)+Cpq(35))*Q%S(9)
        STwo=W(1)+W(2)+W(3)+W(4)+W(5)+W(6)+W(7)
        SPKetC(14)=SPKetC(14)+CTwo*Two
        SPKetS(14)=SPKetS(14)+STwo*Two
      CASE(404)
         CALL CTraX404(Q)
      CASE(405)
         CALL CTraX405(Q)
      CASE(406)
         CALL CTraX406(Q)
      CASE(407)
         CALL CTraX407(Q)
      CASE(500)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)
        SPKetC(0)=SPKetC(0)+COne
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))
        SPKetC(1)=SPKetC(1)+COne
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        CTwo=-(Cpq(2)*Q%C(0))
        STwo=-(Q%C(0)*Spq(2))
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)
        SPKetC(3)=SPKetC(3)+COne
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        CTwo=Cpq(4)*Q%C(0)
        STwo=Q%C(0)*Spq(4)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        CTwo=Cpq(5)*Q%C(0)
        STwo=Q%C(0)*Spq(5)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
!       -------------------------------- Lp=3, Mp = 0-------------------------------------
        COne=-(Cpq(6)*Q%C(0))
        SPKetC(6)=SPKetC(6)+COne
!       -------------------------------- Lp=3, Mp = 1-------------------------------------
        CTwo=-(Cpq(7)*Q%C(0))
        STwo=-(Q%C(0)*Spq(7))
        SPKetC(7)=SPKetC(7)+CTwo*Two
        SPKetS(7)=SPKetS(7)+STwo*Two
!       -------------------------------- Lp=3, Mp = 2-------------------------------------
        CTwo=-(Cpq(8)*Q%C(0))
        STwo=-(Q%C(0)*Spq(8))
        SPKetC(8)=SPKetC(8)+CTwo*Two
        SPKetS(8)=SPKetS(8)+STwo*Two
!       -------------------------------- Lp=3, Mp = 3-------------------------------------
        CTwo=-(Cpq(9)*Q%C(0))
        STwo=-(Q%C(0)*Spq(9))
        SPKetC(9)=SPKetC(9)+CTwo*Two
        SPKetS(9)=SPKetS(9)+STwo*Two
!       -------------------------------- Lp=4, Mp = 0-------------------------------------
        COne=Cpq(10)*Q%C(0)
        SPKetC(10)=SPKetC(10)+COne
!       -------------------------------- Lp=4, Mp = 1-------------------------------------
        CTwo=Cpq(11)*Q%C(0)
        STwo=Q%C(0)*Spq(11)
        SPKetC(11)=SPKetC(11)+CTwo*Two
        SPKetS(11)=SPKetS(11)+STwo*Two
!       -------------------------------- Lp=4, Mp = 2-------------------------------------
        CTwo=Cpq(12)*Q%C(0)
        STwo=Q%C(0)*Spq(12)
        SPKetC(12)=SPKetC(12)+CTwo*Two
        SPKetS(12)=SPKetS(12)+STwo*Two
!       -------------------------------- Lp=4, Mp = 3-------------------------------------
        CTwo=Cpq(13)*Q%C(0)
        STwo=Q%C(0)*Spq(13)
        SPKetC(13)=SPKetC(13)+CTwo*Two
        SPKetS(13)=SPKetS(13)+STwo*Two
!       -------------------------------- Lp=4, Mp = 4-------------------------------------
        CTwo=Cpq(14)*Q%C(0)
        STwo=Q%C(0)*Spq(14)
        SPKetC(14)=SPKetC(14)+CTwo*Two
        SPKetS(14)=SPKetS(14)+STwo*Two
!       -------------------------------- Lp=5, Mp = 0-------------------------------------
        COne=-(Cpq(15)*Q%C(0))
        SPKetC(15)=SPKetC(15)+COne
!       -------------------------------- Lp=5, Mp = 1-------------------------------------
        CTwo=-(Cpq(16)*Q%C(0))
        STwo=-(Q%C(0)*Spq(16))
        SPKetC(16)=SPKetC(16)+CTwo*Two
        SPKetS(16)=SPKetS(16)+STwo*Two
!       -------------------------------- Lp=5, Mp = 2-------------------------------------
        CTwo=-(Cpq(17)*Q%C(0))
        STwo=-(Q%C(0)*Spq(17))
        SPKetC(17)=SPKetC(17)+CTwo*Two
        SPKetS(17)=SPKetS(17)+STwo*Two
!       -------------------------------- Lp=5, Mp = 3-------------------------------------
        CTwo=-(Cpq(18)*Q%C(0))
        STwo=-(Q%C(0)*Spq(18))
        SPKetC(18)=SPKetC(18)+CTwo*Two
        SPKetS(18)=SPKetS(18)+STwo*Two
!       -------------------------------- Lp=5, Mp = 4-------------------------------------
        CTwo=-(Cpq(19)*Q%C(0))
        STwo=-(Q%C(0)*Spq(19))
        SPKetC(19)=SPKetC(19)+CTwo*Two
        SPKetS(19)=SPKetS(19)+STwo*Two
!       -------------------------------- Lp=5, Mp = 5-------------------------------------
        CTwo=-(Cpq(20)*Q%C(0))
        STwo=-(Q%C(0)*Spq(20))
        SPKetC(20)=SPKetC(20)+CTwo*Two
        SPKetS(20)=SPKetS(20)+STwo*Two
      CASE(501)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)
        CTwo=Cpq(2)*Q%C(2)+Spq(2)*Q%S(2)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)
        CTwo=-(Cpq(4)*Q%C(2))-Spq(4)*Q%S(2)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)
        W(2)=Cpq(3)*Q%C(2)-Cpq(5)*Q%C(2)-Spq(5)*Q%S(2)
        CTwo=W(1)+W(2)
        STwo=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)+(Cpq(3)+Cpq(5))*Q%S(2)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)+Cpq(6)*Q%C(1)
        CTwo=Cpq(7)*Q%C(2)+Spq(7)*Q%S(2)
        SPKetC(3)=SPKetC(3)+COne+CTwo*Two
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        CTwo=Cpq(4)*Q%C(0)+Cpq(7)*Q%C(1)+(-Cpq(6)+Cpq(8))*Q%C(2)+Spq(8)*Q%S(2)
        W(1)=Q%C(0)*Spq(4)+Q%C(1)*Spq(7)
        W(2)=Q%C(2)*Spq(8)-Cpq(6)*Q%S(2)-Cpq(8)*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        W(1)=Cpq(5)*Q%C(0)+Cpq(8)*Q%C(1)
        W(2)=(-Cpq(7)+Cpq(9))*Q%C(2)+(Spq(7)+Spq(9))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(5)+Q%C(1)*Spq(8)
        W(2)=Q%C(2)*(-Spq(7)+Spq(9))-(Cpq(7)+Cpq(9))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
!       -------------------------------- Lp=3, Mp = 0-------------------------------------
        COne=-(Cpq(6)*Q%C(0))-Cpq(10)*Q%C(1)
        CTwo=-(Cpq(11)*Q%C(2))-Spq(11)*Q%S(2)
        SPKetC(6)=SPKetC(6)+COne+CTwo*Two
!       -------------------------------- Lp=3, Mp = 1-------------------------------------
        W(1)=-(Cpq(7)*Q%C(0))-Cpq(11)*Q%C(1)
        W(2)=Cpq(10)*Q%C(2)-Cpq(12)*Q%C(2)-Spq(12)*Q%S(2)
        CTwo=W(1)+W(2)
        STwo=-(Q%C(0)*Spq(7))-Q%C(1)*Spq(11)-Q%C(2)*Spq(12)+(Cpq(10)+Cpq(12))*Q%S(2)
        SPKetC(7)=SPKetC(7)+CTwo*Two
        SPKetS(7)=SPKetS(7)+STwo*Two
!       -------------------------------- Lp=3, Mp = 2-------------------------------------
        W(1)=-(Cpq(8)*Q%C(0))-Cpq(12)*Q%C(1)+Cpq(11)*Q%C(2)
        W(2)=-(Cpq(13)*Q%C(2))-Spq(11)*Q%S(2)-Spq(13)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(8))-Q%C(1)*Spq(12)
        W(2)=Q%C(2)*(Spq(11)-Spq(13))+(Cpq(11)+Cpq(13))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(8)=SPKetC(8)+CTwo*Two
        SPKetS(8)=SPKetS(8)+STwo*Two
!       -------------------------------- Lp=3, Mp = 3-------------------------------------
        W(1)=-(Cpq(9)*Q%C(0))-Cpq(13)*Q%C(1)+Cpq(12)*Q%C(2)
        W(2)=-(Cpq(14)*Q%C(2))-Spq(12)*Q%S(2)-Spq(14)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(9))-Q%C(1)*Spq(13)
        W(2)=Q%C(2)*(Spq(12)-Spq(14))+(Cpq(12)+Cpq(14))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(9)=SPKetC(9)+CTwo*Two
        SPKetS(9)=SPKetS(9)+STwo*Two
!       -------------------------------- Lp=4, Mp = 0-------------------------------------
        COne=Cpq(10)*Q%C(0)+Cpq(15)*Q%C(1)
        CTwo=Cpq(16)*Q%C(2)+Spq(16)*Q%S(2)
        SPKetC(10)=SPKetC(10)+COne+CTwo*Two
!       -------------------------------- Lp=4, Mp = 1-------------------------------------
        CTwo=Cpq(11)*Q%C(0)+Cpq(16)*Q%C(1)+(-Cpq(15)+Cpq(17))*Q%C(2)+Spq(17)*Q%S(2)
        W(1)=Q%C(0)*Spq(11)+Q%C(1)*Spq(16)
        W(2)=Q%C(2)*Spq(17)-Cpq(15)*Q%S(2)-Cpq(17)*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(11)=SPKetC(11)+CTwo*Two
        SPKetS(11)=SPKetS(11)+STwo*Two
!       -------------------------------- Lp=4, Mp = 2-------------------------------------
        W(1)=Cpq(12)*Q%C(0)+Cpq(17)*Q%C(1)
        W(2)=(-Cpq(16)+Cpq(18))*Q%C(2)+(Spq(16)+Spq(18))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(12)+Q%C(1)*Spq(17)
        W(2)=Q%C(2)*(-Spq(16)+Spq(18))-(Cpq(16)+Cpq(18))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(12)=SPKetC(12)+CTwo*Two
        SPKetS(12)=SPKetS(12)+STwo*Two
!       -------------------------------- Lp=4, Mp = 3-------------------------------------
        W(1)=Cpq(13)*Q%C(0)+Cpq(18)*Q%C(1)
        W(2)=(-Cpq(17)+Cpq(19))*Q%C(2)+(Spq(17)+Spq(19))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(13)+Q%C(1)*Spq(18)
        W(2)=Q%C(2)*(-Spq(17)+Spq(19))-(Cpq(17)+Cpq(19))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(13)=SPKetC(13)+CTwo*Two
        SPKetS(13)=SPKetS(13)+STwo*Two
!       -------------------------------- Lp=4, Mp = 4-------------------------------------
        W(1)=Cpq(14)*Q%C(0)+Cpq(19)*Q%C(1)
        W(2)=(-Cpq(18)+Cpq(20))*Q%C(2)+(Spq(18)+Spq(20))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(14)+Q%C(1)*Spq(19)
        W(2)=Q%C(2)*(-Spq(18)+Spq(20))-(Cpq(18)+Cpq(20))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(14)=SPKetC(14)+CTwo*Two
        SPKetS(14)=SPKetS(14)+STwo*Two
!       -------------------------------- Lp=5, Mp = 0-------------------------------------
        COne=-(Cpq(15)*Q%C(0))-Cpq(21)*Q%C(1)
        CTwo=-(Cpq(22)*Q%C(2))-Spq(22)*Q%S(2)
        SPKetC(15)=SPKetC(15)+COne+CTwo*Two
!       -------------------------------- Lp=5, Mp = 1-------------------------------------
        W(1)=-(Cpq(16)*Q%C(0))-Cpq(22)*Q%C(1)
        W(2)=Cpq(21)*Q%C(2)-Cpq(23)*Q%C(2)-Spq(23)*Q%S(2)
        CTwo=W(1)+W(2)
        STwo=-(Q%C(0)*Spq(16))-Q%C(1)*Spq(22)-Q%C(2)*Spq(23)+(Cpq(21)+Cpq(23))*Q%S(2)
        SPKetC(16)=SPKetC(16)+CTwo*Two
        SPKetS(16)=SPKetS(16)+STwo*Two
!       -------------------------------- Lp=5, Mp = 2-------------------------------------
        W(1)=-(Cpq(17)*Q%C(0))-Cpq(23)*Q%C(1)+Cpq(22)*Q%C(2)
        W(2)=-(Cpq(24)*Q%C(2))-Spq(22)*Q%S(2)-Spq(24)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(17))-Q%C(1)*Spq(23)
        W(2)=Q%C(2)*(Spq(22)-Spq(24))+(Cpq(22)+Cpq(24))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(17)=SPKetC(17)+CTwo*Two
        SPKetS(17)=SPKetS(17)+STwo*Two
!       -------------------------------- Lp=5, Mp = 3-------------------------------------
        W(1)=-(Cpq(18)*Q%C(0))-Cpq(24)*Q%C(1)+Cpq(23)*Q%C(2)
        W(2)=-(Cpq(25)*Q%C(2))-Spq(23)*Q%S(2)-Spq(25)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(18))-Q%C(1)*Spq(24)
        W(2)=Q%C(2)*(Spq(23)-Spq(25))+(Cpq(23)+Cpq(25))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(18)=SPKetC(18)+CTwo*Two
        SPKetS(18)=SPKetS(18)+STwo*Two
!       -------------------------------- Lp=5, Mp = 4-------------------------------------
        W(1)=-(Cpq(19)*Q%C(0))-Cpq(25)*Q%C(1)+Cpq(24)*Q%C(2)
        W(2)=-(Cpq(26)*Q%C(2))-Spq(24)*Q%S(2)-Spq(26)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(19))-Q%C(1)*Spq(25)
        W(2)=Q%C(2)*(Spq(24)-Spq(26))+(Cpq(24)+Cpq(26))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(19)=SPKetC(19)+CTwo*Two
        SPKetS(19)=SPKetS(19)+STwo*Two
!       -------------------------------- Lp=5, Mp = 5-------------------------------------
        W(1)=-(Cpq(20)*Q%C(0))-Cpq(26)*Q%C(1)+Cpq(25)*Q%C(2)
        W(2)=-(Cpq(27)*Q%C(2))-Spq(25)*Q%S(2)-Spq(27)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(20))-Q%C(1)*Spq(26)
        W(2)=Q%C(2)*(Spq(25)-Spq(27))+(Cpq(25)+Cpq(27))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(20)=SPKetC(20)+CTwo*Two
        SPKetS(20)=SPKetS(20)+STwo*Two
      CASE(502)
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
      CASE(503)
         CALL CTraX503(Q)
      CASE(504)
         CALL CTraX504(Q)
      CASE(505)
         CALL CTraX505(Q)
      CASE(506)
         CALL CTraX506(Q)
      CASE(507)
         CALL CTraX507(Q)
      CASE(600)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)
        SPKetC(0)=SPKetC(0)+COne
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))
        SPKetC(1)=SPKetC(1)+COne
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        CTwo=-(Cpq(2)*Q%C(0))
        STwo=-(Q%C(0)*Spq(2))
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)
        SPKetC(3)=SPKetC(3)+COne
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        CTwo=Cpq(4)*Q%C(0)
        STwo=Q%C(0)*Spq(4)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        CTwo=Cpq(5)*Q%C(0)
        STwo=Q%C(0)*Spq(5)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
!       -------------------------------- Lp=3, Mp = 0-------------------------------------
        COne=-(Cpq(6)*Q%C(0))
        SPKetC(6)=SPKetC(6)+COne
!       -------------------------------- Lp=3, Mp = 1-------------------------------------
        CTwo=-(Cpq(7)*Q%C(0))
        STwo=-(Q%C(0)*Spq(7))
        SPKetC(7)=SPKetC(7)+CTwo*Two
        SPKetS(7)=SPKetS(7)+STwo*Two
!       -------------------------------- Lp=3, Mp = 2-------------------------------------
        CTwo=-(Cpq(8)*Q%C(0))
        STwo=-(Q%C(0)*Spq(8))
        SPKetC(8)=SPKetC(8)+CTwo*Two
        SPKetS(8)=SPKetS(8)+STwo*Two
!       -------------------------------- Lp=3, Mp = 3-------------------------------------
        CTwo=-(Cpq(9)*Q%C(0))
        STwo=-(Q%C(0)*Spq(9))
        SPKetC(9)=SPKetC(9)+CTwo*Two
        SPKetS(9)=SPKetS(9)+STwo*Two
!       -------------------------------- Lp=4, Mp = 0-------------------------------------
        COne=Cpq(10)*Q%C(0)
        SPKetC(10)=SPKetC(10)+COne
!       -------------------------------- Lp=4, Mp = 1-------------------------------------
        CTwo=Cpq(11)*Q%C(0)
        STwo=Q%C(0)*Spq(11)
        SPKetC(11)=SPKetC(11)+CTwo*Two
        SPKetS(11)=SPKetS(11)+STwo*Two
!       -------------------------------- Lp=4, Mp = 2-------------------------------------
        CTwo=Cpq(12)*Q%C(0)
        STwo=Q%C(0)*Spq(12)
        SPKetC(12)=SPKetC(12)+CTwo*Two
        SPKetS(12)=SPKetS(12)+STwo*Two
!       -------------------------------- Lp=4, Mp = 3-------------------------------------
        CTwo=Cpq(13)*Q%C(0)
        STwo=Q%C(0)*Spq(13)
        SPKetC(13)=SPKetC(13)+CTwo*Two
        SPKetS(13)=SPKetS(13)+STwo*Two
!       -------------------------------- Lp=4, Mp = 4-------------------------------------
        CTwo=Cpq(14)*Q%C(0)
        STwo=Q%C(0)*Spq(14)
        SPKetC(14)=SPKetC(14)+CTwo*Two
        SPKetS(14)=SPKetS(14)+STwo*Two
!       -------------------------------- Lp=5, Mp = 0-------------------------------------
        COne=-(Cpq(15)*Q%C(0))
        SPKetC(15)=SPKetC(15)+COne
!       -------------------------------- Lp=5, Mp = 1-------------------------------------
        CTwo=-(Cpq(16)*Q%C(0))
        STwo=-(Q%C(0)*Spq(16))
        SPKetC(16)=SPKetC(16)+CTwo*Two
        SPKetS(16)=SPKetS(16)+STwo*Two
!       -------------------------------- Lp=5, Mp = 2-------------------------------------
        CTwo=-(Cpq(17)*Q%C(0))
        STwo=-(Q%C(0)*Spq(17))
        SPKetC(17)=SPKetC(17)+CTwo*Two
        SPKetS(17)=SPKetS(17)+STwo*Two
!       -------------------------------- Lp=5, Mp = 3-------------------------------------
        CTwo=-(Cpq(18)*Q%C(0))
        STwo=-(Q%C(0)*Spq(18))
        SPKetC(18)=SPKetC(18)+CTwo*Two
        SPKetS(18)=SPKetS(18)+STwo*Two
!       -------------------------------- Lp=5, Mp = 4-------------------------------------
        CTwo=-(Cpq(19)*Q%C(0))
        STwo=-(Q%C(0)*Spq(19))
        SPKetC(19)=SPKetC(19)+CTwo*Two
        SPKetS(19)=SPKetS(19)+STwo*Two
!       -------------------------------- Lp=5, Mp = 5-------------------------------------
        CTwo=-(Cpq(20)*Q%C(0))
        STwo=-(Q%C(0)*Spq(20))
        SPKetC(20)=SPKetC(20)+CTwo*Two
        SPKetS(20)=SPKetS(20)+STwo*Two
!       -------------------------------- Lp=6, Mp = 0-------------------------------------
        COne=Cpq(21)*Q%C(0)
        SPKetC(21)=SPKetC(21)+COne
!       -------------------------------- Lp=6, Mp = 1-------------------------------------
        CTwo=Cpq(22)*Q%C(0)
        STwo=Q%C(0)*Spq(22)
        SPKetC(22)=SPKetC(22)+CTwo*Two
        SPKetS(22)=SPKetS(22)+STwo*Two
!       -------------------------------- Lp=6, Mp = 2-------------------------------------
        CTwo=Cpq(23)*Q%C(0)
        STwo=Q%C(0)*Spq(23)
        SPKetC(23)=SPKetC(23)+CTwo*Two
        SPKetS(23)=SPKetS(23)+STwo*Two
!       -------------------------------- Lp=6, Mp = 3-------------------------------------
        CTwo=Cpq(24)*Q%C(0)
        STwo=Q%C(0)*Spq(24)
        SPKetC(24)=SPKetC(24)+CTwo*Two
        SPKetS(24)=SPKetS(24)+STwo*Two
!       -------------------------------- Lp=6, Mp = 4-------------------------------------
        CTwo=Cpq(25)*Q%C(0)
        STwo=Q%C(0)*Spq(25)
        SPKetC(25)=SPKetC(25)+CTwo*Two
        SPKetS(25)=SPKetS(25)+STwo*Two
!       -------------------------------- Lp=6, Mp = 5-------------------------------------
        CTwo=Cpq(26)*Q%C(0)
        STwo=Q%C(0)*Spq(26)
        SPKetC(26)=SPKetC(26)+CTwo*Two
        SPKetS(26)=SPKetS(26)+STwo*Two
!       -------------------------------- Lp=6, Mp = 6-------------------------------------
        CTwo=Cpq(27)*Q%C(0)
        STwo=Q%C(0)*Spq(27)
        SPKetC(27)=SPKetC(27)+CTwo*Two
        SPKetS(27)=SPKetS(27)+STwo*Two
      CASE(601)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)+Cpq(1)*Q%C(1)
        CTwo=Cpq(2)*Q%C(2)+Spq(2)*Q%S(2)
        SPKetC(0)=SPKetC(0)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))-Cpq(3)*Q%C(1)
        CTwo=-(Cpq(4)*Q%C(2))-Spq(4)*Q%S(2)
        SPKetC(1)=SPKetC(1)+COne+CTwo*Two
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        W(1)=-(Cpq(2)*Q%C(0))-Cpq(4)*Q%C(1)
        W(2)=Cpq(3)*Q%C(2)-Cpq(5)*Q%C(2)-Spq(5)*Q%S(2)
        CTwo=W(1)+W(2)
        STwo=-(Q%C(0)*Spq(2))-Q%C(1)*Spq(4)-Q%C(2)*Spq(5)+(Cpq(3)+Cpq(5))*Q%S(2)
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)+Cpq(6)*Q%C(1)
        CTwo=Cpq(7)*Q%C(2)+Spq(7)*Q%S(2)
        SPKetC(3)=SPKetC(3)+COne+CTwo*Two
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        CTwo=Cpq(4)*Q%C(0)+Cpq(7)*Q%C(1)+(-Cpq(6)+Cpq(8))*Q%C(2)+Spq(8)*Q%S(2)
        W(1)=Q%C(0)*Spq(4)+Q%C(1)*Spq(7)
        W(2)=Q%C(2)*Spq(8)-Cpq(6)*Q%S(2)-Cpq(8)*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        W(1)=Cpq(5)*Q%C(0)+Cpq(8)*Q%C(1)
        W(2)=(-Cpq(7)+Cpq(9))*Q%C(2)+(Spq(7)+Spq(9))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(5)+Q%C(1)*Spq(8)
        W(2)=Q%C(2)*(-Spq(7)+Spq(9))-(Cpq(7)+Cpq(9))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
!       -------------------------------- Lp=3, Mp = 0-------------------------------------
        COne=-(Cpq(6)*Q%C(0))-Cpq(10)*Q%C(1)
        CTwo=-(Cpq(11)*Q%C(2))-Spq(11)*Q%S(2)
        SPKetC(6)=SPKetC(6)+COne+CTwo*Two
!       -------------------------------- Lp=3, Mp = 1-------------------------------------
        W(1)=-(Cpq(7)*Q%C(0))-Cpq(11)*Q%C(1)
        W(2)=Cpq(10)*Q%C(2)-Cpq(12)*Q%C(2)-Spq(12)*Q%S(2)
        CTwo=W(1)+W(2)
        STwo=-(Q%C(0)*Spq(7))-Q%C(1)*Spq(11)-Q%C(2)*Spq(12)+(Cpq(10)+Cpq(12))*Q%S(2)
        SPKetC(7)=SPKetC(7)+CTwo*Two
        SPKetS(7)=SPKetS(7)+STwo*Two
!       -------------------------------- Lp=3, Mp = 2-------------------------------------
        W(1)=-(Cpq(8)*Q%C(0))-Cpq(12)*Q%C(1)+Cpq(11)*Q%C(2)
        W(2)=-(Cpq(13)*Q%C(2))-Spq(11)*Q%S(2)-Spq(13)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(8))-Q%C(1)*Spq(12)
        W(2)=Q%C(2)*(Spq(11)-Spq(13))+(Cpq(11)+Cpq(13))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(8)=SPKetC(8)+CTwo*Two
        SPKetS(8)=SPKetS(8)+STwo*Two
!       -------------------------------- Lp=3, Mp = 3-------------------------------------
        W(1)=-(Cpq(9)*Q%C(0))-Cpq(13)*Q%C(1)+Cpq(12)*Q%C(2)
        W(2)=-(Cpq(14)*Q%C(2))-Spq(12)*Q%S(2)-Spq(14)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(9))-Q%C(1)*Spq(13)
        W(2)=Q%C(2)*(Spq(12)-Spq(14))+(Cpq(12)+Cpq(14))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(9)=SPKetC(9)+CTwo*Two
        SPKetS(9)=SPKetS(9)+STwo*Two
!       -------------------------------- Lp=4, Mp = 0-------------------------------------
        COne=Cpq(10)*Q%C(0)+Cpq(15)*Q%C(1)
        CTwo=Cpq(16)*Q%C(2)+Spq(16)*Q%S(2)
        SPKetC(10)=SPKetC(10)+COne+CTwo*Two
!       -------------------------------- Lp=4, Mp = 1-------------------------------------
        CTwo=Cpq(11)*Q%C(0)+Cpq(16)*Q%C(1)+(-Cpq(15)+Cpq(17))*Q%C(2)+Spq(17)*Q%S(2)
        W(1)=Q%C(0)*Spq(11)+Q%C(1)*Spq(16)
        W(2)=Q%C(2)*Spq(17)-Cpq(15)*Q%S(2)-Cpq(17)*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(11)=SPKetC(11)+CTwo*Two
        SPKetS(11)=SPKetS(11)+STwo*Two
!       -------------------------------- Lp=4, Mp = 2-------------------------------------
        W(1)=Cpq(12)*Q%C(0)+Cpq(17)*Q%C(1)
        W(2)=(-Cpq(16)+Cpq(18))*Q%C(2)+(Spq(16)+Spq(18))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(12)+Q%C(1)*Spq(17)
        W(2)=Q%C(2)*(-Spq(16)+Spq(18))-(Cpq(16)+Cpq(18))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(12)=SPKetC(12)+CTwo*Two
        SPKetS(12)=SPKetS(12)+STwo*Two
!       -------------------------------- Lp=4, Mp = 3-------------------------------------
        W(1)=Cpq(13)*Q%C(0)+Cpq(18)*Q%C(1)
        W(2)=(-Cpq(17)+Cpq(19))*Q%C(2)+(Spq(17)+Spq(19))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(13)+Q%C(1)*Spq(18)
        W(2)=Q%C(2)*(-Spq(17)+Spq(19))-(Cpq(17)+Cpq(19))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(13)=SPKetC(13)+CTwo*Two
        SPKetS(13)=SPKetS(13)+STwo*Two
!       -------------------------------- Lp=4, Mp = 4-------------------------------------
        W(1)=Cpq(14)*Q%C(0)+Cpq(19)*Q%C(1)
        W(2)=(-Cpq(18)+Cpq(20))*Q%C(2)+(Spq(18)+Spq(20))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(14)+Q%C(1)*Spq(19)
        W(2)=Q%C(2)*(-Spq(18)+Spq(20))-(Cpq(18)+Cpq(20))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(14)=SPKetC(14)+CTwo*Two
        SPKetS(14)=SPKetS(14)+STwo*Two
!       -------------------------------- Lp=5, Mp = 0-------------------------------------
        COne=-(Cpq(15)*Q%C(0))-Cpq(21)*Q%C(1)
        CTwo=-(Cpq(22)*Q%C(2))-Spq(22)*Q%S(2)
        SPKetC(15)=SPKetC(15)+COne+CTwo*Two
!       -------------------------------- Lp=5, Mp = 1-------------------------------------
        W(1)=-(Cpq(16)*Q%C(0))-Cpq(22)*Q%C(1)
        W(2)=Cpq(21)*Q%C(2)-Cpq(23)*Q%C(2)-Spq(23)*Q%S(2)
        CTwo=W(1)+W(2)
        STwo=-(Q%C(0)*Spq(16))-Q%C(1)*Spq(22)-Q%C(2)*Spq(23)+(Cpq(21)+Cpq(23))*Q%S(2)
        SPKetC(16)=SPKetC(16)+CTwo*Two
        SPKetS(16)=SPKetS(16)+STwo*Two
!       -------------------------------- Lp=5, Mp = 2-------------------------------------
        W(1)=-(Cpq(17)*Q%C(0))-Cpq(23)*Q%C(1)+Cpq(22)*Q%C(2)
        W(2)=-(Cpq(24)*Q%C(2))-Spq(22)*Q%S(2)-Spq(24)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(17))-Q%C(1)*Spq(23)
        W(2)=Q%C(2)*(Spq(22)-Spq(24))+(Cpq(22)+Cpq(24))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(17)=SPKetC(17)+CTwo*Two
        SPKetS(17)=SPKetS(17)+STwo*Two
!       -------------------------------- Lp=5, Mp = 3-------------------------------------
        W(1)=-(Cpq(18)*Q%C(0))-Cpq(24)*Q%C(1)+Cpq(23)*Q%C(2)
        W(2)=-(Cpq(25)*Q%C(2))-Spq(23)*Q%S(2)-Spq(25)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(18))-Q%C(1)*Spq(24)
        W(2)=Q%C(2)*(Spq(23)-Spq(25))+(Cpq(23)+Cpq(25))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(18)=SPKetC(18)+CTwo*Two
        SPKetS(18)=SPKetS(18)+STwo*Two
!       -------------------------------- Lp=5, Mp = 4-------------------------------------
        W(1)=-(Cpq(19)*Q%C(0))-Cpq(25)*Q%C(1)+Cpq(24)*Q%C(2)
        W(2)=-(Cpq(26)*Q%C(2))-Spq(24)*Q%S(2)-Spq(26)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(19))-Q%C(1)*Spq(25)
        W(2)=Q%C(2)*(Spq(24)-Spq(26))+(Cpq(24)+Cpq(26))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(19)=SPKetC(19)+CTwo*Two
        SPKetS(19)=SPKetS(19)+STwo*Two
!       -------------------------------- Lp=5, Mp = 5-------------------------------------
        W(1)=-(Cpq(20)*Q%C(0))-Cpq(26)*Q%C(1)+Cpq(25)*Q%C(2)
        W(2)=-(Cpq(27)*Q%C(2))-Spq(25)*Q%S(2)-Spq(27)*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=-(Q%C(0)*Spq(20))-Q%C(1)*Spq(26)
        W(2)=Q%C(2)*(Spq(25)-Spq(27))+(Cpq(25)+Cpq(27))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(20)=SPKetC(20)+CTwo*Two
        SPKetS(20)=SPKetS(20)+STwo*Two
!       -------------------------------- Lp=6, Mp = 0-------------------------------------
        COne=Cpq(21)*Q%C(0)+Cpq(28)*Q%C(1)
        CTwo=Cpq(29)*Q%C(2)+Spq(29)*Q%S(2)
        SPKetC(21)=SPKetC(21)+COne+CTwo*Two
!       -------------------------------- Lp=6, Mp = 1-------------------------------------
        CTwo=Cpq(22)*Q%C(0)+Cpq(29)*Q%C(1)+(-Cpq(28)+Cpq(30))*Q%C(2)+Spq(30)*Q%S(2)
        W(1)=Q%C(0)*Spq(22)+Q%C(1)*Spq(29)
        W(2)=Q%C(2)*Spq(30)-Cpq(28)*Q%S(2)-Cpq(30)*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(22)=SPKetC(22)+CTwo*Two
        SPKetS(22)=SPKetS(22)+STwo*Two
!       -------------------------------- Lp=6, Mp = 2-------------------------------------
        W(1)=Cpq(23)*Q%C(0)+Cpq(30)*Q%C(1)
        W(2)=(-Cpq(29)+Cpq(31))*Q%C(2)+(Spq(29)+Spq(31))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(23)+Q%C(1)*Spq(30)
        W(2)=Q%C(2)*(-Spq(29)+Spq(31))-(Cpq(29)+Cpq(31))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(23)=SPKetC(23)+CTwo*Two
        SPKetS(23)=SPKetS(23)+STwo*Two
!       -------------------------------- Lp=6, Mp = 3-------------------------------------
        W(1)=Cpq(24)*Q%C(0)+Cpq(31)*Q%C(1)
        W(2)=(-Cpq(30)+Cpq(32))*Q%C(2)+(Spq(30)+Spq(32))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(24)+Q%C(1)*Spq(31)
        W(2)=Q%C(2)*(-Spq(30)+Spq(32))-(Cpq(30)+Cpq(32))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(24)=SPKetC(24)+CTwo*Two
        SPKetS(24)=SPKetS(24)+STwo*Two
!       -------------------------------- Lp=6, Mp = 4-------------------------------------
        W(1)=Cpq(25)*Q%C(0)+Cpq(32)*Q%C(1)
        W(2)=(-Cpq(31)+Cpq(33))*Q%C(2)+(Spq(31)+Spq(33))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(25)+Q%C(1)*Spq(32)
        W(2)=Q%C(2)*(-Spq(31)+Spq(33))-(Cpq(31)+Cpq(33))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(25)=SPKetC(25)+CTwo*Two
        SPKetS(25)=SPKetS(25)+STwo*Two
!       -------------------------------- Lp=6, Mp = 5-------------------------------------
        W(1)=Cpq(26)*Q%C(0)+Cpq(33)*Q%C(1)
        W(2)=(-Cpq(32)+Cpq(34))*Q%C(2)+(Spq(32)+Spq(34))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(26)+Q%C(1)*Spq(33)
        W(2)=Q%C(2)*(-Spq(32)+Spq(34))-(Cpq(32)+Cpq(34))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(26)=SPKetC(26)+CTwo*Two
        SPKetS(26)=SPKetS(26)+STwo*Two
!       -------------------------------- Lp=6, Mp = 6-------------------------------------
        W(1)=Cpq(27)*Q%C(0)+Cpq(34)*Q%C(1)
        W(2)=(-Cpq(33)+Cpq(35))*Q%C(2)+(Spq(33)+Spq(35))*Q%S(2)
        CTwo=W(1)+W(2)
        W(1)=Q%C(0)*Spq(27)+Q%C(1)*Spq(34)
        W(2)=Q%C(2)*(-Spq(33)+Spq(35))-(Cpq(33)+Cpq(35))*Q%S(2)
        STwo=W(1)+W(2)
        SPKetC(27)=SPKetC(27)+CTwo*Two
        SPKetS(27)=SPKetS(27)+STwo*Two
      CASE(602)
         CALL CTraX602(Q)
      CASE(603)
         CALL CTraX603(Q)
      CASE(604)
         CALL CTraX604(Q)
      CASE(605)
         CALL CTraX605(Q)
      CASE(606)
         CALL CTraX606(Q)
      CASE(607)
         CALL CTraX607(Q)
      CASE(700)
!       -------------------------------- Lp=0, Mp = 0-------------------------------------
        COne=Cpq(0)*Q%C(0)
        SPKetC(0)=SPKetC(0)+COne
!       -------------------------------- Lp=1, Mp = 0-------------------------------------
        COne=-(Cpq(1)*Q%C(0))
        SPKetC(1)=SPKetC(1)+COne
!       -------------------------------- Lp=1, Mp = 1-------------------------------------
        CTwo=-(Cpq(2)*Q%C(0))
        STwo=-(Q%C(0)*Spq(2))
        SPKetC(2)=SPKetC(2)+CTwo*Two
        SPKetS(2)=SPKetS(2)+STwo*Two
!       -------------------------------- Lp=2, Mp = 0-------------------------------------
        COne=Cpq(3)*Q%C(0)
        SPKetC(3)=SPKetC(3)+COne
!       -------------------------------- Lp=2, Mp = 1-------------------------------------
        CTwo=Cpq(4)*Q%C(0)
        STwo=Q%C(0)*Spq(4)
        SPKetC(4)=SPKetC(4)+CTwo*Two
        SPKetS(4)=SPKetS(4)+STwo*Two
!       -------------------------------- Lp=2, Mp = 2-------------------------------------
        CTwo=Cpq(5)*Q%C(0)
        STwo=Q%C(0)*Spq(5)
        SPKetC(5)=SPKetC(5)+CTwo*Two
        SPKetS(5)=SPKetS(5)+STwo*Two
!       -------------------------------- Lp=3, Mp = 0-------------------------------------
        COne=-(Cpq(6)*Q%C(0))
        SPKetC(6)=SPKetC(6)+COne
!       -------------------------------- Lp=3, Mp = 1-------------------------------------
        CTwo=-(Cpq(7)*Q%C(0))
        STwo=-(Q%C(0)*Spq(7))
        SPKetC(7)=SPKetC(7)+CTwo*Two
        SPKetS(7)=SPKetS(7)+STwo*Two
!       -------------------------------- Lp=3, Mp = 2-------------------------------------
        CTwo=-(Cpq(8)*Q%C(0))
        STwo=-(Q%C(0)*Spq(8))
        SPKetC(8)=SPKetC(8)+CTwo*Two
        SPKetS(8)=SPKetS(8)+STwo*Two
!       -------------------------------- Lp=3, Mp = 3-------------------------------------
        CTwo=-(Cpq(9)*Q%C(0))
        STwo=-(Q%C(0)*Spq(9))
        SPKetC(9)=SPKetC(9)+CTwo*Two
        SPKetS(9)=SPKetS(9)+STwo*Two
!       -------------------------------- Lp=4, Mp = 0-------------------------------------
        COne=Cpq(10)*Q%C(0)
        SPKetC(10)=SPKetC(10)+COne
!       -------------------------------- Lp=4, Mp = 1-------------------------------------
        CTwo=Cpq(11)*Q%C(0)
        STwo=Q%C(0)*Spq(11)
        SPKetC(11)=SPKetC(11)+CTwo*Two
        SPKetS(11)=SPKetS(11)+STwo*Two
!       -------------------------------- Lp=4, Mp = 2-------------------------------------
        CTwo=Cpq(12)*Q%C(0)
        STwo=Q%C(0)*Spq(12)
        SPKetC(12)=SPKetC(12)+CTwo*Two
        SPKetS(12)=SPKetS(12)+STwo*Two
!       -------------------------------- Lp=4, Mp = 3-------------------------------------
        CTwo=Cpq(13)*Q%C(0)
        STwo=Q%C(0)*Spq(13)
        SPKetC(13)=SPKetC(13)+CTwo*Two
        SPKetS(13)=SPKetS(13)+STwo*Two
!       -------------------------------- Lp=4, Mp = 4-------------------------------------
        CTwo=Cpq(14)*Q%C(0)
        STwo=Q%C(0)*Spq(14)
        SPKetC(14)=SPKetC(14)+CTwo*Two
        SPKetS(14)=SPKetS(14)+STwo*Two
!       -------------------------------- Lp=5, Mp = 0-------------------------------------
        COne=-(Cpq(15)*Q%C(0))
        SPKetC(15)=SPKetC(15)+COne
!       -------------------------------- Lp=5, Mp = 1-------------------------------------
        CTwo=-(Cpq(16)*Q%C(0))
        STwo=-(Q%C(0)*Spq(16))
        SPKetC(16)=SPKetC(16)+CTwo*Two
        SPKetS(16)=SPKetS(16)+STwo*Two
!       -------------------------------- Lp=5, Mp = 2-------------------------------------
        CTwo=-(Cpq(17)*Q%C(0))
        STwo=-(Q%C(0)*Spq(17))
        SPKetC(17)=SPKetC(17)+CTwo*Two
        SPKetS(17)=SPKetS(17)+STwo*Two
!       -------------------------------- Lp=5, Mp = 3-------------------------------------
        CTwo=-(Cpq(18)*Q%C(0))
        STwo=-(Q%C(0)*Spq(18))
        SPKetC(18)=SPKetC(18)+CTwo*Two
        SPKetS(18)=SPKetS(18)+STwo*Two
!       -------------------------------- Lp=5, Mp = 4-------------------------------------
        CTwo=-(Cpq(19)*Q%C(0))
        STwo=-(Q%C(0)*Spq(19))
        SPKetC(19)=SPKetC(19)+CTwo*Two
        SPKetS(19)=SPKetS(19)+STwo*Two
!       -------------------------------- Lp=5, Mp = 5-------------------------------------
        CTwo=-(Cpq(20)*Q%C(0))
        STwo=-(Q%C(0)*Spq(20))
        SPKetC(20)=SPKetC(20)+CTwo*Two
        SPKetS(20)=SPKetS(20)+STwo*Two
!       -------------------------------- Lp=6, Mp = 0-------------------------------------
        COne=Cpq(21)*Q%C(0)
        SPKetC(21)=SPKetC(21)+COne
!       -------------------------------- Lp=6, Mp = 1-------------------------------------
        CTwo=Cpq(22)*Q%C(0)
        STwo=Q%C(0)*Spq(22)
        SPKetC(22)=SPKetC(22)+CTwo*Two
        SPKetS(22)=SPKetS(22)+STwo*Two
!       -------------------------------- Lp=6, Mp = 2-------------------------------------
        CTwo=Cpq(23)*Q%C(0)
        STwo=Q%C(0)*Spq(23)
        SPKetC(23)=SPKetC(23)+CTwo*Two
        SPKetS(23)=SPKetS(23)+STwo*Two
!       -------------------------------- Lp=6, Mp = 3-------------------------------------
        CTwo=Cpq(24)*Q%C(0)
        STwo=Q%C(0)*Spq(24)
        SPKetC(24)=SPKetC(24)+CTwo*Two
        SPKetS(24)=SPKetS(24)+STwo*Two
!       -------------------------------- Lp=6, Mp = 4-------------------------------------
        CTwo=Cpq(25)*Q%C(0)
        STwo=Q%C(0)*Spq(25)
        SPKetC(25)=SPKetC(25)+CTwo*Two
        SPKetS(25)=SPKetS(25)+STwo*Two
!       -------------------------------- Lp=6, Mp = 5-------------------------------------
        CTwo=Cpq(26)*Q%C(0)
        STwo=Q%C(0)*Spq(26)
        SPKetC(26)=SPKetC(26)+CTwo*Two
        SPKetS(26)=SPKetS(26)+STwo*Two
!       -------------------------------- Lp=6, Mp = 6-------------------------------------
        CTwo=Cpq(27)*Q%C(0)
        STwo=Q%C(0)*Spq(27)
        SPKetC(27)=SPKetC(27)+CTwo*Two
        SPKetS(27)=SPKetS(27)+STwo*Two
!       -------------------------------- Lp=7, Mp = 0-------------------------------------
        COne=-(Cpq(28)*Q%C(0))
        SPKetC(28)=SPKetC(28)+COne
!       -------------------------------- Lp=7, Mp = 1-------------------------------------
        CTwo=-(Cpq(29)*Q%C(0))
        STwo=-(Q%C(0)*Spq(29))
        SPKetC(29)=SPKetC(29)+CTwo*Two
        SPKetS(29)=SPKetS(29)+STwo*Two
!       -------------------------------- Lp=7, Mp = 2-------------------------------------
        CTwo=-(Cpq(30)*Q%C(0))
        STwo=-(Q%C(0)*Spq(30))
        SPKetC(30)=SPKetC(30)+CTwo*Two
        SPKetS(30)=SPKetS(30)+STwo*Two
!       -------------------------------- Lp=7, Mp = 3-------------------------------------
        CTwo=-(Cpq(31)*Q%C(0))
        STwo=-(Q%C(0)*Spq(31))
        SPKetC(31)=SPKetC(31)+CTwo*Two
        SPKetS(31)=SPKetS(31)+STwo*Two
!       -------------------------------- Lp=7, Mp = 4-------------------------------------
        CTwo=-(Cpq(32)*Q%C(0))
        STwo=-(Q%C(0)*Spq(32))
        SPKetC(32)=SPKetC(32)+CTwo*Two
        SPKetS(32)=SPKetS(32)+STwo*Two
!       -------------------------------- Lp=7, Mp = 5-------------------------------------
        CTwo=-(Cpq(33)*Q%C(0))
        STwo=-(Q%C(0)*Spq(33))
        SPKetC(33)=SPKetC(33)+CTwo*Two
        SPKetS(33)=SPKetS(33)+STwo*Two
!       -------------------------------- Lp=7, Mp = 6-------------------------------------
        CTwo=-(Cpq(34)*Q%C(0))
        STwo=-(Q%C(0)*Spq(34))
        SPKetC(34)=SPKetC(34)+CTwo*Two
        SPKetS(34)=SPKetS(34)+STwo*Two
!       -------------------------------- Lp=7, Mp = 7-------------------------------------
        CTwo=-(Cpq(35)*Q%C(0))
        STwo=-(Q%C(0)*Spq(35))
        SPKetC(35)=SPKetC(35)+CTwo*Two
        SPKetS(35)=SPKetS(35)+STwo*Two
      CASE(701)
         CALL CTraX701(Q)
      CASE(702)
         CALL CTraX702(Q)
      CASE(703)
         CALL CTraX703(Q)
      CASE(704)
         CALL CTraX704(Q)
      CASE(705)
         CALL CTraX705(Q)
      CASE(706)
         CALL CTraX706(Q)
      CASE(707)
         CALL CTraX707(Q)
      END SELECT
