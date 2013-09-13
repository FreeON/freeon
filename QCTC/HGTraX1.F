      SUBROUTINE MD3TRR__1_0(PQx,PQy,PQz,AuxR,R)
      IMPLICIT NONE
      REAL*8 AuxR(0:1)
      REAL*8 R(4)
      REAL*8 PQx,PQy,PQz
      R(1)=AuxR(1)
      R(2)=PQx*R(1)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(0)
      RETURN
      END
      SUBROUTINE KTrax__1_0(R,Co,Ket)
      IMPLICIT NONE
      REAL*8 R(4)
      REAL*8 Co(1)
      REAL*8 Ket(4)
      Ket(1)=Ket(1)+R(1)*Co(1)
      Ket(2)=Ket(2)+R(2)*Co(1)
      Ket(3)=Ket(3)+R(3)*Co(1)
      Ket(4)=Ket(4)+R(4)*Co(1)
      RETURN
      END
      SUBROUTINE MD3TRR__1_1(PQx,PQy,PQz,AuxR,R)
      IMPLICIT NONE
      REAL*8 AuxR(0:2)
      REAL*8 R(10)
      REAL*8 PQx,PQy,PQz
      R(1)=AuxR(2)
      R(2)=PQx*R(1)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(1)
      R(5)=PQx*R(2)+R(1)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(6)=PQx*R(3)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(0)
      RETURN
      END
      SUBROUTINE KTrax__1_1(R,Co,Ket)
      IMPLICIT NONE
      REAL*8 R(10)
      REAL*8 Co(4)
      REAL*8 Ket(4)
      Ket(1)=Ket(1)+R(1)*Co(1)
      Ket(2)=Ket(2)+R(2)*Co(1)
      Ket(3)=Ket(3)+R(3)*Co(1)
      Ket(1)=Ket(1)+R(2)*Co(2)
      Ket(1)=Ket(1)+R(3)*Co(3)
      Ket(4)=Ket(4)+R(4)*Co(1)
      Ket(2)=Ket(2)+R(5)*Co(2)
      Ket(3)=Ket(3)+R(6)*Co(2)
      Ket(2)=Ket(2)+R(6)*Co(3)
      Ket(3)=Ket(3)+R(7)*Co(3)
      Ket(1)=Ket(1)+R(4)*Co(4)
      Ket(4)=Ket(4)+R(8)*Co(2)
      Ket(4)=Ket(4)+R(9)*Co(3)
      Ket(2)=Ket(2)+R(8)*Co(4)
      Ket(3)=Ket(3)+R(9)*Co(4)
      Ket(4)=Ket(4)+R(10)*Co(4)
      RETURN
      END
      SUBROUTINE MD3TRR__1_2(PQx,PQy,PQz,AuxR,R)
      IMPLICIT NONE
      REAL*8 AuxR(0:3)
      REAL*8 R(20)
      REAL*8 PQx,PQy,PQz
      R(1)=AuxR(3)
      R(2)=PQx*R(1)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(2)
      R(5)=PQx*R(2)+R(1)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(6)=PQx*R(3)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(1)
      R(11)=PQx*R(5)+0.20D+01*R(2)
      R(14)=PQy*R(7)+0.20D+01*R(3)
      R(20)=PQz*R(10)+0.20D+01*R(4)
      R(5)=PQx*R(2)+R(1)
      R(13)=PQx*R(7)
      R(19)=PQy*R(10)
      R(18)=PQx*R(10)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(12)=PQx*R(6)+R(3)
      R(17)=PQy*R(9)+R(4)
      R(15)=PQx*R(8)+R(4)
      R(6)=PQx*R(3)
      R(16)=PQx*R(9)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(0)
      RETURN
      END
      SUBROUTINE KTrax__1_2(R,Co,Ket)
      IMPLICIT NONE
      REAL*8 R(20)
      REAL*8 Co(10)
      REAL*8 Ket(4)
      Ket(1)=Ket(1)+R(1)*Co(1)
      Ket(2)=Ket(2)+R(2)*Co(1)
      Ket(3)=Ket(3)+R(3)*Co(1)
      Ket(1)=Ket(1)+R(2)*Co(2)
      Ket(1)=Ket(1)+R(3)*Co(3)
      Ket(4)=Ket(4)+R(4)*Co(1)
      Ket(2)=Ket(2)+R(5)*Co(2)
      Ket(3)=Ket(3)+R(6)*Co(2)
      Ket(2)=Ket(2)+R(6)*Co(3)
      Ket(3)=Ket(3)+R(7)*Co(3)
      Ket(1)=Ket(1)+R(4)*Co(4)
      Ket(1)=Ket(1)+R(5)*Co(5)
      Ket(1)=Ket(1)+R(6)*Co(6)
      Ket(1)=Ket(1)+R(7)*Co(7)
      Ket(4)=Ket(4)+R(8)*Co(2)
      Ket(4)=Ket(4)+R(9)*Co(3)
      Ket(2)=Ket(2)+R(8)*Co(4)
      Ket(3)=Ket(3)+R(9)*Co(4)
      Ket(4)=Ket(4)+R(10)*Co(4)
      Ket(2)=Ket(2)+R(11)*Co(5)
      Ket(3)=Ket(3)+R(12)*Co(5)
      Ket(4)=Ket(4)+R(15)*Co(5)
      Ket(2)=Ket(2)+R(12)*Co(6)
      Ket(3)=Ket(3)+R(13)*Co(6)
      Ket(2)=Ket(2)+R(13)*Co(7)
      Ket(3)=Ket(3)+R(14)*Co(7)
      Ket(1)=Ket(1)+R(8)*Co(8)
      Ket(1)=Ket(1)+R(9)*Co(9)
      Ket(1)=Ket(1)+R(10)*Co(10)
      Ket(2)=Ket(2)+R(15)*Co(8)
      Ket(4)=Ket(4)+R(16)*Co(6)
      Ket(4)=Ket(4)+R(17)*Co(7)
      Ket(3)=Ket(3)+R(16)*Co(8)
      Ket(2)=Ket(2)+R(16)*Co(9)
      Ket(3)=Ket(3)+R(17)*Co(9)
      Ket(4)=Ket(4)+R(18)*Co(8)
      Ket(4)=Ket(4)+R(19)*Co(9)
      Ket(2)=Ket(2)+R(18)*Co(10)
      Ket(3)=Ket(3)+R(19)*Co(10)
      Ket(4)=Ket(4)+R(20)*Co(10)
      RETURN
      END
      SUBROUTINE MD3TRR__1_3(PQx,PQy,PQz,AuxR,R)
      IMPLICIT NONE
      REAL*8 AuxR(0:4)
      REAL*8 R(35)
      REAL*8 PQx,PQy,PQz
      R(1)=AuxR(4)
      R(2)=PQx*R(1)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(3)
      R(5)=PQx*R(2)+R(1)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(6)=PQx*R(3)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(2)
      R(11)=PQx*R(5)+0.20D+01*R(2)
      R(14)=PQy*R(7)+0.20D+01*R(3)
      R(20)=PQz*R(10)+0.20D+01*R(4)
      R(5)=PQx*R(2)+R(1)
      R(13)=PQx*R(7)
      R(19)=PQy*R(10)
      R(18)=PQx*R(10)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(12)=PQx*R(6)+R(3)
      R(17)=PQy*R(9)+R(4)
      R(15)=PQx*R(8)+R(4)
      R(6)=PQx*R(3)
      R(16)=PQx*R(9)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(1)
      R(21)=PQx*R(11)+0.30D+01*R(5)
      R(25)=PQy*R(14)+0.30D+01*R(7)
      R(35)=PQz*R(20)+0.30D+01*R(10)
      R(11)=PQx*R(5)+0.20D+01*R(2)
      R(24)=PQx*R(14)
      R(34)=PQy*R(20)
      R(33)=PQx*R(20)
      R(14)=PQy*R(7)+0.20D+01*R(3)
      R(20)=PQz*R(10)+0.20D+01*R(4)
      R(5)=PQx*R(2)+R(1)
      R(23)=PQx*R(13)+R(7)
      R(32)=PQy*R(19)+R(10)
      R(30)=PQx*R(18)+R(10)
      R(13)=PQx*R(7)
      R(31)=PQx*R(19)
      R(19)=PQy*R(10)
      R(18)=PQx*R(10)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(22)=PQx*R(12)+0.20D+01*R(6)
      R(29)=PQy*R(17)+0.20D+01*R(9)
      R(26)=PQx*R(15)+0.20D+01*R(8)
      R(12)=PQx*R(6)+R(3)
      R(28)=PQx*R(17)
      R(17)=PQy*R(9)+R(4)
      R(15)=PQx*R(8)+R(4)
      R(6)=PQx*R(3)
      R(27)=PQx*R(16)+R(9)
      R(16)=PQx*R(9)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(0)
      RETURN
      END
      SUBROUTINE KTrax__1_3(R,Co,Ket)
      IMPLICIT NONE
      REAL*8 R(35)
      REAL*8 Co(20)
      REAL*8 Ket(4)
      Ket(1)=Ket(1)+R(1)*Co(1)
      Ket(2)=Ket(2)+R(2)*Co(1)
      Ket(3)=Ket(3)+R(3)*Co(1)
      Ket(1)=Ket(1)+R(2)*Co(2)
      Ket(1)=Ket(1)+R(3)*Co(3)
      Ket(4)=Ket(4)+R(4)*Co(1)
      Ket(2)=Ket(2)+R(5)*Co(2)
      Ket(3)=Ket(3)+R(6)*Co(2)
      Ket(2)=Ket(2)+R(6)*Co(3)
      Ket(3)=Ket(3)+R(7)*Co(3)
      Ket(1)=Ket(1)+R(4)*Co(4)
      Ket(1)=Ket(1)+R(5)*Co(5)
      Ket(1)=Ket(1)+R(6)*Co(6)
      Ket(1)=Ket(1)+R(7)*Co(7)
      Ket(4)=Ket(4)+R(8)*Co(2)
      Ket(4)=Ket(4)+R(9)*Co(3)
      Ket(2)=Ket(2)+R(8)*Co(4)
      Ket(3)=Ket(3)+R(9)*Co(4)
      Ket(4)=Ket(4)+R(10)*Co(4)
      Ket(2)=Ket(2)+R(11)*Co(5)
      Ket(3)=Ket(3)+R(12)*Co(5)
      Ket(4)=Ket(4)+R(15)*Co(5)
      Ket(2)=Ket(2)+R(12)*Co(6)
      Ket(3)=Ket(3)+R(13)*Co(6)
      Ket(2)=Ket(2)+R(13)*Co(7)
      Ket(3)=Ket(3)+R(14)*Co(7)
      Ket(1)=Ket(1)+R(8)*Co(8)
      Ket(1)=Ket(1)+R(9)*Co(9)
      Ket(1)=Ket(1)+R(10)*Co(10)
      Ket(1)=Ket(1)+R(11)*Co(11)
      Ket(2)=Ket(2)+R(15)*Co(8)
      Ket(1)=Ket(1)+R(12)*Co(12)
      Ket(1)=Ket(1)+R(13)*Co(13)
      Ket(1)=Ket(1)+R(14)*Co(14)
      Ket(1)=Ket(1)+R(15)*Co(15)
      Ket(4)=Ket(4)+R(16)*Co(6)
      Ket(4)=Ket(4)+R(17)*Co(7)
      Ket(3)=Ket(3)+R(16)*Co(8)
      Ket(2)=Ket(2)+R(16)*Co(9)
      Ket(3)=Ket(3)+R(17)*Co(9)
      Ket(4)=Ket(4)+R(18)*Co(8)
      Ket(4)=Ket(4)+R(19)*Co(9)
      Ket(2)=Ket(2)+R(18)*Co(10)
      Ket(3)=Ket(3)+R(19)*Co(10)
      Ket(4)=Ket(4)+R(20)*Co(10)
      Ket(2)=Ket(2)+R(21)*Co(11)
      Ket(3)=Ket(3)+R(22)*Co(11)
      Ket(2)=Ket(2)+R(22)*Co(12)
      Ket(3)=Ket(3)+R(23)*Co(12)
      Ket(2)=Ket(2)+R(23)*Co(13)
      Ket(4)=Ket(4)+R(26)*Co(11)
      Ket(3)=Ket(3)+R(24)*Co(13)
      Ket(4)=Ket(4)+R(27)*Co(12)
      Ket(2)=Ket(2)+R(24)*Co(14)
      Ket(3)=Ket(3)+R(25)*Co(14)
      Ket(2)=Ket(2)+R(26)*Co(15)
      Ket(3)=Ket(3)+R(27)*Co(15)
      Ket(4)=Ket(4)+R(28)*Co(13)
      Ket(4)=Ket(4)+R(29)*Co(14)
      Ket(4)=Ket(4)+R(30)*Co(15)
      Ket(1)=Ket(1)+R(16)*Co(16)
      Ket(1)=Ket(1)+R(17)*Co(17)
      Ket(1)=Ket(1)+R(18)*Co(18)
      Ket(1)=Ket(1)+R(19)*Co(19)
      Ket(1)=Ket(1)+R(20)*Co(20)
      Ket(2)=Ket(2)+R(27)*Co(16)
      Ket(3)=Ket(3)+R(28)*Co(16)
      Ket(2)=Ket(2)+R(28)*Co(17)
      Ket(3)=Ket(3)+R(29)*Co(17)
      Ket(4)=Ket(4)+R(31)*Co(16)
      Ket(2)=Ket(2)+R(30)*Co(18)
      Ket(3)=Ket(3)+R(31)*Co(18)
      Ket(2)=Ket(2)+R(31)*Co(19)
      Ket(4)=Ket(4)+R(32)*Co(17)
      Ket(4)=Ket(4)+R(33)*Co(18)
      Ket(3)=Ket(3)+R(32)*Co(19)
      Ket(4)=Ket(4)+R(34)*Co(19)
      Ket(2)=Ket(2)+R(33)*Co(20)
      Ket(3)=Ket(3)+R(34)*Co(20)
      Ket(4)=Ket(4)+R(35)*Co(20)
      RETURN
      END
      SUBROUTINE MD3TRR__1_4(PQx,PQy,PQz,AuxR,R)
      IMPLICIT NONE
      REAL*8 AuxR(0:5)
      REAL*8 R(56)
      REAL*8 PQx,PQy,PQz
      R(1)=AuxR(5)
      R(2)=PQx*R(1)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(4)
      R(5)=PQx*R(2)+R(1)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(6)=PQx*R(3)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(3)
      R(11)=PQx*R(5)+0.20D+01*R(2)
      R(14)=PQy*R(7)+0.20D+01*R(3)
      R(20)=PQz*R(10)+0.20D+01*R(4)
      R(5)=PQx*R(2)+R(1)
      R(13)=PQx*R(7)
      R(19)=PQy*R(10)
      R(18)=PQx*R(10)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(12)=PQx*R(6)+R(3)
      R(17)=PQy*R(9)+R(4)
      R(15)=PQx*R(8)+R(4)
      R(6)=PQx*R(3)
      R(16)=PQx*R(9)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(2)
      R(21)=PQx*R(11)+0.30D+01*R(5)
      R(25)=PQy*R(14)+0.30D+01*R(7)
      R(35)=PQz*R(20)+0.30D+01*R(10)
      R(11)=PQx*R(5)+0.20D+01*R(2)
      R(24)=PQx*R(14)
      R(34)=PQy*R(20)
      R(33)=PQx*R(20)
      R(14)=PQy*R(7)+0.20D+01*R(3)
      R(20)=PQz*R(10)+0.20D+01*R(4)
      R(5)=PQx*R(2)+R(1)
      R(23)=PQx*R(13)+R(7)
      R(32)=PQy*R(19)+R(10)
      R(30)=PQx*R(18)+R(10)
      R(13)=PQx*R(7)
      R(31)=PQx*R(19)
      R(19)=PQy*R(10)
      R(18)=PQx*R(10)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(22)=PQx*R(12)+0.20D+01*R(6)
      R(29)=PQy*R(17)+0.20D+01*R(9)
      R(26)=PQx*R(15)+0.20D+01*R(8)
      R(12)=PQx*R(6)+R(3)
      R(28)=PQx*R(17)
      R(17)=PQy*R(9)+R(4)
      R(15)=PQx*R(8)+R(4)
      R(6)=PQx*R(3)
      R(27)=PQx*R(16)+R(9)
      R(16)=PQx*R(9)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(1)
      R(36)=PQx*R(21)+0.40D+01*R(11)
      R(41)=PQy*R(25)+0.40D+01*R(14)
      R(56)=PQz*R(35)+0.40D+01*R(20)
      R(21)=PQx*R(11)+0.30D+01*R(5)
      R(40)=PQx*R(25)
      R(55)=PQy*R(35)
      R(54)=PQx*R(35)
      R(25)=PQy*R(14)+0.30D+01*R(7)
      R(35)=PQz*R(20)+0.30D+01*R(10)
      R(11)=PQx*R(5)+0.20D+01*R(2)
      R(39)=PQx*R(24)+R(14)
      R(53)=PQy*R(34)+R(20)
      R(51)=PQx*R(33)+R(20)
      R(24)=PQx*R(14)
      R(52)=PQx*R(34)
      R(34)=PQy*R(20)
      R(33)=PQx*R(20)
      R(14)=PQy*R(7)+0.20D+01*R(3)
      R(20)=PQz*R(10)+0.20D+01*R(4)
      R(5)=PQx*R(2)+R(1)
      R(38)=PQx*R(23)+0.20D+01*R(13)
      R(50)=PQy*R(32)+0.20D+01*R(19)
      R(47)=PQx*R(30)+0.20D+01*R(18)
      R(23)=PQx*R(13)+R(7)
      R(49)=PQx*R(32)
      R(32)=PQy*R(19)+R(10)
      R(30)=PQx*R(18)+R(10)
      R(13)=PQx*R(7)
      R(48)=PQx*R(31)+R(19)
      R(31)=PQx*R(19)
      R(19)=PQy*R(10)
      R(18)=PQx*R(10)
      R(7)=PQy*R(3)+R(1)
      R(10)=PQz*R(4)+R(1)
      R(2)=PQx*R(1)
      R(37)=PQx*R(22)+0.30D+01*R(12)
      R(46)=PQy*R(29)+0.30D+01*R(17)
      R(42)=PQx*R(26)+0.30D+01*R(15)
      R(22)=PQx*R(12)+0.20D+01*R(6)
      R(45)=PQx*R(29)
      R(29)=PQy*R(17)+0.20D+01*R(9)
      R(26)=PQx*R(15)+0.20D+01*R(8)
      R(12)=PQx*R(6)+R(3)
      R(44)=PQx*R(28)+R(17)
      R(28)=PQx*R(17)
      R(17)=PQy*R(9)+R(4)
      R(15)=PQx*R(8)+R(4)
      R(6)=PQx*R(3)
      R(43)=PQx*R(27)+0.20D+01*R(16)
      R(27)=PQx*R(16)+R(9)
      R(16)=PQx*R(9)
      R(9)=PQy*R(4)
      R(8)=PQx*R(4)
      R(3)=PQy*R(1)
      R(4)=PQz*R(1)
      R(1)=AuxR(0)
      RETURN
      END
      SUBROUTINE KTrax__1_4(R,Co,Ket)
      IMPLICIT NONE
      REAL*8 R(56)
      REAL*8 Co(35)
      REAL*8 Ket(4)
      Ket(1)=Ket(1)+R(1)*Co(1)
      Ket(2)=Ket(2)+R(2)*Co(1)
      Ket(3)=Ket(3)+R(3)*Co(1)
      Ket(1)=Ket(1)+R(2)*Co(2)
      Ket(1)=Ket(1)+R(3)*Co(3)
      Ket(4)=Ket(4)+R(4)*Co(1)
      Ket(2)=Ket(2)+R(5)*Co(2)
      Ket(3)=Ket(3)+R(6)*Co(2)
      Ket(2)=Ket(2)+R(6)*Co(3)
      Ket(3)=Ket(3)+R(7)*Co(3)
      Ket(1)=Ket(1)+R(4)*Co(4)
      Ket(1)=Ket(1)+R(5)*Co(5)
      Ket(1)=Ket(1)+R(6)*Co(6)
      Ket(1)=Ket(1)+R(7)*Co(7)
      Ket(4)=Ket(4)+R(8)*Co(2)
      Ket(4)=Ket(4)+R(9)*Co(3)
      Ket(2)=Ket(2)+R(8)*Co(4)
      Ket(3)=Ket(3)+R(9)*Co(4)
      Ket(4)=Ket(4)+R(10)*Co(4)
      Ket(2)=Ket(2)+R(11)*Co(5)
      Ket(3)=Ket(3)+R(12)*Co(5)
      Ket(4)=Ket(4)+R(15)*Co(5)
      Ket(2)=Ket(2)+R(12)*Co(6)
      Ket(3)=Ket(3)+R(13)*Co(6)
      Ket(2)=Ket(2)+R(13)*Co(7)
      Ket(3)=Ket(3)+R(14)*Co(7)
      Ket(1)=Ket(1)+R(8)*Co(8)
      Ket(1)=Ket(1)+R(9)*Co(9)
      Ket(1)=Ket(1)+R(10)*Co(10)
      Ket(1)=Ket(1)+R(11)*Co(11)
      Ket(2)=Ket(2)+R(15)*Co(8)
      Ket(1)=Ket(1)+R(12)*Co(12)
      Ket(1)=Ket(1)+R(13)*Co(13)
      Ket(1)=Ket(1)+R(14)*Co(14)
      Ket(1)=Ket(1)+R(15)*Co(15)
      Ket(4)=Ket(4)+R(16)*Co(6)
      Ket(4)=Ket(4)+R(17)*Co(7)
      Ket(3)=Ket(3)+R(16)*Co(8)
      Ket(2)=Ket(2)+R(16)*Co(9)
      Ket(3)=Ket(3)+R(17)*Co(9)
      Ket(4)=Ket(4)+R(18)*Co(8)
      Ket(4)=Ket(4)+R(19)*Co(9)
      Ket(2)=Ket(2)+R(18)*Co(10)
      Ket(3)=Ket(3)+R(19)*Co(10)
      Ket(4)=Ket(4)+R(20)*Co(10)
      Ket(2)=Ket(2)+R(21)*Co(11)
      Ket(3)=Ket(3)+R(22)*Co(11)
      Ket(2)=Ket(2)+R(22)*Co(12)
      Ket(3)=Ket(3)+R(23)*Co(12)
      Ket(2)=Ket(2)+R(23)*Co(13)
      Ket(4)=Ket(4)+R(26)*Co(11)
      Ket(3)=Ket(3)+R(24)*Co(13)
      Ket(4)=Ket(4)+R(27)*Co(12)
      Ket(2)=Ket(2)+R(24)*Co(14)
      Ket(3)=Ket(3)+R(25)*Co(14)
      Ket(2)=Ket(2)+R(26)*Co(15)
      Ket(3)=Ket(3)+R(27)*Co(15)
      Ket(4)=Ket(4)+R(28)*Co(13)
      Ket(4)=Ket(4)+R(29)*Co(14)
      Ket(4)=Ket(4)+R(30)*Co(15)
      Ket(1)=Ket(1)+R(16)*Co(16)
      Ket(1)=Ket(1)+R(17)*Co(17)
      Ket(1)=Ket(1)+R(18)*Co(18)
      Ket(1)=Ket(1)+R(19)*Co(19)
      Ket(1)=Ket(1)+R(20)*Co(20)
      Ket(1)=Ket(1)+R(21)*Co(21)
      Ket(1)=Ket(1)+R(22)*Co(22)
      Ket(1)=Ket(1)+R(23)*Co(23)
      Ket(2)=Ket(2)+R(27)*Co(16)
      Ket(3)=Ket(3)+R(28)*Co(16)
      Ket(2)=Ket(2)+R(28)*Co(17)
      Ket(3)=Ket(3)+R(29)*Co(17)
      Ket(4)=Ket(4)+R(31)*Co(16)
      Ket(2)=Ket(2)+R(30)*Co(18)
      Ket(3)=Ket(3)+R(31)*Co(18)
      Ket(2)=Ket(2)+R(31)*Co(19)
      Ket(1)=Ket(1)+R(24)*Co(24)
      Ket(1)=Ket(1)+R(25)*Co(25)
      Ket(1)=Ket(1)+R(26)*Co(26)
      Ket(1)=Ket(1)+R(27)*Co(27)
      Ket(1)=Ket(1)+R(28)*Co(28)
      Ket(1)=Ket(1)+R(29)*Co(29)
      Ket(1)=Ket(1)+R(30)*Co(30)
      Ket(1)=Ket(1)+R(31)*Co(31)
      Ket(4)=Ket(4)+R(32)*Co(17)
      Ket(4)=Ket(4)+R(33)*Co(18)
      Ket(3)=Ket(3)+R(32)*Co(19)
      Ket(4)=Ket(4)+R(34)*Co(19)
      Ket(2)=Ket(2)+R(33)*Co(20)
      Ket(3)=Ket(3)+R(34)*Co(20)
      Ket(4)=Ket(4)+R(35)*Co(20)
      Ket(2)=Ket(2)+R(36)*Co(21)
      Ket(3)=Ket(3)+R(37)*Co(21)
      Ket(2)=Ket(2)+R(37)*Co(22)
      Ket(3)=Ket(3)+R(38)*Co(22)
      Ket(2)=Ket(2)+R(38)*Co(23)
      Ket(3)=Ket(3)+R(39)*Co(23)
      Ket(4)=Ket(4)+R(42)*Co(21)
      Ket(4)=Ket(4)+R(43)*Co(22)
      Ket(4)=Ket(4)+R(44)*Co(23)
      Ket(2)=Ket(2)+R(39)*Co(24)
      Ket(3)=Ket(3)+R(40)*Co(24)
      Ket(2)=Ket(2)+R(40)*Co(25)
      Ket(3)=Ket(3)+R(41)*Co(25)
      Ket(2)=Ket(2)+R(42)*Co(26)
      Ket(3)=Ket(3)+R(43)*Co(26)
      Ket(2)=Ket(2)+R(43)*Co(27)
      Ket(4)=Ket(4)+R(45)*Co(24)
      Ket(4)=Ket(4)+R(46)*Co(25)
      Ket(3)=Ket(3)+R(44)*Co(27)
      Ket(4)=Ket(4)+R(47)*Co(26)
      Ket(2)=Ket(2)+R(44)*Co(28)
      Ket(3)=Ket(3)+R(45)*Co(28)
      Ket(2)=Ket(2)+R(45)*Co(29)
      Ket(3)=Ket(3)+R(46)*Co(29)
      Ket(2)=Ket(2)+R(47)*Co(30)
      Ket(4)=Ket(4)+R(48)*Co(27)
      Ket(4)=Ket(4)+R(49)*Co(28)
      Ket(4)=Ket(4)+R(50)*Co(29)
      Ket(3)=Ket(3)+R(48)*Co(30)
      Ket(2)=Ket(2)+R(48)*Co(31)
      Ket(3)=Ket(3)+R(49)*Co(31)
      Ket(4)=Ket(4)+R(51)*Co(30)
      Ket(4)=Ket(4)+R(52)*Co(31)
      Ket(1)=Ket(1)+R(32)*Co(32)
      Ket(1)=Ket(1)+R(33)*Co(33)
      Ket(1)=Ket(1)+R(34)*Co(34)
      Ket(1)=Ket(1)+R(35)*Co(35)
      Ket(2)=Ket(2)+R(49)*Co(32)
      Ket(3)=Ket(3)+R(50)*Co(32)
      Ket(2)=Ket(2)+R(51)*Co(33)
      Ket(4)=Ket(4)+R(53)*Co(32)
      Ket(3)=Ket(3)+R(52)*Co(33)
      Ket(4)=Ket(4)+R(54)*Co(33)
      Ket(2)=Ket(2)+R(52)*Co(34)
      Ket(3)=Ket(3)+R(53)*Co(34)
      Ket(4)=Ket(4)+R(55)*Co(34)
      Ket(2)=Ket(2)+R(54)*Co(35)
      Ket(3)=Ket(3)+R(55)*Co(35)
      Ket(4)=Ket(4)+R(56)*Co(35)
      RETURN
      END
      SUBROUTINE HGTraX10(PQx,PQy,PQz,U,O,Co,Ket)
      IMPLICIT NONE
      REAL*8 PQx,PQy,PQz,U,O,o1,o2,OneOvT
      REAL*8 G(0:1)
      REAL*8 AuxR(0:1)
      REAL*8 R(4)
      REAL*8 Co(1)
      REAL*8 Ket(4)
      INTEGER Mesh,I,J,K
      REAL*8 ET,TwoT,T,T2,T3,T4
      REAL*8 SqrtT,SqrtPi,One,Two
      REAL*8 Switch,Grid,F0Asymp,F1Asymp,F2Asymp,F3Asymp
      REAL*8 F4Asymp,F5Asymp,F6Asymp,F7Asymp,F8Asymp
      REAL*8 F9Asymp,F10Asymp,F11Asymp,F12Asymp
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "Gamma_Asymptotics.Inc"
      INCLUDE "Mesh.Inc"
      REAL*8 F1_0(0:Mesh)
      REAL*8 F1_1(0:Mesh)
      REAL*8 F1_2(0:Mesh)
      REAL*8 F1_3(0:Mesh)
      INCLUDE "Gamma_1.Inc"
      T=O*(PQx*PQx+PQy*PQy+PQz*PQz)
      IF(T.LT.Switch)THEN
         T2=T*T
         T3=T*T2
         J=AINT(T*Grid)
         G(1)=(F1_0(J)+T*F1_1(J)+T2*F1_2(J)+T3*F1_3(J))
         ET=DEXP(-T)
         TwoT=2.0D0*T
         G(0)=TwoT*G(1)+ET
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         AuxR(1)=o1*G(1)
      ELSE
         SqrtT=DSQRT(T)
         OneOvT=One/T
         G(0)=SqrtPi/(Two*SqrtT)
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         G(1)=G(0)*0.5000000000000000D+00*OneOvT
         AuxR(1)=o1*G(1)
      ENDIF
      CALL MD3TRR__1_0(PQx,PQy,PQz,AuxR,R)
      CALL KTraX__1_0(R,Co,Ket)
      RETURN
      END

      SUBROUTINE HGTraX11(PQx,PQy,PQz,U,O,Co,Ket)
      IMPLICIT NONE
      REAL*8 PQx,PQy,PQz,U,O,o1,o2,OneOvT
      REAL*8 G(0:2)
      REAL*8 AuxR(0:2)
      REAL*8 R(10)
      REAL*8 Co(4)
      REAL*8 Ket(4)
      INTEGER Mesh,I,J,K
      REAL*8 ET,TwoT,T,T2,T3,T4
      REAL*8 SqrtT,SqrtPi,One,Two
      REAL*8 Switch,Grid,F0Asymp,F1Asymp,F2Asymp,F3Asymp
      REAL*8 F4Asymp,F5Asymp,F6Asymp,F7Asymp,F8Asymp
      REAL*8 F9Asymp,F10Asymp,F11Asymp,F12Asymp
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "Gamma_Asymptotics.Inc"
      INCLUDE "Mesh.Inc"
      REAL*8 F2_0(0:Mesh)
      REAL*8 F2_1(0:Mesh)
      REAL*8 F2_2(0:Mesh)
      REAL*8 F2_3(0:Mesh)
      INCLUDE "Gamma_2.Inc"
      T=O*(PQx*PQx+PQy*PQy+PQz*PQz)
      IF(T.LT.Switch)THEN
         T2=T*T
         T3=T*T2
         J=AINT(T*Grid)
         G(2)=(F2_0(J)+T*F2_1(J)+T2*F2_2(J)+T3*F2_3(J))
         ET=DEXP(-T)
         TwoT=2.0D0*T
         G(1)=0.3333333333333333D+00*(TwoT*G(2)+ET)
         G(0)=TwoT*G(1)+ET
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         AuxR(1)=o1*G(1)
         o1=o2*o1
         AuxR(2)=o1*G(2)
      ELSE
         SqrtT=DSQRT(T)
         OneOvT=One/T
         G(0)=SqrtPi/(Two*SqrtT)
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         G(1)=G(0)*0.5000000000000000D+00*OneOvT
         AuxR(1)=o1*G(1)
         o1=o2*o1
         G(2)=G(1)*0.1500000000000000D+01*OneOvT
         AuxR(2)=o1*G(2)
      ENDIF
      CALL MD3TRR__1_1(PQx,PQy,PQz,AuxR,R)
      CALL KTraX__1_1(R,Co,Ket)
      RETURN
      END

      SUBROUTINE HGTraX12(PQx,PQy,PQz,U,O,Co,Ket)
      IMPLICIT NONE
      REAL*8 PQx,PQy,PQz,U,O,o1,o2,OneOvT
      REAL*8 G(0:3)
      REAL*8 AuxR(0:3)
      REAL*8 R(20)
      REAL*8 Co(10)
      REAL*8 Ket(4)
      INTEGER Mesh,I,J,K
      REAL*8 ET,TwoT,T,T2,T3,T4
      REAL*8 SqrtT,SqrtPi,One,Two
      REAL*8 Switch,Grid,F0Asymp,F1Asymp,F2Asymp,F3Asymp
      REAL*8 F4Asymp,F5Asymp,F6Asymp,F7Asymp,F8Asymp
      REAL*8 F9Asymp,F10Asymp,F11Asymp,F12Asymp
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "Gamma_Asymptotics.Inc"
      INCLUDE "Mesh.Inc"
      REAL*8 F3_0(0:Mesh)
      REAL*8 F3_1(0:Mesh)
      REAL*8 F3_2(0:Mesh)
      REAL*8 F3_3(0:Mesh)
      INCLUDE "Gamma_3.Inc"
      T=O*(PQx*PQx+PQy*PQy+PQz*PQz)
      IF(T.LT.Switch)THEN
         T2=T*T
         T3=T*T2
         J=AINT(T*Grid)
         G(3)=(F3_0(J)+T*F3_1(J)+T2*F3_2(J)+T3*F3_3(J))
         ET=DEXP(-T)
         TwoT=2.0D0*T
         G(2)=0.2000000000000000D+00*(TwoT*G(3)+ET)
         G(1)=0.3333333333333333D+00*(TwoT*G(2)+ET)
         G(0)=TwoT*G(1)+ET
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         AuxR(1)=o1*G(1)
         o1=o2*o1
         AuxR(2)=o1*G(2)
         o1=o2*o1
         AuxR(3)=o1*G(3)
      ELSE
         SqrtT=DSQRT(T)
         OneOvT=One/T
         G(0)=SqrtPi/(Two*SqrtT)
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         G(1)=G(0)*0.5000000000000000D+00*OneOvT
         AuxR(1)=o1*G(1)
         o1=o2*o1
         G(2)=G(1)*0.1500000000000000D+01*OneOvT
         AuxR(2)=o1*G(2)
         o1=o2*o1
         G(3)=G(2)*0.2500000000000000D+01*OneOvT
         AuxR(3)=o1*G(3)
      ENDIF
      CALL MD3TRR__1_2(PQx,PQy,PQz,AuxR,R)
      CALL KTraX__1_2(R,Co,Ket)
      RETURN
      END

      SUBROUTINE HGTraX13(PQx,PQy,PQz,U,O,Co,Ket)
      IMPLICIT NONE
      REAL*8 PQx,PQy,PQz,U,O,o1,o2,OneOvT
      REAL*8 G(0:4)
      REAL*8 AuxR(0:4)
      REAL*8 R(35)
      REAL*8 Co(20)
      REAL*8 Ket(4)
      INTEGER Mesh,I,J,K
      REAL*8 ET,TwoT,T,T2,T3,T4
      REAL*8 SqrtT,SqrtPi,One,Two
      REAL*8 Switch,Grid,F0Asymp,F1Asymp,F2Asymp,F3Asymp
      REAL*8 F4Asymp,F5Asymp,F6Asymp,F7Asymp,F8Asymp
      REAL*8 F9Asymp,F10Asymp,F11Asymp,F12Asymp
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "Gamma_Asymptotics.Inc"
      INCLUDE "Mesh.Inc"
      REAL*8 F4_0(0:Mesh)
      REAL*8 F4_1(0:Mesh)
      REAL*8 F4_2(0:Mesh)
      REAL*8 F4_3(0:Mesh)
      INCLUDE "Gamma_4.Inc"
      T=O*(PQx*PQx+PQy*PQy+PQz*PQz)
      IF(T.LT.Switch)THEN
         T2=T*T
         T3=T*T2
         J=AINT(T*Grid)
         G(4)=(F4_0(J)+T*F4_1(J)+T2*F4_2(J)+T3*F4_3(J))
         ET=DEXP(-T)
         TwoT=2.0D0*T
         G(3)=0.1428571428571428D+00*(TwoT*G(4)+ET)
         G(2)=0.2000000000000000D+00*(TwoT*G(3)+ET)
         G(1)=0.3333333333333333D+00*(TwoT*G(2)+ET)
         G(0)=TwoT*G(1)+ET
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         AuxR(1)=o1*G(1)
         o1=o2*o1
         AuxR(2)=o1*G(2)
         o1=o2*o1
         AuxR(3)=o1*G(3)
         o1=o2*o1
         AuxR(4)=o1*G(4)
      ELSE
         SqrtT=DSQRT(T)
         OneOvT=One/T
         G(0)=SqrtPi/(Two*SqrtT)
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         G(1)=G(0)*0.5000000000000000D+00*OneOvT
         AuxR(1)=o1*G(1)
         o1=o2*o1
         G(2)=G(1)*0.1500000000000000D+01*OneOvT
         AuxR(2)=o1*G(2)
         o1=o2*o1
         G(3)=G(2)*0.2500000000000000D+01*OneOvT
         AuxR(3)=o1*G(3)
         o1=o2*o1
         G(4)=G(3)*0.3500000000000000D+01*OneOvT
         AuxR(4)=o1*G(4)
      ENDIF
      CALL MD3TRR__1_3(PQx,PQy,PQz,AuxR,R)
      CALL KTraX__1_3(R,Co,Ket)
      RETURN
      END

      SUBROUTINE HGTraX14(PQx,PQy,PQz,U,O,Co,Ket)
      IMPLICIT NONE
      REAL*8 PQx,PQy,PQz,U,O,o1,o2,OneOvT
      REAL*8 G(0:5)
      REAL*8 AuxR(0:5)
      REAL*8 R(56)
      REAL*8 Co(35)
      REAL*8 Ket(4)
      INTEGER Mesh,I,J,K
      REAL*8 ET,TwoT,T,T2,T3,T4
      REAL*8 SqrtT,SqrtPi,One,Two
      REAL*8 Switch,Grid,F0Asymp,F1Asymp,F2Asymp,F3Asymp
      REAL*8 F4Asymp,F5Asymp,F6Asymp,F7Asymp,F8Asymp
      REAL*8 F9Asymp,F10Asymp,F11Asymp,F12Asymp
      PARAMETER(SqrtPi=1.7724538509055160273D0)
      PARAMETER(One=1D0)
      PARAMETER(Two=2D0)
      INCLUDE "Gamma_Asymptotics.Inc"
      INCLUDE "Mesh.Inc"
      REAL*8 F5_0(0:Mesh)
      REAL*8 F5_1(0:Mesh)
      REAL*8 F5_2(0:Mesh)
      REAL*8 F5_3(0:Mesh)
      INCLUDE "Gamma_5.Inc"
      T=O*(PQx*PQx+PQy*PQy+PQz*PQz)
      IF(T.LT.Switch)THEN
         T2=T*T
         T3=T*T2
         J=AINT(T*Grid)
         G(5)=(F5_0(J)+T*F5_1(J)+T2*F5_2(J)+T3*F5_3(J))
         ET=DEXP(-T)
         TwoT=2.0D0*T
         G(4)=0.1111111111111111D+00*(TwoT*G(5)+ET)
         G(3)=0.1428571428571428D+00*(TwoT*G(4)+ET)
         G(2)=0.2000000000000000D+00*(TwoT*G(3)+ET)
         G(1)=0.3333333333333333D+00*(TwoT*G(2)+ET)
         G(0)=TwoT*G(1)+ET
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         AuxR(1)=o1*G(1)
         o1=o2*o1
         AuxR(2)=o1*G(2)
         o1=o2*o1
         AuxR(3)=o1*G(3)
         o1=o2*o1
         AuxR(4)=o1*G(4)
         o1=o2*o1
         AuxR(5)=o1*G(5)
      ELSE
         SqrtT=DSQRT(T)
         OneOvT=One/T
         G(0)=SqrtPi/(Two*SqrtT)
         o1=U
         o2=-2D0*O
         AuxR(0)=o1*G(0)
         o1=o2*o1
         G(1)=G(0)*0.5000000000000000D+00*OneOvT
         AuxR(1)=o1*G(1)
         o1=o2*o1
         G(2)=G(1)*0.1500000000000000D+01*OneOvT
         AuxR(2)=o1*G(2)
         o1=o2*o1
         G(3)=G(2)*0.2500000000000000D+01*OneOvT
         AuxR(3)=o1*G(3)
         o1=o2*o1
         G(4)=G(3)*0.3500000000000000D+01*OneOvT
         AuxR(4)=o1*G(4)
         o1=o2*o1
         G(5)=G(4)*0.4500000000000000D+01*OneOvT
         AuxR(5)=o1*G(5)
      ENDIF
      CALL MD3TRR__1_4(PQx,PQy,PQz,AuxR,R)
      CALL KTraX__1_4(R,Co,Ket)
      RETURN
      END
