      if(lt.eq.0)then 
        rC(0)=zz*rhgtf(1)
        rS(0)=re0
      elseif(lt.eq.1)then 
        rC(0)=zz*rhgtf(1)
        rS(0)=re0
        rC(1)=zz*rhgtf(4)
        rS(1)=re0
        rC(2)=-5.D-1*zz*rhgtf(3)
        rS(2)=-5.D-1*zz*rhgtf(2)
      elseif(lt.eq.2)then 
        rC(0)=zz*rhgtf(1)
        rS(0)=re0
        rC(1)=zz*rhgtf(4)
        rS(1)=re0
        rC(2)=-5.D-1*zz*rhgtf(3)
        rS(2)=-5.D-1*zz*rhgtf(2)
        rC(3)=zz*(-5.D-1*rhgtf(5)-5.D-1*rhgtf(7)+rhgtf(10))
        rS(3)=re0
        rC(4)=-5.D-1*zz*rhgtf(9)
        rS(4)=-5.D-1*zz*rhgtf(8)
        rC(5)=zz*(-2.5D-1*rhgtf(5)+2.5D-1*rhgtf(7))
        rS(5)=2.5D-1*zz*rhgtf(6)
      elseif(lt.eq.3)then 
        rC(0)=zz*rhgtf(1)
        rS(0)=re0
        rC(1)=zz*rhgtf(4)
        rS(1)=re0
        rC(2)=-5.D-1*zz*rhgtf(3)
        rS(2)=-5.D-1*zz*rhgtf(2)
        rC(3)=zz*(-5.D-1*rhgtf(5)-5.D-1*rhgtf(7)+rhgtf(10))
        rS(3)=re0
        rC(4)=-5.D-1*zz*rhgtf(9)
        rS(4)=-5.D-1*zz*rhgtf(8)
        rC(5)=zz*(-2.5D-1*rhgtf(5)+2.5D-1*rhgtf(7))
        rS(5)=2.5D-1*zz*rhgtf(6)
        rC(6)=zz*(-5.D-1*rhgtf(15)-5.D-1*rhgtf(17)+rhgtf(20))
        rS(6)=re0
        rC(7)=zz*(1.25D-1*rhgtf(12)+3.75D-1*rhgtf(14)-5.D-1*rhgtf(19)
     >   )
        rS(7)=zz*(3.75D-1*rhgtf(11)+1.25D-1*rhgtf(13)-5.D-1*rhgtf(18)
     >   )
        rC(8)=zz*(-2.5D-1*rhgtf(15)+2.5D-1*rhgtf(17))
        rS(8)=2.5D-1*zz*rhgtf(16)
        rC(9)=zz*(1.25D-1*rhgtf(12)-1.25D-1*rhgtf(14))
        rS(9)=zz*(1.25D-1*rhgtf(11)-1.25D-1*rhgtf(13))
      elseif(lt.eq.4)then 
        rC(0)=zz*rhgtf(1)
        rS(0)=re0
        rC(1)=zz*rhgtf(4)
        rS(1)=re0
        rC(2)=-5.D-1*zz*rhgtf(3)
        rS(2)=-5.D-1*zz*rhgtf(2)
        rC(3)=zz*(-5.D-1*rhgtf(5)-5.D-1*rhgtf(7)+rhgtf(10))
        rS(3)=re0
        rC(4)=-5.D-1*zz*rhgtf(9)
        rS(4)=-5.D-1*zz*rhgtf(8)
        rC(5)=zz*(-2.5D-1*rhgtf(5)+2.5D-1*rhgtf(7))
        rS(5)=2.5D-1*zz*rhgtf(6)
        rC(6)=zz*(-5.D-1*rhgtf(15)-5.D-1*rhgtf(17)+rhgtf(20))
        rS(6)=re0
        rC(7)=zz*(1.25D-1*rhgtf(12)+3.75D-1*rhgtf(14)-5.D-1*rhgtf(19)
     >   )
        rS(7)=zz*(3.75D-1*rhgtf(11)+1.25D-1*rhgtf(13)-5.D-1*rhgtf(18)
     >   )
        rC(8)=zz*(-2.5D-1*rhgtf(15)+2.5D-1*rhgtf(17))
        rS(8)=2.5D-1*zz*rhgtf(16)
        rC(9)=zz*(1.25D-1*rhgtf(12)-1.25D-1*rhgtf(14))
        rS(9)=zz*(1.25D-1*rhgtf(11)-1.25D-1*rhgtf(13))
        rC(10)=zz*(3.75D-1*rhgtf(21)+1.25D-1*rhgtf(23)+3.75D-1*rhgtf(
     >   25)-5.D-1*rhgtf(30)-5.D-1*rhgtf(32)+rhgtf(35))
        rS(10)=re0
        rC(11)=zz*(1.25D-1*rhgtf(27)+3.75D-1*rhgtf(29)-5.D-1*rhgtf(34
     >   ))
        rS(11)=zz*(3.75D-1*rhgtf(26)+1.25D-1*rhgtf(28)-5.D-1*rhgtf(33
     >   ))
        rC(12)=zz*(2.5D-1*rhgtf(21)-2.5D-1*rhgtf(25)-2.5D-1*rhgtf(30)
     >   +2.5D-1*rhgtf(32))
        rS(12)=zz*(-1.25D-1*rhgtf(22)-1.25D-1*rhgtf(24)+2.5D-1*rhgtf(
     >   31))
        rC(13)=zz*(1.25D-1*rhgtf(27)-1.25D-1*rhgtf(29))
        rS(13)=zz*(1.25D-1*rhgtf(26)-1.25D-1*rhgtf(28))
        rC(14)=zz*(6.25D-2*rhgtf(21)-6.25D-2*rhgtf(23)+6.25D-2*rhgtf(
     >   25))
        rS(14)=zz*(-6.25D-2*rhgtf(22)+6.25D-2*rhgtf(24))
      elseif(lt.eq.5)then 
        rC(0)=zz*rhgtf(1)
        rS(0)=re0
        rC(1)=zz*rhgtf(4)
        rS(1)=re0
        rC(2)=-5.D-1*zz*rhgtf(3)
        rS(2)=-5.D-1*zz*rhgtf(2)
        rC(3)=zz*(-5.D-1*rhgtf(5)-5.D-1*rhgtf(7)+rhgtf(10))
        rS(3)=re0
        rC(4)=-5.D-1*zz*rhgtf(9)
        rS(4)=-5.D-1*zz*rhgtf(8)
        rC(5)=zz*(-2.5D-1*rhgtf(5)+2.5D-1*rhgtf(7))
        rS(5)=2.5D-1*zz*rhgtf(6)
        rC(6)=zz*(-5.D-1*rhgtf(15)-5.D-1*rhgtf(17)+rhgtf(20))
        rS(6)=re0
        rC(7)=zz*(1.25D-1*rhgtf(12)+3.75D-1*rhgtf(14)-5.D-1*rhgtf(19)
     >   )
        rS(7)=zz*(3.75D-1*rhgtf(11)+1.25D-1*rhgtf(13)-5.D-1*rhgtf(18)
     >   )
        rC(8)=zz*(-2.5D-1*rhgtf(15)+2.5D-1*rhgtf(17))
        rS(8)=2.5D-1*zz*rhgtf(16)
        rC(9)=zz*(1.25D-1*rhgtf(12)-1.25D-1*rhgtf(14))
        rS(9)=zz*(1.25D-1*rhgtf(11)-1.25D-1*rhgtf(13))
        rC(10)=zz*(3.75D-1*rhgtf(21)+1.25D-1*rhgtf(23)+3.75D-1*rhgtf(
     >   25)-5.D-1*rhgtf(30)-5.D-1*rhgtf(32)+rhgtf(35))
        rS(10)=re0
        rC(11)=zz*(1.25D-1*rhgtf(27)+3.75D-1*rhgtf(29)-5.D-1*rhgtf(34
     >   ))
        rS(11)=zz*(3.75D-1*rhgtf(26)+1.25D-1*rhgtf(28)-5.D-1*rhgtf(33
     >   ))
        rC(12)=zz*(2.5D-1*rhgtf(21)-2.5D-1*rhgtf(25)-2.5D-1*rhgtf(30)
     >   +2.5D-1*rhgtf(32))
        rS(12)=zz*(-1.25D-1*rhgtf(22)-1.25D-1*rhgtf(24)+2.5D-1*rhgtf(
     >   31))
        rC(13)=zz*(1.25D-1*rhgtf(27)-1.25D-1*rhgtf(29))
        rS(13)=zz*(1.25D-1*rhgtf(26)-1.25D-1*rhgtf(28))
        rC(14)=zz*(6.25D-2*rhgtf(21)-6.25D-2*rhgtf(23)+6.25D-2*rhgtf(
     >   25))
        rS(14)=zz*(-6.25D-2*rhgtf(22)+6.25D-2*rhgtf(24))
        rC(15)=zz*(3.75D-1*rhgtf(42)+1.25D-1*rhgtf(44)+3.75D-1*rhgtf(
     >   46)-5.D-1*rhgtf(51)-5.D-1*rhgtf(53)+rhgtf(56))
        rS(15)=re0
        rC(16)=zz*(-6.25D-2*rhgtf(37)-6.25D-2*rhgtf(39)-3.125D-1*rhgt
     >   f(41)+1.25D-1*rhgtf(48)+3.75D-1*rhgtf(50)-5.D-1*rhgtf(55))
        rS(16)=zz*(-3.125D-1*rhgtf(36)-6.25D-2*rhgtf(38)-6.25D-2*rhgt
     >   f(40)+3.75D-1*rhgtf(47)+1.25D-1*rhgtf(49)-5.D-1*rhgtf(54))
        rC(17)=zz*(2.5D-1*rhgtf(42)-2.5D-1*rhgtf(46)-2.5D-1*rhgtf(51)
     >   +2.5D-1*rhgtf(53))
        rS(17)=zz*(-1.25D-1*rhgtf(43)-1.25D-1*rhgtf(45)+2.5D-1*rhgtf(
     >   52))
        rC(18)=zz*(-9.375D-2*rhgtf(37)-3.125D-2*rhgtf(39)+1.5625D-1*r
     >   hgtf(41)+1.25D-1*rhgtf(48)-1.25D-1*rhgtf(50))
        rS(18)=zz*(-1.5625D-1*rhgtf(36)+3.125D-2*rhgtf(38)+9.375D-2*r
     >   hgtf(40)+1.25D-1*rhgtf(47)-1.25D-1*rhgtf(49))
        rC(19)=zz*(6.25D-2*rhgtf(42)-6.25D-2*rhgtf(44)+6.25D-2*rhgtf(
     >   46))
        rS(19)=zz*(-6.25D-2*rhgtf(43)+6.25D-2*rhgtf(45))
        rC(20)=zz*(-3.125D-2*rhgtf(37)+3.125D-2*rhgtf(39)-3.125D-2*rh
     >   gtf(41))
        rS(20)=zz*(-3.125D-2*rhgtf(36)+3.125D-2*rhgtf(38)-3.125D-2*rh
     >   gtf(40))
      elseif(lt.eq.6)then 
        rC(0)=zz*rhgtf(1)
        rS(0)=re0
        rC(1)=zz*rhgtf(4)
        rS(1)=re0
        rC(2)=-5.D-1*zz*rhgtf(3)
        rS(2)=-5.D-1*zz*rhgtf(2)
        rC(3)=zz*(-5.D-1*rhgtf(5)-5.D-1*rhgtf(7)+rhgtf(10))
        rS(3)=re0
        rC(4)=-5.D-1*zz*rhgtf(9)
        rS(4)=-5.D-1*zz*rhgtf(8)
        rC(5)=zz*(-2.5D-1*rhgtf(5)+2.5D-1*rhgtf(7))
        rS(5)=2.5D-1*zz*rhgtf(6)
        rC(6)=zz*(-5.D-1*rhgtf(15)-5.D-1*rhgtf(17)+rhgtf(20))
        rS(6)=re0
        rC(7)=zz*(1.25D-1*rhgtf(12)+3.75D-1*rhgtf(14)-5.D-1*rhgtf(19)
     >   )
        rS(7)=zz*(3.75D-1*rhgtf(11)+1.25D-1*rhgtf(13)-5.D-1*rhgtf(18)
     >   )
        rC(8)=zz*(-2.5D-1*rhgtf(15)+2.5D-1*rhgtf(17))
        rS(8)=2.5D-1*zz*rhgtf(16)
        rC(9)=zz*(1.25D-1*rhgtf(12)-1.25D-1*rhgtf(14))
        rS(9)=zz*(1.25D-1*rhgtf(11)-1.25D-1*rhgtf(13))
        rC(10)=zz*(3.75D-1*rhgtf(21)+1.25D-1*rhgtf(23)+3.75D-1*rhgtf(
     >   25)-5.D-1*rhgtf(30)-5.D-1*rhgtf(32)+rhgtf(35))
        rS(10)=re0
        rC(11)=zz*(1.25D-1*rhgtf(27)+3.75D-1*rhgtf(29)-5.D-1*rhgtf(34
     >   ))
        rS(11)=zz*(3.75D-1*rhgtf(26)+1.25D-1*rhgtf(28)-5.D-1*rhgtf(33
     >   ))
        rC(12)=zz*(2.5D-1*rhgtf(21)-2.5D-1*rhgtf(25)-2.5D-1*rhgtf(30)
     >   +2.5D-1*rhgtf(32))
        rS(12)=zz*(-1.25D-1*rhgtf(22)-1.25D-1*rhgtf(24)+2.5D-1*rhgtf(
     >   31))
        rC(13)=zz*(1.25D-1*rhgtf(27)-1.25D-1*rhgtf(29))
        rS(13)=zz*(1.25D-1*rhgtf(26)-1.25D-1*rhgtf(28))
        rC(14)=zz*(6.25D-2*rhgtf(21)-6.25D-2*rhgtf(23)+6.25D-2*rhgtf(
     >   25))
        rS(14)=zz*(-6.25D-2*rhgtf(22)+6.25D-2*rhgtf(24))
        rC(15)=zz*(3.75D-1*rhgtf(42)+1.25D-1*rhgtf(44)+3.75D-1*rhgtf(
     >   46)-5.D-1*rhgtf(51)-5.D-1*rhgtf(53)+rhgtf(56))
        rS(15)=re0
        rC(16)=zz*(-6.25D-2*rhgtf(37)-6.25D-2*rhgtf(39)-3.125D-1*rhgt
     >   f(41)+1.25D-1*rhgtf(48)+3.75D-1*rhgtf(50)-5.D-1*rhgtf(55))
        rS(16)=zz*(-3.125D-1*rhgtf(36)-6.25D-2*rhgtf(38)-6.25D-2*rhgt
     >   f(40)+3.75D-1*rhgtf(47)+1.25D-1*rhgtf(49)-5.D-1*rhgtf(54))
        rC(17)=zz*(2.5D-1*rhgtf(42)-2.5D-1*rhgtf(46)-2.5D-1*rhgtf(51)
     >   +2.5D-1*rhgtf(53))
        rS(17)=zz*(-1.25D-1*rhgtf(43)-1.25D-1*rhgtf(45)+2.5D-1*rhgtf(
     >   52))
        rC(18)=zz*(-9.375D-2*rhgtf(37)-3.125D-2*rhgtf(39)+1.5625D-1*r
     >   hgtf(41)+1.25D-1*rhgtf(48)-1.25D-1*rhgtf(50))
        rS(18)=zz*(-1.5625D-1*rhgtf(36)+3.125D-2*rhgtf(38)+9.375D-2*r
     >   hgtf(40)+1.25D-1*rhgtf(47)-1.25D-1*rhgtf(49))
        rC(19)=zz*(6.25D-2*rhgtf(42)-6.25D-2*rhgtf(44)+6.25D-2*rhgtf(
     >   46))
        rS(19)=zz*(-6.25D-2*rhgtf(43)+6.25D-2*rhgtf(45))
        rC(20)=zz*(-3.125D-2*rhgtf(37)+3.125D-2*rhgtf(39)-3.125D-2*rh
     >   gtf(41))
        rS(20)=zz*(-3.125D-2*rhgtf(36)+3.125D-2*rhgtf(38)-3.125D-2*rh
     >   gtf(40))
        rC(21)=zz*(-3.125D-1*rhgtf(57)-6.25D-2*rhgtf(59)-6.25D-2*rhgt
     >   f(61)-3.125D-1*rhgtf(63)+3.75D-1*rhgtf(70)+1.25D-1*rhgtf(72)
     >   +3.75D-1*rhgtf(74)-5.D-1*rhgtf(79)-5.D-1*rhgtf(81)+rhgtf(84)
     >   )
        rS(21)=re0
        rC(22)=zz*(-6.25D-2*rhgtf(65)-6.25D-2*rhgtf(67)-3.125D-1*rhgt
     >   f(69)+1.25D-1*rhgtf(76)+3.75D-1*rhgtf(78)-5.D-1*rhgtf(83))
        rS(22)=zz*(-3.125D-1*rhgtf(64)-6.25D-2*rhgtf(66)-6.25D-2*rhgt
     >   f(68)+3.75D-1*rhgtf(75)+1.25D-1*rhgtf(77)-5.D-1*rhgtf(82))
        rC(23)=zz*(-2.34375D-1*rhgtf(57)-1.5625D-2*rhgtf(59)+1.5625D-
     >   2*rhgtf(61)+2.34375D-1*rhgtf(63)+2.5D-1*rhgtf(70)-2.5D-1*rhg
     >   tf(74)-2.5D-1*rhgtf(79)+2.5D-1*rhgtf(81))
        rS(23)=zz*(7.8125D-2*rhgtf(58)+4.6875D-2*rhgtf(60)+7.8125D-2*
     >   rhgtf(62)-1.25D-1*rhgtf(71)-1.25D-1*rhgtf(73)+2.5D-1*rhgtf(8
     >   0))
        rC(24)=zz*(-9.375D-2*rhgtf(65)-3.125D-2*rhgtf(67)+1.5625D-1*r
     >   hgtf(69)+1.25D-1*rhgtf(76)-1.25D-1*rhgtf(78))
        rS(24)=zz*(-1.5625D-1*rhgtf(64)+3.125D-2*rhgtf(66)+9.375D-2*r
     >   hgtf(68)+1.25D-1*rhgtf(75)-1.25D-1*rhgtf(77))
        rC(25)=zz*(-9.375D-2*rhgtf(57)+3.125D-2*rhgtf(59)+3.125D-2*rh
     >   gtf(61)-9.375D-2*rhgtf(63)+6.25D-2*rhgtf(70)-6.25D-2*rhgtf(7
     >   2)+6.25D-2*rhgtf(74))
        rS(25)=zz*(6.25D-2*rhgtf(58)-6.25D-2*rhgtf(62)-6.25D-2*rhgtf(
     >   71)+6.25D-2*rhgtf(73))
        rC(26)=zz*(-3.125D-2*rhgtf(65)+3.125D-2*rhgtf(67)-3.125D-2*rh
     >   gtf(69))
        rS(26)=zz*(-3.125D-2*rhgtf(64)+3.125D-2*rhgtf(66)-3.125D-2*rh
     >   gtf(68))
        rC(27)=zz*(-1.5625D-2*rhgtf(57)+1.5625D-2*rhgtf(59)-1.5625D-2
     >   *rhgtf(61)+1.5625D-2*rhgtf(63))
        rS(27)=zz*(1.5625D-2*rhgtf(58)-1.5625D-2*rhgtf(60)+1.5625D-2*
     >   rhgtf(62))
      else 
        write(*,*)' lt = ',lt 
        stop ' bad indexing in QCrc ' 
      endif 
