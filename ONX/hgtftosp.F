c                                                                       
      SUBROUTINE HGTFToSP(lt,zz,rc,rs,rhgtf)                          
c                                                                       
      implicit real*8 (a-h,o-z)                                         
c                                                                       
      parameter (re0=0.d0)                                              
c                                                                       
      REAL*8 rc(0:*),rs(0:*),rhgtf(*)                                
c                                                                       
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
      else
        write(*,*)' lt = ',lt
        stop
      endif
      RETURN
      END                                                            
