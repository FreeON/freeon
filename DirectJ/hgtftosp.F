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
        rc(0)=zz*rhgtf(1)
        rs(0)=re0
      elseif(lt.eq.1)then 
        rc(0)=zz*rhgtf(1)
        rs(0)=re0
        rc(1)=zz*rhgtf(4)
        rs(1)=re0
        rc(2)=-5.d-1*zz*rhgtf(3)
        rs(2)=-5.d-1*zz*rhgtf(2)
      elseif(lt.eq.2)then 
        rc(0)=zz*rhgtf(1)
        rs(0)=re0
        rc(1)=zz*rhgtf(4)
        rs(1)=re0
        rc(2)=-5.d-1*zz*rhgtf(3)
        rs(2)=-5.d-1*zz*rhgtf(2)
        rc(3)=zz*(-5.d-1*rhgtf(5)-5.d-1*rhgtf(8)+rhgtf(10))
        rs(3)=re0
        rc(4)=-5.d-1*zz*rhgtf(9)
        rs(4)=-5.d-1*zz*rhgtf(7)
        rc(5)=zz*(-2.5d-1*rhgtf(5)+2.5d-1*rhgtf(8))
        rs(5)=2.5d-1*zz*rhgtf(6)
      elseif(lt.eq.3)then 
        rc(0)=zz*rhgtf(1)
        rs(0)=re0
        rc(1)=zz*rhgtf(4)
        rs(1)=re0
        rc(2)=-5.d-1*zz*rhgtf(3)
        rs(2)=-5.d-1*zz*rhgtf(2)
        rc(3)=zz*(-5.d-1*rhgtf(5)-5.d-1*rhgtf(8)+rhgtf(10))
        rs(3)=re0
        rc(4)=-5.d-1*zz*rhgtf(9)
        rs(4)=-5.d-1*zz*rhgtf(7)
        rc(5)=zz*(-2.5d-1*rhgtf(5)+2.5d-1*rhgtf(8))
        rs(5)=2.5d-1*zz*rhgtf(6)
        rc(6)=zz*(-5.d-1*rhgtf(13)-5.d-1*rhgtf(18)+rhgtf(20))
        rs(6)=re0
        rc(7)=zz*(1.25d-1*rhgtf(12)+3.75d-1*rhgtf(17)-5.d-1*rhgtf(19)
     >   )
        rs(7)=zz*(3.75d-1*rhgtf(11)+1.25d-1*rhgtf(14)-5.d-1*rhgtf(16)
     >   )
        rc(8)=zz*(-2.5d-1*rhgtf(13)+2.5d-1*rhgtf(18))
        rs(8)=2.5d-1*zz*rhgtf(15)
        rc(9)=zz*(1.25d-1*rhgtf(12)-1.25d-1*rhgtf(17))
        rs(9)=zz*(1.25d-1*rhgtf(11)-1.25d-1*rhgtf(14))
      elseif(lt.eq.4)then 
        rc(0)=zz*rhgtf(1)
        rs(0)=re0
        rc(1)=zz*rhgtf(4)
        rs(1)=re0
        rc(2)=-5.d-1*zz*rhgtf(3)
        rs(2)=-5.d-1*zz*rhgtf(2)
        rc(3)=zz*(-5.d-1*rhgtf(5)-5.d-1*rhgtf(8)+rhgtf(10))
        rs(3)=re0
        rc(4)=-5.d-1*zz*rhgtf(9)
        rs(4)=-5.d-1*zz*rhgtf(7)
        rc(5)=zz*(-2.5d-1*rhgtf(5)+2.5d-1*rhgtf(8))
        rs(5)=2.5d-1*zz*rhgtf(6)
        rc(6)=zz*(-5.d-1*rhgtf(13)-5.d-1*rhgtf(18)+rhgtf(20))
        rs(6)=re0
        rc(7)=zz*(1.25d-1*rhgtf(12)+3.75d-1*rhgtf(17)-5.d-1*rhgtf(19)
     >   )
        rs(7)=zz*(3.75d-1*rhgtf(11)+1.25d-1*rhgtf(14)-5.d-1*rhgtf(16)
     >   )
        rc(8)=zz*(-2.5d-1*rhgtf(13)+2.5d-1*rhgtf(18))
        rs(8)=2.5d-1*zz*rhgtf(15)
        rc(9)=zz*(1.25d-1*rhgtf(12)-1.25d-1*rhgtf(17))
        rs(9)=zz*(1.25d-1*rhgtf(11)-1.25d-1*rhgtf(14))
        rc(10)=zz*(3.75d-1*rhgtf(21)+1.25d-1*rhgtf(24)-5.d-1*rhgtf(26
     >   )+3.75d-1*rhgtf(31)-5.d-1*rhgtf(33)+rhgtf(35))
        rs(10)=re0
        rc(11)=zz*(1.25d-1*rhgtf(25)+3.75d-1*rhgtf(32)-5.d-1*rhgtf(34
     >   ))
        rs(11)=zz*(3.75d-1*rhgtf(23)+1.25d-1*rhgtf(28)-5.d-1*rhgtf(30
     >   ))
        rc(12)=zz*(2.5d-1*rhgtf(21)-2.5d-1*rhgtf(26)-2.5d-1*rhgtf(31)
     >   +2.5d-1*rhgtf(33))
        rs(12)=zz*(-1.25d-1*rhgtf(22)-1.25d-1*rhgtf(27)+2.5d-1*rhgtf(
     >   29))
        rc(13)=zz*(1.25d-1*rhgtf(25)-1.25d-1*rhgtf(32))
        rs(13)=zz*(1.25d-1*rhgtf(23)-1.25d-1*rhgtf(28))
        rc(14)=zz*(6.25d-2*rhgtf(21)-6.25d-2*rhgtf(24)+6.25d-2*rhgtf(
     >   31))
        rs(14)=zz*(-6.25d-2*rhgtf(22)+6.25d-2*rhgtf(27))
      else 
        write(*,*)' lt = ',lt 
        stop ' bad indexing in QCrc ' 
      endif 
c                                                                       
      RETURN
      END                                                            
