c
c===================================================================
c
c     DEMO OF LATTICE SUM COMPUTATION WITH EWALD PARTITION OF THE
c     BIPOLAR MULTIPOLE EXPANSION.  
c
c     MATT CHALLACOMBE AND MARTIN HEAD-GORDON NOVEMBER 1996
c
c===================================================================
c
      program lattice_sum_demo
c
      implicit double precision (a-h,o-z)     
c
      parameter(mxell=128)
c
      complex*16 Mpq(0:mxell,-mxell:mxell)
c
      parameter(m=216,n=2)
      dimension rx(m),ry(m),rz(m),z(m)
      integer ell
c
      pi=2.0d0*dacos(0.0d0)
c
c===================================================================
c
c     SET UP THE UNIT CELL
c
c===================================================================
c
c     Note: If delta is to large, it will degrade the accuracy of the
c     conventional ewald method (at least with this simple minded 
c     implimentation).
c
      delta=3.0*1.889725988578923D0
      d2=delta*0.50
      chg=2.0d0
c
      write(*,*)' square unit cell of side ',delta
      write(*,*)' max charge in cell = ',chg
c
      iseed=10
        rx(1)=0.d0
        ry(1)=0.d0
        rz(1)=0.d0
        rx(2)=2.834588982868384D0
        ry(2)=2.834588982868384D0
        rz(2)=2.834588982868384D0
!       rx(1)=delta/4.d0
!       ry(1)=delta/4.d0
!       rz(1)=delta/4.d0
!       rx(2)=3.d0*delta/4.d0
!       ry(2)=3.d0*delta/4.d0
!       rz(2)=3.d0*delta/4.d0
!      do i=1,n
!       rx(i)=ran(iseed)*d2
!       ry(i)=ran(iseed)*d2
!       rz(i)=ran(iseed)*d2
!     enddo
c
        z(1)=2.D0
        z(2)=-2.D0
!     do i=1,n/2
!       z(i)=ran(iseed)*chg
!       z(i+n/2)=-z(i)
!     enddo
c
      do i=1,n
         write(*,*)rx(i),ry(i),rz(i),z(i)
      enddo
c
c===================================================================
c
c     LATTICE SUMS VIA THE EWALD PARTITION OF THE BIPOLAR EXPANSION
c
c===================================================================
c
c     Compute interaction tensors, and write them to fort.22;
c     these values are subsequently included in subroutine 
c     load_interaction_tensors.
c
c        ell=40
c        call compute_interaction_tensors(2*ell,1.0d0,Mpq)
c        do l=4,2*ell,2
c           do mv=0,l,4
c              realm=Mpq(l,mv)*(1.0d0,0.0d0)
c              write(22,88)l,mv,realm
c           enddo
c        enddo
c        stop
c      endif
c
c     Demonstrate convergence of the bipolar expansion with 
c     increasing order.
c
c       do ell=4,40,4
c          call compute_interaction_tensors(2*ell,delta,Mpq)
c          call load_interaction_tensors(2*ell,delta,Mpq)
c          call compute_cell_energy(n,ell,delta,Mpq,rx,ry,rz,z,alphaM)
c          write(*,44)ell,alphaM
c       enddo
c
c     Find the conventional ewald sum.
c
      WRITE(*,*) "Computing Ewald Sum"
      call conventional_ewald_sum(m,n,delta,rx,ry,rz,z,alphaM)
c
      write(*,55)alphaM
c
c     For ammusment, try and do the calculation by summing boxes (interupt). 
c
c      do mbox=10,100,10
c         call direct_sum(n,delta,mbox,rx,ry,rz,z,alphaM)
c         write(*,66)alphaM
c      enddo
c
      stop
c
 44   format('  ell = ',i2,' ewald-multipole sum = ',D23.16)
 55   format('  the correct lattice sum is     ',D23.16)
 66   format('  the direct lattice sum is      ',D23.16)         
 88   format('      Mpq(',i2,',',i2,')=',D26.18)
c
      end
c
      real*8 function ran(jran)
      implicit real*8 (a-h,o-z)
      parameter (ia=16807,im=2147483647,am=1.0d0/im)
      parameter (iq=127773,ir=2836,mask=123459876)
      jran=mod(jran*ia+ic,im)
      ran=float(jran)/float(im)
      return
      end      
c
      subroutine load_interaction_tensors(ell,delta,Mpq)
c
      implicit real*8 (a-h,o-z)
c     
      parameter(mxell=128)    
c
      complex*16 Mpq(0:mxell,-mxell:mxell)
c
      integer ell
c
      do l=0,ell
         do m=-l,l
            Mpq(l,m)=0.0d0
         enddo
      enddo
c
      include "fort.22"
c
      do l=4,ell,2
         dl=delta**(-l-1)
         do m=0,l,4             
            Mpq(l,m)=Mpq(l,m)*dl
            Mpq(l,-m)=Mpq(l,m)
         enddo
      enddo
c
      return
      end   
c
      subroutine compute_interaction_tensors(ell,delta,Mpq)
c
      implicit real*8 (a-h,o-z)
c     
      parameter(mxell=128)
      parameter(beta=1.5d0,beta2=beta*beta)
      parameter(mxbox=150,tiny=1.d-60)
c
      dimension factorial(0:mxell)
      dimension gammas(0:mxell)
c
      real*8 kx,ky,kz,k2,k
c      
      complex*16 tmp(0:mxell,-mxell:mxell),Mpq(0:mxell,-mxell:mxell)
      complex*16 cmplxL,cmplxI,prefactor(0:50)
c
      integer ell
c
      logical true
c
      pi=2.0d0*dacos(0.0d0)
      const=dsqrt(2.0d0*pi/delta)**3
      expnt=1.0d0/(4.0d0*beta**2)
      recip=2.0d0*pi/delta
c
      factorial(0)=1.d0
      do m=1,ell
         factorial(m)=factorial(m-1)*dfloat(m)         
      enddo
c
      cmplxI=(0.0d0,1.0d0)
      cmplxL=(1.0d0,0.0d0)
c
      fr1=0.50d0
      sq2=dsqrt(2.0d0)
      gamma=dsqrt(pi)
      do l=0,ell
         prefactor(l)=const*cmplxL/(gamma*2.0d0**(l-0.5))
         gamma=fr1*gamma
         fr1=fr1+1.0d0
         sq2=sq2*0.50d0
         cmplxL=cmplxL*cmplxI
      enddo     
c
       do l=0,ell
          do m=-l,l
             Mpq(l,m)=0.0
          enddo
       enddo
c
      do mx=-mxbox,mxbox
         write(*,*)' mx = ',mx
         do my=-mxbox,mxbox
            do mz=-mxbox,mxbox
               kx=dfloat(mx)*recip
               ky=dfloat(my)*recip
               kz=dfloat(mz)*recip
               k2=kx**2+ky**2+kz**2
               term=dexp(-expnt*k2)
               if(k2.eq.0.0.or.term.lt.tiny)goto 111
               count=count+1
               k=dsqrt(k2)     
               call calculate_mlms(ell,kx,ky,kz,.true.,tmp)
               do l=0,ell
                  tt=term*k**(l-2)
                  do m=-l,l                     
                     Mpq(l,m)=Mpq(l,m)+tmp(l,m)*tt       
                  enddo
               enddo
 111          continue
            enddo
         enddo
      enddo      
c
       do l=0,ell
          do m=-l,l
             Mpq(l,m)=prefactor(l)*Mpq(l,m)
          enddo
       enddo
c
      do nx=-mxbox,mxbox
         write(*,*)' nx = ',nx
         do ny=-mxbox,mxbox
            do nz=-mxbox,mxbox
               nxa=abs(nx)
               nya=abs(ny)
               nza=abs(nz)
               rnx=dfloat(nx)
               rny=dfloat(ny)
               rnz=dfloat(nz)
               r2=(rnx**2+rny**2+rnz**2)*delta
               true=nxa.le.1.and.nya.le.1.and.nza.le.1
               term=dexp(-beta2*r2)
               true=true.or.term.lt.tiny
               if(true)goto 222
               r=dsqrt(r2)
               call make_gammas(ell,r*beta,gammas)
               call calculate_mlms(ell,rnx,rny,rnz,.false.,tmp)
                do l=0,ell
                  rl=r**(-(l+1))
                  do m=-l,l
                     Mpq(l,m)=Mpq(l,m)+tmp(l,m)*gammas(l)
                  enddo
               enddo
 222           continue
            enddo
         enddo
      enddo
c
      do nx=-1,1
         do ny=-1,1
            do nz=-1,1
               rnx=dfloat(nx)
               rny=dfloat(ny)
               rnz=dfloat(nz)
               r=sqrt(rnx**2+rny**2+rnz**2)*delta
               if(r.eq.0.d0)goto 333
               call make_gammas2(ell,r*beta,gammas)
               call calculate_mlms(ell,rnx,rny,rnz,.false.,tmp)
               do l=0,ell
                  rl=r**(-(l+1))
                  do m=-l,l
                     Mpq(l,m)=Mpq(l,m)-tmp(l,m)*gammas(l)
                  enddo
               enddo
 333           continue
            enddo
         enddo
      enddo
c
      return
      end   
c
      subroutine compute_cell_energy
     >           (n,ell,delta,Mpq,rx,ry,rz,z,alphaM)
c
      implicit double precision (a-h,o-z)     
c
      dimension rx(n),ry(n),rz(n),z(n)
c
      integer ell
      parameter (mxell=128)
c
      complex*16 Qp(0:mxell,-mxell:mxell),tmp(0:mxell,-mxell:mxell),
     >           Mpq(0:mxell,-mxell:mxell)
c
       do l=0,ell
          do m=-l,l
             Qp(l,m)=0.0
          enddo
       enddo
c
       do i=1,n
          x =rx(i)
          y =ry(i)
          zz=rz(i)
          q=z(i)
          do l=0,ell
             do m=-l,l
                tmp(l,m)=0.0
             enddo
          enddo
          call calculate_olms(ell,x,y,zz,tmp)
          do l=0,ell
             do m=-l,l
                Qp(l,m)=Qp(l,m)+q*tmp(l,m)
             enddo
          enddo
       enddo
c
      e_if=0.d0
c
      do lp=0,ell
         do lq=0,ell
            do mp=-lp,lp
                do mq=-lq,lq
                  e_if=e_if+(-1)**(lp+mp+mq)
     >            *Qp(lp,mp)*Mpq(lp+lq,-mp-mq)*Qp(lq,mq)         
                enddo
             enddo
          enddo
       enddo
c

      e_nf=0.d0
c      
      do nx=-1,1
         do ny=-1,1
            do nz=-1,1
               rnx=nx*delta
               rny=ny*delta
               rnz=nz*delta
               do i=1,n
                  do j=1,n
                     x =rnx+rx(i)-rx(j)
                     y =rny+ry(i)-ry(j)
                     zz=rnz+rz(i)-rz(j)
                     r=sqrt(x**2+y**2+zz**2)
                     if(r.ne.0.0)then
                       e_nf=e_nf+z(i)*z(j)/r
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      e_if=e_if/2.0d0
      e_nf=e_nf/2.0d0
c      
       write(*,*)' image field contribution = ',e_if
       write(*,*)' near  filed contribution = ',e_nf
c
      alphaM=e_if+e_nf
c
      return
      end
c
c=================================================================
c
c     LATTICE SUMS VIA THE CONVENTIONAL EWALD PARTITION
c
c=================================================================
c
      subroutine conventional_ewald_sum
     >           (m,n,delta,rx,ry,rz,z,alphaM)
c
      implicit double precision (a-h,o-z)     
c
      parameter(beta=1.5d0,beta2=beta*beta)
      parameter(mxbox=600,tiny=1.d-60)
c
      dimension rx(m),ry(m),rz(m),z(m)
c
      pi=2.0d0*dacos(0.0d0)

      const=4.0d0*pi/delta**3
      expnt=1.0d0/(4.0d0*beta**2)
      recip=2.0d0*pi/delta
      recip2=recip*recip
c
      recip_term=0.0d0
c
      do nx=-mxbox,mxbox         
         icount=0
         do ny=-mxbox,mxbox
            do nz=-mxbox,mxbox
               rnx=float(nx)*recip
               rny=float(ny)*recip
               rnz=float(nz)*recip
               rk2n=rnx**2+rny**2+rnz**2
               term=dexp(-expnt*rk2n)
               if(term.lt.tiny.or.rk2n.eq.0.0d0)goto 111
               icount=icount+1
               term=term/rk2n
               do i=1,n
                  do j=1,n
                    x =rx(i)-rx(j)
                    y =ry(i)-ry(j)
                    zz=rz(i)-rz(j)
                    dot=x*rnx+y*rny+zz*rnz
                    phase=dcos(dot)
                    recip_term=recip_term+term*phase*z(i)*z(j)
                  enddo
               enddo
 111           continue
            enddo
         enddo
         write(*,*)' nx ',nx,' count ',icount
      enddo      
c
      recip_term=0.5d0*const*recip_term
c
      self_term=0.0d0
      do i=1,n
         self_term=self_term+z(i)*z(i)
      enddo
c
      recip_term=recip_term-self_term*beta/dsqrt(pi)
c
      realspc_term=0.0d0
c
      do nx=-mxbox,mxbox
         icount=0
         do ny=-mxbox,mxbox
            do nz=-mxbox,mxbox
               rnx=nx*delta
               rny=ny*delta
               rnz=nz*delta
               rn2=rnx**2+rny**2+rnzz**2
               term=dexp(-rn2*beta2)
               if(term.lt.tiny)goto 222
               icount=icount+1
               rn=dsqrt(rn2)
               do i=1,n
                  do j=1,n
                    x =rnx+rx(i)-rx(j)
                    y =rny+ry(i)-ry(j)
                    zz=rnz+rz(i)-rz(j)
                    r=sqrt(x**2+y**2+zz**2)
                    if(r.ne.0.0d0)then
                       realspc_term=realspc_term
     >                             +z(i)*z(j)*derfc(beta*r)/r
                    endif
                  enddo
               enddo
 222           continue
            enddo
         enddo
         write(*,*)' nx ',nx,' count ',icount
      enddo
c
      realspc_term=realspc_term*0.5d0
c
      dx=0.d0
      dy=0.d0
      dz=0.d0
      do i=1,n
         dx=dx+z(i)*rx(i)
         dy=dy+z(i)*ry(i)
         dz=dz+z(i)*rz(i)
      enddo
      dp=dx*dx+dy*dy+dz*dz
c      dipole_correction=2.0d0*pi*dp/(3.0d0*delta**3)
c
c      write(*,*)' real space = ',realspc_term
c      write(*,*)' recp space = ',recip_term
c
      alphaM=recip_term+realspc_term
c
      return
      end
c
      subroutine direct_sum(n,delta,l,rx,ry,rz,z,alphaM)
c
      implicit double precision (a-h,o-z)     
c
      dimension rx(n),ry(n),rz(n),z(n)
c
      vr=0.0
      do nx=-l,l
         do ny=-l,l
            do nz=-l,l
               rnx=nx*delta
               rny=ny*delta
               rnz=nz*delta
               do i=1,n
                  do j=1,n
                    x =rnx+rx(i)-rx(j)
                    y =rny+ry(i)-ry(j)
                    zz=rnz+rz(i)-rz(j)
                    r=sqrt(x**2+y**2+zz**2)
                    if(r.ne.0.0)then
                       vr=vr+z(i)*z(j)/r
                    endif
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      alphaM=vr/2.0d0
c
      return
      end
c
c````````````````````````````````````````````````````````````````````````````
c
      subroutine make_gammas(l,x,gammas)
c
      parameter(mxl=40)
c
      implicit double precision (a-h,o-z)           
      dimension gammas(0:l)
c
      pi=2.0d0*dacos(0.d0)
      fr1=0.50d0
      x1=x
      x2=x*x 
      ex2=dexp(-x2)
      gamma1=dsqrt(pi)*derfc(x)
      gamma2=dsqrt(pi)
      gammas(0)=gamma1/gamma2
      do n=1,l
         gamma1=fr1*gamma1+x1*ex2
         gamma2=fr1*gamma2
         x1=x1*x2
         fr1=fr1+1.0d0
         gammas(n)=gamma1/gamma2
      enddo
c
      return
      end
c
      subroutine make_gammas2(l,x,gammas)
c
      parameter(mxl=40)
c
      implicit double precision (a-h,o-z)           
      dimension gammas(0:l)
c
      pi=2.0d0*dacos(0.d0)
      fr1=0.50d0
      x1=x
      x2=x*x 
      ex2=dexp(-x2)
      gamma1=dsqrt(pi)*derfc(x)
      gamma2=dsqrt(pi)
      gammas(0)=(gamma2-gamma1)/gamma2
      do n=1,l
         gamma1=fr1*gamma1+x1*ex2
         gamma2=fr1*gamma2
         x1=x1*x2
         fr1=fr1+1.0d0
         gammas(n)=(gamma2-gamma1)/gamma2
      enddo
c
      return
      end
c
c====================================================================
c
      subroutine calculate_olms(ell,x,y,z,sh)
c
      implicit real*8 (a-h,o-z)
      integer ell
      parameter(mxell=128)
      dimension p(0:mxell,0:mxell),cp(-mxell:mxell),sp(-mxell:mxell)
      dimension factorial(0:2*mxell)
      complex*16 sh(0:mxell,-mxell:mxell)
c
      r=dsqrt(x*x+y*y+z*z)
      ct=z/r
      phi=datan2(x,y)
c
      cp(0)=1.0
      sp(0)=0.0
      cp(1)=dcos(phi)
      sp(1)=dsin(phi)
      do m=2,ell+2
         cp(m)=2.*cp(1)*cp(m-1)-cp(m-2)
         sp(m)=2.*cp(1)*sp(m-1)-sp(m-2)
      enddo
c
      factorial(0)=1.0d0
      do l=1,2*ell
         factorial(l)=factorial(l-1)*dfloat(l)
      enddo
c
      sq=dsqrt(1.d0-ct*ct)
      rs=1.d0
      sgn=1.0d0
      dblfact=1.d0
      twotimes=1.d0
      do m=0,ell
         cofact=sgn*dblfact/factorial(2*m)
         p(m,m)=cofact*rs
         p(m+1,m)=ct*p(m,m)
         rs=rs*sq
         dblfact=dblfact*twotimes
         twotimes=twotimes+2.0d0
         sgn=-sgn
      enddo
c
      do m=0,ell
         do l=m+2,ell         
            p(l,m)=(ct*dfloat(2*l-1)*p(l-1,m)
     >            -p(l-2,m))/dfloat((l+m)*(l-m))
         enddo
      enddo
c
      rr=1.d0
      do l=0,ell
         do m=0,l
            sh(l, m)=rr*p(l,m)*((1.0d0,0.0d0)*cp(m)
     >             +(0.0d0,1.0d0)*sp(m))
          enddo
         do m=1,l
            phs=(-1)**m
          sh(l,-m)=rr*phs*p(l,m)*((1.0d0,0.0d0)*cp(m)
     >           +(0.,-1.)*sp(m))
         enddo
         rr=rr*r
      enddo
c
c      do l=0,ell
c         do m=-l,l
c            rre=(1.,0.)*sh(l,m)
c            rim=(0.,-1.)*sh(l,m)
c            write(*,77)l,m,rre,rim
c         enddo
c      enddo
c
 77   format('one l = ',i2,' m = ',i2,'(',F12.8,',',F12.8,')')
c
      return
      end
c
      subroutine calculate_mlms(ell,x,y,z,nor,sh)
c
      implicit real*8 (a-h,o-z)
c
      integer ell
      parameter(mxell=128)
      dimension p(0:mxell,0:mxell),cp(-mxell:mxell),sp(-mxell:mxell)
      complex*16 sh(0:mxell,-mxell:mxell)
      logical nor
c
      if(ell.gt.mxell)then
        write(*,*)' ell = ',ell,' mxell = ',mxell
        stop ' in calc_sh2 '
      endif
c
      r=dsqrt(x*x+y*y+z*z)
      oneoverr=1.0d0/r
      ct=z*oneoverr
      phi=datan2(x,y)
c
      cp(0)=1.0d0
      sp(0)=0.0d0
      cp(1)=dcos(phi)
      sp(1)=dsin(phi)
      do m=2,ell+2
         cp(m)=2.0d0*cp(1)*cp(m-1)-cp(m-2)
         sp(m)=2.0d0*cp(1)*sp(m-1)-sp(m-2)
      enddo
c
      sq=sqrt(1.-ct*ct)
      rs=1.d0
      dblfact=1.d0
      twotimes=1.d0
      sgn=1.0d0
      do m=0,ell
         p(m,m)=sgn*dblfact*rs
         p(m+1,m)=ct*dfloat(2*m+1)*p(m,m)
         rs=rs*sq
         dblfact=dblfact*twotimes
         twotimes=twotimes+2.0d0
         sgn=-sgn
      enddo
c
      do m=0,ell
         do l=m+2,ell         
            p(l,m)=ct*dfloat(2*l-1)*p(l-1,m)
     >            -dfloat((l+m-1)*(l-m-1))*p(l-2,m)
         enddo
      enddo
c    
      if(nor)then
        do l=0,ell
           do m=0,l
              sh(l, m)=p(l,m)*((1.d0,0.d0)*cp(m)+ (0.d0,1.d0)*sp(m) )
           enddo
           do m=1,l
              phs=(-1)**m
              sh(l,-m)=phs*p(l,m)*((1.d0,0.d0)*cp(m)
     >                +(0.d0,-1.d0)*sp(m) )
           enddo
        enddo
      else
        rr=1.0d0
        do l=0,ell
           rr=rr*oneoverr
           do m=0,l
              sh(l, m)=rr*p(l,m)*((1.d0,0.d0)*cp(m)+ (0.d0,1.d0)*sp(m) )
           enddo
           do m=1,l
              phs=(-1)**m
              sh(l,-m)=rr*phs*p(l,m)*((1.d0,0.d0)*cp(m)
     >                +(0.d0,-1.d0)*sp(m) )
           enddo
        enddo
      endif
c
c      do l=0,ell
c         do m=-l,l
c            rre=(1.d0,0.d0)*sh(l,m)
c            rim=(0.d0,-1.d0)*sh(l,m)
c            write(*,77)l,m,rre,rim
c         enddo
c      enddo
c
 77   format('two l = ',i2,' m = ',i2,'(',F14.8,',',F14.8,')')
c
      return
      end











