c:X_LDAsubrstart

c    Generated: Wed Jan 29 09:08:46 GMT 2003

      subroutine uks_x_lda
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     P.A.M. Dirac
c     Proceedings of the Cambridge Philosophical Society, 26 (1930) 376
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt,i
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      t1 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhob
      elseif(rhob.lt.tol) then
      rho = rhoa
      t1 = rhoa**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      t1 = rhoa**(1.D0/3.D0)
      t4 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa-0.9305257363491D0*t4*rhob
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      t1 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhob
      vrhoa(i) = 0.D0
      vrhob(i) = -0.12407009817988D1*t1
      elseif(rhob.lt.tol) then
      rho = rhoa
      t1 = rhoa**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa
      vrhoa(i) = -0.12407009817988D1*t1
      vrhob(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      t1 = rhoa**(1.D0/3.D0)
      t4 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa-0.9305257363491D0*t4*rhob
      vrhoa(i) = -0.12407009817988D1*t1
      vrhob(i) = -0.12407009817988D1*t4
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      t1 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhob
      vrhoa(i) = 0.D0
      vrhob(i) = -0.12407009817988D1*t1
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t5 = t1**2
      v2rhob2(i) = -0.4135669939329333D0/t5
      elseif(rhob.lt.tol) then
      rho = rhoa
      t1 = rhoa**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa
      vrhoa(i) = -0.12407009817988D1*t1
      vrhob(i) = 0.D0
      t5 = t1**2
      v2rhoa2(i) = -0.4135669939329333D0/t5
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      t1 = rhoa**(1.D0/3.D0)
      t4 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa-0.9305257363491D0*t4*rhob
      vrhoa(i) = -0.12407009817988D1*t1
      vrhob(i) = -0.12407009817988D1*t4
      t9 = t1**2
      v2rhoa2(i) = -0.4135669939329333D0/t9
      t12 = t4**2
      v2rhob2(i) = -0.4135669939329333D0/t12
      v2rhoab(i) = 0.D0
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_x_lda
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     P.A.M. Dirac
c     Proceedings of the Cambridge Philosophical Society, 26 (1930) 376
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt,i
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = rho**(1.D0/3.D0)
      zk(i) = -0.7385587663820224D0*t1*rho
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = rho**(1.D0/3.D0)
      zk(i) = -0.7385587663820224D0*t1*rho
      vrhoa(i) = -0.9847450218426965D0*t1
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = rho**(1.D0/3.D0)
      zk(i) = -0.7385587663820224D0*t1*rho
      vrhoa(i) = -0.9847450218426965D0*t1
      t5 = t1**2
      v2rhoa2(i) = -0.6564966812284644D0/t5
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end

c:X_LDAsubrend
c:XC_PW91subrstart

c    Generated: Mon Oct 25 09:31:43 BST 2004

      subroutine uks_xc_pw91
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson,
c     M.R. Pederson, D.J. Singh, C. Fiolhais
c     Atoms, molecules, solids and surfaces:
c     Applications of the generalized gradient approximation
c     for exchange and correlation
c     Phys. Rev. B 46 (1992) 6671--6687
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt,i
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(tol,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = dsqrt(sigmabb)
      t6 = t4/t3
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhob**2
      t12 = t2**2
      t14 = 1/t12/t11
      t17 = dexp(-0.1645530784602056D1*sigmabb*t14)
      t25 = sigmabb**2
      t26 = t11**2
      t36 = 1/rhob
      t37 = t36**(1.D0/3.D0)
      t40 = t36**(1.D0/6.D0)
      t43 = dsqrt(t36)
      t45 = t37**2
      t51 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t40
     #+0.3844746237447211D1*t37+0.1644733775567609D1*t43
     #+0.2405871291288192D0*t45))
      t52 = (1.D0+0.1274696188700087D0*t37)*t51
      t55 = 1/t2/t11
      t59 = dexp(0.2513789108010112D1*t52)
      t60 = t59-1.D0
      t61 = 1/t60
      t64 = 1/t12/t26
      t71 = t60**2
      t81 = dlog(1.D0+0.2697586091519874D1*(0.1007494971260294D0
     #*sigmabb*t55+0.2738174287780083D-1*t61*t25*t64)/(1.D0
     #+0.2717804421747983D0*t61*sigmabb*t55+0.7386460874872889D-1/t71
     #*t25*t64))
      zk(i) = -0.9305257363491D0*t3*(1.D0+t10+0.1645530784602056D-1*
     #(0.2743D0-0.1508D0*t17)*sigmabb*t14)/(1.D0+t10
     #+0.1083108625229223D-5*t25/t2/t26/rhob)+rhob*(-0.3109D-1*t52
     #+0.1236778371778789D-1*t81+0.7937005259840997D0*(0.1D-2*
     #(0.2568D1+0.1443307452126544D2*t37+0.2843543831490386D-2*t45)/
     #(1.D0+0.5411317332115466D1*t37+0.1816419932959077D0*t45
     #+0.1763993811759022D-1*t36)-0.1853571428571429D-2)*sigmabb*t55
     #*t17)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(tol,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = dsqrt(sigmaaa)
      t6 = t4/t3
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhoa**2
      t12 = t2**2
      t14 = 1/t12/t11
      t17 = dexp(-0.1645530784602056D1*sigmaaa*t14)
      t25 = sigmaaa**2
      t26 = t11**2
      t36 = 1/rhoa
      t37 = t36**(1.D0/3.D0)
      t40 = t36**(1.D0/6.D0)
      t43 = dsqrt(t36)
      t45 = t37**2
      t51 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t40
     #+0.3844746237447211D1*t37+0.1644733775567609D1*t43
     #+0.2405871291288192D0*t45))
      t52 = (1.D0+0.1274696188700087D0*t37)*t51
      t55 = 1/t2/t11
      t59 = dexp(0.2513789108010112D1*t52)
      t60 = t59-1.D0
      t61 = 1/t60
      t64 = 1/t12/t26
      t71 = t60**2
      t81 = dlog(1.D0+0.2697586091519874D1*(0.1007494971260294D0
     #*sigmaaa*t55+0.2738174287780083D-1*t61*t25*t64)/(1.D0
     #+0.2717804421747983D0*t61*sigmaaa*t55+0.7386460874872889D-1/t71
     #*t25*t64))
      zk(i) = -0.9305257363491D0*t3*(1.D0+t10+0.1645530784602056D-1*
     #(0.2743D0-0.1508D0*t17)*sigmaaa*t14)/(1.D0+t10
     #+0.1083108625229223D-5*t25/t2/t26/rhoa)+rhoa*(-0.3109D-1*t52
     #+0.1236778371778789D-1*t81+0.7937005259840997D0*(0.1D-2*
     #(0.2568D1+0.1443307452126544D2*t37+0.2843543831490386D-2*t45)/
     #(1.D0+0.5411317332115466D1*t37+0.1816419932959077D0*t45
     #+0.1763993811759022D-1*t36)-0.1853571428571429D-2)*sigmaaa*t55
     #*t17)
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(tol,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(tol,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = dsqrt(sigmaaa)
      t8 = t6/t5
      t10 = dlog(0.1000005877780776D1*t8+dsqrt(1+0.10000117555961D1*t8
     #**2))
      t12 = 0.2520026100493014D-1*t8*t10
      t13 = rhoa**2
      t14 = t4**2
      t16 = 1/t14/t13
      t19 = dexp(-0.1645530784602056D1*sigmaaa*t16)
      t27 = sigmaaa**2
      t28 = t13**2
      t38 = rhob**(1.D0/3.D0)
      t39 = t38*rhob
      t40 = dsqrt(sigmabb)
      t42 = t40/t39
      t44 = dlog(0.1000005877780776D1*t42+dsqrt(1+0.10000117555961D1
     #*t42**2))
      t46 = 0.2520026100493014D-1*t42*t44
      t47 = rhob**2
      t48 = t38**2
      t50 = 1/t48/t47
      t53 = dexp(-0.1645530784602056D1*sigmabb*t50)
      t61 = sigmabb**2
      t62 = t47**2
      t72 = 1/rho
      t73 = t72**(1.D0/3.D0)
      t76 = t72**(1.D0/6.D0)
      t79 = dsqrt(t72)
      t81 = t73**2
      t87 = dlog(1.D0+0.160818243221511D2/(0.598255043577108D1*t76
     #+0.2225569421150687D1*t73+0.8004286349993634D0*t79
     #+0.1897004325747559D0*t81))
      t89 = 0.62182D-1*(1.D0+0.1325688999052018D0*t73)*t87
      t100 = dlog(1.D0+0.2960857464321668D2/(0.8157414703487641D1*t76
     #+0.2247591863577616D1*t73+0.4300972471276643D0*t79
     #+0.1911512595127338D0*t81))
      t103 = rhoa-1.D0*rhob
      t104 = t103*t72
      t105 = 1.D0+t104
      t106 = t105**(1.D0/3.D0)
      t109 = 1.D0-1.D0*t104
      t110 = t109**(1.D0/3.D0)
      t112 = t106*t105+t110*t109-2.D0
      t113 = t103**2
      t114 = t113**2
      t115 = rho**2
      t116 = t115**2
      t118 = t114/t116
      t123 = 0.3799575D-1*(1.D0+0.6901399211255825D-1*t73)*t100*t112*
     #(1.D0-1.D0*t118)
      t134 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t76
     #+0.3844746237447211D1*t73+0.1644733775567609D1*t79
     #+0.2405871291288192D0*t81))
      t140 = 0.1923661050931536D1*(-0.3109D-1*(1.D0
     #+0.1274696188700087D0*t73)*t134+t89)*t112*t118
      t141 = t106**2
      t143 = t110**2
      t145 = 0.5D0*t141+0.5D0*t143
      t146 = t145**2
      t147 = t146*t145
      t148 = 1/t146
      t150 = rho**(1.D0/3.D0)
      t152 = 1/t150/t115
      t159 = dexp(-0.4042761511756372D2*(-t89+t123+t140)/t147)
      t160 = t159-1.D0
      t161 = 1/t160
      t162 = sigma**2
      t164 = t146**2
      t166 = t150**2
      t169 = 1/t164/t166/t116
      t177 = t160**2
      t187 = dlog(1.D0+0.2697586091519874D1*(0.634682060977037D-1
     #*sigma*t148*t152+0.1086645186223595D-1*t161*t162*t169)/(1.D0
     #+0.1712109500228824D0*t161*sigma*t148*t152+0.2931318940773793D-1
     #/t177*t162*t169))
      t208 = dexp(-0.261211729852336D1*t146/t166/t115*sigma)
      zk(i) = -0.9305257363491D0*t5*(1.D0+t12+0.1645530784602056D-1*
     #(0.2743D0-0.1508D0*t19)*sigmaaa*t16)/(1.D0+t12
     #+0.1083108625229223D-5*t27/t4/t28/rhoa)-0.9305257363491D0*t39*
     #(1.D0+t46+0.1645530784602056D-1*(0.2743D0-0.1508D0*t53)*sigmabb
     #*t50)/(1.D0+t46+0.1083108625229223D-5*t61/t38/t62/rhob)+rho*(
     #-t89+t123+t140+0.2473556743557577D-1*t147*t187+0.1D1*(0.1D-2*
     #(0.2568D1+0.1443307452126544D2*t73+0.2843543831490386D-2*t81)/
     #(1.D0+0.5411317332115466D1*t73+0.1816419932959077D0*t81
     #+0.1763993811759022D-1*t72)-0.1853571428571429D-2)*t145*sigma
     #*t152*t208)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(tol,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = dsqrt(sigmabb)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhob**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigmabb*t14
      t17 = dexp(-0.1645530784602056D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigmabb
      t23 = 1.D0+t10+0.1645530784602056D-1*t20*t14
      t24 = t3*t23
      t25 = sigmabb**2
      t26 = t11**2
      t27 = t26*rhob
      t29 = 1/t2/t27
      t32 = 1.D0+t10+0.1083108625229223D-5*t25*t29
      t33 = 1/t32
      t36 = 1/rhob
      t37 = t36**(1.D0/3.D0)
      t39 = 1.D0+0.1274696188700087D0*t37
      t40 = t36**(1.D0/6.D0)
      t43 = dsqrt(t36)
      t45 = t37**2
      t47 = 0.1112037486309468D2*t40+0.3844746237447211D1*t37
     #+0.1644733775567609D1*t43+0.2405871291288192D0*t45
      t50 = 1.D0+0.321646831778707D2/t47
      t51 = dlog(t50)
      t52 = t39*t51
      t53 = 0.3109D-1*t52
      t55 = 1/t2/t11
      t59 = dexp(0.2513789108010112D1*t52)
      t60 = t59-1.D0
      t61 = 1/t60
      t62 = t61*t25
      t64 = 1/t12/t26
      t67 = 0.1007494971260294D0*sigmabb*t55+0.2738174287780083D-1*t62
     #*t64
      t68 = t61*sigmabb
      t71 = t60**2
      t72 = 1/t71
      t73 = t72*t25
      t76 = 1.D0+0.2717804421747983D0*t68*t55+0.7386460874872889D-1
     #*t73*t64
      t77 = 1/t76
      t80 = 1.D0+0.2697586091519874D1*t67*t77
      t81 = dlog(t80)
      t82 = 0.1236778371778789D-1*t81
      t85 = 0.2568D1+0.1443307452126544D2*t37+0.2843543831490386D-2*t45
      t89 = 1.D0+0.5411317332115466D1*t37+0.1816419932959077D0*t45
     #+0.1763993811759022D-1*t36
      t90 = 1/t89
      t93 = 0.1D-2*t85*t90-0.1853571428571429D-2
      t94 = t93*sigmabb
      t95 = t55*t17
      t97 = 0.7937005259840997D0*t94*t95
      zk(i) = -0.9305257363491D0*t24*t33+rhob*(-t53+t82+t97)
      vrhoa(i) = 0.D0
      t105 = 0.3360034800657352D-1*t4*t55*t8
      t106 = t11*rhob
      t108 = 1/t12/t106
      t112 = dsqrt(1.D0+0.10000117555961D1*t15)
      t113 = 1/t112
      t115 = 0.3360054550205309D-1*sigmabb*t108*t113
      t116 = t26*t11
      t119 = t25/t2/t116
      t128 = t32**2
      t129 = 1/t128
      t136 = 1/t11
      t137 = 1/t45*t136
      t139 = 0.1321010150222857D-2*t137*t51
      t140 = t47**2
      t143 = t40**2
      t144 = t143**2
      t154 = 1/t37*t136
      t160 = 0.1D1*t39/t140*(-0.1853395810515781D1/t144/t40*t136
     #-0.128158207914907D1*t137-0.8223668877838045D0/t43*t136
     #-0.1603914194192128D0*t154)/t50
      t162 = 1/t2/t106
      t165 = t139+t160
      t167 = t64*t165*t59
      t171 = 1/t12/t27
      t177 = t76**2
      t179 = t67/t177
      t180 = t72*sigmabb
      t198 = 1/t80
      t206 = t89**2
      s1 = -0.12407009817988D1*t2*t23*t33-0.9305257363491D0*t3*(-t105
     #-t115-0.1088885204563779D-1*t119*t17-0.4388082092272149D-1*t20
     #*t108)*t33+0.9305257363491D0*t24*t129*(-t105-t115
     #-0.5776579334555855D-5*t119)
      vrhob(i) = s1-t53+t82+t97+rhob*(t139+t160+0.1236778371778789D-1*
     #(0.2697586091519874D1*(-0.2350821599607352D0*sigmabb*t162
     #+0.2213957124623647D1*t73*t167-0.1277814667630705D0*t62*t171)
     #*t77-0.2697586091519874D1*t179*(0.2197487022544806D2*t180*t55
     #*t165*t59-0.6341543650745294D0*t68*t162+0.1194467989321217D2/t71
     #/t60*t25*t167-0.3447015074940681D0*t73*t171))*t198
     #+0.7937005259840997D0*(0.1D-2*(-0.4811024840421814D1*t137
     #-0.1895695887660258D-2*t154)*t90-0.1D-2*t85/t206*(
     #-0.1803772444038489D1*t137-0.1210946621972718D0*t154
     #-0.1763993811759022D-1*t136))*sigmabb*t95-0.1851967893962899D1
     #*t94*t162*t17+0.3482823064697813D1*t93*t25/t116*t17)
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t232 = 0.1260013050246507D-1/t4*t5*t8
      t234 = 0.1260020456326991D-1*t14*t113
      vsigmabb(i) = -0.9305257363491D0*t3*(t232+t234
     #+0.408331951711417D-2*t29*t17*sigmabb+0.1645530784602056D-1*t19
     #*t14)*t33+0.9305257363491D0*t24*t129*(t232+t234
     #+0.2166217250458446D-5*sigmabb*t29)+rhob*(0.1236778371778789D-1*
     #(0.2697586091519874D1*(0.1007494971260294D0*t55
     #+0.5476348575560166D-1*t68*t64)*t77-0.2697586091519874D1*t179*
     #(0.2717804421747983D0*t61*t55+0.1477292174974578D0*t180*t64))
     #*t198+0.7937005259840997D0*t93*t55*t17-0.130605864926168D1*t94
     #/t27*t17)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(tol,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = dsqrt(sigmaaa)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhoa**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigmaaa*t14
      t17 = dexp(-0.1645530784602056D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigmaaa
      t23 = 1.D0+t10+0.1645530784602056D-1*t20*t14
      t24 = t3*t23
      t25 = sigmaaa**2
      t26 = t11**2
      t27 = t26*rhoa
      t29 = 1/t2/t27
      t32 = 1.D0+t10+0.1083108625229223D-5*t25*t29
      t33 = 1/t32
      t36 = 1/rhoa
      t37 = t36**(1.D0/3.D0)
      t39 = 1.D0+0.1274696188700087D0*t37
      t40 = t36**(1.D0/6.D0)
      t43 = dsqrt(t36)
      t45 = t37**2
      t47 = 0.1112037486309468D2*t40+0.3844746237447211D1*t37
     #+0.1644733775567609D1*t43+0.2405871291288192D0*t45
      t50 = 1.D0+0.321646831778707D2/t47
      t51 = dlog(t50)
      t52 = t39*t51
      t53 = 0.3109D-1*t52
      t55 = 1/t2/t11
      t59 = dexp(0.2513789108010112D1*t52)
      t60 = t59-1.D0
      t61 = 1/t60
      t62 = t61*t25
      t64 = 1/t12/t26
      t67 = 0.1007494971260294D0*sigmaaa*t55+0.2738174287780083D-1*t62
     #*t64
      t68 = t61*sigmaaa
      t71 = t60**2
      t72 = 1/t71
      t73 = t72*t25
      t76 = 1.D0+0.2717804421747983D0*t68*t55+0.7386460874872889D-1
     #*t73*t64
      t77 = 1/t76
      t80 = 1.D0+0.2697586091519874D1*t67*t77
      t81 = dlog(t80)
      t82 = 0.1236778371778789D-1*t81
      t85 = 0.2568D1+0.1443307452126544D2*t37+0.2843543831490386D-2*t45
      t89 = 1.D0+0.5411317332115466D1*t37+0.1816419932959077D0*t45
     #+0.1763993811759022D-1*t36
      t90 = 1/t89
      t93 = 0.1D-2*t85*t90-0.1853571428571429D-2
      t94 = t93*sigmaaa
      t95 = t55*t17
      t97 = 0.7937005259840997D0*t94*t95
      zk(i) = -0.9305257363491D0*t24*t33+rhoa*(-t53+t82+t97)
      t105 = 0.3360034800657352D-1*t4*t55*t8
      t106 = t11*rhoa
      t108 = 1/t12/t106
      t112 = dsqrt(1.D0+0.10000117555961D1*t15)
      t113 = 1/t112
      t115 = 0.3360054550205309D-1*sigmaaa*t108*t113
      t116 = t26*t11
      t119 = t25/t2/t116
      t128 = t32**2
      t129 = 1/t128
      t136 = 1/t11
      t137 = 1/t45*t136
      t139 = 0.1321010150222857D-2*t137*t51
      t140 = t47**2
      t143 = t40**2
      t144 = t143**2
      t154 = 1/t37*t136
      t160 = 0.1D1*t39/t140*(-0.1853395810515781D1/t144/t40*t136
     #-0.128158207914907D1*t137-0.8223668877838045D0/t43*t136
     #-0.1603914194192128D0*t154)/t50
      t162 = 1/t2/t106
      t165 = t139+t160
      t167 = t64*t165*t59
      t171 = 1/t12/t27
      t177 = t76**2
      t179 = t67/t177
      t180 = t72*sigmaaa
      t198 = 1/t80
      t206 = t89**2
      s1 = -0.12407009817988D1*t2*t23*t33-0.9305257363491D0*t3*(-t105
     #-t115-0.1088885204563779D-1*t119*t17-0.4388082092272149D-1*t20
     #*t108)*t33+0.9305257363491D0*t24*t129*(-t105-t115
     #-0.5776579334555855D-5*t119)
      vrhoa(i) = s1-t53+t82+t97+rhoa*(t139+t160+0.1236778371778789D-1*
     #(0.2697586091519874D1*(-0.2350821599607352D0*sigmaaa*t162
     #+0.2213957124623647D1*t73*t167-0.1277814667630705D0*t62*t171)
     #*t77-0.2697586091519874D1*t179*(0.2197487022544806D2*t180*t55
     #*t165*t59-0.6341543650745294D0*t68*t162+0.1194467989321217D2/t71
     #/t60*t25*t167-0.3447015074940681D0*t73*t171))*t198
     #+0.7937005259840997D0*(0.1D-2*(-0.4811024840421814D1*t137
     #-0.1895695887660258D-2*t154)*t90-0.1D-2*t85/t206*(
     #-0.1803772444038489D1*t137-0.1210946621972718D0*t154
     #-0.1763993811759022D-1*t136))*sigmaaa*t95-0.1851967893962899D1
     #*t94*t162*t17+0.3482823064697813D1*t93*t25/t116*t17)
      vrhob(i) = 0.D0
      t232 = 0.1260013050246507D-1/t4*t5*t8
      t234 = 0.1260020456326991D-1*t14*t113
      vsigmaaa(i) = -0.9305257363491D0*t3*(t232+t234
     #+0.408331951711417D-2*t29*t17*sigmaaa+0.1645530784602056D-1*t19
     #*t14)*t33+0.9305257363491D0*t24*t129*(t232+t234
     #+0.2166217250458446D-5*sigmaaa*t29)+rhoa*(0.1236778371778789D-1*
     #(0.2697586091519874D1*(0.1007494971260294D0*t55
     #+0.5476348575560166D-1*t68*t64)*t77-0.2697586091519874D1*t179*
     #(0.2717804421747983D0*t61*t55+0.1477292174974578D0*t180*t64))
     #*t198+0.7937005259840997D0*t93*t55*t17-0.130605864926168D1*t94
     #/t27*t17)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(tol,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(tol,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = dsqrt(sigmaaa)
      t7 = 1/t5
      t8 = t6*t7
      t10 = dlog(0.1000005877780776D1*t8+dsqrt(1+0.10000117555961D1*t8
     #**2))
      t12 = 0.2520026100493014D-1*t8*t10
      t13 = rhoa**2
      t14 = t4**2
      t16 = 1/t14/t13
      t17 = sigmaaa*t16
      t19 = dexp(-0.1645530784602056D1*t17)
      t21 = 0.2743D0-0.1508D0*t19
      t22 = t21*sigmaaa
      t25 = 1.D0+t12+0.1645530784602056D-1*t22*t16
      t26 = t5*t25
      t27 = sigmaaa**2
      t28 = t13**2
      t31 = 1/t4/t28/rhoa
      t34 = 1.D0+t12+0.1083108625229223D-5*t27*t31
      t35 = 1/t34
      t38 = rhob**(1.D0/3.D0)
      t39 = t38*rhob
      t40 = dsqrt(sigmabb)
      t41 = 1/t39
      t42 = t40*t41
      t44 = dlog(0.1000005877780776D1*t42+dsqrt(1+0.10000117555961D1
     #*t42**2))
      t46 = 0.2520026100493014D-1*t42*t44
      t47 = rhob**2
      t48 = t38**2
      t50 = 1/t48/t47
      t51 = sigmabb*t50
      t53 = dexp(-0.1645530784602056D1*t51)
      t55 = 0.2743D0-0.1508D0*t53
      t56 = t55*sigmabb
      t59 = 1.D0+t46+0.1645530784602056D-1*t56*t50
      t60 = t39*t59
      t61 = sigmabb**2
      t62 = t47**2
      t65 = 1/t38/t62/rhob
      t68 = 1.D0+t46+0.1083108625229223D-5*t61*t65
      t69 = 1/t68
      t72 = 1/rho
      t73 = t72**(1.D0/3.D0)
      t75 = 1.D0+0.1325688999052018D0*t73
      t76 = t72**(1.D0/6.D0)
      t79 = dsqrt(t72)
      t81 = t73**2
      t83 = 0.598255043577108D1*t76+0.2225569421150687D1*t73
     #+0.8004286349993634D0*t79+0.1897004325747559D0*t81
      t86 = 1.D0+0.160818243221511D2/t83
      t87 = dlog(t86)
      t89 = 0.62182D-1*t75*t87
      t91 = 1.D0+0.6901399211255825D-1*t73
      t96 = 0.8157414703487641D1*t76+0.2247591863577616D1*t73
     #+0.4300972471276643D0*t79+0.1911512595127338D0*t81
      t99 = 1.D0+0.2960857464321668D2/t96
      t100 = dlog(t99)
      t101 = t91*t100
      t103 = rhoa-1.D0*rhob
      t104 = t103*t72
      t105 = 1.D0+t104
      t106 = t105**(1.D0/3.D0)
      t109 = 1.D0-1.D0*t104
      t110 = t109**(1.D0/3.D0)
      t112 = t106*t105+t110*t109-2.D0
      t113 = t103**2
      t114 = t113**2
      t115 = rho**2
      t116 = t115**2
      t117 = 1/t116
      t118 = t114*t117
      t120 = 1.D0-1.D0*t118
      t123 = 0.3799575D-1*t101*t112*t120
      t125 = 1.D0+0.1274696188700087D0*t73
      t130 = 0.1112037486309468D2*t76+0.3844746237447211D1*t73
     #+0.1644733775567609D1*t79+0.2405871291288192D0*t81
      t133 = 1.D0+0.321646831778707D2/t130
      t134 = dlog(t133)
      t137 = -0.3109D-1*t125*t134+t89
      t138 = t137*t112
      t140 = 0.1923661050931536D1*t138*t118
      t141 = t106**2
      t143 = t110**2
      t145 = 0.5D0*t141+0.5D0*t143
      t146 = t145**2
      t147 = t146*t145
      t148 = 1/t146
      t149 = sigma*t148
      t150 = rho**(1.D0/3.D0)
      t152 = 1/t150/t115
      t155 = -t89+t123+t140
      t156 = 1/t147
      t159 = dexp(-0.4042761511756372D2*t155*t156)
      t160 = t159-1.D0
      t161 = 1/t160
      t162 = sigma**2
      t163 = t161*t162
      t164 = t146**2
      t165 = 1/t164
      t166 = t150**2
      t168 = 1/t166/t116
      t169 = t165*t168
      t172 = 0.634682060977037D-1*t149*t152+0.1086645186223595D-1*t163
     #*t169
      t173 = t161*sigma
      t174 = t148*t152
      t177 = t160**2
      t178 = 1/t177
      t179 = t178*t162
      t182 = 1.D0+0.1712109500228824D0*t173*t174+0.2931318940773793D-1
     #*t179*t169
      t183 = 1/t182
      t186 = 1.D0+0.2697586091519874D1*t172*t183
      t187 = dlog(t186)
      t189 = 0.2473556743557577D-1*t147*t187
      t192 = 0.2568D1+0.1443307452126544D2*t73+0.2843543831490386D-2*t81
      t196 = 1.D0+0.5411317332115466D1*t73+0.1816419932959077D0*t81
     #+0.1763993811759022D-1*t72
      t197 = 1/t196
      t200 = 0.1D-2*t192*t197-0.1853571428571429D-2
      t201 = t200*t145
      t204 = 1/t166/t115
      t208 = dexp(-0.261211729852336D1*t146*t204*sigma)
      t209 = sigma*t152*t208
      t211 = 0.1D1*t201*t209
      zk(i) = -0.9305257363491D0*t26*t35-0.9305257363491D0*t60*t69+rho
     #*(-t89+t123+t140+t189+t211)
      t221 = 0.3360034800657352D-1*t6/t4/t13*t10
      t224 = 1/t14/t13/rhoa
      t228 = dsqrt(1.D0+0.10000117555961D1*t17)
      t229 = 1/t228
      t231 = 0.3360054550205309D-1*sigmaaa*t224*t229
      t235 = t27/t4/t28/t13
      t244 = t34**2
      t245 = 1/t244
      t255 = 0.1333333333333333D1*t106*t72-0.1333333333333333D1*t110*t72
      t258 = 0.3799575D-1*t101*t255*t120
      t259 = t113*t103
      t262 = t101*t112*t259*t117
      t263 = 0.151983D0*t262
      t266 = 0.1923661050931536D1*t137*t255*t118
      t268 = t138*t259*t117
      t269 = 0.7694644203726145D1*t268
      t270 = t146*t187
      t271 = 1/t106
      t274 = 1/t110
      t277 = 0.3333333333333333D0*t271*t72-0.3333333333333333D0*t274*t72
      t280 = sigma*t156
      t284 = t179*t165
      t288 = t155*t165
      t291 = -0.4042761511756372D2*(t258-t263+t266+t269)*t156
     #+0.1212828453526912D3*t288*t277
      t293 = t168*t291*t159
      t298 = 1/t164/t145*t168
      t299 = t298*t277
      t305 = t182**2
      t307 = t172/t305
      t308 = t178*sigma
      t309 = t308*t148
      t314 = t156*t152
      t321 = 1/t177/t160*t162*t165
      t331 = 1/t186
      t338 = t200*t146*t162
      t339 = t116*rho
      t340 = 1/t339
      t348 = 1/t115
      t349 = 1/t81*t348
      t350 = t349*t87
      t351 = 0.2747799777968419D-2*t350
      t352 = t83**2
      t355 = t76**2
      t356 = t355**2
      t359 = 1/t356/t76*t348
      t363 = 1/t79*t348
      t366 = 1/t73*t348
      t371 = t75/t352*(-0.99709173929518D0*t359-0.7418564737168958D0
     #*t349-0.4002143174996817D0*t363-0.1264669550498372D0*t366)/t86
      t372 = 0.1D1*t371
      t376 = 0.8740794636035784D-3*t349*t100*t112*t120
      t377 = t96**2
      t390 = 0.1125D1*t91/t377*(-0.135956911724794D1*t359
     #-0.7491972878592054D0*t349-0.2150486235638321D0*t363
     #-0.1274341730084892D0*t366)/t99*t112*t120
      t397 = -0.1333333333333333D1*t106*t103*t348+0.1333333333333333D1
     #*t110*t103*t348
      t400 = 0.3799575D-1*t101*t397*t120
      t404 = 0.151983D0*t101*t112*t114*t340
      t407 = t130**2
      t424 = 0.1923661050931536D1*(0.1321010150222857D-2*t349*t134
     #+0.1D1*t125/t407*(-0.1853395810515781D1*t359-0.128158207914907D1
     #*t349-0.8223668877838045D0*t363-0.1603914194192128D0*t366)/t133
     #-0.2747799777968419D-2*t350-0.1D1*t371)*t112*t118
      t427 = 0.1923661050931536D1*t137*t397*t118
      t430 = 0.7694644203726145D1*t138*t114*t340
      t437 = -0.3333333333333333D0*t271*t103*t348+0.3333333333333333D0
     #*t274*t103*t348
      t443 = t115*rho
      t445 = 1/t150/t443
      t453 = -0.4042761511756372D2*(t351+t372-t376-t390+t400+t404+t424
     #+t427-t430)*t156+0.1212828453526912D3*t288*t437
      t455 = t168*t453*t159
      t458 = t298*t437
      t463 = t165/t166/t339
      t497 = t196**2
      s1 = t351+t372-t376-t390+t400+t404+t424
      s2 = s1+t427-t430+0.7420670230672731D-1*t270*t437
      s3 = s2+0.2473556743557577D-1*t147*(0.2697586091519874D1*(
     #-0.1269364121954074D0*t280*t152*t437-0.148092480894642D0*t149
     #*t445-0.1086645186223595D-1*t284*t455-0.4346580744894379D-1*t163
     #*t458-0.5071010869043442D-1*t163*t463)*t183-0.2697586091519874D1
     #*t307*(-0.1712109500228824D0*t309*t152*t453*t159
     #-0.3424219000457648D0*t173*t314*t437-0.3994922167200589D0*t173
     #*t148*t445-0.5862637881547585D-1*t321*t455-0.1172527576309517D0
     #*t179*t458-0.136794883902777D0*t179*t463))*t331+0.1D1*(0.1D-2*(
     #-0.4811024840421814D1*t349-0.1895695887660258D-2*t366)*t197
     #-0.1D-2*t192/t497*(-0.1803772444038489D1*t349
     #-0.1210946621972718D0*t366-0.1763993811759022D-1*t348))*t145*t209
      t532 = s3+0.1D1*t200*t437*t209-0.2333333333333333D1*t201*sigma
     #*t445*t208+0.1D1*t201*sigma*t152*(-0.522423459704672D1*t145*t204
     #*sigma*t437+0.6965646129395627D1*t146/t166/t443*sigma)*t208
      t533 = rho*t532
      vrhoa(i) = -0.12407009817988D1*t4*t25*t35-0.9305257363491D0*t5*(
     #-t221-t231-0.1088885204563779D-1*t235*t19-0.4388082092272149D-1
     #*t22*t224)*t35+0.9305257363491D0*t26*t245*(-t221-t231
     #-0.5776579334555855D-5*t235)+rho*(t258-t263+t266+t269
     #+0.7420670230672731D-1*t270*t277+0.2473556743557577D-1*t147*
     #(0.2697586091519874D1*(-0.1269364121954074D0*t280*t152*t277
     #-0.1086645186223595D-1*t284*t293-0.4346580744894379D-1*t163*t299
     #)*t183-0.2697586091519874D1*t307*(-0.1712109500228824D0*t309
     #*t152*t291*t159-0.3424219000457648D0*t173*t314*t277
     #-0.5862637881547585D-1*t321*t293-0.1172527576309517D0*t179*t299)
     #)*t331+0.1D1*t200*t277*t209-0.522423459704672D1*t338*t340*t277
     #*t208)-t89+t123+t140+t189+t211+t533
      t541 = 0.3360034800657352D-1*t40/t38/t47*t44
      t544 = 1/t48/t47/rhob
      t548 = dsqrt(1.D0+0.10000117555961D1*t51)
      t549 = 1/t548
      t551 = 0.3360054550205309D-1*sigmabb*t544*t549
      t555 = t61/t38/t62/t47
      t564 = t68**2
      t565 = 1/t564
      t571 = -t255
      t574 = 0.3799575D-1*t101*t571*t120
      t575 = 0.151983D0*t262
      t578 = 0.1923661050931536D1*t137*t571*t118
      t579 = 0.7694644203726145D1*t268
      t580 = -t277
      t591 = -0.4042761511756372D2*(t574+t575+t578-t579)*t156
     #+0.1212828453526912D3*t288*t580
      t593 = t168*t591*t159
      t596 = t298*t580
      vrhob(i) = -0.12407009817988D1*t38*t59*t69-0.9305257363491D0*t39
     #*(-t541-t551-0.1088885204563779D-1*t555*t53
     #-0.4388082092272149D-1*t56*t544)*t69+0.9305257363491D0*t60*t565*
     #(-t541-t551-0.5776579334555855D-5*t555)+rho*(t574+t575+t578-t579
     #+0.7420670230672731D-1*t270*t580+0.2473556743557577D-1*t147*
     #(0.2697586091519874D1*(-0.1269364121954074D0*t280*t152*t580
     #-0.1086645186223595D-1*t284*t593-0.4346580744894379D-1*t163*t596
     #)*t183-0.2697586091519874D1*t307*(-0.1712109500228824D0*t309
     #*t152*t591*t159-0.3424219000457648D0*t173*t314*t580
     #-0.5862637881547585D-1*t321*t593-0.1172527576309517D0*t179*t596)
     #)*t331+0.1D1*t200*t580*t209-0.522423459704672D1*t338*t340*t580
     #*t208)-t89+t123+t140+t189+t211+t533
      t632 = 0.1260013050246507D-1/t6*t7*t10
      t634 = 0.1260020456326991D-1*t16*t229
      t677 = rho*(0.2473556743557577D-1*t147*(0.2697586091519874D1*
     #(0.634682060977037D-1*t174+0.2173290372447189D-1*t173*t169)*t183
     #-0.2697586091519874D1*t307*(0.1712109500228824D0*t161*t148*t152
     #+0.5862637881547585D-1*t308*t169))*t331+0.1D1*t201*t152*t208
     #-0.261211729852336D1*t200*t147*sigma*t340*t208)
      vsigmaaa(i) = -0.9305257363491D0*t5*(t632+t634
     #+0.408331951711417D-2*t31*t19*sigmaaa+0.1645530784602056D-1*t21
     #*t16)*t35+0.9305257363491D0*t26*t245*(t632+t634
     #+0.2166217250458446D-5*sigmaaa*t31)+t677
      vsigmaab(i) = 2.D0*t677
      t681 = 0.1260013050246507D-1/t40*t41*t44
      t683 = 0.1260020456326991D-1*t50*t549
      vsigmabb(i) = -0.9305257363491D0*t39*(t681+t683
     #+0.408331951711417D-2*t65*t53*sigmabb+0.1645530784602056D-1*t55
     #*t50)*t69+0.9305257363491D0*t60*t565*(t681+t683
     #+0.2166217250458446D-5*sigmabb*t65)+t677
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(tol,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = dsqrt(sigmabb)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhob**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigmabb*t14
      t17 = dexp(-0.1645530784602056D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigmabb
      t23 = 1.D0+t10+0.1645530784602056D-1*t20*t14
      t24 = t3*t23
      t25 = sigmabb**2
      t26 = t11**2
      t27 = t26*rhob
      t29 = 1/t2/t27
      t32 = 1.D0+t10+0.1083108625229223D-5*t25*t29
      t33 = 1/t32
      t36 = 1/rhob
      t37 = t36**(1.D0/3.D0)
      t39 = 1.D0+0.1274696188700087D0*t37
      t40 = t36**(1.D0/6.D0)
      t43 = dsqrt(t36)
      t45 = t37**2
      t47 = 0.1112037486309468D2*t40+0.3844746237447211D1*t37
     #+0.1644733775567609D1*t43+0.2405871291288192D0*t45
      t50 = 1.D0+0.321646831778707D2/t47
      t51 = dlog(t50)
      t52 = t39*t51
      t53 = 0.3109D-1*t52
      t55 = 1/t2/t11
      t59 = dexp(0.2513789108010112D1*t52)
      t60 = t59-1.D0
      t61 = 1/t60
      t62 = t61*t25
      t64 = 1/t12/t26
      t67 = 0.1007494971260294D0*sigmabb*t55+0.2738174287780083D-1*t62
     #*t64
      t68 = t61*sigmabb
      t71 = t60**2
      t72 = 1/t71
      t73 = t72*t25
      t76 = 1.D0+0.2717804421747983D0*t68*t55+0.7386460874872889D-1
     #*t73*t64
      t77 = 1/t76
      t80 = 1.D0+0.2697586091519874D1*t67*t77
      t81 = dlog(t80)
      t82 = 0.1236778371778789D-1*t81
      t85 = 0.2568D1+0.1443307452126544D2*t37+0.2843543831490386D-2*t45
      t89 = 1.D0+0.5411317332115466D1*t37+0.1816419932959077D0*t45
     #+0.1763993811759022D-1*t36
      t90 = 1/t89
      t93 = 0.1D-2*t85*t90-0.1853571428571429D-2
      t94 = t93*sigmabb
      t95 = t55*t17
      t97 = 0.7937005259840997D0*t94*t95
      zk(i) = -0.9305257363491D0*t24*t33+rhob*(-t53+t82+t97)
      vrhoa(i) = 0.D0
      t100 = t2*t23
      t105 = 0.3360034800657352D-1*t4*t55*t8
      t106 = t11*rhob
      t108 = 1/t12/t106
      t111 = 1.D0+0.10000117555961D1*t15
      t112 = dsqrt(t111)
      t113 = 1/t112
      t115 = 0.3360054550205309D-1*sigmabb*t108*t113
      t116 = t26*t11
      t119 = t25/t2/t116
      t124 = -t105-t115-0.1088885204563779D-1*t119*t17
     #-0.4388082092272149D-1*t20*t108
      t125 = t3*t124
      t128 = t32**2
      t129 = 1/t128
      t131 = -t105-t115-0.5776579334555855D-5*t119
      t132 = t129*t131
      t135 = 1/t45
      t136 = 1/t11
      t137 = t135*t136
      t138 = t137*t51
      t139 = 0.1321010150222857D-2*t138
      t140 = t47**2
      t141 = 1/t140
      t142 = t39*t141
      t143 = t40**2
      t144 = t143**2
      t145 = t144*t40
      t146 = 1/t145
      t150 = 1/t43
      t153 = 1/t37
      t154 = t153*t136
      t156 = -0.1853395810515781D1*t146*t136-0.128158207914907D1*t137
     #-0.8223668877838045D0*t150*t136-0.1603914194192128D0*t154
      t157 = 1/t50
      t159 = t142*t156*t157
      t160 = 0.1D1*t159
      t162 = 1/t2/t106
      t165 = t139+t160
      t167 = t64*t165*t59
      t171 = 1/t12/t27
      t174 = -0.2350821599607352D0*sigmabb*t162+0.2213957124623647D1
     #*t73*t167-0.1277814667630705D0*t62*t171
      t177 = t76**2
      t178 = 1/t177
      t179 = t67*t178
      t180 = t72*sigmabb
      t188 = 1/t71/t60
      t189 = t188*t25
      t194 = 0.2197487022544806D2*t180*t55*t165*t59
     #-0.6341543650745294D0*t68*t162+0.1194467989321217D2*t189*t167
     #-0.3447015074940681D0*t73*t171
      t197 = 0.2697586091519874D1*t174*t77-0.2697586091519874D1*t179
     #*t194
      t198 = 1/t80
      t199 = t197*t198
      t203 = -0.4811024840421814D1*t137-0.1895695887660258D-2*t154
      t206 = t89**2
      t207 = 1/t206
      t208 = t85*t207
      t212 = -0.1803772444038489D1*t137-0.1210946621972718D0*t154
     #-0.1763993811759022D-1*t136
      t215 = 0.1D-2*t203*t90-0.1D-2*t208*t212
      t216 = t215*sigmabb
      t217 = t216*t95
      t219 = t162*t17
      t220 = t94*t219
      t222 = t93*t25
      t224 = 1/t116*t17
      t225 = t222*t224
      vrhob(i) = -0.12407009817988D1*t100*t33-0.9305257363491D0*t125
     #*t33+0.9305257363491D0*t24*t132-t53+t82+t97+rhob*(t139+t160
     #+0.1236778371778789D-1*t199+0.7937005259840997D0*t217
     #-0.1851967893962899D1*t220+0.3482823064697813D1*t225)
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t232 = 0.1260013050246507D-1/t4*t5*t8
      t234 = 0.1260020456326991D-1*t14*t113
      t235 = t29*t17
      t241 = t3*(t232+t234+0.408331951711417D-2*t235*sigmabb
     #+0.1645530784602056D-1*t19*t14)
      t246 = t232+t234+0.2166217250458446D-5*sigmabb*t29
      t247 = t129*t246
      t253 = 0.1007494971260294D0*t55+0.5476348575560166D-1*t68*t64
      t260 = 0.2717804421747983D0*t61*t55+0.1477292174974578D0*t180*t64
      t263 = 0.2697586091519874D1*t253*t77-0.2697586091519874D1*t179
     #*t260
      t269 = 1/t27
      vsigmabb(i) = -0.9305257363491D0*t241*t33+0.9305257363491D0*t24
     #*t247+rhob*(0.1236778371778789D-1*t263*t198+0.7937005259840997D0
     #*t93*t55*t17-0.130605864926168D1*t94*t269*t17)
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t286 = 0.7840081201533821D-1*t4*t162*t8
      t289 = 0.1680027275102654D0*sigmabb*t64*t113
      t290 = t26*t106
      t293 = t25/t2/t290
      t295 = 1/t112/t111
      t297 = 0.4480125399532632D-1*t293*t295
      t300 = t25*sigmabb
      t301 = t26**2
      t316 = 1/t128/t32
      t317 = t131**2
      t334 = dlog(1.D0+0.160818243221511D2/(0.598255043577108D1*t40
     #+0.2225569421150687D1*t37+0.8004286349993634D0*t43
     #+0.1897004325747559D0*t45))
      t346 = t156**2
      t349 = 0.2D1*t39/t140/t47*t346*t157
      t350 = 1/t106
      t351 = t135*t350
      t353 = 0.1D-22*t351*t334
      t356 = 1/t26
      t357 = 1/t45/t36*t356
      t359 = 0.8806734334819047D-3*t357*t51
      t361 = 0.2642020300445714D-2*t351*t51
      t365 = 0.8497974591333914D-1*t137*t141*t156*t157
      t382 = 1/t37/t36*t356
      t384 = t153*t350
      t389 = 0.1D1*t142*(-0.1544496508763151D1/t145/t36*t356
     #+0.3706791621031562D1*t146*t350-0.854388052766047D0*t357
     #+0.2563164158298141D1*t351-0.4111834438919023D0/t43/t36*t356
     #+0.1644733775567609D1*t150*t350-0.5346380647307093D-1*t382
     #+0.3207828388384256D0*t384)*t157
      t390 = t140**2
      t393 = t50**2
      t397 = 0.321646831778707D2*t39/t390*t346/t393
      t411 = t212**2
      t427 = 1/t2/t26
      t430 = t165**2
      t431 = t64*t430
      t433 = dexp(0.5027578216020225D1*t52)
      t434 = t431*t433
      t438 = t171*t165*t59
      t441 = -t353+t359-t361-t365-t349+t389+t397
      t443 = t64*t441*t59
      t446 = t431*t59
      t450 = 1/t12/t116
      t461 = t67/t177/t76
      t462 = t194**2
      t466 = t55*t430
      t483 = t71**2
      t506 = t197**2
      t507 = t80**2
      t508 = 1/t507
      s1 = -t353+t359-t361-t365-t349+t389+t397
      s3 = s1+0.7937005259840997D0*(0.1D-2*(-0.3207349893614542D1*t357
     #+0.9622049680843627D1*t351-0.6318986292200858D-3*t382
     #+0.3791391775320515D-2*t384)*t90-0.2D-2*t203*t207*t212+0.2D-2
     #*t85/t206/t89*t411-0.1D-2*t208*(-0.1202514962692326D1*t357
     #+0.3607544888076978D1*t351-0.4036488739909061D-1*t382
     #+0.2421893243945437D0*t384+0.3527987623518044D-1*t350))*sigmabb
     #*t95
      s2 = s3+0.1236778371778789D-1*(0.2697586091519874D1*
     #(0.7836071998691172D0*sigmabb*t427+0.3580200260842915D3*t189
     #*t434-0.2066359982982071D2*t73*t438+0.2213957124623647D1*t73
     #*t443-0.1790100130421458D3*t73*t446+0.7240949783240664D0*t62
     #*t450)*t77-0.5395172183039748D1*t174*t178*t194
     #+0.5395172183039748D1*t461*t462-0.2697586091519874D1*t179*
     #(0.35535663829313D4*t188*sigmabb*t466*t433-0.1025493943854243D3
     #*t180*t162*t165*t59+0.2197487022544806D2*t180*t55*t441*t59
     #-0.177678319146565D4*t180*t466*t59+0.2113847883581765D1*t68*t427
     #+0.2897369528551702D4/t483*t25*t434-0.1114836790033136D3*t189
     #*t438+0.1194467989321217D2*t189*t443-0.9657898428505673D3*t189
     #*t446+0.1953308542466386D1*t73*t450))*t198-0.2902352553914844D2
     #*t222/t290*t17
      t526 = s2-0.1236778371778789D-1*t506*t508-0.3703935787925799D1
     #*t216*t219+0.6965646129395627D1*t215*t25*t224
     #+0.6173226313209665D1*t94*t427*t17+0.1528291352075288D2*t93*t300
     #/t12/t301/rhob*t17
      v2rhob2(i) = -0.4135669939329333D0/t12*t23*t33
     #-0.24814019635976D1*t2*t124*t33+0.24814019635976D1*t100*t132
     #-0.9305257363491D0*t3*(t286+t289-t297+0.9799966841074009D-1*t293
     #*t17-0.4778117666686413D-1*t300/t301/t11*t17
     #+0.1608963433833121D0*t20*t64)*t33+0.18610514726982D1*t125*t132
     #-0.18610514726982D1*t24*t316*t317+0.9305257363491D0*t24*t129*
     #(t286+t289-t297+0.3658500245218708D-4*t293)-0.1D-22*t137*t334
     #+0.2642020300445714D-2*t138+0.2D1*t159+0.2473556743557577D-1
     #*t199+0.1587401051968199D1*t217-0.3703935787925799D1*t220
     #+0.6965646129395627D1*t225+rhob*t526
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t532 = 0.6300065251232535D-2/t4/sigmabb*t5*t8
      t536 = 0.6300102281634954D-2/sigmabb*t14*t113
      t538 = 0.6300176343092764D-2*t29*t295
      t550 = t246**2
      t565 = t260**2
      t574 = t263**2
      v2sigmabb2(i) = -0.9305257363491D0*t3*(-t532+t536-t538
     #-0.6719227968777768D-2/t301*t17*sigmabb+0.8166639034228341D-2
     #*t235)*t33+0.18610514726982D1*t241*t247-0.18610514726982D1*t24
     #*t316*t550+0.9305257363491D0*t24*t129*(-t532+t536-t538
     #+0.2166217250458446D-5*t29)+rhob*(0.1236778371778789D-1*
     #(0.1477292174974578D0*t61*t64*t77-0.5395172183039748D1*t253*t178
     #*t260+0.5395172183039748D1*t461*t565-0.3985122824322565D0*t179
     #*t72*t64)*t198-0.1236778371778789D-1*t574*t508
     #-0.261211729852336D1*t93*t269*t17+0.2149159713855873D1*t94/t12
     #/t290*t17)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(tol,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = dsqrt(sigmaaa)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhoa**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigmaaa*t14
      t17 = dexp(-0.1645530784602056D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigmaaa
      t23 = 1.D0+t10+0.1645530784602056D-1*t20*t14
      t24 = t3*t23
      t25 = sigmaaa**2
      t26 = t11**2
      t27 = t26*rhoa
      t29 = 1/t2/t27
      t32 = 1.D0+t10+0.1083108625229223D-5*t25*t29
      t33 = 1/t32
      t36 = 1/rhoa
      t37 = t36**(1.D0/3.D0)
      t39 = 1.D0+0.1274696188700087D0*t37
      t40 = t36**(1.D0/6.D0)
      t43 = dsqrt(t36)
      t45 = t37**2
      t47 = 0.1112037486309468D2*t40+0.3844746237447211D1*t37
     #+0.1644733775567609D1*t43+0.2405871291288192D0*t45
      t50 = 1.D0+0.321646831778707D2/t47
      t51 = dlog(t50)
      t52 = t39*t51
      t53 = 0.3109D-1*t52
      t55 = 1/t2/t11
      t59 = dexp(0.2513789108010112D1*t52)
      t60 = t59-1.D0
      t61 = 1/t60
      t62 = t61*t25
      t64 = 1/t12/t26
      t67 = 0.1007494971260294D0*sigmaaa*t55+0.2738174287780083D-1*t62
     #*t64
      t68 = t61*sigmaaa
      t71 = t60**2
      t72 = 1/t71
      t73 = t72*t25
      t76 = 1.D0+0.2717804421747983D0*t68*t55+0.7386460874872889D-1
     #*t73*t64
      t77 = 1/t76
      t80 = 1.D0+0.2697586091519874D1*t67*t77
      t81 = dlog(t80)
      t82 = 0.1236778371778789D-1*t81
      t85 = 0.2568D1+0.1443307452126544D2*t37+0.2843543831490386D-2*t45
      t89 = 1.D0+0.5411317332115466D1*t37+0.1816419932959077D0*t45
     #+0.1763993811759022D-1*t36
      t90 = 1/t89
      t93 = 0.1D-2*t85*t90-0.1853571428571429D-2
      t94 = t93*sigmaaa
      t95 = t55*t17
      t97 = 0.7937005259840997D0*t94*t95
      zk(i) = -0.9305257363491D0*t24*t33+rhoa*(-t53+t82+t97)
      t100 = t2*t23
      t105 = 0.3360034800657352D-1*t4*t55*t8
      t106 = t11*rhoa
      t108 = 1/t12/t106
      t111 = 1.D0+0.10000117555961D1*t15
      t112 = dsqrt(t111)
      t113 = 1/t112
      t115 = 0.3360054550205309D-1*sigmaaa*t108*t113
      t116 = t26*t11
      t119 = t25/t2/t116
      t124 = -t105-t115-0.1088885204563779D-1*t119*t17
     #-0.4388082092272149D-1*t20*t108
      t125 = t3*t124
      t128 = t32**2
      t129 = 1/t128
      t131 = -t105-t115-0.5776579334555855D-5*t119
      t132 = t129*t131
      t135 = 1/t45
      t136 = 1/t11
      t137 = t135*t136
      t138 = t137*t51
      t139 = 0.1321010150222857D-2*t138
      t140 = t47**2
      t141 = 1/t140
      t142 = t39*t141
      t143 = t40**2
      t144 = t143**2
      t145 = t144*t40
      t146 = 1/t145
      t150 = 1/t43
      t153 = 1/t37
      t154 = t153*t136
      t156 = -0.1853395810515781D1*t146*t136-0.128158207914907D1*t137
     #-0.8223668877838045D0*t150*t136-0.1603914194192128D0*t154
      t157 = 1/t50
      t159 = t142*t156*t157
      t160 = 0.1D1*t159
      t162 = 1/t2/t106
      t165 = t139+t160
      t167 = t64*t165*t59
      t171 = 1/t12/t27
      t174 = -0.2350821599607352D0*sigmaaa*t162+0.2213957124623647D1
     #*t73*t167-0.1277814667630705D0*t62*t171
      t177 = t76**2
      t178 = 1/t177
      t179 = t67*t178
      t180 = t72*sigmaaa
      t188 = 1/t71/t60
      t189 = t188*t25
      t194 = 0.2197487022544806D2*t180*t55*t165*t59
     #-0.6341543650745294D0*t68*t162+0.1194467989321217D2*t189*t167
     #-0.3447015074940681D0*t73*t171
      t197 = 0.2697586091519874D1*t174*t77-0.2697586091519874D1*t179
     #*t194
      t198 = 1/t80
      t199 = t197*t198
      t203 = -0.4811024840421814D1*t137-0.1895695887660258D-2*t154
      t206 = t89**2
      t207 = 1/t206
      t208 = t85*t207
      t212 = -0.1803772444038489D1*t137-0.1210946621972718D0*t154
     #-0.1763993811759022D-1*t136
      t215 = 0.1D-2*t203*t90-0.1D-2*t208*t212
      t216 = t215*sigmaaa
      t217 = t216*t95
      t219 = t162*t17
      t220 = t94*t219
      t222 = t93*t25
      t224 = 1/t116*t17
      t225 = t222*t224
      vrhoa(i) = -0.12407009817988D1*t100*t33-0.9305257363491D0*t125
     #*t33+0.9305257363491D0*t24*t132-t53+t82+t97+rhoa*(t139+t160
     #+0.1236778371778789D-1*t199+0.7937005259840997D0*t217
     #-0.1851967893962899D1*t220+0.3482823064697813D1*t225)
      vrhob(i) = 0.D0
      t232 = 0.1260013050246507D-1/t4*t5*t8
      t234 = 0.1260020456326991D-1*t14*t113
      t235 = t29*t17
      t241 = t3*(t232+t234+0.408331951711417D-2*t235*sigmaaa
     #+0.1645530784602056D-1*t19*t14)
      t246 = t232+t234+0.2166217250458446D-5*sigmaaa*t29
      t247 = t129*t246
      t253 = 0.1007494971260294D0*t55+0.5476348575560166D-1*t68*t64
      t260 = 0.2717804421747983D0*t61*t55+0.1477292174974578D0*t180*t64
      t263 = 0.2697586091519874D1*t253*t77-0.2697586091519874D1*t179
     #*t260
      t269 = 1/t27
      vsigmaaa(i) = -0.9305257363491D0*t241*t33+0.9305257363491D0*t24
     #*t247+rhoa*(0.1236778371778789D-1*t263*t198+0.7937005259840997D0
     #*t93*t55*t17-0.130605864926168D1*t94*t269*t17)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t286 = 0.7840081201533821D-1*t4*t162*t8
      t289 = 0.1680027275102654D0*sigmaaa*t64*t113
      t290 = t26*t106
      t293 = t25/t2/t290
      t295 = 1/t112/t111
      t297 = 0.4480125399532632D-1*t293*t295
      t300 = t25*sigmaaa
      t301 = t26**2
      t316 = 1/t128/t32
      t317 = t131**2
      t334 = dlog(1.D0+0.160818243221511D2/(0.598255043577108D1*t40
     #+0.2225569421150687D1*t37+0.8004286349993634D0*t43
     #+0.1897004325747559D0*t45))
      t346 = t156**2
      t349 = 0.2D1*t39/t140/t47*t346*t157
      t350 = 1/t106
      t351 = t135*t350
      t353 = 0.1D-22*t351*t334
      t356 = 1/t26
      t357 = 1/t45/t36*t356
      t359 = 0.8806734334819047D-3*t357*t51
      t361 = 0.2642020300445714D-2*t351*t51
      t365 = 0.8497974591333914D-1*t137*t141*t156*t157
      t382 = 1/t37/t36*t356
      t384 = t153*t350
      t389 = 0.1D1*t142*(-0.1544496508763151D1/t145/t36*t356
     #+0.3706791621031562D1*t146*t350-0.854388052766047D0*t357
     #+0.2563164158298141D1*t351-0.4111834438919023D0/t43/t36*t356
     #+0.1644733775567609D1*t150*t350-0.5346380647307093D-1*t382
     #+0.3207828388384256D0*t384)*t157
      t390 = t140**2
      t393 = t50**2
      t397 = 0.321646831778707D2*t39/t390*t346/t393
      t411 = t212**2
      t427 = 1/t2/t26
      t430 = t165**2
      t431 = t64*t430
      t433 = dexp(0.5027578216020225D1*t52)
      t434 = t431*t433
      t438 = t171*t165*t59
      t441 = -t353+t359-t361-t365-t349+t389+t397
      t443 = t64*t441*t59
      t446 = t431*t59
      t450 = 1/t12/t116
      t461 = t67/t177/t76
      t462 = t194**2
      t466 = t55*t430
      t483 = t71**2
      t506 = t197**2
      t507 = t80**2
      t508 = 1/t507
      s1 = -t353+t359-t361-t365-t349+t389+t397
      s3 = s1+0.7937005259840997D0*(0.1D-2*(-0.3207349893614542D1*t357
     #+0.9622049680843627D1*t351-0.6318986292200858D-3*t382
     #+0.3791391775320515D-2*t384)*t90-0.2D-2*t203*t207*t212+0.2D-2
     #*t85/t206/t89*t411-0.1D-2*t208*(-0.1202514962692326D1*t357
     #+0.3607544888076978D1*t351-0.4036488739909061D-1*t382
     #+0.2421893243945437D0*t384+0.3527987623518044D-1*t350))*sigmaaa
     #*t95
      s2 = s3+0.1236778371778789D-1*(0.2697586091519874D1*
     #(0.7836071998691172D0*sigmaaa*t427+0.3580200260842915D3*t189
     #*t434-0.2066359982982071D2*t73*t438+0.2213957124623647D1*t73
     #*t443-0.1790100130421458D3*t73*t446+0.7240949783240664D0*t62
     #*t450)*t77-0.5395172183039748D1*t174*t178*t194
     #+0.5395172183039748D1*t461*t462-0.2697586091519874D1*t179*
     #(0.35535663829313D4*t188*sigmaaa*t466*t433-0.1025493943854243D3
     #*t180*t162*t165*t59+0.2197487022544806D2*t180*t55*t441*t59
     #-0.177678319146565D4*t180*t466*t59+0.2113847883581765D1*t68*t427
     #+0.2897369528551702D4/t483*t25*t434-0.1114836790033136D3*t189
     #*t438+0.1194467989321217D2*t189*t443-0.9657898428505673D3*t189
     #*t446+0.1953308542466386D1*t73*t450))*t198-0.2902352553914844D2
     #*t222/t290*t17
      t526 = s2-0.1236778371778789D-1*t506*t508-0.3703935787925799D1
     #*t216*t219+0.6965646129395627D1*t215*t25*t224
     #+0.6173226313209665D1*t94*t427*t17+0.1528291352075288D2*t93*t300
     #/t12/t301/rhoa*t17
      v2rhoa2(i) = -0.4135669939329333D0/t12*t23*t33
     #-0.24814019635976D1*t2*t124*t33+0.24814019635976D1*t100*t132
     #-0.9305257363491D0*t3*(t286+t289-t297+0.9799966841074009D-1*t293
     #*t17-0.4778117666686413D-1*t300/t301/t11*t17
     #+0.1608963433833121D0*t20*t64)*t33+0.18610514726982D1*t125*t132
     #-0.18610514726982D1*t24*t316*t317+0.9305257363491D0*t24*t129*
     #(t286+t289-t297+0.3658500245218708D-4*t293)-0.1D-22*t137*t334
     #+0.2642020300445714D-2*t138+0.2D1*t159+0.2473556743557577D-1
     #*t199+0.1587401051968199D1*t217-0.3703935787925799D1*t220
     #+0.6965646129395627D1*t225+rhoa*t526
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t532 = 0.6300065251232535D-2/t4/sigmaaa*t5*t8
      t536 = 0.6300102281634954D-2/sigmaaa*t14*t113
      t538 = 0.6300176343092764D-2*t29*t295
      t550 = t246**2
      t565 = t260**2
      t574 = t263**2
      v2sigmaaa2(i) = -0.9305257363491D0*t3*(-t532+t536-t538
     #-0.6719227968777768D-2/t301*t17*sigmaaa+0.8166639034228341D-2
     #*t235)*t33+0.18610514726982D1*t241*t247-0.18610514726982D1*t24
     #*t316*t550+0.9305257363491D0*t24*t129*(-t532+t536-t538
     #+0.2166217250458446D-5*t29)+rhoa*(0.1236778371778789D-1*
     #(0.1477292174974578D0*t61*t64*t77-0.5395172183039748D1*t253*t178
     #*t260+0.5395172183039748D1*t461*t565-0.3985122824322565D0*t179
     #*t72*t64)*t198-0.1236778371778789D-1*t574*t508
     #-0.261211729852336D1*t93*t269*t17+0.2149159713855873D1*t94/t12
     #/t290*t17)
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(tol,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(tol,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = dsqrt(sigmaaa)
      t7 = 1/t5
      t8 = t6*t7
      t10 = dlog(0.1000005877780776D1*t8+dsqrt(1+0.10000117555961D1*t8
     #**2))
      t12 = 0.2520026100493014D-1*t8*t10
      t13 = rhoa**2
      t14 = t4**2
      t16 = 1/t14/t13
      t17 = sigmaaa*t16
      t19 = dexp(-0.1645530784602056D1*t17)
      t21 = 0.2743D0-0.1508D0*t19
      t22 = t21*sigmaaa
      t25 = 1.D0+t12+0.1645530784602056D-1*t22*t16
      t26 = t5*t25
      t27 = sigmaaa**2
      t28 = t13**2
      t31 = 1/t4/t28/rhoa
      t34 = 1.D0+t12+0.1083108625229223D-5*t27*t31
      t35 = 1/t34
      t38 = rhob**(1.D0/3.D0)
      t39 = t38*rhob
      t40 = dsqrt(sigmabb)
      t41 = 1/t39
      t42 = t40*t41
      t44 = dlog(0.1000005877780776D1*t42+dsqrt(1+0.10000117555961D1
     #*t42**2))
      t46 = 0.2520026100493014D-1*t42*t44
      t47 = rhob**2
      t48 = t38**2
      t50 = 1/t48/t47
      t51 = sigmabb*t50
      t53 = dexp(-0.1645530784602056D1*t51)
      t55 = 0.2743D0-0.1508D0*t53
      t56 = t55*sigmabb
      t59 = 1.D0+t46+0.1645530784602056D-1*t56*t50
      t60 = t39*t59
      t61 = sigmabb**2
      t62 = t47**2
      t65 = 1/t38/t62/rhob
      t68 = 1.D0+t46+0.1083108625229223D-5*t61*t65
      t69 = 1/t68
      t72 = 1/rho
      t73 = t72**(1.D0/3.D0)
      t75 = 1.D0+0.1325688999052018D0*t73
      t76 = t72**(1.D0/6.D0)
      t79 = dsqrt(t72)
      t81 = t73**2
      t83 = 0.598255043577108D1*t76+0.2225569421150687D1*t73
     #+0.8004286349993634D0*t79+0.1897004325747559D0*t81
      t86 = 1.D0+0.160818243221511D2/t83
      t87 = dlog(t86)
      t89 = 0.62182D-1*t75*t87
      t91 = 1.D0+0.6901399211255825D-1*t73
      t96 = 0.8157414703487641D1*t76+0.2247591863577616D1*t73
     #+0.4300972471276643D0*t79+0.1911512595127338D0*t81
      t99 = 1.D0+0.2960857464321668D2/t96
      t100 = dlog(t99)
      t101 = t91*t100
      t103 = rhoa-1.D0*rhob
      t104 = t103*t72
      t105 = 1.D0+t104
      t106 = t105**(1.D0/3.D0)
      t107 = t106*t105
      t109 = 1.D0-1.D0*t104
      t110 = t109**(1.D0/3.D0)
      t111 = t110*t109
      t112 = t107+t111-2.D0
      t113 = t103**2
      t114 = t113**2
      t115 = rho**2
      t116 = t115**2
      t117 = 1/t116
      t118 = t114*t117
      t120 = 1.D0-1.D0*t118
      t121 = t112*t120
      t123 = 0.3799575D-1*t101*t121
      t125 = 1.D0+0.1274696188700087D0*t73
      t130 = 0.1112037486309468D2*t76+0.3844746237447211D1*t73
     #+0.1644733775567609D1*t79+0.2405871291288192D0*t81
      t133 = 1.D0+0.321646831778707D2/t130
      t134 = dlog(t133)
      t137 = -0.3109D-1*t125*t134+t89
      t138 = t137*t112
      t140 = 0.1923661050931536D1*t138*t118
      t141 = t106**2
      t143 = t110**2
      t145 = 0.5D0*t141+0.5D0*t143
      t146 = t145**2
      t147 = t146*t145
      t148 = 1/t146
      t149 = sigma*t148
      t150 = rho**(1.D0/3.D0)
      t152 = 1/t150/t115
      t155 = -t89+t123+t140
      t156 = 1/t147
      t157 = t155*t156
      t159 = dexp(-0.4042761511756372D2*t157)
      t160 = t159-1.D0
      t161 = 1/t160
      t162 = sigma**2
      t163 = t161*t162
      t164 = t146**2
      t165 = 1/t164
      t166 = t150**2
      t168 = 1/t166/t116
      t169 = t165*t168
      t172 = 0.634682060977037D-1*t149*t152+0.1086645186223595D-1*t163
     #*t169
      t173 = t161*sigma
      t174 = t148*t152
      t177 = t160**2
      t178 = 1/t177
      t179 = t178*t162
      t182 = 1.D0+0.1712109500228824D0*t173*t174+0.2931318940773793D-1
     #*t179*t169
      t183 = 1/t182
      t186 = 1.D0+0.2697586091519874D1*t172*t183
      t187 = dlog(t186)
      t189 = 0.2473556743557577D-1*t147*t187
      t192 = 0.2568D1+0.1443307452126544D2*t73+0.2843543831490386D-2*t81
      t196 = 1.D0+0.5411317332115466D1*t73+0.1816419932959077D0*t81
     #+0.1763993811759022D-1*t72
      t197 = 1/t196
      t200 = 0.1D-2*t192*t197-0.1853571428571429D-2
      t201 = t200*t145
      t204 = 1/t166/t115
      t208 = dexp(-0.261211729852336D1*t146*t204*sigma)
      t209 = sigma*t152*t208
      t211 = 0.1D1*t201*t209
      zk(i) = -0.9305257363491D0*t26*t35-0.9305257363491D0*t60*t69+rho
     #*(-t89+t123+t140+t189+t211)
      t214 = t4*t25
      t218 = 1/t4/t13
      t221 = 0.3360034800657352D-1*t6*t218*t10
      t222 = t13*rhoa
      t224 = 1/t14/t222
      t227 = 1.D0+0.10000117555961D1*t17
      t228 = dsqrt(t227)
      t229 = 1/t228
      t231 = 0.3360054550205309D-1*sigmaaa*t224*t229
      t234 = 1/t4/t28/t13
      t235 = t27*t234
      t240 = -t221-t231-0.1088885204563779D-1*t235*t19
     #-0.4388082092272149D-1*t22*t224
      t241 = t5*t240
      t244 = t34**2
      t245 = 1/t244
      t247 = -t221-t231-0.5776579334555855D-5*t235
      t248 = t245*t247
      t255 = 0.1333333333333333D1*t106*t72-0.1333333333333333D1*t110*t72
      t257 = t101*t255*t120
      t258 = 0.3799575D-1*t257
      t259 = t113*t103
      t260 = t112*t259
      t262 = t101*t260*t117
      t263 = 0.151983D0*t262
      t264 = t137*t255
      t265 = t264*t118
      t266 = 0.1923661050931536D1*t265
      t267 = t259*t117
      t268 = t138*t267
      t269 = 0.7694644203726145D1*t268
      t270 = t146*t187
      t271 = 1/t106
      t274 = 1/t110
      t277 = 0.3333333333333333D0*t271*t72-0.3333333333333333D0*t274*t72
      t278 = t270*t277
      t279 = 0.7420670230672731D-1*t278
      t280 = sigma*t156
      t281 = t152*t277
      t284 = t179*t165
      t285 = t258-t263+t266+t269
      t288 = t155*t165
      t291 = -0.4042761511756372D2*t285*t156+0.1212828453526912D3*t288
     #*t277
      t292 = t168*t291
      t293 = t292*t159
      t296 = t164*t145
      t297 = 1/t296
      t298 = t297*t168
      t299 = t298*t277
      t302 = -0.1269364121954074D0*t280*t281-0.1086645186223595D-1
     #*t284*t293-0.4346580744894379D-1*t163*t299
      t305 = t182**2
      t306 = 1/t305
      t307 = t172*t306
      t308 = t178*sigma
      t309 = t308*t148
      t310 = t152*t291
      t311 = t310*t159
      t314 = t156*t152
      t315 = t314*t277
      t319 = 1/t177/t160
      t320 = t319*t162
      t321 = t320*t165
      t326 = -0.1712109500228824D0*t309*t311-0.3424219000457648D0*t173
     #*t315-0.5862637881547585D-1*t321*t293-0.1172527576309517D0*t179
     #*t299
      t329 = 0.2697586091519874D1*t302*t183-0.2697586091519874D1*t307
     #*t326
      t330 = t147*t329
      t331 = 1/t186
      t332 = t330*t331
      t333 = 0.2473556743557577D-1*t332
      t334 = t200*t277
      t335 = t334*t209
      t336 = 0.1D1*t335
      t338 = t200*t146*t162
      t339 = t116*rho
      t340 = 1/t339
      t342 = t340*t277*t208
      t343 = t338*t342
      t344 = 0.522423459704672D1*t343
      t347 = 1/t81
      t348 = 1/t115
      t349 = t347*t348
      t350 = t349*t87
      t351 = 0.2747799777968419D-2*t350
      t352 = t83**2
      t353 = 1/t352
      t354 = t75*t353
      t355 = t76**2
      t356 = t355**2
      t357 = t356*t76
      t358 = 1/t357
      t359 = t358*t348
      t362 = 1/t79
      t363 = t362*t348
      t365 = 1/t73
      t366 = t365*t348
      t368 = -0.99709173929518D0*t359-0.7418564737168958D0*t349
     #-0.4002143174996817D0*t363-0.1264669550498372D0*t366
      t369 = 1/t86
      t371 = t354*t368*t369
      t372 = 0.1D1*t371
      t373 = t100*t112
      t374 = t373*t120
      t375 = t349*t374
      t376 = 0.8740794636035784D-3*t375
      t377 = t96**2
      t378 = 1/t377
      t379 = t91*t378
      t384 = -0.135956911724794D1*t359-0.7491972878592054D0*t349
     #-0.2150486235638321D0*t363-0.1274341730084892D0*t366
      t385 = t379*t384
      t386 = 1/t99
      t387 = t386*t112
      t388 = t387*t120
      t389 = t385*t388
      t390 = 0.1125D1*t389
      t391 = t106*t103
      t394 = t110*t103
      t397 = -0.1333333333333333D1*t391*t348+0.1333333333333333D1*t394
     #*t348
      t399 = t101*t397*t120
      t400 = 0.3799575D-1*t399
      t401 = t112*t114
      t403 = t101*t401*t340
      t404 = 0.151983D0*t403
      t407 = t130**2
      t408 = 1/t407
      t409 = t125*t408
      t414 = -0.1853395810515781D1*t359-0.128158207914907D1*t349
     #-0.8223668877838045D0*t363-0.1603914194192128D0*t366
      t415 = 1/t133
      t421 = 0.1321010150222857D-2*t349*t134+0.1D1*t409*t414*t415
     #-0.2747799777968419D-2*t350-0.1D1*t371
      t422 = t421*t112
      t423 = t422*t118
      t424 = 0.1923661050931536D1*t423
      t425 = t137*t397
      t426 = t425*t118
      t427 = 0.1923661050931536D1*t426
      t428 = t114*t340
      t429 = t138*t428
      t430 = 0.7694644203726145D1*t429
      t431 = t271*t103
      t434 = t274*t103
      t437 = -0.3333333333333333D0*t431*t348+0.3333333333333333D0*t434
     #*t348
      t438 = t270*t437
      t440 = t152*t437
      t443 = t115*rho
      t445 = 1/t150/t443
      t448 = t351+t372-t376-t390+t400+t404+t424+t427-t430
      t453 = -0.4042761511756372D2*t448*t156+0.1212828453526912D3*t288
     #*t437
      t454 = t168*t453
      t455 = t454*t159
      t458 = t298*t437
      t462 = 1/t166/t339
      t463 = t165*t462
      t466 = -0.1269364121954074D0*t280*t440-0.148092480894642D0*t149
     #*t445-0.1086645186223595D-1*t284*t455-0.4346580744894379D-1*t163
     #*t458-0.5071010869043442D-1*t163*t463
      t469 = t152*t453
      t470 = t469*t159
      t473 = t314*t437
      t476 = t148*t445
      t485 = -0.1712109500228824D0*t309*t470-0.3424219000457648D0*t173
     #*t473-0.3994922167200589D0*t173*t476-0.5862637881547585D-1*t321
     #*t455-0.1172527576309517D0*t179*t458-0.136794883902777D0*t179*t463
      t488 = 0.2697586091519874D1*t466*t183-0.2697586091519874D1*t307
     #*t485
      t489 = t147*t488
      t490 = t489*t331
      t494 = -0.4811024840421814D1*t349-0.1895695887660258D-2*t366
      t497 = t196**2
      t498 = 1/t497
      t499 = t192*t498
      t503 = -0.1803772444038489D1*t349-0.1210946621972718D0*t366
     #-0.1763993811759022D-1*t348
      t506 = 0.1D-2*t494*t197-0.1D-2*t499*t503
      t507 = t506*t145
      t508 = t507*t209
      t510 = t200*t437
      t511 = t510*t209
      t514 = sigma*t445*t208
      t515 = t201*t514
      t517 = t201*sigma
      t518 = t145*t204
      t519 = sigma*t437
      t523 = 1/t166/t443
      t524 = t146*t523
      t527 = -0.522423459704672D1*t518*t519+0.6965646129395627D1*t524
     #*sigma
      t529 = t152*t527*t208
      t530 = t517*t529
      t532 = t351+t372-t376-t390+t400+t404+t424+t427-t430
     #+0.7420670230672731D-1*t438+0.2473556743557577D-1*t490+0.1D1
     #*t508+0.1D1*t511-0.2333333333333333D1*t515+0.1D1*t530
      t533 = rho*t532
      vrhoa(i) = -0.12407009817988D1*t214*t35-0.9305257363491D0*t241
     #*t35+0.9305257363491D0*t26*t248+rho*(t258-t263+t266+t269+t279
     #+t333+t336-t344)-t89+t123+t140+t189+t211+t533
      t534 = t38*t59
      t538 = 1/t38/t47
      t541 = 0.3360034800657352D-1*t40*t538*t44
      t542 = t47*rhob
      t544 = 1/t48/t542
      t547 = 1.D0+0.10000117555961D1*t51
      t548 = dsqrt(t547)
      t549 = 1/t548
      t551 = 0.3360054550205309D-1*sigmabb*t544*t549
      t554 = 1/t38/t62/t47
      t555 = t61*t554
      t560 = -t541-t551-0.1088885204563779D-1*t555*t53
     #-0.4388082092272149D-1*t56*t544
      t561 = t39*t560
      t564 = t68**2
      t565 = 1/t564
      t567 = -t541-t551-0.5776579334555855D-5*t555
      t568 = t565*t567
      t571 = -t255
      t573 = t101*t571*t120
      t574 = 0.3799575D-1*t573
      t575 = 0.151983D0*t262
      t576 = t137*t571
      t577 = t576*t118
      t578 = 0.1923661050931536D1*t577
      t579 = 0.7694644203726145D1*t268
      t580 = -t277
      t581 = t270*t580
      t582 = 0.7420670230672731D-1*t581
      t583 = t152*t580
      t586 = t574+t575+t578-t579
      t591 = -0.4042761511756372D2*t586*t156+0.1212828453526912D3*t288
     #*t580
      t592 = t168*t591
      t593 = t592*t159
      t596 = t298*t580
      t599 = -0.1269364121954074D0*t280*t583-0.1086645186223595D-1
     #*t284*t593-0.4346580744894379D-1*t163*t596
      t602 = t152*t591
      t603 = t602*t159
      t606 = t314*t580
      t613 = -0.1712109500228824D0*t309*t603-0.3424219000457648D0*t173
     #*t606-0.5862637881547585D-1*t321*t593-0.1172527576309517D0*t179
     #*t596
      t616 = 0.2697586091519874D1*t599*t183-0.2697586091519874D1*t307
     #*t613
      t617 = t147*t616
      t618 = t617*t331
      t619 = 0.2473556743557577D-1*t618
      t620 = t200*t580
      t621 = t620*t209
      t622 = 0.1D1*t621
      t624 = t340*t580*t208
      t625 = t338*t624
      t626 = 0.522423459704672D1*t625
      vrhob(i) = -0.12407009817988D1*t534*t69-0.9305257363491D0*t561
     #*t69+0.9305257363491D0*t60*t568+rho*(t574+t575+t578-t579+t582
     #+t619+t622-t626)-t89+t123+t140+t189+t211+t533
      t629 = 1/t6
      t632 = 0.1260013050246507D-1*t629*t7*t10
      t634 = 0.1260020456326991D-1*t16*t229
      t635 = t31*t19
      t640 = t632+t634+0.408331951711417D-2*t635*sigmaaa
     #+0.1645530784602056D-1*t21*t16
      t641 = t5*t640
      t646 = t632+t634+0.2166217250458446D-5*sigmaaa*t31
      t647 = t245*t646
      t653 = 0.634682060977037D-1*t174+0.2173290372447189D-1*t173*t169
      t656 = t161*t148
      t661 = 0.1712109500228824D0*t656*t152+0.5862637881547585D-1*t308
     #*t169
      t664 = 0.2697586091519874D1*t653*t183-0.2697586091519874D1*t307
     #*t661
      t666 = t147*t664*t331
      t667 = 0.2473556743557577D-1*t666
      t668 = t152*t208
      t669 = t201*t668
      t670 = 0.1D1*t669
      t671 = t200*t147
      t673 = sigma*t340*t208
      t674 = t671*t673
      t675 = 0.261211729852336D1*t674
      t677 = rho*(t667+t670-t675)
      vsigmaaa(i) = -0.9305257363491D0*t641*t35+0.9305257363491D0*t26
     #*t647+t677
      vsigmaab(i) = 2.D0*t677
      t678 = 1/t40
      t681 = 0.1260013050246507D-1*t678*t41*t44
      t683 = 0.1260020456326991D-1*t50*t549
      t684 = t65*t53
      t689 = t681+t683+0.408331951711417D-2*t684*sigmabb
     #+0.1645530784602056D-1*t55*t50
      t690 = t39*t689
      t695 = t681+t683+0.2166217250458446D-5*sigmabb*t65
      t696 = t565*t695
      vsigmabb(i) = -0.9305257363491D0*t690*t69+0.9305257363491D0*t60
     #*t696+t677
      t705 = 0.2D1*t371
      t706 = 0.759915D-1*t399
      t707 = 0.3847322101863073D1*t423
      t708 = 0.3847322101863073D1*t426
      t709 = 0.1538928840745229D2*t429
      t711 = 0.5495599555936838D-2*t350
      t712 = 0.1484134046134546D0*t438
      t714 = 0.303966D0*t403
      t715 = 0.1748158927207157D-2*t375
      t716 = 0.759915D-1*t257+0.3847322101863073D1*t265
     #+0.1538928840745229D2*t268+0.1484134046134546D0*t278
     #+0.4947113487115154D-1*t332-0.303966D0*t262+t705+t706+t707+t708
     #-t709+0.2D1*t335+t711+t712-0.1044846919409344D2*t343+t714-t715
      t717 = 0.225D1*t389
      t718 = t334*sigma
      t720 = 0.1D1*t718*t529
      t721 = 1/t107
      t723 = 1/t443
      t728 = 1/t111
      t734 = 0.1111111111111111D0*t721*t103*t723-0.3333333333333333D0
     #*t271*t348+0.1111111111111111D0*t728*t103*t723
     #+0.3333333333333333D0*t274*t348
      t736 = 0.7420670230672731D-1*t270*t734
      t737 = t116*t115
      t738 = 1/t737
      t741 = t338*t738*t277*t208
      t743 = sigma*t165
      t744 = t440*t277
      t746 = 0.3808092365862222D0*t743*t744
      t749 = 0.296184961789284D0*t280*t445*t277
      t752 = 0.1269364121954074D0*t280*t152*t734
      t754 = dexp(-0.8085523023512745D2*t157)
      t755 = t754*t291
      t756 = t454*t755
      t758 = 0.2173290372447189D-1*t321*t756
      t759 = t179*t297
      t760 = t168*t437
      t761 = t291*t159
      t762 = t760*t761
      t764 = 0.4346580744894379D-1*t759*t762
      t766 = t462*t291*t159
      t768 = 0.5071010869043442D-1*t284*t766
      t772 = 0.8740794636035784D-3*t349*t100*t255*t120
      t776 = 0.1125D1*t385*t386*t255*t120
      t777 = 1/t141
      t783 = 1/t143
      t789 = -0.4444444444444444D0*t777*t103*t723-0.1333333333333333D1
     #*t106*t348-0.4444444444444444D0*t783*t103*t723
     #+0.1333333333333333D1*t110*t348
      t792 = 0.3799575D-1*t101*t789*t120
      t796 = 0.151983D0*t101*t255*t114*t340
      t799 = t347*t738*t373*t259
      t800 = 0.3496317854414314D-2*t799
      t802 = t385*t387*t267
      t803 = 0.45D1*t802
      t806 = t101*t397*t259*t117
      t807 = 0.151983D0*t806
      t809 = t101*t260*t340
      t810 = 0.607932D0*t809
      t813 = 0.1923661050931536D1*t421*t255*t118
      t816 = 0.1923661050931536D1*t137*t789*t118
      t818 = 0.7694644203726145D1*t264*t428
      t819 = t422*t267
      t820 = 0.7694644203726145D1*t819
      t821 = t425*t267
      t822 = 0.7694644203726145D1*t821
      t824 = t138*t259*t340
      t825 = 0.3077857681490458D2*t824
      t826 = -t772-t776+t792+t796+t800+t803-t807+t810+t813+t816-t818
     #+t820+t822-t825
      t829 = t285*t165
      t831 = 0.1212828453526912D3*t829*t437
      t832 = t448*t165
      t834 = 0.1212828453526912D3*t832*t277
      t835 = t155*t297
      t836 = t437*t277
      t838 = 0.4851313814107647D3*t835*t836
      t840 = 0.1212828453526912D3*t288*t734
      t841 = -0.4042761511756372D2*t826*t156+t831+t834-t838+t840
      t843 = t168*t841*t159
      t846 = t454*t761
      t848 = 0.1086645186223595D-1*t284*t846
      t849 = t159*t277
      t850 = t454*t849
      t852 = 0.4346580744894379D-1*t759*t850
      t854 = 1/t164/t146
      t855 = t163*t854
      t856 = t760*t277
      t858 = 0.2173290372447189D0*t855*t856
      t859 = t297*t462
      t860 = t859*t277
      t862 = 0.2028404347617377D0*t163*t860
      t863 = t298*t734
      t865 = 0.4346580744894379D-1*t163*t863
      t866 = t746+t749-t752+t758+t764+t768-0.1086645186223595D-1*t284
     #*t843-t848+t852+t858+t862-t865
      t869 = t302*t306
      t871 = 0.2697586091519874D1*t869*t485
      t872 = t466*t306
      t874 = 0.2697586091519874D1*t872*t326
      t877 = t172/t305/t182
      t880 = 0.5395172183039748D1*t877*t485*t326
      t881 = t319*sigma
      t882 = t881*t148
      t885 = 0.3424219000457648D0*t882*t469*t755
      t886 = t308*t156
      t889 = 0.3424219000457648D0*t886*t440*t761
      t893 = 0.3994922167200589D0*t309*t445*t291*t159
      t900 = 0.1712109500228824D0*t309*t469*t761
      t903 = 0.3424219000457648D0*t886*t469*t849
      t904 = t173*t165
      t906 = 0.1027265700137294D1*t904*t744
      t907 = t156*t445
      t910 = 0.7989844334401178D0*t173*t907*t277
      t913 = 0.3424219000457648D0*t173*t314*t734
      t914 = t177**2
      t917 = 1/t914*t162*t165
      t919 = 0.1758791364464276D0*t917*t756
      t920 = t320*t297
      t922 = 0.2345055152619034D0*t920*t762
      t924 = 0.273589767805554D0*t321*t766
      t928 = 0.5862637881547585D-1*t321*t846
      t930 = 0.2345055152619034D0*t920*t850
      t931 = t179*t854
      t933 = 0.5862637881547585D0*t931*t856
      t935 = 0.547179535611108D0*t179*t860
      t937 = 0.1172527576309517D0*t179*t863
      t938 = t885+t889+t893-0.1712109500228824D0*t309*t152*t841*t159
     #-t900+t903+t906+t910-t913+t919+t922+t924-0.5862637881547585D-1
     #*t321*t843-t928+t930+t933+t935-t937
      t949 = t186**2
      t950 = 1/t949
      t953 = 0.2473556743557577D-1*t489*t950*t329
      t954 = t510*t162
      t955 = t340*t145
      t956 = t277*t208
      t958 = t954*t955*t956
      t961 = t506*t146*t162
      t963 = 0.522423459704672D1*t961*t342
      t964 = t146*t329
      t965 = t331*t437
      t967 = 0.7420670230672731D-1*t964*t965
      t970 = 0.1D1*t200*t734*t209
      t971 = t146*t488
      t972 = t331*t277
      t974 = 0.7420670230672731D-1*t971*t972
      t975 = t720+t736+0.261211729852336D2*t741+0.2473556743557577D-1
     #*t147*(0.2697586091519874D1*t866*t183-t871-t874+t880
     #-0.2697586091519874D1*t307*t938)*t331-0.522423459704672D1*t338
     #*t340*t734*t208-t953-0.1044846919409344D2*t958-t963+t967+t970
     #+t974+t792+t820+t813
      t976 = t145*t187
      t978 = 0.1484134046134546D0*t976*t836
      t981 = 0.1D1*t506*t277*t209
      t982 = t340*t527
      t985 = 0.522423459704672D1*t338*t982*t956
      t987 = 0.2333333333333333D1*t334*t514
      t988 = t816+t822-t818-t825+t978+t803+t800-t776-t772+t796-t807
     #+t810+t981-t985-t987
      t990 = rho*(t975+t988)
      t991 = 0.4947113487115154D-1*t490
      t993 = 0.45D1*t802
      t994 = -t772+t800-t776+t993+t792-t807+t796+t810+t813+t820+t816
     #+t822-t818-t825
      t997 = -0.4042761511756372D2*t994*t156+t834+t831-t838+t840
      t999 = t168*t997*t159
      t1002 = t746-t752+t749+t758+t852-0.1086645186223595D-1*t284*t999
     #-t848+t764+t858-t865+t768+t862
      t1011 = t885+t903-0.1712109500228824D0*t309*t152*t997*t159-t900
     #+t889+t906-t913+t893+t910+t919+t930-0.5862637881547585D-1*t321
     #*t999-t928+t922+t933-t937+t924+t935
      t1019 = t720+t736+0.1218988072644235D2*t741
     #+0.2473556743557577D-1*t147*(0.2697586091519874D1*t1002*t183
     #-t874-t871+t880-0.2697586091519874D1*t307*t1011)*t331-t953
     #-0.522423459704672D1*t958-t963+t967+t970+t974+t792+t820+t813+t816
      t1026 = t145*t523
      t1035 = t822-t818-t825+t978+t993+t800-t776-t772+t796-t807+t810
     #+t981+0.1D1*t517*t152*(-0.522423459704672D1*t277*t204*t519
     #-0.522423459704672D1*t518*sigma*t734+0.1393129225879125D2*t1026
     #*sigma*t277)*t208-t985-t987
      t1037 = rho*(t1019+t1035)
      t1038 = 0.2D1*t508
      t1039 = 0.2D1*t511
      t1040 = 0.4666666666666667D1*t515
      t1041 = 0.2D1*t530
      t1052 = -0.1111111111111111D0*t721*t113*t117
     #+0.6666666666666667D0*t431*t723-0.1111111111111111D0*t728*t113
     #*t117-0.6666666666666667D0*t434*t723
      t1055 = 0.1D1*t200*t1052*t209
      t1058 = 1/t81/t72*t117
      t1059 = t1058*t87
      t1060 = 0.1831866518645613D-2*t1059
      t1063 = 0.3847322101863073D1*t421*t397*t118
      t1064 = t352**2
      t1067 = t368**2
      t1068 = t86**2
      t1071 = t75/t1064*t1067/t1068
      t1072 = 0.160818243221511D2*t1071
      t1073 = t347*t723
      t1074 = t1073*t87
      t1075 = 0.5495599555936838D-2*t1074
      t1079 = t384**2
      t1082 = 0.225D1*t91/t377/t96*t1079*t388
      t1093 = 0.4444444444444444D0*t777*t113*t117+0.2666666666666667D1
     #*t391*t723+0.4444444444444444D0*t783*t113*t117
     #-0.2666666666666667D1*t394*t723
      t1096 = 0.1923661050931536D1*t137*t1093*t118
      t1098 = 0.1538928840745229D2*t425*t428
      t1110 = t414**2
      t1116 = 1/t357/t72*t117
      t1118 = t358*t723
      t1124 = 1/t79/t72*t117
      t1126 = t362*t723
      t1130 = 1/t73/t72*t117
      t1132 = t365*t723
      t1138 = t407**2
      t1141 = t133**2
      t1150 = t349*t353*t368*t369
      t1156 = t75/t352/t83*t1067*t369
      t1168 = t354*(-0.8309097827459833D0*t1116+0.199418347859036D1
     #*t1118-0.4945709824779306D0*t1058+0.1483712947433792D1*t1073
     #-0.2001071587498409D0*t1124+0.8004286349993634D0*t1126
     #-0.4215565168327908D-1*t1130+0.2529339100996745D0*t1132)*t369
      t1171 = 0.8806734334819047D-3*t1058*t134-0.2642020300445714D-2
     #*t1073*t134-0.8497974591333914D-1*t349*t408*t414*t415-0.2D1*t125
     #/t407/t130*t1110*t415+0.1D1*t409*(-0.1544496508763151D1*t1116
     #+0.3706791621031562D1*t1118-0.854388052766047D0*t1058
     #+0.2563164158298141D1*t1073-0.4111834438919023D0*t1124
     #+0.1644733775567609D1*t1126-0.5346380647307093D-1*t1130
     #+0.3207828388384256D0*t1132)*t415+0.321646831778707D2*t125/t1138
     #*t1110/t1141-0.1831866518645613D-2*t1059+0.5495599555936838D-2
     #*t1074+0.8837926660346786D-1*t1150+0.2D1*t1156-0.1D1*t1168
     #-0.160818243221511D2*t1071
      t1174 = 0.1923661050931536D1*t1171*t112*t118
      t1177 = 0.3799575D-1*t101*t1093*t120
      t1179 = 0.1538928840745229D2*t422*t428
      t1180 = t159*t437
      t1181 = t454*t1180
      t1183 = 0.8693161489788758D-1*t759*t1181
      t1185 = t462*t453*t159
      t1187 = 0.1014202173808688D0*t284*t1185
      t1191 = 0.225D1*t385*t386*t397*t120
      t1194 = 0.9D1*t385*t387*t428
      t1195 = t1060+t1063+t1072-t1075+t1082+t1096-t1098+t1174+t1177
     #-t1179-t1191-t1194
      t1196 = 0.2D1*t1156
      t1200 = 0.303966D0*t101*t397*t114*t340
      t1212 = 0.1125D1*t379*(-0.1132974264373283D1*t1116
     #+0.271913823449588D1*t1118-0.4994648585728036D0*t1058
     #+0.1498394575718411D1*t1073-0.1075243117819161D0*t1124
     #+0.4300972471276643D0*t1126-0.4247805766949639D-1*t1130
     #+0.2548683460169784D0*t1132)*t388
      t1213 = 0.8837926660346786D-1*t1150
      t1216 = 0.3847322101863073D2*t138*t114*t738
      t1217 = t116*t443
      t1222 = 0.6992635708828627D-2*t347/t1217*t373*t114
      t1223 = t377**2
      t1227 = t99**2
      t1232 = 0.3330964647361876D2*t91/t1223*t1079/t1227*t112*t120
      t1234 = 0.5827196424023856D-3*t1058*t374
      t1238 = 0.1748158927207157D-2*t349*t100*t397*t120
      t1240 = 0.1748158927207157D-2*t1073*t374
      t1245 = 0.5176049408441869D-1*t349*t378*t384*t386*t121
      t1248 = 0.759915D0*t101*t401*t738
      t1249 = 0.1D1*t1168
      t1250 = -t1196+t1200-t1212-t1213+t1216-t1222-t1232-t1234-t1238
     #+t1240+t1245-t1248+t1249
      t1256 = t437**2
      t1261 = -0.4042761511756372D2*(t1195+t1250)*t156
     #+0.2425656907053823D3*t832*t437-0.4851313814107647D3*t835*t1256
     #+0.1212828453526912D3*t288*t1052
      t1263 = t168*t1261*t159
      t1265 = 0.1086645186223595D-1*t284*t1263
      t1268 = 0.3808092365862222D0*t743*t152*t1256
      t1271 = 0.5923699235785679D0*t280*t445*t437
      t1274 = 0.1269364121954074D0*t280*t152*t1052
      t1276 = 1/t150/t116
      t1278 = 0.4936416029821399D0*t149*t1276
      t1279 = t453**2
      t1280 = t168*t1279
      t1281 = t1280*t754
      t1283 = 0.2173290372447189D-1*t321*t1281
      t1284 = t1280*t159
      t1286 = 0.1086645186223595D-1*t284*t1284
      t1289 = t165/t166/t737
      t1290 = t163*t1289
      t1292 = t854*t168
      t1293 = t1292*t1256
      t1295 = 0.2173290372447189D0*t163*t1293
      t1296 = t859*t437
      t1298 = 0.4056808695234754D0*t163*t1296
      t1299 = t298*t1052
      t1301 = 0.4346580744894379D-1*t163*t1299
      t1302 = t1183+t1187-t1265+t1268+t1271-t1274+t1278+t1283-t1286
     #+0.2873572825791284D0*t1290+t1295+t1298-t1301
      t1306 = 0.5395172183039748D1*t872*t485
      t1307 = t485**2
      t1309 = 0.5395172183039748D1*t877*t1307
      t1329 = t165*t152
      t1337 = t152*t1279
      t1359 = 0.5862637881547585D0*t179*t1293-0.1712109500228824D0
     #*t309*t152*t1261*t159+0.6848438000915295D0*t886*t469*t1180
     #+0.1597968866880236D1*t173*t907*t437-0.3424219000457648D0*t173
     #*t314*t1052+0.7989844334401178D0*t309*t445*t453*t159
     #+0.1027265700137294D1*t173*t1329*t1256-0.5862637881547585D-1
     #*t321*t1284+0.7751710087824029D0*t179*t1289+0.3424219000457648D0
     #*t882*t1337*t754-0.1712109500228824D0*t309*t1337*t159
     #+0.1331640722400196D1*t173*t148*t1276+0.1094359071222216D1*t179
     #*t1296+0.4690110305238068D0*t920*t1181+0.547179535611108D0*t321
     #*t1185-0.1172527576309517D0*t179*t1299-0.5862637881547585D-1
     #*t321*t1263+0.1758791364464276D0*t917*t1281
      t1361 = 0.2697586091519874D1*t307*t1359
      t1368 = 0.2D1*t507*sigma*t529
      t1369 = t488**2
      t1372 = 0.2473556743557577D-1*t147*t1369*t950
      t1374 = 0.7420670230672731D-1*t270*t1052
      t1375 = t1055+t1060+t1063+t1072-t1075+t1082+t1096-t1098+t1174
     #+t1177-t1179+0.2473556743557577D-1*t147*(0.2697586091519874D1
     #*t1302*t183-t1306+t1309-t1361)*t331-t1191+t1368-t1194-t1196
     #+t1200-t1372+t1374-t1212
      t1377 = 0.4666666666666667D1*t510*t514
      t1381 = 0.4666666666666667D1*t517*t445*t527*t208
      t1385 = 0.7777777777777778D1*t201*sigma*t1276*t208
      t1386 = t510*sigma
      t1388 = 0.2D1*t1386*t529
      t1390 = 0.4666666666666667D1*t507*t514
      t1393 = 0.2D1*t506*t437*t209
      t1394 = t527**2
      t1398 = 0.1D1*t517*t152*t1394*t208
      t1412 = t503**2
      t1426 = 0.1D1*(0.1D-2*(-0.3207349893614542D1*t1058
     #+0.9622049680843627D1*t1073-0.6318986292200858D-3*t1130
     #+0.3791391775320515D-2*t1132)*t197-0.2D-2*t494*t498*t503+0.2D-2
     #*t192/t497/t196*t1412-0.1D-2*t499*(-0.1202514962692326D1*t1058
     #+0.3607544888076978D1*t1073-0.4036488739909061D-1*t1130
     #+0.2421893243945437D0*t1132+0.3527987623518044D-1*t723))*t145*t209
      t1442 = 0.1D1*t517*t152*(-0.522423459704672D1*t1256*t204*sigma
     #+0.2786258451758251D2*t1026*t519-0.522423459704672D1*t518*sigma
     #*t1052-0.2554070247445063D2*t146*t168*sigma)*t208
      t1448 = t1216-t1222-t1232-t1234-t1238+t1240+t1245
     #+0.1484134046134546D0*t971*t965-t1248+t1249+0.1484134046134546D0
     #*t976*t1256
      t1451 = rho*(t1375-t1377-t1381-t1213+t1385+t1388-t1390+t1393
     #+t1398+t1426+t1442+t1448)
      t1463 = 0.7840081201533821D-1*t6/t4/t222*t10
      t1465 = 1/t14/t28
      t1468 = 0.1680027275102654D0*sigmaaa*t1465*t229
      t1472 = t27/t4/t28/t222
      t1474 = 1/t228/t227
      t1476 = 0.4480125399532632D-1*t1472*t1474
      t1480 = t28**2
      t1502 = 1/t244/t34
      t1503 = t247**2
      t1511 = 0.4444444444444444D0*t777*t348+0.4444444444444444D0*t783
     #*t348
      t1514 = 0.3799575D-1*t101*t1511*t120
      t1517 = t101*t255*t259*t117
      t1518 = 0.303966D0*t1517
      t1521 = t101*t112*t113*t117
      t1522 = 0.455949D0*t1521
      t1525 = 0.1923661050931536D1*t137*t1511*t118
      t1526 = t264*t267
      t1527 = 0.1538928840745229D2*t1526
      t1529 = t138*t113*t117
      t1530 = 0.2308393261117844D2*t1529
      t1531 = t277**2
      t1540 = -0.1111111111111111D0*t721*t348-0.1111111111111111D0
     #*t728*t348
      t1542 = 0.7420670230672731D-1*t270*t1540
      t1548 = 0.1269364121954074D0*t280*t152*t1540
      t1549 = t291**2
      t1550 = t168*t1549
      t1551 = t1550*t754
      t1554 = t292*t849
      t1565 = 0.1212828453526912D3*t288*t1540
      t1566 = -0.4042761511756372D2*(t1514-t1518-t1522+t1525+t1527
     #+t1530)*t156+0.2425656907053823D3*t829*t277-0.4851313814107647D3
     #*t835*t1531+t1565
      t1568 = t168*t1566*t159
      t1571 = t1550*t159
      t1574 = t1292*t1531
      t1577 = t298*t1540
      t1579 = 0.4346580744894379D-1*t163*t1577
      t1585 = t326**2
      t1588 = t152*t1549
      t1607 = 0.3424219000457648D0*t173*t314*t1540
      t1619 = 0.1172527576309517D0*t179*t1577
      t1620 = 0.3424219000457648D0*t882*t1588*t754
     #+0.6848438000915295D0*t886*t310*t849-0.1712109500228824D0*t309
     #*t152*t1566*t159-0.1712109500228824D0*t309*t1588*t159
     #+0.1027265700137294D1*t173*t1329*t1531-t1607
     #+0.1758791364464276D0*t917*t1551+0.4690110305238068D0*t920*t1554
     #-0.5862637881547585D-1*t321*t1568-0.5862637881547585D-1*t321
     #*t1571+0.5862637881547585D0*t179*t1574-t1619
      t1627 = t329**2
      t1633 = 0.1D1*t200*t1540*t209
      t1636 = t955*t208
      t1642 = 0.522423459704672D1*t338*t340*t1540*t208
      t1644 = t671*t162*sigma
      t1646 = 1/t166/t1217
      t1651 = t1514-t1518-t1522+t1525+t1527+t1530+0.1484134046134546D0
     #*t976*t1531+0.1484134046134546D0*t964*t972+t1542
     #+0.2473556743557577D-1*t147*(0.2697586091519874D1*
     #(0.3808092365862222D0*t743*t152*t1531-t1548
     #+0.2173290372447189D-1*t321*t1551+0.8693161489788758D-1*t759
     #*t1554-0.1086645186223595D-1*t284*t1568-0.1086645186223595D-1
     #*t284*t1571+0.2173290372447189D0*t163*t1574-t1579)*t183
     #-0.5395172183039748D1*t869*t326+0.5395172183039748D1*t877*t1585
     #-0.2697586091519874D1*t307*t1620)*t331-0.2473556743557577D-1
     #*t147*t1627*t950+t1633-0.1567270379114016D2*t200*t1531*t162
     #*t1636-t1642+0.272926271249799D2*t1644*t1646*t1531*t208
      t1653 = -t717+t990+t991+t1037+t1038+t1039-t1040+t1041+t1451
     #-0.4135669939329333D0/t14*t25*t35-0.24814019635976D1*t4*t240*t35
     #-0.9305257363491D0*t5*(t1463+t1468-t1476+0.9799966841074009D-1
     #*t1472*t19-0.4778117666686413D-1*t27*sigmaaa/t1480/t13*t19
     #+0.1608963433833121D0*t22*t1465)*t35+0.18610514726982D1*t241
     #*t248+0.24814019635976D1*t214*t248+0.9305257363491D0*t26*t245*
     #(t1463+t1468-t1476+0.3658500245218708D-4*t1472)
     #-0.18610514726982D1*t26*t1502*t1503+rho*t1651
      v2rhoa2(i) = t716+t1653
      t1656 = t950*t616
      t1658 = 0.2473556743557577D-1*t489*t1656
      t1661 = 0.1923661050931536D1*t421*t571*t118
      t1662 = -t789
      t1665 = 0.1923661050931536D1*t137*t1662*t118
      t1667 = 0.7694644203726145D1*t576*t428
      t1668 = t146*t616
      t1670 = 0.7420670230672731D-1*t1668*t965
      t1674 = 0.151983D0*t101*t571*t114*t340
      t1675 = -t734
      t1677 = 0.7420670230672731D-1*t270*t1675
      t1680 = 0.1D1*t506*t580*t209
      t1682 = 0.522423459704672D1*t961*t624
      t1683 = t437*t580
      t1685 = 0.1484134046134546D0*t976*t1683
      t1688 = 0.1D1*t200*t1675*t209
      t1695 = t338*t738*t580*t208
      t1698 = 0.2333333333333333D1*t620*t514
      t1699 = -t1658+t1661+t1665-t1667+t1670+t1674+t1677+t1680-t1682
     #+t1685+t1688-0.522423459704672D1*t338*t340*t1675*t208
     #+0.261211729852336D2*t1695-t1698
      t1700 = t580*t208
      t1701 = t955*t1700
      t1702 = t954*t1701
      t1704 = t620*sigma
      t1706 = 0.1D1*t1704*t529
      t1709 = 0.522423459704672D1*t338*t982*t1700
      t1710 = t331*t580
      t1712 = 0.7420670230672731D-1*t971*t1710
      t1715 = 0.3799575D-1*t101*t1662*t120
      t1716 = t440*t580
      t1725 = t754*t591
      t1726 = t454*t1725
      t1729 = t591*t159
      t1730 = t760*t1729
      t1734 = t462*t591*t159
      t1740 = 0.8740794636035784D-3*t349*t100*t571*t120
      t1744 = 0.1125D1*t385*t386*t571*t120
      t1745 = 0.3496317854414314D-2*t799
      t1746 = 0.45D1*t802
      t1747 = 0.151983D0*t806
      t1748 = 0.607932D0*t809
      t1749 = 0.7694644203726145D1*t819
      t1750 = 0.7694644203726145D1*t821
      t1751 = 0.3077857681490458D2*t824
      t1752 = -t1740-t1744+t1715+t1674-t1745-t1746+t1747-t1748+t1661
     #+t1665-t1667-t1749-t1750+t1751
      t1755 = t586*t165
      t1764 = -0.4042761511756372D2*t1752*t156+0.1212828453526912D3
     #*t1755*t437+0.1212828453526912D3*t832*t580-0.4851313814107647D3
     #*t835*t1683+0.1212828453526912D3*t288*t1675
      t1766 = t168*t1764*t159
      t1769 = t454*t1729
      t1772 = t159*t580
      t1773 = t454*t1772
      t1776 = t760*t580
      t1779 = t859*t580
      t1782 = t298*t1675
      t1785 = 0.3808092365862222D0*t743*t1716+0.296184961789284D0*t280
     #*t445*t580-0.1269364121954074D0*t280*t152*t1675
     #+0.2173290372447189D-1*t321*t1726+0.4346580744894379D-1*t759
     #*t1730+0.5071010869043442D-1*t284*t1734-0.1086645186223595D-1
     #*t284*t1766-0.1086645186223595D-1*t284*t1769
     #+0.4346580744894379D-1*t759*t1773+0.2173290372447189D0*t855
     #*t1776+0.2028404347617377D0*t163*t1779-0.4346580744894379D-1
     #*t163*t1782
      t1788 = t599*t306
      t1842 = 0.3424219000457648D0*t882*t469*t1725
     #+0.3424219000457648D0*t886*t440*t1729+0.3994922167200589D0*t309
     #*t445*t591*t159-0.1712109500228824D0*t309*t152*t1764*t159
     #-0.1712109500228824D0*t309*t469*t1729+0.3424219000457648D0*t886
     #*t469*t1772+0.1027265700137294D1*t904*t1716+0.7989844334401178D0
     #*t173*t907*t580-0.3424219000457648D0*t173*t314*t1675
     #+0.1758791364464276D0*t917*t1726+0.2345055152619034D0*t920*t1730
     #+0.273589767805554D0*t321*t1734-0.5862637881547585D-1*t321*t1766
     #-0.5862637881547585D-1*t321*t1769+0.2345055152619034D0*t920
     #*t1773+0.5862637881547585D0*t931*t1776+0.547179535611108D0*t179
     #*t1779-0.1172527576309517D0*t179*t1782
      t1848 = 0.2473556743557577D-1*t147*(0.2697586091519874D1*t1785
     #*t183-0.2697586091519874D1*t1788*t485-0.2697586091519874D1*t872
     #*t613+0.5395172183039748D1*t877*t485*t613-0.2697586091519874D1
     #*t307*t1842)*t331
      t1849 = -0.1044846919409344D2*t1702+t1706-t1709+t1712+t1715
     #+t1848-t1744-t1740-t1749-t1750+t1751-t1746-t1745+t1747-t1748
      t1851 = rho*(t1699+t1849)
      t1854 = t101*t571*t259*t117
      t1855 = 0.303966D0*t1854
      t1856 = t576*t267
      t1857 = 0.1538928840745229D2*t1856
      t1858 = t580**2
      t1866 = t591**2
      t1867 = t168*t1866
      t1868 = t1867*t754
      t1871 = t592*t1772
      t1881 = -0.4042761511756372D2*(t1514+t1855-t1522+t1525-t1857
     #+t1530)*t156+0.2425656907053823D3*t1755*t580
     #-0.4851313814107647D3*t835*t1858+t1565
      t1883 = t168*t1881*t159
      t1886 = t1867*t159
      t1889 = t1292*t1858
      t1897 = t613**2
      t1900 = t152*t1866
      t1927 = 0.3424219000457648D0*t882*t1900*t754
     #+0.6848438000915295D0*t886*t602*t1772-0.1712109500228824D0*t309
     #*t152*t1881*t159-0.1712109500228824D0*t309*t1900*t159
     #+0.1027265700137294D1*t173*t1329*t1858-t1607
     #+0.1758791364464276D0*t917*t1868+0.4690110305238068D0*t920*t1871
     #-0.5862637881547585D-1*t321*t1883-0.5862637881547585D-1*t321
     #*t1886+0.5862637881547585D0*t179*t1889-t1619
      t1934 = t616**2
      t1946 = t1514+t1855-t1522+t1525-t1857+t1530+0.1484134046134546D0
     #*t976*t1858+0.1484134046134546D0*t1668*t1710+t1542
     #+0.2473556743557577D-1*t147*(0.2697586091519874D1*
     #(0.3808092365862222D0*t743*t152*t1858-t1548
     #+0.2173290372447189D-1*t321*t1868+0.8693161489788758D-1*t759
     #*t1871-0.1086645186223595D-1*t284*t1883-0.1086645186223595D-1
     #*t284*t1886+0.2173290372447189D0*t163*t1889-t1579)*t183
     #-0.5395172183039748D1*t1788*t613+0.5395172183039748D1*t877*t1897
     #-0.2697586091519874D1*t307*t1927)*t331-0.2473556743557577D-1
     #*t147*t1934*t950+t1633-0.1567270379114016D2*t200*t1858*t162
     #*t1636-t1642+0.272926271249799D2*t1644*t1646*t1858*t208
      t1951 = 1/t564/t68
      t1952 = t567**2
      t1960 = 0.7840081201533821D-1*t40/t38/t542*t44
      t1962 = 1/t48/t62
      t1965 = 0.1680027275102654D0*sigmabb*t1962*t549
      t1969 = t61/t38/t62/t542
      t1971 = 1/t548/t547
      t1973 = 0.4480125399532632D-1*t1969*t1971
      t1979 = -0.1538928840745229D2*t268+0.303966D0*t262+t705+t706
     #+t707+t708-t709+t711+t712+t1851+t714-t715-t717+rho*t1946
     #+0.24814019635976D1*t534*t568-0.18610514726982D1*t60*t1951*t1952
     #+0.9305257363491D0*t60*t565*(t1960+t1965-t1973
     #+0.3658500245218708D-4*t1969)
      t1999 = -t1658+t1661+t1665-t1667+t1670+t1674+t1677+t1680-t1682
     #+t1685+t1688+0.1D1*t517*t152*(-0.522423459704672D1*t580*t204
     #*t519-0.522423459704672D1*t518*sigma*t1675+0.1393129225879125D2
     #*t1026*sigma*t580)*t208+0.1218988072644235D2*t1695-t1698
      t2001 = -0.522423459704672D1*t1702+t1706+t1848-t1709+t1712+t1715
     #-t1744-t1740-t1749-t1750+t1751-t1746-t1745+t1747-t1748
      t2003 = rho*(t1999+t2001)
      t2012 = t62**2
      t2030 = 0.18610514726982D1*t561*t568+0.759915D-1*t573
     #+0.3847322101863073D1*t577+t2003+t991+0.1484134046134546D0*t581
     #+t1038+t1039-t1040+t1041+0.4947113487115154D-1*t618
     #-0.24814019635976D1*t38*t560*t69-0.9305257363491D0*t39*(t1960
     #+t1965-t1973+0.9799966841074009D-1*t1969*t53
     #-0.4778117666686413D-1*t61*sigmabb/t2012/t47*t53
     #+0.1608963433833121D0*t56*t1962)*t69-0.4135669939329333D0/t48
     #*t59*t69+0.2D1*t621-0.1044846919409344D2*t625+t1451
      v2rhob2(i) = t1979+t2030
      t2032 = t258+t266+t279+t333+t705+t706+t707+t708-t709+t336+t711
     #+t712-t344+0.5D0*t1851+t714-t715
      t2034 = -t1511
      t2037 = 0.3799575D-1*t101*t2034*t120
      t2038 = 0.151983D0*t1517
      t2039 = 0.151983D0*t1854
      t2040 = 0.455949D0*t1521
      t2043 = 0.1923661050931536D1*t137*t2034*t118
      t2044 = 0.7694644203726145D1*t1526
      t2045 = 0.7694644203726145D1*t1856
      t2046 = 0.2308393261117844D2*t1529
      t2047 = t277*t580
      t2052 = -t1540
      t2057 = t281*t580
      t2063 = t292*t1725
      t2066 = t292*t1772
      t2080 = -0.4042761511756372D2*(t2037+t2038-t2039+t2040+t2043
     #-t2044+t2045-t2046)*t156+0.1212828453526912D3*t829*t580
     #+0.1212828453526912D3*t1755*t277-0.4851313814107647D3*t835*t2047
     #+0.1212828453526912D3*t288*t2052
      t2082 = t168*t2080*t159
      t2085 = t292*t1729
      t2088 = t168*t277
      t2089 = t2088*t1729
      t2092 = t2088*t580
      t2095 = t298*t2052
      t2143 = 0.3424219000457648D0*t882*t310*t1725
     #+0.3424219000457648D0*t886*t310*t1772-0.1712109500228824D0*t309
     #*t152*t2080*t159-0.1712109500228824D0*t309*t310*t1729
     #+0.3424219000457648D0*t886*t281*t1729+0.1027265700137294D1*t904
     #*t2057-0.3424219000457648D0*t173*t314*t2052+0.1758791364464276D0
     #*t917*t2063+0.2345055152619034D0*t920*t2066
     #-0.5862637881547585D-1*t321*t2082-0.5862637881547585D-1*t321
     #*t2085+0.2345055152619034D0*t920*t2089+0.5862637881547585D0*t931
     #*t2092-0.1172527576309517D0*t179*t2095
      t2162 = t1646*t277
      s1 = t2037+t2038-t2039+t2040+t2043-t2044+t2045-t2046
     #+0.1484134046134546D0*t976*t2047
      s2 = s1+0.7420670230672731D-1*t1668*t972+0.7420670230672731D-1
     #*t270*t2052+0.7420670230672731D-1*t964*t1710
     #+0.2473556743557577D-1*t147*(0.2697586091519874D1*
     #(0.3808092365862222D0*t743*t2057-0.1269364121954074D0*t280*t152
     #*t2052+0.2173290372447189D-1*t321*t2063+0.4346580744894379D-1
     #*t759*t2066-0.1086645186223595D-1*t284*t2082
     #-0.1086645186223595D-1*t284*t2085+0.4346580744894379D-1*t759
     #*t2089+0.2173290372447189D0*t855*t2092-0.4346580744894379D-1
     #*t163*t2095)*t183-0.2697586091519874D1*t869*t613
     #-0.2697586091519874D1*t1788*t326+0.5395172183039748D1*t877*t326
     #*t613-0.2697586091519874D1*t307*t2143)*t331
      t2166 = s2-0.2473556743557577D-1*t330*t1656+0.1D1*t200*t2052
     #*t209-0.1567270379114016D2*t334*t162*t1701-0.522423459704672D1
     #*t338*t340*t2052*t208+0.272926271249799D2*t1644*t2162*t1700
      t2170 = t1055+t1060+t1063+t1072-t1075+t1082+t1096-t1098+t1174
     #+t1177-t1179-t1191+t1368-t1194-t1196+t1200-t1372+t1374-t1212-t1377
      t2172 = t1278+t1183+t1187-t1265+t1268+t1271-t1274+t1283+t1295
     #+t1298-t1301-t1286+0.2873572825791284D0*t1290
      t2183 = -t717+0.5D0*t990+rho*t2166+t574+t578+0.5D0*t2003+t991
     #+t582+0.5D0*t1037+t1038+t1039-t1040+t1041+t619+t622-t626+rho*
     #(t2170-t1381-t1213+t1385+t1388-t1390+0.2473556743557577D-1*t147*
     #(0.2697586091519874D1*t2172*t183-t1306+t1309-t1361)*t331+t1393
     #+t1398+t1426+t1442+t1448)
      v2rhoab(i) = t2032+t2183
      t2191 = 0.1680017400328676D-1*t629*t218*t10
      t2193 = 0.5040081825307963D-1*t224*t229
      t2194 = sigmaaa*t234
      t2196 = 0.1680047024824737D-1*t2194*t1474
      t2224 = t146*t664
      t2228 = t308*t165
      t2238 = t653*t306
      t2244 = t178*t148
      t2247 = t161*t156
      t2250 = t881*t165
      t2262 = t950*t664
      t2268 = t340*t146*t208
      t2272 = t200*t164*t162
      t2277 = rho*(0.7420670230672731D-1*t2224*t972
     #+0.2473556743557577D-1*t147*(0.2697586091519874D1*(
     #-0.1269364121954074D0*t315-0.2173290372447189D-1*t2228*t293
     #-0.8693161489788758D-1*t173*t299)*t183-0.2697586091519874D1*t869
     #*t661-0.2697586091519874D1*t2238*t326+0.5395172183039748D1*t877
     #*t326*t661-0.2697586091519874D1*t307*(-0.1712109500228824D0
     #*t2244*t311-0.3424219000457648D0*t2247*t281-0.1172527576309517D0
     #*t2250*t293-0.2345055152619034D0*t308*t299))*t331
     #-0.2473556743557577D-1*t330*t2262+0.1D1*t334*t668
     #-0.130605864926168D2*t718*t2268+0.1364631356248995D2*t2272*t2162
     #*t208)
      s1 = 0.7420670230672731D-1*t2224*t965+0.2473556743557577D-1*t147
     #*(0.2697586091519874D1*(-0.1269364121954074D0*t473
     #-0.148092480894642D0*t476-0.2173290372447189D-1*t2228*t455
     #-0.8693161489788758D-1*t173*t458-0.1014202173808688D0*t173*t463)
     #*t183-0.2697586091519874D1*t872*t661-0.2697586091519874D1*t2238
     #*t485+0.5395172183039748D1*t877*t485*t661-0.2697586091519874D1
     #*t307*(-0.1712109500228824D0*t2244*t470-0.3424219000457648D0
     #*t2247*t440-0.3994922167200589D0*t656*t445-0.1172527576309517D0
     #*t2250*t455-0.2345055152619034D0*t308*t458-0.273589767805554D0
     #*t308*t463))*t331-0.2473556743557577D-1*t489*t2262+0.1D1*t507
     #*t668-0.261211729852336D1*t506*t147*t673+0.1D1*t510*t668
      t2349 = s1-0.261211729852336D1*t1386*t2268-0.2333333333333333D1
     #*t201*t445*t208+0.6094940363221173D1*t671*sigma*t738*t208+0.1D1
     #*t201*t529+0.1D1*t517*t152*(-0.522423459704672D1*t518*t437
     #+0.6965646129395627D1*t524)*t208-0.261211729852336D1*t671*sigma
     #*t982*t208
      t2350 = rho*t2349
      v2rhoasigmaaa(i) = -0.12407009817988D1*t4*t640*t35
     #+0.12407009817988D1*t214*t647-0.9305257363491D0*t5*(-t2191-t2193
     #+t2196-0.3266655613691336D-1*t234*t19*sigmaaa
     #+0.1791794125007405D-1*t27/t1480/rhoa*t19-0.4388082092272149D-1
     #*t21*t224)*t35+0.9305257363491D0*t241*t647+0.9305257363491D0
     #*t641*t248-0.18610514726982D1*t26*t1502*t247*t646
     #+0.9305257363491D0*t26*t245*(-t2191-t2193+t2196
     #-0.1155315866911171D-4*t2194)+t2277+t667+t670-t675+t2350
      t2352 = 0.4947113487115154D-1*t666
      t2353 = 0.2D1*t669
      t2354 = 0.522423459704672D1*t674
      t2355 = 2.D0*t2350
      v2rhoasigmaab(i) = 2.D0*t2277+t2352+t2353-t2354+t2355
      v2rhoasigmabb(i) = t2277+t667+t670-t675+t2350
      t2399 = rho*(0.7420670230672731D-1*t2224*t1710
     #+0.2473556743557577D-1*t147*(0.2697586091519874D1*(
     #-0.1269364121954074D0*t606-0.2173290372447189D-1*t2228*t593
     #-0.8693161489788758D-1*t173*t596)*t183-0.2697586091519874D1
     #*t1788*t661-0.2697586091519874D1*t2238*t613+0.5395172183039748D1
     #*t877*t613*t661-0.2697586091519874D1*t307*(-0.1712109500228824D0
     #*t2244*t603-0.3424219000457648D0*t2247*t583-0.1172527576309517D0
     #*t2250*t593-0.2345055152619034D0*t308*t596))*t331
     #-0.2473556743557577D-1*t617*t2262+0.1D1*t620*t668
     #-0.130605864926168D2*t1704*t2268+0.1364631356248995D2*t2272
     #*t1646*t580*t208)
      v2rhobsigmaaa(i) = t2399+t667+t670-t675+t2350
      v2rhobsigmaab(i) = 2.D0*t2399+t2352+t2353-t2354+t2355
      t2408 = 0.1680017400328676D-1*t678*t538*t44
      t2410 = 0.5040081825307963D-1*t544*t549
      t2411 = sigmabb*t554
      t2413 = 0.1680047024824737D-1*t2411*t1971
      v2rhobsigmabb(i) = -0.12407009817988D1*t38*t689*t69
     #+0.12407009817988D1*t534*t696-0.9305257363491D0*t39*(-t2408
     #-t2410+t2413-0.3266655613691336D-1*t554*t53*sigmabb
     #+0.1791794125007405D-1*t61/t2012/rhob*t53-0.4388082092272149D-1
     #*t55*t544)*t69+0.9305257363491D0*t561*t696+0.9305257363491D0
     #*t690*t568-0.18610514726982D1*t60*t1951*t567*t695
     #+0.9305257363491D0*t60*t565*(-t2408-t2410+t2413
     #-0.1155315866911171D-4*t2411)+t2399+t667+t670-t675+t2350
      t2445 = 0.6300065251232535D-2/t6/sigmaaa*t7*t10
      t2449 = 0.6300102281634954D-2/sigmaaa*t16*t229
      t2451 = 0.6300176343092764D-2*t31*t1474
      t2463 = t646**2
      t2478 = t661**2
      t2489 = t664**2
      t2502 = rho*(0.2473556743557577D-1*t147*(0.5862637881547585D-1
     #*t161*t165*t168*t183-0.5395172183039748D1*t2238*t661
     #+0.5395172183039748D1*t877*t2478-0.1581497040888031D0*t307*t178
     #*t165*t168)*t331-0.2473556743557577D-1*t147*t2489*t950
     #-0.522423459704672D1*t671*t340*t208+0.6823156781244976D1*t200
     #*t296*sigma*t1646*t208)
      v2sigmaaa2(i) = -0.9305257363491D0*t5*(-t2445+t2449-t2451
     #-0.6719227968777768D-2/t1480*t19*sigmaaa+0.8166639034228341D-2
     #*t635)*t35+0.18610514726982D1*t641*t647-0.18610514726982D1*t26
     #*t1502*t2463+0.9305257363491D0*t26*t245*(-t2445+t2449-t2451
     #+0.2166217250458446D-5*t31)+t2502
      v2sigmaaaab(i) = 2.D0*t2502
      v2sigmaaabb(i) = t2502
      v2sigmaab2(i) = 4.D0*v2sigmaaabb(i)
      v2sigmaabbb(i) = v2sigmaaaab(i)
      t2507 = 0.6300065251232535D-2/t40/sigmabb*t41*t44
      t2511 = 0.6300102281634954D-2/sigmabb*t50*t549
      t2513 = 0.6300176343092764D-2*t65*t1971
      t2525 = t695**2
      v2sigmabb2(i) = -0.9305257363491D0*t39*(-t2507+t2511-t2513
     #-0.6719227968777768D-2/t2012*t53*sigmabb+0.8166639034228341D-2
     #*t684)*t69+0.18610514726982D1*t690*t696-0.18610514726982D1*t60
     #*t1951*t2525+0.9305257363491D0*t60*t565*(-t2507+t2511-t2513
     #+0.2166217250458446D-5*t65)+v2sigmaaabb(i)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2rhoasigmaab(i) = 0.0d0
      v2rhoasigmabb(i) = 0.0d0
      v2rhobsigmaaa(i) = 0.0d0
      v2rhobsigmaab(i) = 0.0d0
      v2rhobsigmabb(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      v2sigmaab2(i) = 0.0d0
      v2sigmabb2(i) = 0.0d0
      v2sigmaaaab(i) = 0.0d0
      v2sigmaaabb(i) = 0.0d0
      v2sigmaabbb(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_xc_pw91
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson,
c     M.R. Pederson, D.J. Singh, C. Fiolhais
c     Atoms, molecules, solids and surfaces:
c     Applications of the generalized gradient approximation
c     for exchange and correlation
c     Phys. Rev. B 46 (1992) 6671--6687
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt,i
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(tol,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = dsqrt(sigma)
      t6 = t4/t3
      t8 = dlog(0.1259928455434599D1*t6+dsqrt(1+0.1587419712813814D1
     #*t6**2))
      t10 = 0.3175033930295641D-1*t6*t8
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t17 = dexp(-0.261211729852336D1*t14*sigma)
      t25 = sigma**2
      t26 = t11**2
      t36 = 1/rho
      t37 = t36**(1.D0/3.D0)
      t40 = t36**(1.D0/6.D0)
      t43 = dsqrt(t36)
      t45 = t37**2
      t51 = dlog(1.D0+0.160818243221511D2/(0.598255043577108D1*t40
     #+0.2225569421150687D1*t37+0.8004286349993634D0*t43
     #+0.1897004325747559D0*t45))
      t52 = (1.D0+0.1325688999052018D0*t37)*t51
      t55 = 1/t2/t11
      t59 = dexp(0.2513869963240348D1*t52)
      t60 = t59-1.D0
      t61 = 1/t60
      t64 = 1/t12/t26
      t71 = t60**2
      t81 = dlog(1.D0+0.2697586091519874D1*(0.634682060977037D-1*sigma
     #*t55+0.1086645186223595D-1*t61*t25*t64)/(1.D0
     #+0.1712109500228824D0*t61*sigma*t55+0.2931318940773793D-1/t71
     #*t25*t64))
      zk(i) = -0.7385587663820224D0*t3*(1.D0+t10+0.261211729852336D-1*
     #(0.2743D0-0.1508D0*t17)*sigma*t14)/(1.D0+t10
     #+0.272926271249799D-5*t25/t2/t26/rho)+rho*(-0.62182D-1*t52
     #+0.2473556743557577D-1*t81+0.1D1*(0.1D-2*(0.2568D1
     #+0.1443307452126544D2*t37+0.2843543831490386D-2*t45)/(1.D0
     #+0.5411317332115466D1*t37+0.1816419932959077D0*t45
     #+0.1763993811759022D-1*t36)-0.1853571428571429D-2)*sigma*t55*t17)
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(tol,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = dsqrt(sigma)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1259928455434599D1*t6+dsqrt(1+0.1587419712813814D1
     #*t6**2))
      t10 = 0.3175033930295641D-1*t6*t8
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = t14*sigma
      t17 = dexp(-0.261211729852336D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigma
      t23 = 1.D0+t10+0.261211729852336D-1*t20*t14
      t24 = t3*t23
      t25 = sigma**2
      t26 = t11**2
      t27 = t26*rho
      t29 = 1/t2/t27
      t32 = 1.D0+t10+0.272926271249799D-5*t25*t29
      t33 = 1/t32
      t36 = 1/rho
      t37 = t36**(1.D0/3.D0)
      t39 = 1.D0+0.1325688999052018D0*t37
      t40 = t36**(1.D0/6.D0)
      t43 = dsqrt(t36)
      t45 = t37**2
      t47 = 0.598255043577108D1*t40+0.2225569421150687D1*t37
     #+0.8004286349993634D0*t43+0.1897004325747559D0*t45
      t50 = 1.D0+0.160818243221511D2/t47
      t51 = dlog(t50)
      t52 = t39*t51
      t53 = 0.62182D-1*t52
      t55 = 1/t2/t11
      t59 = dexp(0.2513869963240348D1*t52)
      t60 = t59-1.D0
      t61 = 1/t60
      t62 = t61*t25
      t64 = 1/t12/t26
      t67 = 0.634682060977037D-1*sigma*t55+0.1086645186223595D-1*t62*t64
      t68 = t61*sigma
      t71 = t60**2
      t72 = 1/t71
      t73 = t72*t25
      t76 = 1.D0+0.1712109500228824D0*t68*t55+0.2931318940773793D-1
     #*t73*t64
      t77 = 1/t76
      t80 = 1.D0+0.2697586091519874D1*t67*t77
      t81 = dlog(t80)
      t82 = 0.2473556743557577D-1*t81
      t85 = 0.2568D1+0.1443307452126544D2*t37+0.2843543831490386D-2*t45
      t89 = 1.D0+0.5411317332115466D1*t37+0.1816419932959077D0*t45
     #+0.1763993811759022D-1*t36
      t90 = 1/t89
      t93 = 0.1D-2*t85*t90-0.1853571428571429D-2
      t94 = t93*sigma
      t95 = t55*t17
      t97 = 0.1D1*t94*t95
      zk(i) = -0.7385587663820224D0*t24*t33+rho*(-t53+t82+t97)
      t105 = 0.8466757147455043D-1*t4*t55*t8
      t106 = t11*rho
      t108 = 1/t12/t106
      t112 = dsqrt(1.D0+0.1587419712813815D1*t15)
      t113 = 1/t112
      t115 = 0.1066750825533289D0*sigma*t108*t113
      t116 = t26*t11
      t119 = t25/t2/t116
      t128 = t32**2
      t129 = 1/t128
      t136 = 1/t11
      t137 = 1/t45*t136
      t138 = t137*t51
      t140 = t47**2
      t143 = t40**2
      t144 = t143**2
      t154 = 1/t37*t136
      t159 = t39/t140*(-0.99709173929518D0/t144/t40*t136
     #-0.7418564737168958D0*t137-0.4002143174996817D0/t43*t136
     #-0.1264669550498372D0*t154)/t50
      t162 = 1/t2/t106
      t167 = -0.1110869918438343D0*t138-0.4042761511756372D2*t159
      t169 = t64*t167*t59
      t173 = 1/t12/t27
      t179 = t76**2
      t181 = t67/t179
      t182 = t72*sigma
      t200 = 1/t80
      t208 = t89**2
      s1 = -0.9847450218426965D0*t2*t23*t33-0.3692793831910112D0*t3*(
     #-t105-t115-0.5487637560595959D-1*t119*t17-0.1393129225879125D0
     #*t20*t108)*t33+0.3692793831910112D0*t24*t129*(-t105-t115
     #-0.2911213559997856D-4*t119)
      vrhoa(i) = s1-t53+t82+t97+rho*(0.2747799777968419D-2*t138+0.1D1
     #*t159+0.2473556743557577D-1*(0.2697586091519874D1*(
     #-0.148092480894642D0*sigma*t162-0.1086645186223595D-1*t73*t169
     #-0.5071010869043442D-1*t62*t173)*t77-0.2697586091519874D1*t181*(
     #-0.1712109500228824D0*t182*t55*t167*t59-0.3994922167200589D0*t68
     #*t162-0.5862637881547585D-1/t71/t60*t25*t169-0.136794883902777D0
     #*t73*t173))*t200+0.1D1*(0.1D-2*(-0.4811024840421814D1*t137
     #-0.1895695887660258D-2*t154)*t90-0.1D-2*t85/t208*(
     #-0.1803772444038489D1*t137-0.1210946621972718D0*t154
     #-0.1763993811759022D-1*t136))*sigma*t95-0.2333333333333333D1*t94
     #*t162*t17+0.6965646129395627D1*t93*t25/t116*t17)
      t234 = 0.6350067860591282D-1/t4*t5*t8
      t236 = 0.8000631191499664D-1*t14*t113
      vsigmaaa(i) = -0.7385587663820224D0*t3*(t234+t236
     #+0.411572817044697D-1*t29*t17*sigma+0.1044846919409344D0*t19*t14
     #)*t33+0.7385587663820224D0*t24*t129*(t234+t236
     #+0.2183410169998392D-4*sigma*t29)+4.D0*rho*
     #(0.2473556743557577D-1*(0.2697586091519874D1*
     #(0.634682060977037D-1*t55+0.2173290372447189D-1*t68*t64)*t77
     #-0.2697586091519874D1*t181*(0.1712109500228824D0*t61*t55
     #+0.5862637881547585D-1*t182*t64))*t200+0.1D1*t93*t55*t17
     #-0.261211729852336D1*t94/t27*t17)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(tol,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = dsqrt(sigma)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1259928455434599D1*t6+dsqrt(1+0.1587419712813814D1
     #*t6**2))
      t10 = 0.3175033930295641D-1*t6*t8
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = t14*sigma
      t17 = dexp(-0.261211729852336D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigma
      t23 = 1.D0+t10+0.261211729852336D-1*t20*t14
      t24 = t3*t23
      t25 = sigma**2
      t26 = t11**2
      t27 = t26*rho
      t29 = 1/t2/t27
      t32 = 1.D0+t10+0.272926271249799D-5*t25*t29
      t33 = 1/t32
      t36 = 1/rho
      t37 = t36**(1.D0/3.D0)
      t39 = 1.D0+0.1325688999052018D0*t37
      t40 = t36**(1.D0/6.D0)
      t43 = dsqrt(t36)
      t45 = t37**2
      t47 = 0.598255043577108D1*t40+0.2225569421150687D1*t37
     #+0.8004286349993634D0*t43+0.1897004325747559D0*t45
      t50 = 1.D0+0.160818243221511D2/t47
      t51 = dlog(t50)
      t52 = t39*t51
      t53 = 0.62182D-1*t52
      t55 = 1/t2/t11
      t59 = dexp(0.2513869963240348D1*t52)
      t60 = t59-1.D0
      t61 = 1/t60
      t62 = t61*t25
      t64 = 1/t12/t26
      t67 = 0.634682060977037D-1*sigma*t55+0.1086645186223595D-1*t62*t64
      t68 = t61*sigma
      t71 = t60**2
      t72 = 1/t71
      t73 = t72*t25
      t76 = 1.D0+0.1712109500228824D0*t68*t55+0.2931318940773793D-1
     #*t73*t64
      t77 = 1/t76
      t80 = 1.D0+0.2697586091519874D1*t67*t77
      t81 = dlog(t80)
      t82 = 0.2473556743557577D-1*t81
      t85 = 0.2568D1+0.1443307452126544D2*t37+0.2843543831490386D-2*t45
      t89 = 1.D0+0.5411317332115466D1*t37+0.1816419932959077D0*t45
     #+0.1763993811759022D-1*t36
      t90 = 1/t89
      t93 = 0.1D-2*t85*t90-0.1853571428571429D-2
      t94 = t93*sigma
      t95 = t55*t17
      t97 = 0.1D1*t94*t95
      zk(i) = -0.7385587663820224D0*t24*t33+rho*(-t53+t82+t97)
      t100 = t2*t23
      t105 = 0.8466757147455043D-1*t4*t55*t8
      t106 = t11*rho
      t108 = 1/t12/t106
      t111 = 1.D0+0.1587419712813815D1*t15
      t112 = dsqrt(t111)
      t113 = 1/t112
      t115 = 0.1066750825533289D0*sigma*t108*t113
      t116 = t26*t11
      t118 = 1/t2/t116
      t119 = t25*t118
      t124 = -t105-t115-0.5487637560595959D-1*t119*t17
     #-0.1393129225879125D0*t20*t108
      t125 = t3*t124
      t128 = t32**2
      t129 = 1/t128
      t131 = -t105-t115-0.2911213559997856D-4*t119
      t132 = t129*t131
      t135 = 1/t45
      t136 = 1/t11
      t137 = t135*t136
      t138 = t137*t51
      t140 = t47**2
      t141 = 1/t140
      t142 = t39*t141
      t143 = t40**2
      t144 = t143**2
      t145 = t144*t40
      t146 = 1/t145
      t150 = 1/t43
      t153 = 1/t37
      t154 = t153*t136
      t156 = -0.99709173929518D0*t146*t136-0.7418564737168958D0*t137
     #-0.4002143174996817D0*t150*t136-0.1264669550498372D0*t154
      t157 = 1/t50
      t159 = t142*t156*t157
      t162 = 1/t2/t106
      t167 = -0.1110869918438343D0*t138-0.4042761511756372D2*t159
      t169 = t64*t167*t59
      t173 = 1/t12/t27
      t176 = -0.148092480894642D0*sigma*t162-0.1086645186223595D-1*t73
     #*t169-0.5071010869043442D-1*t62*t173
      t179 = t76**2
      t180 = 1/t179
      t181 = t67*t180
      t182 = t72*sigma
      t190 = 1/t71/t60
      t191 = t190*t25
      t196 = -0.1712109500228824D0*t182*t55*t167*t59
     #-0.3994922167200589D0*t68*t162-0.5862637881547585D-1*t191*t169
     #-0.136794883902777D0*t73*t173
      t199 = 0.2697586091519874D1*t176*t77-0.2697586091519874D1*t181
     #*t196
      t200 = 1/t80
      t201 = t199*t200
      t205 = -0.4811024840421814D1*t137-0.1895695887660258D-2*t154
      t208 = t89**2
      t209 = 1/t208
      t210 = t85*t209
      t214 = -0.1803772444038489D1*t137-0.1210946621972718D0*t154
     #-0.1763993811759022D-1*t136
      t217 = 0.1D-2*t205*t90-0.1D-2*t210*t214
      t218 = t217*sigma
      t219 = t218*t95
      t221 = t162*t17
      t222 = t94*t221
      t224 = t93*t25
      t226 = 1/t116*t17
      t227 = t224*t226
      vrhoa(i) = -0.9847450218426965D0*t100*t33-0.3692793831910112D0
     #*t125*t33+0.3692793831910112D0*t24*t132-t53+t82+t97+rho*
     #(0.2747799777968419D-2*t138+0.1D1*t159+0.2473556743557577D-1
     #*t201+0.1D1*t219-0.2333333333333333D1*t222+0.6965646129395627D1
     #*t227)
      t231 = 1/t4
      t234 = 0.6350067860591282D-1*t231*t5*t8
      t236 = 0.8000631191499664D-1*t14*t113
      t237 = t29*t17
      t242 = t234+t236+0.411572817044697D-1*t237*sigma
     #+0.1044846919409344D0*t19*t14
      t243 = t3*t242
      t248 = t234+t236+0.2183410169998392D-4*sigma*t29
      t249 = t129*t248
      t255 = 0.634682060977037D-1*t55+0.2173290372447189D-1*t68*t64
      t262 = 0.1712109500228824D0*t61*t55+0.5862637881547585D-1*t182*t64
      t265 = 0.2697586091519874D1*t255*t77-0.2697586091519874D1*t181
     #*t262
      t266 = t265*t200
      t269 = t93*t55*t17
      t271 = 1/t27
      t272 = t271*t17
      t273 = t94*t272
      vsigmaaa(i) = -0.7385587663820224D0*t243*t33
     #+0.7385587663820224D0*t24*t249+4.D0*rho*(0.2473556743557577D-1
     #*t266+0.1D1*t269-0.261211729852336D1*t273)
      t284 = t26*t106
      t287 = t224/t284*t17
      t290 = 1/t2/t26
      t293 = t93*t290*sigma*t17
      t297 = 1/t26
      t298 = 1/t45/t36*t297
      t300 = 1/t106
      t301 = t135*t300
      t305 = 1/t37/t36*t297
      t307 = t153*t300
      t318 = t214**2
      t333 = t199**2
      t334 = t80**2
      t335 = 1/t334
      t356 = t142*(-0.8309097827459833D0/t145/t36*t297
     #+0.199418347859036D1*t146*t300-0.4945709824779306D0*t298
     #+0.1483712947433792D1*t301-0.2001071587498409D0/t43/t36*t297
     #+0.8004286349993634D0*t150*t300-0.4215565168327908D-1*t305
     #+0.2529339100996745D0*t307)*t157
      t358 = t140**2
      t361 = t156**2
      t362 = t50**2
      t365 = t39/t358*t361/t362
      t367 = t301*t51
      t371 = t137*t141*t156*t157
      t373 = t298*t51
      t379 = t39/t140/t47*t361*t157
      t386 = t25*sigma
      t388 = t26**2
      t389 = t388*rho
      t395 = sigma*t290
      t398 = t173*t167*t59
      t407 = -0.4042761511756372D2*t356-0.6501498040842D3*t365
     #+0.2221739836876686D0*t367+0.3572962974617552D1*t371
     #-0.740579945625562D-1*t373+0.8085523023512745D2*t379
      t409 = t64*t407*t59
      t412 = t167**2
      t413 = t64*t412
      t415 = dexp(0.5027739926480695D1*t52)
      t416 = t413*t415
      t419 = t413*t59
      t423 = 1/t12/t116
      t424 = t62*t423
      t429 = t176*t180
      t434 = t67/t179/t76
      t435 = t196**2
      t448 = t71**2
      t453 = t190*sigma
      t454 = t55*t412
      t461 = t68*t290
      t463 = t73*t423
      s1 = -0.5804705107829689D2*t287+0.7777777777777778D1*t293+0.1D1*
     #(0.1D-2*(-0.3207349893614542D1*t298+0.9622049680843627D1*t301
     #-0.6318986292200858D-3*t305+0.3791391775320515D-2*t307)*t90
     #-0.2D-2*t205*t209*t214+0.2D-2*t85/t208/t89*t318-0.1D-2*t210*(
     #-0.1202514962692326D1*t298+0.3607544888076978D1*t301
     #-0.4036488739909061D-1*t305+0.2421893243945437D0*t307
     #+0.3527987623518044D-1*t300))*sigma*t95-0.2473556743557577D-1
     #*t333*t335+0.1D1*t356+0.160818243221511D2*t365
     #-0.5495599555936838D-2*t367
      s2 = s1-0.8837926660346786D-1*t371+0.1831866518645613D-2*t373
     #-0.2D1*t379
      s3 = s2+0.1393129225879125D2*t217*t25*t226
      t475 = s3-0.4666666666666667D1*t218*t221+0.4852022599996427D2
     #*t93*t386/t12/t389*t17+0.2473556743557577D-1*
     #(0.2697586091519874D1*(0.4936416029821399D0*t395
     #+0.1014202173808688D0*t73*t398-0.1086645186223595D-1*t73*t409
     #+0.2173290372447189D-1*t191*t416-0.1086645186223595D-1*t73*t419
     #+0.2873572825791284D0*t424)*t77-0.5395172183039748D1*t429*t196
     #+0.5395172183039748D1*t434*t435-0.2697586091519874D1*t181*
     #(0.547179535611108D0*t191*t398-0.1712109500228824D0*t182*t55
     #*t407*t59+0.7989844334401178D0*t182*t162*t167*t59
     #+0.1758791364464276D0/t448*t25*t416+0.3424219000457648D0*t453
     #*t454*t415-0.1712109500228824D0*t182*t454*t59
     #+0.1331640722400196D1*t461+0.7751710087824029D0*t463
     #-0.5862637881547585D-1*t191*t409-0.5862637881547585D-1*t191*t419
     #))*t200
      t488 = dlog(1.D0+0.2960857464321668D2/(0.8157414703487641D1*t40
     #+0.2247591863577616D1*t37+0.4300972471276643D0*t43
     #+0.1911512595127338D0*t45))
      t490 = (1.D0+0.6901399211255825D-1*t37)*t488*t136
      t492 = t81*t136
      t498 = -0.1365402272980597D1*t490+0.1675913308826898D1*t52*t136
      t500 = t64*t498*t59
      t528 = 1/t128/t32
      t529 = t131**2
      t535 = 0.395115333547902D0*t4*t162*t8
      t538 = 0.1066750825533289D1*sigma*t64*t113
      t541 = t25/t2/t284
      t543 = 1/t112/t111
      t545 = 0.4515683437631874D0*t541*t543
      t572 = -t498
      t574 = t64*t572*t59
      s1 = 0.4D1*t159+0.1099119911187368D-1*t138-0.6564966812284644D0
     #/t12*t23*t33+2.D0*rho*t475+rho*(0.33774D-1*t490
     #-0.1649037829038385D-1*t492+0.2473556743557577D-1*
     #(0.2697586091519874D1*(0.2820809159897942D-1*t395
     #-0.1086645186223595D-1*t73*t500+0.9659068321987509D-2*t424)*t77
     #-0.2697586091519874D1*t181*(-0.1712109500228824D0*t182*t55*t498
     #*t59+0.760937555657255D-1*t461-0.5862637881547585D-1*t191*t500
     #+0.2605616836243371D-1*t463))*t200-0.2222222222222222D0*t293
     #+0.1160941021565938D1*t287)+0.7385587663820224D0*t125*t132
     #-0.7385587663820224D0*t24*t528*t529+0.3692793831910112D0*t24
     #*t129*(t535+t538-t545+0.3687537175997285D-3*t541)
      v2rhoa2(i) = s1-0.3692793831910112D0*t3*(t535+t538-t545
     #+0.9877747609072727D0*t541*t17-0.764498826669826D0*t386/t388/t11
     #*t17+0.1021628098978025D1*t20*t64)*t33-0.1969490043685393D1*t2
     #*t124*t33+0.1969490043685393D1*t100*t132+rho*(-0.33774D-1*t490
     #+0.1649037829038385D-1*t492+0.2473556743557577D-1*
     #(0.2697586091519874D1*(-0.2820809159897942D-1*t395
     #-0.1086645186223595D-1*t73*t574-0.9659068321987509D-2*t424)*t77
     #-0.2697586091519874D1*t181*(-0.1712109500228824D0*t182*t55*t572
     #*t59-0.760937555657255D-1*t461-0.5862637881547585D-1*t191*t574
     #-0.2605616836243371D-1*t463))*t200+0.2222222222222222D0*t293
     #-0.1160941021565938D1*t287)+0.9894226974230308D-1*t201+0.4D1
     #*t219-0.9333333333333333D1*t222+0.2786258451758251D2*t227
      t610 = 0.1693351429491009D0*t231*t55*t8
      t612 = 0.6400504953199731D0*t108*t113
      t613 = sigma*t118
      t615 = 0.3386762578223905D0*t613*t543
      t655 = t255*t180
      s1 = -0.9847450218426965D0*t2*t242*t33+0.9847450218426965D0*t100
     #*t249-0.3692793831910112D0*t3*(-t610-t612+t615
     #-0.6585165072715151D0*t118*t17*sigma+0.5733741200023695D0*t25
     #/t389*t17-0.5572516903516501D0*t19*t108)*t33
     #+0.3692793831910112D0*t125*t249+0.3692793831910112D0*t243*t132
      s2 = s1-0.7385587663820224D0*t24*t528*t131*t248
     #+0.3692793831910112D0*t24*t129*(-t610-t612+t615
     #-0.2328970847998285D-3*t613)
      v2rhoasigmaaa(i) = s2+0.9894226974230308D-1*t266+0.4D1*t269
     #-0.1044846919409344D2*t273+4.D0*rho*(0.2473556743557577D-1*
     #(0.2697586091519874D1*(-0.148092480894642D0*t162
     #-0.2173290372447189D-1*t182*t169-0.1014202173808688D0*t68*t173)
     #*t77-0.2697586091519874D1*t429*t262-0.2697586091519874D1*t655
     #*t196+0.5395172183039748D1*t434*t196*t262-0.2697586091519874D1
     #*t181*(-0.1712109500228824D0*t72*t55*t167*t59
     #-0.3994922167200589D0*t61*t162-0.1172527576309517D0*t453*t169
     #-0.273589767805554D0*t182*t173))*t200-0.2473556743557577D-1*t199
     #*t335*t265+0.1D1*t217*t55*t17-0.261211729852336D1*t218*t272
     #-0.2333333333333333D1*t93*t162*t17+0.2002623262201243D2*t94*t226
     #-0.181950847499866D2*t224/t12/t388*t17)
      t702 = 0.1270013572118256D0/t4/sigma*t5*t8
      t706 = 0.1600126238299933D0/sigma*t14*t113
      t708 = 0.2540071933667929D0*t29*t543
      t720 = t248**2
      t734 = t262**2
      t743 = t265**2
      v2sigmaaa2(i) = -0.7385587663820224D0*t3*(-t702+t706-t708
     #-0.4300305900017772D0/t388*t17*sigma+0.3292582536357576D0*t237)
     #*t33+0.1477117532764045D1*t243*t249-0.1477117532764045D1*t24
     #*t528*t720+0.7385587663820224D0*t24*t129*(-t702+t706-t708
     #+0.8733640679993569D-4*t29)+16.D0*rho*(0.2473556743557577D-1*
     #(0.5862637881547585D-1*t61*t64*t77-0.5395172183039748D1*t655
     #*t262+0.5395172183039748D1*t434*t734-0.1581497040888031D0*t181
     #*t72*t64)*t200-0.2473556743557577D-1*t743*t335
     #-0.522423459704672D1*t93*t271*t17+0.6823156781244976D1*t94/t12
     #/t284*t17)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end

c:XC_PW91subrend
c:XC_B3LYPsubrstart

c    Generated: Thu Feb 13 08:52:41 GMT 2003

      subroutine uks_xc_b3lyp
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     P.J. Stephens, F.J. Devlin, C.F. Chabalowski, M.J. Frisch
c     Ab initio calculation of vibrational absorption and circular 
c     dichroism spectra using density functional force fields
c     J. Phys. Chem. 98 (1994) 11623-11627
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt,i
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(tol,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t5 = 1/t3
      t7 = dsqrt(sigmabb)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      t16 = 1/rhob
      t17 = t16**(1.D0/3.D0)
      t19 = t16**(1.D0/6.D0)
      t22 = 1/(0.6203504908994D0*t17+0.1584942278842832D2*t19
     #+0.101578D3)
      t25 = dlog(0.6203504908994D0*t17*t22)
      t31 = datan(0.1171685277708993D1/(0.1575246635799487D1*t19
     #+0.201231D2))
      t35 = (0.7876233178997433D0*t19+0.743294D0)**2
      t37 = dlog(t35*t22)
      zk(i) = -0.74442058907928D0*t3-0.3024D-2*t5*sigmabb/(1.D0
     #+0.252D-1*t8*t9)+0.19D0*rhob*(0.1554535D-1*t25
     #+0.6188180297906063D0*t31+0.2667310007273315D-2*t37)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(tol,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t5 = 1/t3
      t7 = dsqrt(sigmaaa)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      t16 = 1/rhoa
      t17 = t16**(1.D0/3.D0)
      t19 = t16**(1.D0/6.D0)
      t22 = 1/(0.6203504908994D0*t17+0.1584942278842832D2*t19
     #+0.101578D3)
      t25 = dlog(0.6203504908994D0*t17*t22)
      t31 = datan(0.1171685277708993D1/(0.1575246635799487D1*t19
     #+0.201231D2))
      t35 = (0.7876233178997433D0*t19+0.743294D0)**2
      t37 = dlog(t35*t22)
      zk(i) = -0.74442058907928D0*t3-0.3024D-2*t5*sigmaaa/(1.D0
     #+0.252D-1*t8*t9)+0.19D0*rhoa*(0.1554535D-1*t25
     #+0.6188180297906063D0*t31+0.2667310007273315D-2*t37)
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(tol,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(tol,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t7 = 1/t5
      t9 = dsqrt(sigmaaa)
      t10 = t9*t7
      t11 = dlog(t10+dsqrt(1+t10**2))
      t18 = rhob**(1.D0/3.D0)
      t19 = t18*rhob
      t21 = 1/t19
      t23 = dsqrt(sigmabb)
      t24 = t23*t21
      t25 = dlog(t24+dsqrt(1+t24**2))
      t32 = rho**(1.D0/3.D0)
      t33 = 1/t32
      t36 = 1/(1.D0+0.349D0*t33)
      t38 = 1/rho
      t39 = rhob*t38
      t42 = 0.2533D0*t33
      t43 = dexp(-t42)
      t45 = rho**2
      t47 = t32**2
      t51 = rhoa**2
      t52 = t4**2
      t55 = rhob**2
      t56 = t18**2
      t60 = t33*t36
      t82 = 0.6666666666666667D0*t45
      t93 = t38**(1.D0/3.D0)
      t94 = 0.6203504908994D0*t93
      t95 = t38**(1.D0/6.D0)
      t98 = 1/(t94+0.1029581201158544D2*t95+0.427198D2)
      t101 = dlog(0.6203504908994D0*t93*t98)
      t103 = 0.1575246635799487D1*t95
      t107 = datan(0.448998886412873D-1/(t103+0.13072D2))
      t109 = 0.7876233178997433D0*t95
      t111 = (t109+0.409286D0)**2
      t113 = dlog(t111*t98)
      t117 = 1/(t94+0.1584942278842832D2*t95+0.101578D3)
      t120 = dlog(0.6203504908994D0*t93*t117)
      t125 = datan(0.1171685277708993D1/(t103+0.201231D2))
      t128 = (t109+0.743294D0)**2
      t130 = dlog(t128*t117)
      t138 = (rhoa-1.D0*rhob)*t38
      t139 = 1.D0+t138
      t140 = t139**(1.D0/3.D0)
      t143 = 1.D0-1.D0*t138
      t144 = t143**(1.D0/3.D0)
      s1 = -0.74442058907928D0*t5-0.3024D-2*t7*sigmaaa/(1.D0+0.252D-1
     #*t10*t11)-0.74442058907928D0*t19
      s2 = s1-0.3024D-2*t21*sigmabb/(1.D0+0.252D-1*t24*t25)
      zk(i) = s2-0.1593432D0*t36*rhoa*t39-0.52583256D-2*t43*t36/t47
     #/t45/rho*(rhoa*rhob*(0.3646239897876478D2*t52*t51
     #+0.3646239897876478D2*t56*t55+(0.2611111111111111D1
     #-0.9850555555555556D-1*t33-0.1357222222222222D0*t60)*sigma-1.D0*
     #(0.25D1-0.1407222222222222D-1*t33-0.1938888888888889D-1*t60)*
     #(sigmaaa+sigmabb)-0.1111111111111111D0*(t42+0.349D0*t60-11.D0)*
     #(rhoa*t38*sigmaaa+t39*sigmabb))-0.6666666666666667D0*t45*sigma+
     #(t82-1.D0*t51)*sigmabb+(t82-1.D0*t55)*sigmaaa)+0.19D0*rho*
     #(0.310907D-1*t101+0.205219729378375D2*t107+0.4431373767749538D-2
     #*t113+0.1923661050931536D1*(0.1554535D-1*t120
     #+0.6188180297906063D0*t125+0.2667310007273315D-2*t130
     #-0.310907D-1*t101-0.205219729378375D2*t107-0.4431373767749538D-2
     #*t113)*(t140*t139+t144*t143-2.D0))
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(tol,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t5 = 1/t3
      t6 = t5*sigmabb
      t7 = dsqrt(sigmabb)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      t12 = 1.D0+0.252D-1*t8*t9
      t13 = 1/t12
      t16 = 1/rhob
      t17 = t16**(1.D0/3.D0)
      t19 = t16**(1.D0/6.D0)
      t21 = 0.6203504908994D0*t17+0.1584942278842832D2*t19+0.101578D3
      t22 = 1/t21
      t25 = dlog(0.6203504908994D0*t17*t22)
      t28 = 0.1575246635799487D1*t19+0.201231D2
      t31 = datan(0.1171685277708993D1/t28)
      t34 = 0.7876233178997433D0*t19+0.743294D0
      t35 = t34**2
      t37 = dlog(t35*t22)
      zk(i) = -0.74442058907928D0*t3-0.3024D-2*t6*t13+0.19D0*rhob*
     #(0.1554535D-1*t25+0.6188180297906063D0*t31+0.2667310007273315D-2
     #*t37)
      vrhoa(i) = 0.D0
      t43 = rhob**2
      t45 = 1/t2/t43
      t49 = t12**2
      t50 = 1/t49
      t55 = t2**2
      t60 = 1/t55/t43
      t63 = dsqrt(1.D0+sigmabb*t60)
      t64 = 1/t63
      t74 = t17**2
      t75 = 1/t74
      t77 = 1/t43
      t80 = t21**2
      t81 = 1/t80
      t85 = t19**2
      t86 = t85**2
      t88 = 1/t86/t19
      t89 = t88*t77
      t91 = -0.2067834969664667D0*t75*t77-0.2641570464738054D1*t89
      t99 = t28**2
      t100 = 1/t99
      vrhob(i) = -0.99256078543904D0*t2+0.4032D-2*t45*sigmabb*t13
     #+0.3024D-2*t6*t50*(-0.336D-1*t7*t45*t9-0.336D-1*sigmabb/t55/t43
     #/rhob*t64)+0.29536165D-2*t25+0.1175754256602152D0*t31
     #+0.5067889013819299D-3*t37+0.19D0*rhob*(0.2505897912236993D-1*(
     #-0.2067834969664667D0*t75*t22*t77-0.6203504908994D0*t17*t81*t91)
     #/t17*t21+0.1903580477513215D0*t100*t88*t77/(1.D0+0.137284639D1
     #*t100)+0.2667310007273315D-2*(-0.2625411059665811D0*t34*t22*t89
     #-1.D0*t35*t81*t91)/t35*t21)
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = -0.3024D-2*t5*t13+0.3024D-2*t6*t50*(0.126D-1/t7*t5
     #*t9+0.126D-1*t60*t64)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(tol,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t5 = 1/t3
      t6 = t5*sigmaaa
      t7 = dsqrt(sigmaaa)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      t12 = 1.D0+0.252D-1*t8*t9
      t13 = 1/t12
      t16 = 1/rhoa
      t17 = t16**(1.D0/3.D0)
      t19 = t16**(1.D0/6.D0)
      t21 = 0.6203504908994D0*t17+0.1584942278842832D2*t19+0.101578D3
      t22 = 1/t21
      t25 = dlog(0.6203504908994D0*t17*t22)
      t28 = 0.1575246635799487D1*t19+0.201231D2
      t31 = datan(0.1171685277708993D1/t28)
      t34 = 0.7876233178997433D0*t19+0.743294D0
      t35 = t34**2
      t37 = dlog(t35*t22)
      zk(i) = -0.74442058907928D0*t3-0.3024D-2*t6*t13+0.19D0*rhoa*
     #(0.1554535D-1*t25+0.6188180297906063D0*t31+0.2667310007273315D-2
     #*t37)
      t43 = rhoa**2
      t45 = 1/t2/t43
      t49 = t12**2
      t50 = 1/t49
      t55 = t2**2
      t60 = 1/t55/t43
      t63 = dsqrt(1.D0+sigmaaa*t60)
      t64 = 1/t63
      t74 = t17**2
      t75 = 1/t74
      t77 = 1/t43
      t80 = t21**2
      t81 = 1/t80
      t85 = t19**2
      t86 = t85**2
      t88 = 1/t86/t19
      t89 = t88*t77
      t91 = -0.2067834969664667D0*t75*t77-0.2641570464738054D1*t89
      t99 = t28**2
      t100 = 1/t99
      vrhoa(i) = -0.99256078543904D0*t2+0.4032D-2*t45*sigmaaa*t13
     #+0.3024D-2*t6*t50*(-0.336D-1*t7*t45*t9-0.336D-1*sigmaaa/t55/t43
     #/rhoa*t64)+0.29536165D-2*t25+0.1175754256602152D0*t31
     #+0.5067889013819299D-3*t37+0.19D0*rhoa*(0.2505897912236993D-1*(
     #-0.2067834969664667D0*t75*t22*t77-0.6203504908994D0*t17*t81*t91)
     #/t17*t21+0.1903580477513215D0*t100*t88*t77/(1.D0+0.137284639D1
     #*t100)+0.2667310007273315D-2*(-0.2625411059665811D0*t34*t22*t89
     #-1.D0*t35*t81*t91)/t35*t21)
      vrhob(i) = 0.D0
      vsigmaaa(i) = -0.3024D-2*t5*t13+0.3024D-2*t6*t50*(0.126D-1/t7*t5
     #*t9+0.126D-1*t60*t64)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(tol,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(tol,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t7 = 1/t5
      t8 = t7*sigmaaa
      t9 = dsqrt(sigmaaa)
      t10 = t9*t7
      t11 = dlog(t10+dsqrt(1+t10**2))
      t14 = 1.D0+0.252D-1*t10*t11
      t15 = 1/t14
      t18 = rhob**(1.D0/3.D0)
      t19 = t18*rhob
      t21 = 1/t19
      t22 = t21*sigmabb
      t23 = dsqrt(sigmabb)
      t24 = t23*t21
      t25 = dlog(t24+dsqrt(1+t24**2))
      t28 = 1.D0+0.252D-1*t24*t25
      t29 = 1/t28
      t32 = rho**(1.D0/3.D0)
      t33 = 1/t32
      t35 = 1.D0+0.349D0*t33
      t36 = 1/t35
      t37 = t36*rhoa
      t38 = 1/rho
      t39 = rhob*t38
      t42 = 0.2533D0*t33
      t43 = dexp(-t42)
      t44 = t43*t36
      t45 = rho**2
      t47 = t32**2
      t49 = 1/t47/t45/rho
      t50 = rhoa*rhob
      t51 = rhoa**2
      t52 = t4**2
      t53 = t52*t51
      t55 = rhob**2
      t56 = t18**2
      t57 = t56*t55
      t60 = t33*t36
      t62 = 0.2611111111111111D1-0.9850555555555556D-1*t33
     #-0.1357222222222222D0*t60
      t67 = sigmaaa+sigmabb
      t71 = t42+0.349D0*t60-11.D0
      t75 = rhoa*t38*sigmaaa+t39*sigmabb
      t78 = 0.3646239897876478D2*t53+0.3646239897876478D2*t57+t62
     #*sigma-1.D0*(0.25D1-0.1407222222222222D-1*t33
     #-0.1938888888888889D-1*t60)*t67-0.1111111111111111D0*t71*t75
      t82 = 0.6666666666666667D0*t45
      t83 = 1.D0*t51
      t86 = 1.D0*t55
      t89 = t50*t78-0.6666666666666667D0*t45*sigma+(t82-t83)*sigmabb+
     #(t82-t86)*sigmaaa
      t93 = t38**(1.D0/3.D0)
      t94 = 0.6203504908994D0*t93
      t95 = t38**(1.D0/6.D0)
      t97 = t94+0.1029581201158544D2*t95+0.427198D2
      t98 = 1/t97
      t101 = dlog(0.6203504908994D0*t93*t98)
      t103 = 0.1575246635799487D1*t95
      t104 = t103+0.13072D2
      t107 = datan(0.448998886412873D-1/t104)
      t109 = 0.7876233178997433D0*t95
      t110 = t109+0.409286D0
      t111 = t110**2
      t113 = dlog(t111*t98)
      t116 = t94+0.1584942278842832D2*t95+0.101578D3
      t117 = 1/t116
      t120 = dlog(0.6203504908994D0*t93*t117)
      t122 = t103+0.201231D2
      t125 = datan(0.1171685277708993D1/t122)
      t127 = t109+0.743294D0
      t128 = t127**2
      t130 = dlog(t128*t117)
      t135 = 0.1554535D-1*t120+0.6188180297906063D0*t125
     #+0.2667310007273315D-2*t130-0.310907D-1*t101-0.205219729378375D2
     #*t107-0.4431373767749538D-2*t113
      t137 = rhoa-1.D0*rhob
      t138 = t137*t38
      t139 = 1.D0+t138
      t140 = t139**(1.D0/3.D0)
      t143 = 1.D0-1.D0*t138
      t144 = t143**(1.D0/3.D0)
      t146 = t140*t139+t144*t143-2.D0
      t147 = t135*t146
      zk(i) = -0.74442058907928D0*t5-0.3024D-2*t8*t15
     #-0.74442058907928D0*t19-0.3024D-2*t22*t29-0.1593432D0*t37*t39
     #-0.52583256D-2*t44*t49*t89+0.19D0*rho*(0.310907D-1*t101
     #+0.205219729378375D2*t107+0.4431373767749538D-2*t113
     #+0.1923661050931536D1*t147)
      t154 = 1/t4/t51
      t158 = t14**2
      t159 = 1/t158
      t167 = 1/t53
      t170 = dsqrt(1.D0+sigmaaa*t167)
      t171 = 1/t170
      t184 = t71*t38
      t195 = rho*t135
      t200 = 0.1333333333333333D1*t140*t38-0.1333333333333333D1*t144*t38
      t203 = t35**2
      t204 = 1/t203
      t210 = 0.185369256D-1*t204*rhoa*rhob/t32/t45
      t211 = 1/t45
      t212 = rhob*t211
      t214 = 0.1593432D0*t37*t212
      t215 = t45**2
      t217 = 1/t215/rho
      t221 = 0.44397795816D-3*t217*t43*t36*t89
      t225 = 0.6117185448D-3*t43*t204*t217*t89
      t230 = 0.192805272D-1*t44/t47/t215*t89
      t232 = 1/t32/rho
      t234 = t232*t36
      t238 = 1/t47/rho*t204
      t273 = 0.52583256D-2*t44*t49*(t50*((0.3283518518518519D-1*t232
     #+0.4524074074074074D-1*t234-0.1578901851851852D-1*t238)*sigma
     #-1.D0*(0.4690740740740741D-2*t232+0.6462962962962963D-2*t234
     #-0.2255574074074074D-2*t238)*t67-0.1111111111111111D0*(
     #-0.8443333333333333D-1*t232-0.1163333333333333D0*t234
     #+0.4060033333333333D-1*t238)*t75-0.1111111111111111D0*t71*(-1.D0
     #*rhoa*t211*sigmaaa-1.D0*t212*sigmabb))-0.1333333333333333D1*rho
     #*sigma+0.1333333333333333D1*rho*sigmabb+0.1333333333333333D1*rho
     #*sigmaaa)
      t274 = 0.5907233D-2*t101
      t275 = 0.3899174858189126D1*t107
      t276 = 0.8419610158724123D-3*t113
      t277 = 0.3654955996769919D0*t147
      t278 = t93**2
      t279 = 1/t278
      t283 = t97**2
      t284 = 1/t283
      t287 = 0.2067834969664667D0*t279*t211
      t288 = t95**2
      t289 = t288**2
      t291 = 1/t289/t95
      t292 = t291*t211
      t294 = -t287-0.1715968668597574D1*t292
      t298 = 1/t93
      t300 = (-0.2067834969664667D0*t279*t98*t211-0.6203504908994D0
     #*t93*t284*t294)*t298*t97
      t302 = t104**2
      t303 = 1/t302
      t309 = t303*t291*t211/(1.D0+0.2016D-2*t303)
      t320 = (-0.2625411059665811D0*t110*t98*t292-1.D0*t111*t284*t294)
     #/t111*t97
      t325 = t116**2
      t326 = 1/t325
      t329 = -t287-0.2641570464738054D1*t292
      t336 = t122**2
      t337 = 1/t336
      t373 = 0.19D0*rho*(0.5011795824473985D-1*t300
     #+0.2419143800947354D0*t309+0.4431373767749538D-2*t320
     #+0.1923661050931536D1*(0.2505897912236993D-1*(
     #-0.2067834969664667D0*t279*t117*t211-0.6203504908994D0*t93*t326
     #*t329)*t298*t116+0.1903580477513215D0*t337*t291*t211/(1.D0
     #+0.137284639D1*t337)+0.2667310007273315D-2*(
     #-0.2625411059665811D0*t127*t117*t292-1.D0*t128*t326*t329)/t128
     #*t116-0.5011795824473985D-1*t300-0.2419143800947354D0*t309
     #-0.4431373767749538D-2*t320)*t146+0.1923661050931536D1*t135*(
     #-0.1333333333333333D1*t140*t137*t211+0.1333333333333333D1*t144
     #*t137*t211))
      vrhoa(i) = -0.99256078543904D0*t4+0.4032D-2*t154*sigmaaa*t15
     #+0.3024D-2*t8*t159*(-0.336D-1*t9*t154*t11-0.336D-1*sigmaaa/t52
     #/t51/rhoa*t171)-0.1593432D0*t36*rhob*t38-0.52583256D-2*t44*t49*
     #(rhob*t78+t50*(0.9723306394337274D2*t52*rhoa
     #-0.1111111111111111D0*t184*sigmaaa)-2.D0*rhoa*sigmabb)
     #+0.3654955996769919D0*t195*t200-t210+t214-t221-t225+t230-t273
     #+t274+t275+t276+t277+t373
      t376 = 1/t18/t55
      t380 = t28**2
      t381 = 1/t380
      t389 = 1/t57
      t392 = dsqrt(1.D0+sigmabb*t389)
      t393 = 1/t392
      vrhob(i) = -0.99256078543904D0*t18+0.4032D-2*t376*sigmabb*t29
     #+0.3024D-2*t22*t381*(-0.336D-1*t23*t376*t25-0.336D-1*sigmabb/t56
     #/t55/rhob*t393)-0.1593432D0*t37*t38-0.52583256D-2*t44*t49*(rhoa
     #*t78+t50*(0.9723306394337274D2*t56*rhob-0.1111111111111111D0
     #*t184*sigmabb)-2.D0*rhob*sigmaaa)-0.3654955996769919D0*t195*t200
     #-t210+t214-t221-t225+t230-t273+t274+t275+t276+t277+t373
      t430 = 0.1407222222222222D-1*t33
      t431 = 0.1938888888888889D-1*t60
      t445 = t44*t49*(t50*t62-0.6666666666666667D0*t45)
      t446 = 0.52583256D-2*t445
      vsigmaaa(i) = -0.3024D-2*t7*t15+0.3024D-2*t8*t159*(0.126D-1/t9
     #*t7*t11+0.126D-1*t167*t171)-0.52583256D-2*t44*t49*(t50*(-0.25D1
     #+t430+t431-0.1111111111111111D0*t71*rhoa*t38)+t82-t86)-t446
      vsigmaab(i) = -0.105166512D-1*t445
      vsigmabb(i) = -0.3024D-2*t21*t29+0.3024D-2*t22*t381*(0.126D-1
     #/t23*t21*t25+0.126D-1*t389*t393)-0.52583256D-2*t44*t49*(t50*(
     #-0.25D1+t430+t431-0.1111111111111111D0*t71*rhob*t38)+t82-t83)-t446
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(tol,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t5 = 1/t3
      t6 = t5*sigmabb
      t7 = dsqrt(sigmabb)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      t12 = 1.D0+0.252D-1*t8*t9
      t13 = 1/t12
      t16 = 1/rhob
      t17 = t16**(1.D0/3.D0)
      t19 = t16**(1.D0/6.D0)
      t21 = 0.6203504908994D0*t17+0.1584942278842832D2*t19+0.101578D3
      t22 = 1/t21
      t25 = dlog(0.6203504908994D0*t17*t22)
      t28 = 0.1575246635799487D1*t19+0.201231D2
      t31 = datan(0.1171685277708993D1/t28)
      t34 = 0.7876233178997433D0*t19+0.743294D0
      t35 = t34**2
      t37 = dlog(t35*t22)
      zk(i) = -0.74442058907928D0*t3-0.3024D-2*t6*t13+0.19D0*rhob*
     #(0.1554535D-1*t25+0.6188180297906063D0*t31+0.2667310007273315D-2
     #*t37)
      vrhoa(i) = 0.D0
      t43 = rhob**2
      t45 = 1/t2/t43
      t46 = t45*sigmabb
      t49 = t12**2
      t50 = 1/t49
      t54 = t43*rhob
      t55 = t2**2
      t60 = 1/t55/t43
      t62 = 1.D0+sigmabb*t60
      t63 = dsqrt(t62)
      t64 = 1/t63
      t67 = -0.336D-1*t7*t45*t9-0.336D-1*sigmabb/t55/t54*t64
      t68 = t50*t67
      t74 = t17**2
      t75 = 1/t74
      t76 = t75*t22
      t77 = 1/t43
      t80 = t21**2
      t81 = 1/t80
      t82 = t17*t81
      t85 = t19**2
      t86 = t85**2
      t87 = t86*t19
      t88 = 1/t87
      t89 = t88*t77
      t91 = -0.2067834969664667D0*t75*t77-0.2641570464738054D1*t89
      t94 = -0.2067834969664667D0*t76*t77-0.6203504908994D0*t82*t91
      t95 = 1/t17
      t96 = t94*t95
      t97 = t96*t21
      t99 = t28**2
      t100 = 1/t99
      t101 = t100*t88
      t103 = 1.D0+0.137284639D1*t100
      t104 = 1/t103
      t106 = t101*t77*t104
      t108 = t34*t22
      t111 = t35*t81
      t114 = -0.2625411059665811D0*t108*t89-1.D0*t111*t91
      t115 = 1/t35
      t116 = t114*t115
      t117 = t116*t21
      vrhob(i) = -0.99256078543904D0*t2+0.4032D-2*t46*t13+0.3024D-2*t6
     #*t68+0.29536165D-2*t25+0.1175754256602152D0*t31
     #+0.5067889013819299D-3*t37+0.19D0*rhob*(0.2505897912236993D-1
     #*t97+0.1903580477513215D0*t106+0.2667310007273315D-2*t117)
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t130 = 0.126D-1/t7*t5*t9+0.126D-1*t60*t64
      vsigmabb(i) = -0.3024D-2*t5*t13+0.3024D-2*t6*t50*t130
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t137 = 1/t2/t54
      t144 = 1/t49/t12
      t145 = t67**2
      t152 = t43**2
      t158 = sigmabb**2
      t164 = 1/t63/t62
      t175 = 1/t74/t16
      t177 = 1/t152
      t178 = t175*t22*t177
      t184 = 1/t54
      t188 = 1/t80/t21
      t190 = t91**2
      t198 = 1/t87/t16
      t199 = t198*t177
      t201 = t88*t184
      t203 = -0.1378556646443111D0*t175*t177+0.4135669939329333D0*t75
     #*t184-0.2201308720615045D1*t199+0.5283140929476108D1*t201
      t221 = t177*t104
      t230 = t99**2
      t234 = t103**2
      s1 = -0.3308535951463467D0/t55-0.9408D-2*t137*sigmabb*t13
     #-0.8064D-2*t46*t68-0.6048D-2*t6*t144*t145
      s2 = s1+0.3024D-2*t6*t50*(0.784D-1*t7*t137*t9+0.168D0*sigmabb
     #/t55/t152*t64-0.448D-1*t158/t2/t152/t54*t164)
     #+0.9522412066500572D-2*t97
      s3 = s2+0.7233605814550218D-1*t106
      s4 = s3
      s5 = 0.101357780276386D-2*t117+0.19D0*rhob*
     #(0.2505897912236993D-1*(-0.1378556646443111D0*t178
     #+0.4135669939329333D0*t75*t81*t77*t91+0.4135669939329333D0*t76
     #*t184+0.12407009817988D1*t17*t188*t190-0.6203504908994D0*t82
     #*t203)*t95*t21+0.8352993040789975D-2*t94/t17/t16*t21*t77
     #+0.2505897912236993D-1*t96*t91+0.9995362477254242D-1/t99/t28
     #*t175*t221+0.1586317064594346D0*t100*t198*t221
     #-0.3807160955026431D0*t101*t184*t104-0.1372209729363994D0/t230
     #/t28*t175*t177/t234+0.2667310007273315D-2*(0.3446391616107778D-1
     #*t178+0.5250822119331622D0*t34*t81*t89*t91-0.2187842549721509D0
     #*t108*t199+0.5250822119331622D0*t108*t201+2.D0*t35*t188*t190
     #-1.D0*t111*t203)*t115*t21+0.7002785192652656D-3*t114/t35/t34*t21
     #*t88*t77+0.2667310007273315D-2*t116*t91)
      v2rhob2(i) = s4+s5
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t272 = t130**2
      v2sigmabb2(i) = 0.6048D-2*t5*t50*t130-0.6048D-2*t6*t144*t272
     #+0.3024D-2*t6*t50*(-0.63D-2/t7/sigmabb*t5*t9+0.63D-2/sigmabb*t60
     #*t64-0.63D-2/t2/t152/rhob*t164)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(tol,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t5 = 1/t3
      t6 = t5*sigmaaa
      t7 = dsqrt(sigmaaa)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      t12 = 1.D0+0.252D-1*t8*t9
      t13 = 1/t12
      t16 = 1/rhoa
      t17 = t16**(1.D0/3.D0)
      t19 = t16**(1.D0/6.D0)
      t21 = 0.6203504908994D0*t17+0.1584942278842832D2*t19+0.101578D3
      t22 = 1/t21
      t25 = dlog(0.6203504908994D0*t17*t22)
      t28 = 0.1575246635799487D1*t19+0.201231D2
      t31 = datan(0.1171685277708993D1/t28)
      t34 = 0.7876233178997433D0*t19+0.743294D0
      t35 = t34**2
      t37 = dlog(t35*t22)
      zk(i) = -0.74442058907928D0*t3-0.3024D-2*t6*t13+0.19D0*rhoa*
     #(0.1554535D-1*t25+0.6188180297906063D0*t31+0.2667310007273315D-2
     #*t37)
      t43 = rhoa**2
      t45 = 1/t2/t43
      t46 = t45*sigmaaa
      t49 = t12**2
      t50 = 1/t49
      t54 = t43*rhoa
      t55 = t2**2
      t60 = 1/t55/t43
      t62 = 1.D0+sigmaaa*t60
      t63 = dsqrt(t62)
      t64 = 1/t63
      t67 = -0.336D-1*t7*t45*t9-0.336D-1*sigmaaa/t55/t54*t64
      t68 = t50*t67
      t74 = t17**2
      t75 = 1/t74
      t76 = t75*t22
      t77 = 1/t43
      t80 = t21**2
      t81 = 1/t80
      t82 = t17*t81
      t85 = t19**2
      t86 = t85**2
      t87 = t86*t19
      t88 = 1/t87
      t89 = t88*t77
      t91 = -0.2067834969664667D0*t75*t77-0.2641570464738054D1*t89
      t94 = -0.2067834969664667D0*t76*t77-0.6203504908994D0*t82*t91
      t95 = 1/t17
      t96 = t94*t95
      t97 = t96*t21
      t99 = t28**2
      t100 = 1/t99
      t101 = t100*t88
      t103 = 1.D0+0.137284639D1*t100
      t104 = 1/t103
      t106 = t101*t77*t104
      t108 = t34*t22
      t111 = t35*t81
      t114 = -0.2625411059665811D0*t108*t89-1.D0*t111*t91
      t115 = 1/t35
      t116 = t114*t115
      t117 = t116*t21
      vrhoa(i) = -0.99256078543904D0*t2+0.4032D-2*t46*t13+0.3024D-2*t6
     #*t68+0.29536165D-2*t25+0.1175754256602152D0*t31
     #+0.5067889013819299D-3*t37+0.19D0*rhoa*(0.2505897912236993D-1
     #*t97+0.1903580477513215D0*t106+0.2667310007273315D-2*t117)
      vrhob(i) = 0.D0
      t130 = 0.126D-1/t7*t5*t9+0.126D-1*t60*t64
      vsigmaaa(i) = -0.3024D-2*t5*t13+0.3024D-2*t6*t50*t130
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t137 = 1/t2/t54
      t144 = 1/t49/t12
      t145 = t67**2
      t152 = t43**2
      t158 = sigmaaa**2
      t164 = 1/t63/t62
      t175 = 1/t74/t16
      t177 = 1/t152
      t178 = t175*t22*t177
      t184 = 1/t54
      t188 = 1/t80/t21
      t190 = t91**2
      t198 = 1/t87/t16
      t199 = t198*t177
      t201 = t88*t184
      t203 = -0.1378556646443111D0*t175*t177+0.4135669939329333D0*t75
     #*t184-0.2201308720615045D1*t199+0.5283140929476108D1*t201
      t221 = t177*t104
      t230 = t99**2
      t234 = t103**2
      s1 = -0.3308535951463467D0/t55-0.9408D-2*t137*sigmaaa*t13
     #-0.8064D-2*t46*t68-0.6048D-2*t6*t144*t145
      s2 = s1+0.3024D-2*t6*t50*(0.784D-1*t7*t137*t9+0.168D0*sigmaaa
     #/t55/t152*t64-0.448D-1*t158/t2/t152/t54*t164)
     #+0.9522412066500572D-2*t97
      s3 = s2+0.7233605814550218D-1*t106
      s4 = s3
      s5 = 0.101357780276386D-2*t117+0.19D0*rhoa*
     #(0.2505897912236993D-1*(-0.1378556646443111D0*t178
     #+0.4135669939329333D0*t75*t81*t77*t91+0.4135669939329333D0*t76
     #*t184+0.12407009817988D1*t17*t188*t190-0.6203504908994D0*t82
     #*t203)*t95*t21+0.8352993040789975D-2*t94/t17/t16*t21*t77
     #+0.2505897912236993D-1*t96*t91+0.9995362477254242D-1/t99/t28
     #*t175*t221+0.1586317064594346D0*t100*t198*t221
     #-0.3807160955026431D0*t101*t184*t104-0.1372209729363994D0/t230
     #/t28*t175*t177/t234+0.2667310007273315D-2*(0.3446391616107778D-1
     #*t178+0.5250822119331622D0*t34*t81*t89*t91-0.2187842549721509D0
     #*t108*t199+0.5250822119331622D0*t108*t201+2.D0*t35*t188*t190
     #-1.D0*t111*t203)*t115*t21+0.7002785192652656D-3*t114/t35/t34*t21
     #*t88*t77+0.2667310007273315D-2*t116*t91)
      v2rhoa2(i) = s4+s5
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t272 = t130**2
      v2sigmaaa2(i) = 0.6048D-2*t5*t50*t130-0.6048D-2*t6*t144*t272
     #+0.3024D-2*t6*t50*(-0.63D-2/t7/sigmaaa*t5*t9+0.63D-2/sigmaaa*t60
     #*t64-0.63D-2/t2/t152/rhoa*t164)
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(tol,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(tol,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t7 = 1/t5
      t8 = t7*sigmaaa
      t9 = dsqrt(sigmaaa)
      t10 = t9*t7
      t11 = dlog(t10+dsqrt(1+t10**2))
      t14 = 1.D0+0.252D-1*t10*t11
      t15 = 1/t14
      t18 = rhob**(1.D0/3.D0)
      t19 = t18*rhob
      t21 = 1/t19
      t22 = t21*sigmabb
      t23 = dsqrt(sigmabb)
      t24 = t23*t21
      t25 = dlog(t24+dsqrt(1+t24**2))
      t28 = 1.D0+0.252D-1*t24*t25
      t29 = 1/t28
      t32 = rho**(1.D0/3.D0)
      t33 = 1/t32
      t35 = 1.D0+0.349D0*t33
      t36 = 1/t35
      t37 = t36*rhoa
      t38 = 1/rho
      t39 = rhob*t38
      t42 = 0.2533D0*t33
      t43 = dexp(-t42)
      t44 = t43*t36
      t45 = rho**2
      t46 = t45*rho
      t47 = t32**2
      t49 = 1/t47/t46
      t50 = rhoa*rhob
      t51 = rhoa**2
      t52 = t4**2
      t53 = t52*t51
      t54 = 0.3646239897876478D2*t53
      t55 = rhob**2
      t56 = t18**2
      t57 = t56*t55
      t58 = 0.3646239897876478D2*t57
      t60 = t33*t36
      t62 = 0.2611111111111111D1-0.9850555555555556D-1*t33
     #-0.1357222222222222D0*t60
      t63 = t62*sigma
      t67 = sigmaaa+sigmabb
      t69 = 1.D0*(0.25D1-0.1407222222222222D-1*t33
     #-0.1938888888888889D-1*t60)*t67
      t71 = t42+0.349D0*t60-11.D0
      t75 = rhoa*t38*sigmaaa+t39*sigmabb
      t77 = 0.1111111111111111D0*t71*t75
      t78 = t54+t58+t63-t69-t77
      t82 = 0.6666666666666667D0*t45
      t83 = 1.D0*t51
      t86 = 1.D0*t55
      t89 = t50*t78-0.6666666666666667D0*t45*sigma+(t82-t83)*sigmabb+
     #(t82-t86)*sigmaaa
      t93 = t38**(1.D0/3.D0)
      t94 = 0.6203504908994D0*t93
      t95 = t38**(1.D0/6.D0)
      t97 = t94+0.1029581201158544D2*t95+0.427198D2
      t98 = 1/t97
      t101 = dlog(0.6203504908994D0*t93*t98)
      t103 = 0.1575246635799487D1*t95
      t104 = t103+0.13072D2
      t107 = datan(0.448998886412873D-1/t104)
      t109 = 0.7876233178997433D0*t95
      t110 = t109+0.409286D0
      t111 = t110**2
      t113 = dlog(t111*t98)
      t116 = t94+0.1584942278842832D2*t95+0.101578D3
      t117 = 1/t116
      t120 = dlog(0.6203504908994D0*t93*t117)
      t122 = t103+0.201231D2
      t125 = datan(0.1171685277708993D1/t122)
      t127 = t109+0.743294D0
      t128 = t127**2
      t130 = dlog(t128*t117)
      t135 = 0.1554535D-1*t120+0.6188180297906063D0*t125
     #+0.2667310007273315D-2*t130-0.310907D-1*t101-0.205219729378375D2
     #*t107-0.4431373767749538D-2*t113
      t137 = rhoa-1.D0*rhob
      t138 = t137*t38
      t139 = 1.D0+t138
      t140 = t139**(1.D0/3.D0)
      t143 = 1.D0-1.D0*t138
      t144 = t143**(1.D0/3.D0)
      t146 = t140*t139+t144*t143-2.D0
      t147 = t135*t146
      zk(i) = -0.74442058907928D0*t5-0.3024D-2*t8*t15
     #-0.74442058907928D0*t19-0.3024D-2*t22*t29-0.1593432D0*t37*t39
     #-0.52583256D-2*t44*t49*t89+0.19D0*rho*(0.310907D-1*t101
     #+0.205219729378375D2*t107+0.4431373767749538D-2*t113
     #+0.1923661050931536D1*t147)
      t154 = 1/t4/t51
      t155 = t154*sigmaaa
      t158 = t14**2
      t159 = 1/t158
      t163 = t51*rhoa
      t165 = 1/t52/t163
      t167 = 1/t53
      t169 = 1.D0+sigmaaa*t167
      t170 = dsqrt(t169)
      t171 = 1/t170
      t174 = -0.336D-1*t9*t154*t11-0.336D-1*sigmaaa*t165*t171
      t175 = t159*t174
      t178 = t36*rhob
      t182 = t52*rhoa
      t184 = t71*t38
      t187 = 0.9723306394337274D2*t182-0.1111111111111111D0*t184*sigmaaa
      t191 = rhob*t78+t50*t187-2.D0*rhoa*sigmabb
      t195 = rho*t135
      t200 = 0.1333333333333333D1*t140*t38-0.1333333333333333D1*t144*t38
      t203 = t35**2
      t204 = 1/t203
      t205 = t204*rhoa
      t207 = 1/t32/t45
      t210 = 0.185369256D-1*t205*rhob*t207
      t211 = 1/t45
      t212 = rhob*t211
      t214 = 0.1593432D0*t37*t212
      t215 = t45**2
      t216 = t215*rho
      t217 = 1/t216
      t218 = t217*t43
      t219 = t36*t89
      t221 = 0.44397795816D-3*t218*t219
      t222 = t43*t204
      t225 = 0.6117185448D-3*t222*t217*t89
      t227 = 1/t47/t215
      t230 = 0.192805272D-1*t44*t227*t89
      t232 = 1/t32/rho
      t234 = t232*t36
      t238 = 1/t47/rho*t204
      t240 = 0.3283518518518519D-1*t232+0.4524074074074074D-1*t234
     #-0.1578901851851852D-1*t238
      t251 = -0.8443333333333333D-1*t232-0.1163333333333333D0*t234
     #+0.4060033333333333D-1*t238
      t259 = -1.D0*rhoa*t211*sigmaaa-1.D0*t212*sigmabb
      t262 = t240*sigma-1.D0*(0.4690740740740741D-2*t232
     #+0.6462962962962963D-2*t234-0.2255574074074074D-2*t238)*t67
     #-0.1111111111111111D0*t251*t75-0.1111111111111111D0*t71*t259
      t270 = t50*t262-0.1333333333333333D1*rho*sigma
     #+0.1333333333333333D1*rho*sigmabb+0.1333333333333333D1*rho*sigmaaa
      t273 = 0.52583256D-2*t44*t49*t270
      t274 = 0.5907233D-2*t101
      t275 = 0.3899174858189126D1*t107
      t276 = 0.8419610158724123D-3*t113
      t277 = 0.3654955996769919D0*t147
      t278 = t93**2
      t279 = 1/t278
      t280 = t279*t98
      t283 = t97**2
      t284 = 1/t283
      t285 = t93*t284
      t287 = 0.2067834969664667D0*t279*t211
      t288 = t95**2
      t289 = t288**2
      t290 = t289*t95
      t291 = 1/t290
      t292 = t291*t211
      t294 = -t287-0.1715968668597574D1*t292
      t297 = -0.2067834969664667D0*t280*t211-0.6203504908994D0*t285*t294
      t298 = 1/t93
      t299 = t297*t298
      t300 = t299*t97
      t302 = t104**2
      t303 = 1/t302
      t304 = t303*t291
      t306 = 1.D0+0.2016D-2*t303
      t307 = 1/t306
      t309 = t304*t211*t307
      t311 = t110*t98
      t314 = t111*t284
      t317 = -0.2625411059665811D0*t311*t292-1.D0*t314*t294
      t318 = 1/t111
      t319 = t317*t318
      t320 = t319*t97
      t322 = t279*t117
      t325 = t116**2
      t326 = 1/t325
      t327 = t93*t326
      t329 = -t287-0.2641570464738054D1*t292
      t332 = -0.2067834969664667D0*t322*t211-0.6203504908994D0*t327*t329
      t333 = t332*t298
      t336 = t122**2
      t337 = 1/t336
      t338 = t337*t291
      t340 = 1.D0+0.137284639D1*t337
      t341 = 1/t340
      t345 = t127*t117
      t348 = t128*t326
      t351 = -0.2625411059665811D0*t345*t292-1.D0*t348*t329
      t352 = 1/t128
      t353 = t351*t352
      t359 = 0.2505897912236993D-1*t333*t116+0.1903580477513215D0*t338
     #*t211*t341+0.2667310007273315D-2*t353*t116-0.5011795824473985D-1
     #*t300-0.2419143800947354D0*t309-0.4431373767749538D-2*t320
      t360 = t359*t146
      t362 = t140*t137
      t365 = t144*t137
      t368 = -0.1333333333333333D1*t362*t211+0.1333333333333333D1*t365
     #*t211
      t369 = t135*t368
      t373 = 0.19D0*rho*(0.5011795824473985D-1*t300
     #+0.2419143800947354D0*t309+0.4431373767749538D-2*t320
     #+0.1923661050931536D1*t360+0.1923661050931536D1*t369)
      vrhoa(i) = -0.99256078543904D0*t4+0.4032D-2*t155*t15+0.3024D-2
     #*t8*t175-0.1593432D0*t178*t38-0.52583256D-2*t44*t49*t191
     #+0.3654955996769919D0*t195*t200-t210+t214-t221-t225+t230-t273
     #+t274+t275+t276+t277+t373
      t376 = 1/t18/t55
      t377 = t376*sigmabb
      t380 = t28**2
      t381 = 1/t380
      t385 = t55*rhob
      t387 = 1/t56/t385
      t389 = 1/t57
      t391 = 1.D0+sigmabb*t389
      t392 = dsqrt(t391)
      t393 = 1/t392
      t396 = -0.336D-1*t23*t376*t25-0.336D-1*sigmabb*t387*t393
      t397 = t381*t396
      t403 = t56*rhob
      t407 = 0.9723306394337274D2*t403-0.1111111111111111D0*t184*sigmabb
      t411 = rhoa*t78+t50*t407-2.D0*rhob*sigmaaa
      t415 = -t200
      vrhob(i) = -0.99256078543904D0*t18+0.4032D-2*t377*t29+0.3024D-2
     #*t22*t397-0.1593432D0*t37*t38-0.52583256D-2*t44*t49*t411
     #+0.3654955996769919D0*t195*t415-t210+t214-t221-t225+t230-t273
     #+t274+t275+t276+t277+t373
      t420 = 1/t9
      t426 = 0.126D-1*t420*t7*t11+0.126D-1*t167*t171
      t427 = t159*t426
      t430 = 0.1407222222222222D-1*t33
      t431 = 0.1938888888888889D-1*t60
      t432 = t71*rhoa
      t435 = -0.25D1+t430+t431-0.1111111111111111D0*t432*t38
      t437 = t50*t435+t82-t86
      t443 = t50*t62-0.6666666666666667D0*t45
      t445 = t44*t49*t443
      t446 = 0.52583256D-2*t445
      vsigmaaa(i) = -0.3024D-2*t7*t15+0.3024D-2*t8*t427-0.52583256D-2
     #*t44*t49*t437-t446
      vsigmaab(i) = -0.105166512D-1*t445
      t450 = 1/t23
      t456 = 0.126D-1*t450*t21*t25+0.126D-1*t389*t393
      t457 = t381*t456
      t460 = t71*rhob
      t463 = -0.25D1+t430+t431-0.1111111111111111D0*t460*t38
      t465 = t50*t463+t82-t83
      vsigmabb(i) = -0.3024D-2*t21*t29+0.3024D-2*t22*t457
     #-0.52583256D-2*t44*t49*t465-t446
      t469 = 0.9192746443599945D-1*t309
      t470 = 0.1904482413300114D-1*t300
      t471 = 0.1683922031744825D-2*t320
      t472 = 0.7309911993539838D0*t360
      t473 = 0.7309911993539838D0*t369
      t478 = 0.899757936D-1*t44/t47/t216*t89
      t483 = 0.61789752D-1*t205*rhob/t32/t46
      t484 = t215*t45
      t485 = 1/t484
      t488 = 0.53015607216D-2*t222*t485*t89
      t490 = 1/t32/t484
      t491 = t490*t43
      t494 = 0.10329887159856D-3*t491*t204*t89
      t497 = 0.385610544D-1*t44*t227*t270
      t500 = 0.12234370896D-2*t222*t217*t270
      t501 = 1/t46
      t502 = rhob*t501
      t503 = t37*t502
      t504 = 0.3186864D0*t503
      t506 = t207*t36
      t510 = 1/t47/t45*t204
      t513 = 1/t203/t35
      t514 = t501*t513
      t549 = t44*t49*(t50*((-0.4378024691358025D-1*t207
     #-0.6032098765432099D-1*t506+0.3157803703703704D-1*t510
     #-0.3673578308641975D-2*t514)*sigma-1.D0*(-0.6254320987654321D-2
     #*t207-0.8617283950617284D-2*t506+0.4511148148148148D-2*t510
     #-0.5247969012345679D-3*t514)*t67-0.1111111111111111D0*
     #(0.1125777777777778D0*t207+0.1551111111111111D0*t506
     #-0.8120066666666667D-1*t510+0.9446344222222222D-2*t514)*t75
     #-0.2222222222222222D0*t251*t259-0.1111111111111111D0*t71*(2.D0
     #*rhoa*t501*sigmaaa+2.D0*t502*sigmabb))-0.1333333333333333D1
     #*sigma+0.1333333333333333D1*sigmabb+0.1333333333333333D1*sigmaaa)
      t550 = 0.52583256D-2*t549
      t553 = 0.88795591632D-3*t218*t36*t270
      t555 = 0.37486538933976D-4*t491*t219
      t559 = 1/t158/t14
      t560 = t174**2
      t565 = 1/t4/t163
      t569 = t51**2
      t575 = sigmaaa**2
      t581 = 1/t170/t169
      t588 = t469+t470+t471+t472+t473-t478+t483+t488-t494+t497-t500
     #-t504-t550-t553-t555-0.8064D-2*t155*t175-0.6048D-2*t8*t559*t560
     #+0.3024D-2*t8*t159*(0.784D-1*t9*t565*t11+0.168D0*sigmaaa/t52
     #/t569*t171-0.448D-1*t575/t4/t569/t163*t581)
      t589 = rho*t359
      t590 = t589*t200
      t593 = 1/t278/t38
      t595 = 1/t215
      t596 = t593*t98*t595
      t605 = 1/t283/t97
      t607 = t294**2
      t611 = 0.1378556646443111D0*t593*t595
      t613 = 0.4135669939329333D0*t279*t501
      t615 = 1/t290/t38
      t616 = t615*t595
      t618 = t291*t501
      t620 = -t611+t613-0.1429973890497978D1*t616+0.3431937337195148D1
     #*t618
      t625 = (-0.1378556646443111D0*t596+0.4135669939329333D0*t279
     #*t284*t211*t294+0.4135669939329333D0*t280*t501
     #+0.12407009817988D1*t93*t605*t607-0.6203504908994D0*t285*t620)
     #*t298*t97
      t628 = 1/t93/t38
      t631 = t297*t628*t97*t211
      t633 = t299*t294
      t638 = t595*t307
      t639 = 1/t302/t104*t593*t638
      t642 = t303*t615*t638
      t645 = t304*t501*t307
      t647 = t302**2
      t651 = t306**2
      t654 = 1/t647/t104*t593*t595/t651
      t672 = (0.3446391616107778D-1*t596+0.5250822119331622D0*t110
     #*t284*t292*t294-0.2187842549721509D0*t311*t616
     #+0.5250822119331622D0*t311*t618+2.D0*t111*t605*t607-1.D0*t314
     #*t620)*t318*t97
      t679 = t317/t111/t110*t97*t291*t211
      t681 = t319*t294
      t684 = t593*t117*t595
      t693 = 1/t325/t116
      t695 = t329**2
      t700 = -t611+t613-0.2201308720615045D1*t616+0.5283140929476108D1
     #*t618
      t716 = t595*t341
      t725 = t336**2
      t729 = t340**2
      s1 = 0.2505897912236993D-1*(-0.1378556646443111D0*t684
     #+0.4135669939329333D0*t279*t326*t211*t329+0.4135669939329333D0
     #*t322*t501+0.12407009817988D1*t93*t693*t695-0.6203504908994D0
     #*t327*t700)*t298*t116+0.8352993040789975D-2*t332*t628*t116*t211
     #+0.2505897912236993D-1*t333*t329+0.9995362477254242D-1/t336/t122
     #*t593*t716+0.1586317064594346D0*t337*t615*t716
     #-0.3807160955026431D0*t338*t501*t341-0.1372209729363994D0/t725
     #/t122*t593*t595/t729+0.2667310007273315D-2*
     #(0.3446391616107778D-1*t684+0.5250822119331622D0*t127*t326*t292
     #*t329-0.2187842549721509D0*t345*t616+0.5250822119331622D0*t345
     #*t618+2.D0*t128*t693*t695-1.D0*t348*t700)*t352*t116
     #+0.7002785192652656D-3*t351/t128/t127*t116*t291*t211
     #+0.2667310007273315D-2*t353*t329
      t771 = s1-0.5011795824473985D-1*t625-0.1670598608157995D-1*t631
     #-0.5011795824473985D-1*t633-0.1270249377985834D0*t639
     #-0.2015953167456128D0*t642+0.4838287601894708D0*t645
     #+0.2560822746019441D-3*t654-0.4431373767749538D-2*t672
     #-0.1163417769936259D-2*t679-0.4431373767749538D-2*t681
      t776 = t140**2
      t777 = 1/t776
      t778 = t137**2
      t784 = t144**2
      t785 = 1/t784
      t794 = 0.5011795824473985D-1*t625+0.1670598608157995D-1*t631
     #+0.5011795824473985D-1*t633+0.1270249377985834D0*t639
     #+0.2015953167456128D0*t642-0.4838287601894708D0*t645
     #-0.2560822746019441D-3*t654+0.4431373767749538D-2*t672
     #+0.1163417769936259D-2*t679+0.4431373767749538D-2*t681
     #+0.1923661050931536D1*t771*t146+0.3847322101863073D1*t359*t368
     #+0.1923661050931536D1*t135*(0.4444444444444444D0*t777*t778*t595
     #+0.2666666666666667D1*t362*t501+0.4444444444444444D0*t785*t778
     #*t595-0.2666666666666667D1*t365*t501)
      t795 = rho*t794
      t796 = 0.19D0*t795
      t814 = -0.4444444444444444D0*t777*t137*t501-0.1333333333333333D1
     #*t140*t211-0.4444444444444444D0*t785*t137*t501
     #+0.1333333333333333D1*t144*t211
      t818 = rho*(0.1923661050931536D1*t359*t200+0.1923661050931536D1
     #*t135*t814)
      t823 = 0.1423265147568D-3*t43*t513*t490*t89
      t835 = 0.384780897072D-2*t485*t43*t219
      t837 = t218*t36*t191
      t840 = t222*t217*t191
      t843 = t44*t227*t191
      t846 = t251*t38
      t849 = t71*t211
      t856 = t44*t49*(rhob*t262+t50*(-0.1111111111111111D0*t846
     #*sigmaaa+0.1111111111111111D0*t849*sigmaaa))
      t859 = rhob*t49
      t861 = 0.43129246896D-2*t513*rhoa*t859
      t862 = t135*t200
      t868 = 0.4444444444444444D0*t777*t211+0.4444444444444444D0*t785
     #*t211
      t870 = 0.3654955996769919D0*t195*t868
      t872 = t204*rhob*t207
      t874 = t178*t211
      t876 = t195*t814
      t878 = 0.3654955996769919D0*t590+t796-0.9408D-2*t565*sigmaaa*t15
     #-0.3308535951463467D0/t52+0.19D0*t818-t823-0.52583256D-2*t44*t49
     #*(2.D0*rhob*t187+0.1620551065722879D3*t182*rhob-2.D0*sigmabb)
     #+t835-0.88795591632D-3*t837-0.12234370896D-2*t840+0.385610544D-1
     #*t843-0.105166512D-1*t856-t861+0.7309911993539838D0*t862+t870
     #-0.370738512D-1*t872+0.3186864D0*t874+0.3654955996769919D0*t876
      v2rhoa2(i) = t588+t878
      t879 = t135*t415
      t884 = 1/t380/t28
      t885 = t396**2
      t890 = 1/t18/t385
      t894 = t55**2
      t900 = sigmabb**2
      t906 = 1/t392/t391
      t923 = t218*t36*t411
      t926 = t222*t217*t411
      t929 = t44*t227*t411
      t940 = t44*t49*(rhoa*t262+t50*(-0.1111111111111111D0*t846
     #*sigmabb+0.1111111111111111D0*t849*sigmabb))
      t942 = t469+t470+t471+t472+t473+0.7309911993539838D0*t879
     #-0.8064D-2*t377*t397-0.6048D-2*t22*t884*t885+0.3024D-2*t22*t381*
     #(0.784D-1*t23*t890*t25+0.168D0*sigmabb/t56/t894*t393-0.448D-1
     #*t900/t18/t894/t385*t906)-0.52583256D-2*t44*t49*(2.D0*rhoa*t407
     #+0.1620551065722879D3*rhoa*t403-2.D0*sigmaaa)-0.88795591632D-3
     #*t923-0.12234370896D-2*t926+0.385610544D-1*t929-0.105166512D-1
     #*t940-t478+t483+t488-t494
      t950 = -t814
      t954 = rho*(0.1923661050931536D1*t359*t415+0.1923661050931536D1
     #*t135*t950)
      t956 = t205*t207
      t958 = t37*t211
      t960 = t195*t950
      t962 = t589*t415
      t964 = t497-t500-t504-t550-t553-t555-0.9408D-2*t890*sigmabb*t29
     #-0.3308535951463467D0/t56+0.19D0*t954-0.370738512D-1*t956
     #+0.3186864D0*t958+0.3654955996769919D0*t960+0.3654955996769919D0
     #*t962+t796-t823+t835-t861+t870
      v2rhob2(i) = t942+t964
      t974 = -t478+t483+t488-t494+t497-t500-0.3186864D0*t503
     #-0.52583256D-2*t549-t553-t555+0.182747799838496D0*t590
      t1003 = -0.6117185448D-3*t840+0.192805272D-1*t843-0.52583256D-2
     #*t856-t861+0.3654955996769919D0*t862-0.185369256D-1*t872
     #+0.1593432D0*t874+0.182747799838496D0*t876-0.52583256D-2*t44*t49
     #*(t54+t58+t63-t69-t77+rhob*t407+rhoa*t187)-0.1593432D0*t36*t38
     #-0.3654955996769919D0*t195*t868
      v2rhoab(i) = t469+t470+t471+t472+t473+0.3654955996769919D0*t879
     #-0.44397795816D-3*t923-0.6117185448D-3*t926+0.192805272D-1*t929
     #-0.52583256D-2*t940+t974+0.95D-1*t954-0.185369256D-1*t956
     #+0.1593432D0*t958+0.182747799838496D0*t960+0.182747799838496D0
     #*t962+0.19D0*t795+0.95D-1*t818-t823+t835-0.44397795816D-3*t837
     #+t1003
      t1009 = t7*t159
      t1033 = 0.1111111111111111D0*t50*t184
      t1040 = 0.44397795816D-3*t218*t36*t437
      t1043 = 0.6117185448D-3*t222*t217*t437
      t1046 = 0.192805272D-1*t44*t227*t437
      t1047 = 0.4690740740740741D-2*t232
      t1048 = 0.6462962962962963D-2*t234
      t1049 = 0.2255574074074074D-2*t238
      t1057 = 0.1333333333333333D1*rho
      t1061 = 0.52583256D-2*t44*t49*(t50*(-t1047-t1048+t1049
     #-0.1111111111111111D0*t251*rhoa*t38+0.1111111111111111D0*t432
     #*t211)+t1057)
      t1063 = t44*t859*t62
      t1064 = 0.52583256D-2*t1063
      t1066 = t218*t36*t443
      t1067 = 0.44397795816D-3*t1066
      t1069 = t222*t217*t443
      t1070 = 0.6117185448D-3*t1069
      t1072 = t44*t227*t443
      t1073 = 0.192805272D-1*t1072
      t1078 = t44*t49*(t50*t240-0.1333333333333333D1*rho)
      t1079 = 0.52583256D-2*t1078
      v2rhoasigmaaa(i) = 0.4032D-2*t154*t15-0.4032D-2*t155*t427
     #+0.3024D-2*t1009*t174-0.6048D-2*t8*t559*t174*t426+0.3024D-2*t8
     #*t159*(-0.168D-1*t420*t154*t11-0.504D-1*t165*t171+0.168D-1
     #*sigmaaa/t4/t569/t51*t581)-0.52583256D-2*t44*t49*(rhob*t435
     #-t1033)-t1040-t1043+t1046-t1061-t1064-t1067-t1070+t1073-t1079
      t1081 = 0.88795591632D-3*t1066
      t1082 = 0.12234370896D-2*t1069
      t1083 = 0.385610544D-1*t1072
      t1084 = 0.105166512D-1*t1078
      v2rhoasigmaab(i) = -0.105166512D-1*t1063-t1081-t1082+t1083-t1084
      t1093 = 0.44397795816D-3*t218*t36*t465
      t1096 = 0.6117185448D-3*t222*t217*t465
      t1099 = 0.192805272D-1*t44*t227*t465
      t1110 = 0.52583256D-2*t44*t49*(t50*(-t1047-t1048+t1049
     #-0.1111111111111111D0*t251*rhob*t38+0.1111111111111111D0*t460
     #*t211)+t1057)
      v2rhoasigmabb(i) = -0.52583256D-2*t44*t49*(rhob*t463-2.D0*rhoa)
     #-t1093-t1096+t1099-t1110-t1064-t1067-t1070+t1073-t1079
      t1119 = t44*t49*rhoa*t62
      t1120 = 0.52583256D-2*t1119
      v2rhobsigmaaa(i) = -0.52583256D-2*t44*t49*(rhoa*t435-2.D0*rhob)
     #-t1040-t1043+t1046-t1061-t1120-t1067-t1070+t1073-t1079
      v2rhobsigmaab(i) = -0.105166512D-1*t1119-t1081-t1082+t1083-t1084
      t1126 = t21*t381
      v2rhobsigmabb(i) = 0.4032D-2*t376*t29-0.4032D-2*t377*t457
     #+0.3024D-2*t1126*t396-0.6048D-2*t22*t884*t396*t456+0.3024D-2*t22
     #*t381*(-0.168D-1*t450*t376*t25-0.504D-1*t387*t393+0.168D-1
     #*sigmabb/t18/t894/t55*t906)-0.52583256D-2*t44*t49*(rhoa*t463
     #-t1033)-t1093-t1096+t1099-t1110-t1120-t1067-t1070+t1073-t1079
      t1155 = t426**2
      v2sigmaaa2(i) = 0.6048D-2*t1009*t426-0.6048D-2*t8*t559*t1155
     #+0.3024D-2*t8*t159*(-0.63D-2/t9/sigmaaa*t7*t11+0.63D-2/sigmaaa
     #*t167*t171-0.63D-2/t4/t569/rhoa*t581)
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t1179 = t456**2
      v2sigmabb2(i) = 0.6048D-2*t1126*t456-0.6048D-2*t22*t884*t1179
     #+0.3024D-2*t22*t381*(-0.63D-2/t23/sigmabb*t21*t25+0.63D-2
     #/sigmabb*t389*t393-0.63D-2/t18/t894/rhob*t906)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2rhoasigmaab(i) = 0.0d0
      v2rhoasigmabb(i) = 0.0d0
      v2rhobsigmaaa(i) = 0.0d0
      v2rhobsigmaab(i) = 0.0d0
      v2rhobsigmabb(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      v2sigmaab2(i) = 0.0d0
      v2sigmabb2(i) = 0.0d0
      v2sigmaaaab(i) = 0.0d0
      v2sigmaaabb(i) = 0.0d0
      v2sigmaabbb(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_xc_b3lyp
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     P.J. Stephens, F.J. Devlin, C.F. Chabalowski, M.J. Frisch
c     Ab initio calculation of vibrational absorption and circular 
c     dichroism spectra using density functional force fields
c     J. Phys. Chem. 98 (1994) 11623-11627
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt,i
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(tol,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t5 = 1/t3
      t7 = dsqrt(sigma)
      t8 = t7*t5
      t10 = dlog(0.1259921049894873D1*t8+dsqrt(1+0.1587401051968199D1
     #*t8**2))
      t17 = 1/t2
      t20 = 1/(1.D0+0.349D0*t17)
      t23 = 0.2533D0*t17
      t24 = dexp(-t23)
      t26 = rho**2
      t28 = t2**2
      t34 = t17*t20
      t56 = 1/rho
      t57 = t56**(1.D0/3.D0)
      t59 = t56**(1.D0/6.D0)
      t62 = 1/(0.6203504908994D0*t57+0.1029581201158544D2*t59
     #+0.427198D2)
      t65 = dlog(0.6203504908994D0*t57*t62)
      t71 = datan(0.448998886412873D-1/(0.1575246635799487D1*t59
     #+0.13072D2))
      t75 = (0.7876233178997433D0*t59+0.409286D0)**2
      t77 = dlog(t75*t62)
      zk(i) = -0.5908470131056179D0*t3-0.3810001254882096D-2*t5*sigma/
     #(1.D0+0.317500104573508D-1*t8*t10)-0.398358D-1*t20*rho
     #-0.52583256D-2*t24*t20/t28/t26/rho*(0.25D0*t26*
     #(0.1148493600075277D2*t28*t26+(0.2611111111111111D1
     #-0.9850555555555556D-1*t17-0.1357222222222222D0*t34)*sigma-0.5D0
     #*(0.25D1-0.1407222222222222D-1*t17-0.1938888888888889D-1*t34)
     #*sigma-0.2777777777777778D-1*(t23+0.349D0*t34-11.D0)*sigma)
     #-0.4583333333333333D0*t26*sigma)+0.19D0*rho*(0.310907D-1*t65
     #+0.205219729378375D2*t71+0.4431373767749538D-2*t77)
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(tol,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t5 = 1/t3
      t6 = t5*sigma
      t7 = dsqrt(sigma)
      t8 = t7*t5
      t10 = dlog(0.1259921049894873D1*t8+dsqrt(1+0.1587401051968199D1
     #*t8**2))
      t13 = 1.D0+0.317500104573508D-1*t8*t10
      t14 = 1/t13
      t17 = 1/t2
      t19 = 1.D0+0.349D0*t17
      t20 = 1/t19
      t23 = 0.2533D0*t17
      t24 = dexp(-t23)
      t25 = t24*t20
      t26 = rho**2
      t28 = t2**2
      t30 = 1/t28/t26/rho
      t31 = t28*t26
      t34 = t17*t20
      t36 = 0.2611111111111111D1-0.9850555555555556D-1*t17
     #-0.1357222222222222D0*t34
      t44 = t23+0.349D0*t34-11.D0
      t47 = 0.1148493600075277D2*t31+t36*sigma-0.5D0*(0.25D1
     #-0.1407222222222222D-1*t17-0.1938888888888889D-1*t34)*sigma
     #-0.2777777777777778D-1*t44*sigma
      t52 = 0.25D0*t26*t47-0.4583333333333333D0*t26*sigma
      t56 = 1/rho
      t57 = t56**(1.D0/3.D0)
      t59 = t56**(1.D0/6.D0)
      t61 = 0.6203504908994D0*t57+0.1029581201158544D2*t59+0.427198D2
      t62 = 1/t61
      t65 = dlog(0.6203504908994D0*t57*t62)
      t68 = 0.1575246635799487D1*t59+0.13072D2
      t71 = datan(0.448998886412873D-1/t68)
      t74 = 0.7876233178997433D0*t59+0.409286D0
      t75 = t74**2
      t77 = dlog(t75*t62)
      zk(i) = -0.5908470131056179D0*t3-0.3810001254882096D-2*t6*t14
     #-0.398358D-1*t20*rho-0.52583256D-2*t25*t30*t52+0.19D0*rho*
     #(0.310907D-1*t65+0.205219729378375D2*t71+0.4431373767749538D-2
     #*t77)
      t84 = 1/t2/t26
      t88 = t13**2
      t89 = 1/t88
      t94 = 1/t31
      t98 = dsqrt(1.D0+0.1587401051968199D1*sigma*t94)
      t99 = 1/t98
      t109 = t28*rho
      t112 = t44*t56*sigma
      t117 = rho*sigma
      t123 = t19**2
      t124 = 1/t123
      t127 = t26**2
      t129 = 1/t127/rho
      t144 = t5*t20
      t146 = 1/t109
      t147 = t146*t124
      t175 = t57**2
      t176 = 1/t175
      t178 = 1/t26
      t181 = t61**2
      t182 = 1/t181
      t186 = t59**2
      t187 = t186**2
      t189 = 1/t187/t59
      t190 = t189*t178
      t192 = -0.2067834969664667D0*t176*t178-0.1715968668597574D1*t190
      t200 = t68**2
      t201 = 1/t200
      s1 = -0.7877960174741572D0*t2+0.5080001673176129D-2*t84*sigma
     #*t14+0.1905000627441048D-2*t6*t89*(-0.8466669455293548D-1*t7*t84
     #*t10-0.106673350692263D0*sigma*t30*t99)-0.398358D-1*t20
     #-0.52583256D-2*t25*t30*(0.5D0*rho*t47+0.25D0*t26*
     #(0.3062649600200738D2*t109-0.2777777777777778D-1*t112)-0.25D0
     #*t117)-0.46342314D-2*t124*t17-0.44397795816D-3*t129*t24*t20*t52
      s2 = s1-0.6117185448D-3*t24*t124*t129*t52+0.192805272D-1*t25/t28
     #/t127*t52-0.52583256D-2*t25*t30*(0.25D0*t26*(
     #(0.3283518518518519D-1*t5+0.4524074074074074D-1*t144
     #-0.1578901851851852D-1*t147)*sigma-0.5D0*(0.4690740740740741D-2
     #*t5+0.6462962962962963D-2*t144-0.2255574074074074D-2*t147)*sigma
     #-0.2777777777777778D-1*(-0.8443333333333333D-1*t5
     #-0.1163333333333333D0*t144+0.4060033333333333D-1*t147)*sigma
     #+0.2777777777777778D-1*t112)-0.6666666666666667D0*t117)
      vrhoa(i) = s2+0.5907233D-2*t65+0.3899174858189126D1*t71
     #+0.8419610158724123D-3*t77+0.19D0*rho*(0.5011795824473985D-1*(
     #-0.2067834969664667D0*t176*t62*t178-0.6203504908994D0*t57*t182
     #*t192)/t57*t61+0.2419143800947354D0*t201*t189*t178/(1.D0
     #+0.2016D-2*t201)+0.4431373767749538D-2*(-0.2625411059665811D0
     #*t74*t62*t190-1.D0*t75*t182*t192)/t75*t61)
      vsigmaaa(i) = -0.1524000501952839D-1*t5*t14
     #+0.3810001254882096D-2*t6*t89*(0.6350002091470161D-1/t7*t5*t10
     #+0.8000501301919725D-1*t94*t99)+0.5842584D-3*t25*t146
     #-0.210333024D-1*t25*t30*(0.25D0*t26*t36-0.6666666666666667D0*t26)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(tol,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t5 = 1/t3
      t6 = t5*sigma
      t7 = dsqrt(sigma)
      t8 = t7*t5
      t10 = dlog(0.1259921049894873D1*t8+dsqrt(1+0.1587401051968199D1
     #*t8**2))
      t13 = 1.D0+0.317500104573508D-1*t8*t10
      t14 = 1/t13
      t17 = 1/t2
      t19 = 1.D0+0.349D0*t17
      t20 = 1/t19
      t23 = 0.2533D0*t17
      t24 = dexp(-t23)
      t25 = t24*t20
      t26 = rho**2
      t27 = t26*rho
      t28 = t2**2
      t30 = 1/t28/t27
      t31 = t28*t26
      t32 = 0.1148493600075277D2*t31
      t34 = t17*t20
      t36 = 0.2611111111111111D1-0.9850555555555556D-1*t17
     #-0.1357222222222222D0*t34
      t37 = t36*sigma
      t42 = 0.5D0*(0.25D1-0.1407222222222222D-1*t17
     #-0.1938888888888889D-1*t34)*sigma
      t44 = t23+0.349D0*t34-11.D0
      t46 = 0.2777777777777778D-1*t44*sigma
      t47 = t32+t37-t42-t46
      t52 = 0.25D0*t26*t47-0.4583333333333333D0*t26*sigma
      t56 = 1/rho
      t57 = t56**(1.D0/3.D0)
      t59 = t56**(1.D0/6.D0)
      t61 = 0.6203504908994D0*t57+0.1029581201158544D2*t59+0.427198D2
      t62 = 1/t61
      t65 = dlog(0.6203504908994D0*t57*t62)
      t68 = 0.1575246635799487D1*t59+0.13072D2
      t71 = datan(0.448998886412873D-1/t68)
      t74 = 0.7876233178997433D0*t59+0.409286D0
      t75 = t74**2
      t77 = dlog(t75*t62)
      zk(i) = -0.5908470131056179D0*t3-0.3810001254882096D-2*t6*t14
     #-0.398358D-1*t20*rho-0.52583256D-2*t25*t30*t52+0.19D0*rho*
     #(0.310907D-1*t65+0.205219729378375D2*t71+0.4431373767749538D-2
     #*t77)
      t84 = 1/t2/t26
      t85 = t84*sigma
      t88 = t13**2
      t89 = 1/t88
      t94 = 1/t31
      t97 = 1.D0+0.1587401051968199D1*sigma*t94
      t98 = dsqrt(t97)
      t99 = 1/t98
      t102 = -0.8466669455293548D-1*t7*t84*t10-0.106673350692263D0
     #*sigma*t30*t99
      t103 = t89*t102
      t109 = t28*rho
      t111 = t44*t56
      t112 = t111*sigma
      t114 = 0.3062649600200738D2*t109-0.2777777777777778D-1*t112
      t117 = rho*sigma
      t119 = 0.5D0*rho*t47+0.25D0*t26*t114-0.25D0*t117
      t123 = t19**2
      t124 = 1/t123
      t127 = t26**2
      t128 = t127*rho
      t129 = 1/t128
      t130 = t129*t24
      t131 = t20*t52
      t134 = t24*t124
      t139 = 1/t28/t127
      t144 = t5*t20
      t146 = 1/t109
      t147 = t146*t124
      t149 = 0.3283518518518519D-1*t5+0.4524074074074074D-1*t144
     #-0.1578901851851852D-1*t147
      t160 = -0.8443333333333333D-1*t5-0.1163333333333333D0*t144
     #+0.4060033333333333D-1*t147
      t164 = t149*sigma-0.5D0*(0.4690740740740741D-2*t5
     #+0.6462962962962963D-2*t144-0.2255574074074074D-2*t147)*sigma
     #-0.2777777777777778D-1*t160*sigma+0.2777777777777778D-1*t112
      t168 = 0.25D0*t26*t164-0.6666666666666667D0*t117
      t175 = t57**2
      t176 = 1/t175
      t177 = t176*t62
      t178 = 1/t26
      t181 = t61**2
      t182 = 1/t181
      t183 = t57*t182
      t186 = t59**2
      t187 = t186**2
      t188 = t187*t59
      t189 = 1/t188
      t190 = t189*t178
      t192 = -0.2067834969664667D0*t176*t178-0.1715968668597574D1*t190
      t195 = -0.2067834969664667D0*t177*t178-0.6203504908994D0*t183*t192
      t196 = 1/t57
      t197 = t195*t196
      t198 = t197*t61
      t200 = t68**2
      t201 = 1/t200
      t202 = t201*t189
      t204 = 1.D0+0.2016D-2*t201
      t205 = 1/t204
      t207 = t202*t178*t205
      t209 = t74*t62
      t212 = t75*t182
      t215 = -0.2625411059665811D0*t209*t190-1.D0*t212*t192
      t216 = 1/t75
      t217 = t215*t216
      t218 = t217*t61
      vrhoa(i) = -0.7877960174741572D0*t2+0.5080001673176129D-2*t85
     #*t14+0.1905000627441048D-2*t6*t103-0.398358D-1*t20-0.52583256D-2
     #*t25*t30*t119-0.46342314D-2*t124*t17-0.44397795816D-3*t130*t131
     #-0.6117185448D-3*t134*t129*t52+0.192805272D-1*t25*t139*t52
     #-0.52583256D-2*t25*t30*t168+0.5907233D-2*t65
     #+0.3899174858189126D1*t71+0.8419610158724123D-3*t77+0.19D0*rho*
     #(0.5011795824473985D-1*t198+0.2419143800947354D0*t207
     #+0.4431373767749538D-2*t218)
      t225 = 1/t7
      t231 = 0.6350002091470161D-1*t225*t5*t10+0.8000501301919725D-1
     #*t94*t99
      t232 = t89*t231
      t240 = 0.25D0*t26*t36-0.6666666666666667D0*t26
      vsigmaaa(i) = -0.1524000501952839D-1*t5*t14
     #+0.3810001254882096D-2*t6*t232+0.5842584D-3*t25*t146
     #-0.210333024D-1*t25*t30*t240
      t254 = 1/t2/t27
      t258 = rho*t114
      t269 = 1/t123/t19
      t271 = t127*t26
      t273 = 1/t2/t271
      t277 = 1/t271
      t287 = t84*t20
      t289 = t94*t124
      t291 = 1/t27
      t292 = t291*t269
      t311 = t160*t56*sigma
      t314 = t44*t178*sigma
      s1 = -0.177591183264D-2*t130*t20*t119-0.24468741792D-2*t134*t129
     #*t119+0.771221088D-1*t25*t139*t119-0.2370667447482193D-1*t254
     #*sigma*t14-0.52583256D-2*t25*t30*(t258+0.2552208000167282D2*t31
     #-0.5D0*sigma)+0.771221088D-1*t25*t139*t168
      s2 = s1-0.2846530295136D-3*t24*t269*t273*t52+0.106031214432D-1
     #*t134*t277*t52-0.1799515872D0*t25/t28/t128*t52
      t342 = s2-0.105166512D-1*t25*t30*(0.25D0*t26*((
     #-0.4378024691358025D-1*t84-0.6032098765432099D-1*t287
     #+0.3157803703703704D-1*t289-0.3673578308641975D-2*t292)*sigma
     #-0.5D0*(-0.6254320987654321D-2*t84-0.8617283950617284D-2*t287
     #+0.4511148148148148D-2*t289-0.5247969012345679D-3*t292)*sigma
     #-0.2777777777777778D-1*(0.1125777777777778D0*t84
     #+0.1551111111111111D0*t287-0.8120066666666667D-1*t289
     #+0.9446344222222222D-2*t292)*sigma+0.5555555555555556D-1*t311
     #-0.5555555555555556D-1*t314)-0.6666666666666667D0*sigma)
     #+0.769561794144D-2*t277*t24*t131-0.52583256D-2*t25*t30*(t32+t37
     #-t42-t46+t258)-0.210333024D-1*t25*t30*(0.5D0*rho*t164+0.25D0*t26
     #*(-0.2777777777777778D-1*t311+0.2777777777777778D-1*t314))
      t346 = t273*t24
      t363 = sigma**2
      t369 = 1/t98/t97
      t377 = 1/t88/t13
      t378 = t102**2
      t392 = 1/t175/t56
      t394 = 1/t127
      t395 = t392*t62*t394
      t404 = 1/t181/t61
      t406 = t192**2
      t414 = 1/t188/t56
      t415 = t414*t394
      t417 = t189*t291
      t419 = -0.1378556646443111D0*t392*t394+0.4135669939329333D0*t176
     #*t291-0.1429973890497978D1*t415+0.3431937337195148D1*t417
      t437 = t394*t205
      t446 = t200**2
      t450 = t204**2
      s1 = -0.177591183264D-2*t130*t20*t168-0.74973077867952D-4*t346
     #*t131-0.20659774319712D-3*t346*t124*t52+0.1838549288719989D0
     #*t207+0.3808964826600229D-1*t198+0.3367844063489649D-2*t218
     #-0.1016000334635226D-1*t85*t103
      s2 = s1+0.1905000627441048D-2*t6*t89*(0.3951112412470322D0*t7
     #*t254*t10+0.106673350692263D1*sigma*t139*t99
     #-0.4515557042823225D0*t363/t2/t127/t27*t369)
     #-0.3810001254882096D-2*t6*t377*t378-0.24468741792D-2*t134*t129
     #*t168
      s3 = s2-0.5251973449827715D0/t28
      s4 = s3-0.61789752D-2*t124*t5
      s5 = s4
      s7 = -0.21564623448D-2*t269*t146
      s8 = 0.38D0*rho*(0.5011795824473985D-1*(-0.1378556646443111D0
     #*t395+0.4135669939329333D0*t176*t182*t178*t192
     #+0.4135669939329333D0*t177*t291+0.12407009817988D1*t57*t404*t406
     #-0.6203504908994D0*t183*t419)*t196*t61+0.1670598608157995D-1
     #*t195/t57/t56*t61*t178+0.5011795824473985D-1*t197*t192
     #+0.1270249377985834D0/t200/t68*t392*t437+0.2015953167456128D0
     #*t201*t414*t437-0.4838287601894708D0*t202*t291*t205
     #-0.2560822746019441D-3/t446/t68*t392*t394/t450
     #+0.4431373767749538D-2*(0.3446391616107778D-1*t395
     #+0.5250822119331622D0*t74*t182*t190*t192-0.2187842549721509D0
     #*t209*t415+0.5250822119331622D0*t209*t417+2.D0*t75*t404*t406
     #-1.D0*t212*t419)*t216*t61+0.1163417769936259D-2*t215/t75/t74*t61
     #*t189*t178+0.4431373767749538D-2*t217*t192)
      s6 = s7+s8
      t485 = s5+s6
      v2rhoa2(i) = t342+t485
      t490 = t5*t89
      s1 = 0.2032000669270451D-1*t84*t14-0.5080001673176129D-2*t85
     #*t232+0.7620002509764193D-2*t490*t102-0.3810001254882096D-2*t6
     #*t377*t102*t231+0.1905000627441048D-2*t6*t89*(
     #-0.169333389105871D0*t225*t84*t10-0.640040104153578D0*t30*t99
     #+0.3386667782117419D0*sigma*t273*t369)-0.52583256D-2*t25*t30*(
     #-0.9444444444444444D0*rho-0.2777777777777778D-1*rho*t44)
     #+0.4933088424D-4*t291*t24*t20
      v2rhoasigmaaa(i) = s1+0.679687272D-4*t134*t291+0.80822412D-2*t25
     #*t94-0.105166512D-1*t25*t30*(0.25D0*t26*(-0.3D-22*t144
     #+0.5555555555555556D-1*t111)+0.1333333333333333D1*rho)
     #-0.105166512D-1*t25*t94*t36-0.177591183264D-2*t130*t20*t240
     #-0.24468741792D-2*t134*t129*t240+0.771221088D-1*t25*t139*t240
     #-0.210333024D-1*t25*t30*(0.25D0*t26*t149-0.1333333333333333D1*rho)
      t554 = t231**2
      v2sigmaaa2(i) = 0.3048001003905677D-1*t490*t231
     #-0.7620002509764193D-2*t6*t377*t554+0.3810001254882096D-2*t6*t89
     #*(-0.1270000418294032D0/t7/sigma*t5*t10+0.1600100260383945D0
     #/sigma*t94*t99-0.2540000836588064D0/t2/t128*t369)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end

c:XC_B3LYPsubrend
