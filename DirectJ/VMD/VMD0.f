        SUBROUTINE VMD0(NQ,Px,Py,Pz,Omega,Upq,Tol,                       
     $          rQx,rQy,rQz,rQt,RhoCo,HGKet)
       IMPLICIT REAL*8(a-h,o-z)
       IMPLICIT INTEGER(i-n)
       REAL*8 rQx(*)
       REAL*8 rQy(*)
       REAL*8 rQz(*)
       REAL*8 rQt(*)
       REAL*8 RhoCo(*)
       REAL*8 HGKet(*)
       INCLUDE "Mesh.Inc"
       REAL*8 F0_0(0:Mesh),F0_1(0:Mesh),F0_2(0:Mesh),F0_3(0:Mesh)        
       INCLUDE "Gamma_0.Inc"                                                            
       INCLUDE "Gamma_Asymptotics.Inc"
       IF(rQt(1).LT.Tol)THEN
          RETURN
       ELSEIF(rQt(NQ).GT.Tol)THEN
          N=NQ
       ELSE
          K=0
          L=NQ
          N=NQ/2
  102     CONTINUE
          IF(rQt(N).LT.Tol)THEN
             L=N
          ELSE
             K=N
          ENDIF
          N=(K+L)/2
          IF((K.LT.N).AND.(N.LT.L))GOTO 102
       ENDIF
       DO 100 iq=1,N
          QPx=rQx(iq)-Px
          QPy=rQy(iq)-Py
          QPz=rQz(iq)-Pz
          T=Omega*(QPx*QPx+QPy*QPy+QPz*QPz)
          o1=Upq
          o2=-2.0D0*Omega
          IF(T.LT.Switch)THEN
             j=AINT(T*Grid)                                              
             T2=T*T                                                      
             T3=T2*T                                                     
             o1=upq                                                      
             Aux0=o1*(F0_0(j)+T*F0_1(j)+T2*F0_2(j)+T3*F0_3(j))           
          ELSE                                                           
             t1=1.0D0/T                                                  
             t1x2=Upq*DSQRT(T1)                                          
             Aux0=F0Asymp*t1x2                                           
          ENDIF                                                          
          r1=Aux0                                                        
          HGKet(1)=HGKet(1)+r1*RhoCo(iq)                                 
!          WRITE(*,33)omega,rQz(iq),RhoCo(iq),T,r1*RhoCo(iq)
!33        format('Z = ',D12.6,' Q = ',D14.6,' Co = ',D14.6,
!     > ' T = ',D14.6,' Ket = ',D14.6)

  100  CONTINUE
       RETURN
       END
  
