        SUBROUTINE VMD11(NQ,Px,Py,Pz,Omega,Upq,Tol,                      
     $          rQx,rQy,rQz,rQt,RhoCo,HGKet)
       IMPLICIT REAL*8(a-h,o-z)
       IMPLICIT INTEGER(i-n)
       REAL*8 rQx(*)
       REAL*8 rQy(*)
       REAL*8 rQz(*)
       REAL*8 rQt(*)
       REAL*8 RhoCo(*)
       REAL*8 HGKet(*)
       PARAMETER (Mesh=968)
       REAL*8 F0_0(0:Mesh),F0_1(0:Mesh),F0_2(0:Mesh),F0_3(0:Mesh)        
       REAL*8 F1_0(0:Mesh),F1_1(0:Mesh),F1_2(0:Mesh),F1_3(0:Mesh)        
       REAL*8 F2_0(0:Mesh),F2_1(0:Mesh),F2_2(0:Mesh),F2_3(0:Mesh)        
       INCLUDE "Gamma_0.Inc"                                                            
       INCLUDE "Gamma_1.Inc"                                                            
       INCLUDE "Gamma_2.Inc"                                                            
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
          jadd=(iq-1)*4                                                        
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
             o1=o2*o1                                                    
             Aux1=o1*(F1_0(j)+T*F1_1(j)+T2*F1_2(j)+T3*F1_3(j))           
             o1=o2*o1                                                    
             Aux2=o1*(F2_0(j)+T*F2_1(j)+T2*F2_2(j)+T3*F2_3(j))           
          ELSE                                                           
             t1=1.0D0/T                                                  
             t1x2=Upq*DSQRT(T1)                                          
             t1=o2*t1                                                    
             Aux0=F0Asymp*t1x2                                           
             t1x2=t1x2*t1                                                
             Aux1=F1Asymp*t1x2                                           
             t1x2=t1x2*t1                                                
             Aux2=F2Asymp*t1x2                                           
          ENDIF                                                          
          r1=Aux2                                                        
          r2=QPx*r1                                                      
          r3=QPy*r1                                                      
          r4=QPz*r1                                                      
          r1=Aux1                                                        
          r5=QPx*r2+r1                                                   
          r7=QPy*r3+r1                                                   
          r10=QPz*r4+r1                                                  
          r2=QPx*r1                                                      
          r6=QPx*r3                                                      
          r9=QPy*r4                                                      
          r8=QPx*r4                                                      
          r3=QPy*r1                                                      
          r4=QPz*r1                                                      
          r1=Aux0                                                        
          HGKet(1)=HGKet(1)+r1*RhoCo(jadd+1)                             
          HGKet(1)=HGKet(1)+r2*RhoCo(jadd+2)                             
          HGKet(1)=HGKet(1)+r3*RhoCo(jadd+3)                             
          HGKet(1)=HGKet(1)+r4*RhoCo(jadd+4)                             
          HGKet(2)=HGKet(2)+r2*RhoCo(jadd+1)                             
          HGKet(2)=HGKet(2)+r5*RhoCo(jadd+2)                             
          HGKet(2)=HGKet(2)+r6*RhoCo(jadd+3)                             
          HGKet(2)=HGKet(2)+r8*RhoCo(jadd+4)                             
          HGKet(3)=HGKet(3)+r3*RhoCo(jadd+1)                             
          HGKet(3)=HGKet(3)+r6*RhoCo(jadd+2)                             
          HGKet(3)=HGKet(3)+r7*RhoCo(jadd+3)                             
          HGKet(3)=HGKet(3)+r9*RhoCo(jadd+4)                             
          HGKet(4)=HGKet(4)+r4*RhoCo(jadd+1)                             
          HGKet(4)=HGKet(4)+r8*RhoCo(jadd+2)                             
          HGKet(4)=HGKet(4)+r9*RhoCo(jadd+3)                             
          HGKet(4)=HGKet(4)+r10*RhoCo(jadd+4)                            
  100  CONTINUE
       RETURN
       END
  
