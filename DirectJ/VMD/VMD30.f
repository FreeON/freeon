        SUBROUTINE VMD30(NQ,Px,Py,Pz,Omega,Upq,Tol,                      
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
       REAL*8 F3_0(0:Mesh),F3_1(0:Mesh),F3_2(0:Mesh),F3_3(0:Mesh)        
       INCLUDE "Gamma_0.Inc"                                                            
       INCLUDE "Gamma_1.Inc"                                                            
       INCLUDE "Gamma_2.Inc"                                                            
       INCLUDE "Gamma_3.Inc"                                                            
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
             o1=o2*o1                                                    
             Aux1=o1*(F1_0(j)+T*F1_1(j)+T2*F1_2(j)+T3*F1_3(j))           
             o1=o2*o1                                                    
             Aux2=o1*(F2_0(j)+T*F2_1(j)+T2*F2_2(j)+T3*F2_3(j))           
             o1=o2*o1                                                    
             Aux3=o1*(F3_0(j)+T*F3_1(j)+T2*F3_2(j)+T3*F3_3(j))           
          ELSE                                                           
             t1=1.0D0/T                                                  
             t1x2=Upq*DSQRT(T1)                                          
             t1=o2*t1                                                    
             Aux0=F0Asymp*t1x2                                           
             t1x2=t1x2*t1                                                
             Aux1=F1Asymp*t1x2                                           
             t1x2=t1x2*t1                                                
             Aux2=F2Asymp*t1x2                                           
             t1x2=t1x2*t1                                                
             Aux3=F3Asymp*t1x2                                           
          ENDIF                                                          
          r1=Aux3                                                        
          r2=QPx*r1                                                      
          r3=QPy*r1                                                      
          r4=QPz*r1                                                      
          r1=Aux2                                                        
          r5=QPx*r2+r1                                                   
          r7=QPy*r3+r1                                                   
          r10=QPz*r4+r1                                                  
          r2=QPx*r1                                                      
          r6=QPx*r3                                                      
          r9=QPy*r4                                                      
          r8=QPx*r4                                                      
          r3=QPy*r1                                                      
          r4=QPz*r1                                                      
          r1=Aux1                                                        
          r11=QPx*r5+.20000D+01*r2                                       
          r14=QPy*r7+.20000D+01*r3                                       
          r20=QPz*r10+.20000D+01*r4                                      
          r5=QPx*r2+r1                                                   
          r13=QPx*r7                                                     
          r19=QPy*r10                                                    
          r18=QPx*r10                                                    
          r7=QPy*r3+r1                                                   
          r10=QPz*r4+r1                                                  
          r2=QPx*r1                                                      
          r12=QPx*r6+r3                                                  
          r17=QPy*r9+r4                                                  
          r15=QPx*r8+r4                                                  
          r6=QPx*r3                                                      
          r16=QPx*r9                                                     
          r9=QPy*r4                                                      
          r8=QPx*r4                                                      
          r3=QPy*r1                                                      
          r4=QPz*r1                                                      
          r1=Aux0                                                        
          HGKet(1)=HGKet(1)+r1*RhoCo(iq)                                 
          HGKet(2)=HGKet(2)+r2*RhoCo(iq)                                 
          HGKet(5)=HGKet(5)+r5*RhoCo(iq)                                 
          HGKet(11)=HGKet(11)+r11*RhoCo(iq)                              
          HGKet(3)=HGKet(3)+r3*RhoCo(iq)                                 
          HGKet(6)=HGKet(6)+r6*RhoCo(iq)                                 
          HGKet(12)=HGKet(12)+r12*RhoCo(iq)                              
          HGKet(7)=HGKet(7)+r7*RhoCo(iq)                                 
          HGKet(13)=HGKet(13)+r13*RhoCo(iq)                              
          HGKet(14)=HGKet(14)+r14*RhoCo(iq)                              
          HGKet(4)=HGKet(4)+r4*RhoCo(iq)                                 
          HGKet(8)=HGKet(8)+r8*RhoCo(iq)                                 
          HGKet(15)=HGKet(15)+r15*RhoCo(iq)                              
          HGKet(9)=HGKet(9)+r9*RhoCo(iq)                                 
          HGKet(16)=HGKet(16)+r16*RhoCo(iq)                              
          HGKet(17)=HGKet(17)+r17*RhoCo(iq)                              
          HGKet(10)=HGKet(10)+r10*RhoCo(iq)                              
          HGKet(18)=HGKet(18)+r18*RhoCo(iq)                              
          HGKet(19)=HGKet(19)+r19*RhoCo(iq)                              
          HGKet(20)=HGKet(20)+r20*RhoCo(iq)                              
  100  CONTINUE
       RETURN
       END
  
