        SUBROUTINE VMD40(NQ,Px,Py,Pz,Omega,Upq,Tol,                      
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
       REAL*8 F4_0(0:Mesh),F4_1(0:Mesh),F4_2(0:Mesh),F4_3(0:Mesh)        
       INCLUDE "Gamma_4.Inc"                                                            
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
             G4=(F4_0(j)+T*F4_1(j)+T2*F4_2(j)+T3*F4_3(j))                
             ET=DEXP(-T)
             TwoT=2.0D0*T
             G3=0.142857142857143D+00*(TwoT*G4+ET)                       
             G2=0.200000000000000D+00*(TwoT*G3+ET)                       
             G1=0.333333333333333D+00*(TwoT*G2+ET)                       
             G0=TwoT*G1+ET                                               
             Aux0=o1*G0                                                  
             o1=o2*o1
             Aux1=o1*G1                                                  
             o1=o2*o1
             Aux2=o1*G2                                                  
             o1=o2*o1
             Aux3=o1*G3                                                  
             o1=o2*o1
             Aux4=o1*G4                                                  
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
             t1x2=t1x2*t1                                                
             Aux4=F4Asymp*t1x2                                           
          ENDIF                                                          
          r1=Aux4                                                        
          r2=QPx*r1                                                      
          r3=QPy*r1                                                      
          r4=QPz*r1                                                      
          r1=Aux3                                                        
          r5=QPx*r2+r1                                                   
          r7=QPy*r3+r1                                                   
          r10=QPz*r4+r1                                                  
          r2=QPx*r1                                                      
          r6=QPx*r3                                                      
          r9=QPy*r4                                                      
          r8=QPx*r4                                                      
          r3=QPy*r1                                                      
          r4=QPz*r1                                                      
          r1=Aux2                                                        
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
          r1=Aux1                                                        
          r21=QPx*r11+.30000D+01*r5                                      
          r25=QPy*r14+.30000D+01*r7                                      
          r35=QPz*r20+.30000D+01*r10                                     
          r11=QPx*r5+.20000D+01*r2                                       
          r24=QPx*r14                                                    
          r34=QPy*r20                                                    
          r33=QPx*r20                                                    
          r14=QPy*r7+.20000D+01*r3                                       
          r20=QPz*r10+.20000D+01*r4                                      
          r5=QPx*r2+r1                                                   
          r23=QPx*r13+r7                                                 
          r32=QPy*r19+r10                                                
          r30=QPx*r18+r10                                                
          r13=QPx*r7                                                     
          r31=QPx*r19                                                    
          r19=QPy*r10                                                    
          r18=QPx*r10                                                    
          r7=QPy*r3+r1                                                   
          r10=QPz*r4+r1                                                  
          r2=QPx*r1                                                      
          r22=QPx*r12+.20000D+01*r6                                      
          r29=QPy*r17+.20000D+01*r9                                      
          r26=QPx*r15+.20000D+01*r8                                      
          r12=QPx*r6+r3                                                  
          r28=QPx*r17                                                    
          r17=QPy*r9+r4                                                  
          r15=QPx*r8+r4                                                  
          r6=QPx*r3                                                      
          r27=QPx*r16+r9                                                 
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
          HGKet(21)=HGKet(21)+r21*RhoCo(iq)                              
          HGKet(3)=HGKet(3)+r3*RhoCo(iq)                                 
          HGKet(6)=HGKet(6)+r6*RhoCo(iq)                                 
          HGKet(12)=HGKet(12)+r12*RhoCo(iq)                              
          HGKet(22)=HGKet(22)+r22*RhoCo(iq)                              
          HGKet(7)=HGKet(7)+r7*RhoCo(iq)                                 
          HGKet(13)=HGKet(13)+r13*RhoCo(iq)                              
          HGKet(23)=HGKet(23)+r23*RhoCo(iq)                              
          HGKet(14)=HGKet(14)+r14*RhoCo(iq)                              
          HGKet(24)=HGKet(24)+r24*RhoCo(iq)                              
          HGKet(25)=HGKet(25)+r25*RhoCo(iq)                              
          HGKet(4)=HGKet(4)+r4*RhoCo(iq)                                 
          HGKet(8)=HGKet(8)+r8*RhoCo(iq)                                 
          HGKet(15)=HGKet(15)+r15*RhoCo(iq)                              
          HGKet(26)=HGKet(26)+r26*RhoCo(iq)                              
          HGKet(9)=HGKet(9)+r9*RhoCo(iq)                                 
          HGKet(16)=HGKet(16)+r16*RhoCo(iq)                              
          HGKet(27)=HGKet(27)+r27*RhoCo(iq)                              
          HGKet(17)=HGKet(17)+r17*RhoCo(iq)                              
          HGKet(28)=HGKet(28)+r28*RhoCo(iq)                              
          HGKet(29)=HGKet(29)+r29*RhoCo(iq)                              
          HGKet(10)=HGKet(10)+r10*RhoCo(iq)                              
          HGKet(18)=HGKet(18)+r18*RhoCo(iq)                              
          HGKet(30)=HGKet(30)+r30*RhoCo(iq)                              
          HGKet(19)=HGKet(19)+r19*RhoCo(iq)                              
          HGKet(31)=HGKet(31)+r31*RhoCo(iq)                              
          HGKet(32)=HGKet(32)+r32*RhoCo(iq)                              
          HGKet(20)=HGKet(20)+r20*RhoCo(iq)                              
          HGKet(33)=HGKet(33)+r33*RhoCo(iq)                              
          HGKet(34)=HGKet(34)+r34*RhoCo(iq)                              
          HGKet(35)=HGKet(35)+r35*RhoCo(iq)                              
  100  CONTINUE
       RETURN
       END
  
