      SUBROUTINE AuxIGen0(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:0,*)                                
      REAL*8 ET,TwoT,T,T2,T3,T4,RTE,RPE,Upq,TwoPi5x2                    
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(SqrtPi=1.7724538509055160273D0)                         
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)                       
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      INCLUDE "Gamma_Asymptotics.Inc"                                   
      INCLUDE "Mesh.Inc"                                                
      REAL*8 F0_0(0:Mesh)                                               
      REAL*8 F0_1(0:Mesh)                                               
      REAL*8 F0_2(0:Mesh)                                               
      REAL*8 F0_3(0:Mesh)                                               
      INCLUDE "Gamma_0.Inc"                                             
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=TwoPi5x2/(RTE*SQRT(RPE))                                   
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         IF(T.LT.Switch)THEN                                            
            T2=T*T                                                      
            T3=T*T2                                                     
            J=AINT(T*Grid)                                              
            G0=(F0_0(J)+T*F0_1(J)+T2*F0_2(J)+T3*F0_3(J))                
            o1=Upq                                                      
            AuxR(0,I)=o1*G0                                             
         ELSE                                                           
            SqrtT=DSQRT(T)                                              
            OneOvT=One/T                                                
            G0=SqrtPi/(Two*SqrtT)                                       
            o1=Upq                                                      
            AuxR(0,I)=o1*G0                                             
         ENDIF                                                          
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen1(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:1,*)                                
      REAL*8 ET,TwoT,T,T2,T3,T4,RTE,RPE,Upq,TwoPi5x2                    
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(SqrtPi=1.7724538509055160273D0)                         
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)                       
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      INCLUDE "Gamma_Asymptotics.Inc"                                   
      INCLUDE "Mesh.Inc"                                                
      REAL*8 F1_0(0:Mesh)                                               
      REAL*8 F1_1(0:Mesh)                                               
      REAL*8 F1_2(0:Mesh)                                               
      REAL*8 F1_3(0:Mesh)                                               
      INCLUDE "Gamma_1.Inc"                                             
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=TwoPi5x2/(RTE*SQRT(RPE))                                   
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         IF(T.LT.Switch)THEN                                            
            T2=T*T                                                      
            T3=T*T2                                                     
            J=AINT(T*Grid)                                              
            G1=(F1_0(J)+T*F1_1(J)+T2*F1_2(J)+T3*F1_3(J))                
            ET=DEXP(-T)                                                 
            TwoT=2.0D0*T                                                
            G0=TwoT*G1+ET                                               
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            AuxR(1,I)=o1*G1                                             
         ELSE                                                           
            SqrtT=DSQRT(T)                                              
            OneOvT=One/T                                                
            G0=SqrtPi/(Two*SqrtT)                                       
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            G1=G0*0.5000000000000000D+00*OneOvT                         
            AuxR(1,I)=o1*G1                                             
         ENDIF                                                          
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen2(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:2,*)                                
      REAL*8 ET,TwoT,T,T2,T3,T4,RTE,RPE,Upq,TwoPi5x2                    
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(SqrtPi=1.7724538509055160273D0)                         
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)                       
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      INCLUDE "Gamma_Asymptotics.Inc"                                   
      INCLUDE "Mesh.Inc"                                                
      REAL*8 F2_0(0:Mesh)                                               
      REAL*8 F2_1(0:Mesh)                                               
      REAL*8 F2_2(0:Mesh)                                               
      REAL*8 F2_3(0:Mesh)                                               
      INCLUDE "Gamma_2.Inc"                                             
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=TwoPi5x2/(RTE*SQRT(RPE))                                   
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         IF(T.LT.Switch)THEN                                            
            T2=T*T                                                      
            T3=T*T2                                                     
            J=AINT(T*Grid)                                              
            G2=(F2_0(J)+T*F2_1(J)+T2*F2_2(J)+T3*F2_3(J))                
            ET=DEXP(-T)                                                 
            TwoT=2.0D0*T                                                
            G1=0.3333333333333333D+00*(TwoT*G2+ET)                      
            G0=TwoT*G1+ET                                               
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            AuxR(2,I)=o1*G2                                             
         ELSE                                                           
            SqrtT=DSQRT(T)                                              
            OneOvT=One/T                                                
            G0=SqrtPi/(Two*SqrtT)                                       
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            G1=G0*0.5000000000000000D+00*OneOvT                         
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            G2=G1*0.1500000000000000D+01*OneOvT                         
            AuxR(2,I)=o1*G2                                             
         ENDIF                                                          
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen3(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:3,*)                                
      REAL*8 ET,TwoT,T,T2,T3,T4,RTE,RPE,Upq,TwoPi5x2                    
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(SqrtPi=1.7724538509055160273D0)                         
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)                       
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      INCLUDE "Gamma_Asymptotics.Inc"                                   
      INCLUDE "Mesh.Inc"                                                
      REAL*8 F3_0(0:Mesh)                                               
      REAL*8 F3_1(0:Mesh)                                               
      REAL*8 F3_2(0:Mesh)                                               
      REAL*8 F3_3(0:Mesh)                                               
      INCLUDE "Gamma_3.Inc"                                             
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=TwoPi5x2/(RTE*SQRT(RPE))                                   
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         IF(T.LT.Switch)THEN                                            
            T2=T*T                                                      
            T3=T*T2                                                     
            J=AINT(T*Grid)                                              
            G3=(F3_0(J)+T*F3_1(J)+T2*F3_2(J)+T3*F3_3(J))                
            ET=DEXP(-T)                                                 
            TwoT=2.0D0*T                                                
            G2=0.2000000000000000D+00*(TwoT*G3+ET)                      
            G1=0.3333333333333333D+00*(TwoT*G2+ET)                      
            G0=TwoT*G1+ET                                               
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            AuxR(3,I)=o1*G3                                             
         ELSE                                                           
            SqrtT=DSQRT(T)                                              
            OneOvT=One/T                                                
            G0=SqrtPi/(Two*SqrtT)                                       
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            G1=G0*0.5000000000000000D+00*OneOvT                         
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            G2=G1*0.1500000000000000D+01*OneOvT                         
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            G3=G2*0.2500000000000000D+01*OneOvT                         
            AuxR(3,I)=o1*G3                                             
         ENDIF                                                          
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen4(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:4,*)                                
      REAL*8 ET,TwoT,T,T2,T3,T4,RTE,RPE,Upq,TwoPi5x2                    
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(SqrtPi=1.7724538509055160273D0)                         
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)                       
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      INCLUDE "Gamma_Asymptotics.Inc"                                   
      INCLUDE "Mesh.Inc"                                                
      REAL*8 F4_0(0:Mesh)                                               
      REAL*8 F4_1(0:Mesh)                                               
      REAL*8 F4_2(0:Mesh)                                               
      REAL*8 F4_3(0:Mesh)                                               
      INCLUDE "Gamma_4.Inc"                                             
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=TwoPi5x2/(RTE*SQRT(RPE))                                   
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         IF(T.LT.Switch)THEN                                            
            T2=T*T                                                      
            T3=T*T2                                                     
            J=AINT(T*Grid)                                              
            G4=(F4_0(J)+T*F4_1(J)+T2*F4_2(J)+T3*F4_3(J))                
            ET=DEXP(-T)                                                 
            TwoT=2.0D0*T                                                
            G3=0.1428571428571428D+00*(TwoT*G4+ET)                      
            G2=0.2000000000000000D+00*(TwoT*G3+ET)                      
            G1=0.3333333333333333D+00*(TwoT*G2+ET)                      
            G0=TwoT*G1+ET                                               
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            AuxR(4,I)=o1*G4                                             
         ELSE                                                           
            SqrtT=DSQRT(T)                                              
            OneOvT=One/T                                                
            G0=SqrtPi/(Two*SqrtT)                                       
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            G1=G0*0.5000000000000000D+00*OneOvT                         
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            G2=G1*0.1500000000000000D+01*OneOvT                         
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            G3=G2*0.2500000000000000D+01*OneOvT                         
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            G4=G3*0.3500000000000000D+01*OneOvT                         
            AuxR(4,I)=o1*G4                                             
         ENDIF                                                          
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen5(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:5,*)                                
      REAL*8 ET,TwoT,T,T2,T3,T4,RTE,RPE,Upq,TwoPi5x2                    
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(SqrtPi=1.7724538509055160273D0)                         
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)                       
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      INCLUDE "Gamma_Asymptotics.Inc"                                   
      INCLUDE "Mesh.Inc"                                                
      REAL*8 F5_0(0:Mesh)                                               
      REAL*8 F5_1(0:Mesh)                                               
      REAL*8 F5_2(0:Mesh)                                               
      REAL*8 F5_3(0:Mesh)                                               
      INCLUDE "Gamma_5.Inc"                                             
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=TwoPi5x2/(RTE*SQRT(RPE))                                   
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         IF(T.LT.Switch)THEN                                            
            T2=T*T                                                      
            T3=T*T2                                                     
            J=AINT(T*Grid)                                              
            G5=(F5_0(J)+T*F5_1(J)+T2*F5_2(J)+T3*F5_3(J))                
            ET=DEXP(-T)                                                 
            TwoT=2.0D0*T                                                
            G4=0.1111111111111111D+00*(TwoT*G5+ET)                      
            G3=0.1428571428571428D+00*(TwoT*G4+ET)                      
            G2=0.2000000000000000D+00*(TwoT*G3+ET)                      
            G1=0.3333333333333333D+00*(TwoT*G2+ET)                      
            G0=TwoT*G1+ET                                               
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            AuxR(5,I)=o1*G5                                             
         ELSE                                                           
            SqrtT=DSQRT(T)                                              
            OneOvT=One/T                                                
            G0=SqrtPi/(Two*SqrtT)                                       
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            G1=G0*0.5000000000000000D+00*OneOvT                         
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            G2=G1*0.1500000000000000D+01*OneOvT                         
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            G3=G2*0.2500000000000000D+01*OneOvT                         
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            G4=G3*0.3500000000000000D+01*OneOvT                         
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            G5=G4*0.4500000000000000D+01*OneOvT                         
            AuxR(5,I)=o1*G5                                             
         ENDIF                                                          
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen6(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:6,*)                                
      REAL*8 ET,TwoT,T,T2,T3,T4,RTE,RPE,Upq,TwoPi5x2                    
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(SqrtPi=1.7724538509055160273D0)                         
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)                       
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      INCLUDE "Gamma_Asymptotics.Inc"                                   
      INCLUDE "Mesh.Inc"                                                
      REAL*8 F6_0(0:Mesh)                                               
      REAL*8 F6_1(0:Mesh)                                               
      REAL*8 F6_2(0:Mesh)                                               
      REAL*8 F6_3(0:Mesh)                                               
      INCLUDE "Gamma_6.Inc"                                             
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=TwoPi5x2/(RTE*SQRT(RPE))                                   
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         IF(T.LT.Switch)THEN                                            
            T2=T*T                                                      
            T3=T*T2                                                     
            J=AINT(T*Grid)                                              
            G6=(F6_0(J)+T*F6_1(J)+T2*F6_2(J)+T3*F6_3(J))                
            ET=DEXP(-T)                                                 
            TwoT=2.0D0*T                                                
            G5=0.9090909090909091D-01*(TwoT*G6+ET)                      
            G4=0.1111111111111111D+00*(TwoT*G5+ET)                      
            G3=0.1428571428571428D+00*(TwoT*G4+ET)                      
            G2=0.2000000000000000D+00*(TwoT*G3+ET)                      
            G1=0.3333333333333333D+00*(TwoT*G2+ET)                      
            G0=TwoT*G1+ET                                               
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            AuxR(5,I)=o1*G5                                             
            o1=o2*o1                                                    
            AuxR(6,I)=o1*G6                                             
         ELSE                                                           
            SqrtT=DSQRT(T)                                              
            OneOvT=One/T                                                
            G0=SqrtPi/(Two*SqrtT)                                       
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            G1=G0*0.5000000000000000D+00*OneOvT                         
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            G2=G1*0.1500000000000000D+01*OneOvT                         
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            G3=G2*0.2500000000000000D+01*OneOvT                         
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            G4=G3*0.3500000000000000D+01*OneOvT                         
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            G5=G4*0.4500000000000000D+01*OneOvT                         
            AuxR(5,I)=o1*G5                                             
            o1=o2*o1                                                    
            G6=G5*0.5500000000000000D+01*OneOvT                         
            AuxR(6,I)=o1*G6                                             
         ENDIF                                                          
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen7(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:7,*)                                
      REAL*8 ET,TwoT,T,T2,T3,T4,RTE,RPE,Upq,TwoPi5x2                    
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(SqrtPi=1.7724538509055160273D0)                         
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)                       
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      INCLUDE "Gamma_Asymptotics.Inc"                                   
      INCLUDE "Mesh.Inc"                                                
      REAL*8 F7_0(0:Mesh)                                               
      REAL*8 F7_1(0:Mesh)                                               
      REAL*8 F7_2(0:Mesh)                                               
      REAL*8 F7_3(0:Mesh)                                               
      INCLUDE "Gamma_7.Inc"                                             
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=TwoPi5x2/(RTE*SQRT(RPE))                                   
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         IF(T.LT.Switch)THEN                                            
            T2=T*T                                                      
            T3=T*T2                                                     
            J=AINT(T*Grid)                                              
            G7=(F7_0(J)+T*F7_1(J)+T2*F7_2(J)+T3*F7_3(J))                
            ET=DEXP(-T)                                                 
            TwoT=2.0D0*T                                                
            G6=0.7692307692307693D-01*(TwoT*G7+ET)                      
            G5=0.9090909090909091D-01*(TwoT*G6+ET)                      
            G4=0.1111111111111111D+00*(TwoT*G5+ET)                      
            G3=0.1428571428571428D+00*(TwoT*G4+ET)                      
            G2=0.2000000000000000D+00*(TwoT*G3+ET)                      
            G1=0.3333333333333333D+00*(TwoT*G2+ET)                      
            G0=TwoT*G1+ET                                               
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            AuxR(5,I)=o1*G5                                             
            o1=o2*o1                                                    
            AuxR(6,I)=o1*G6                                             
            o1=o2*o1                                                    
            AuxR(7,I)=o1*G7                                             
         ELSE                                                           
            SqrtT=DSQRT(T)                                              
            OneOvT=One/T                                                
            G0=SqrtPi/(Two*SqrtT)                                       
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            G1=G0*0.5000000000000000D+00*OneOvT                         
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            G2=G1*0.1500000000000000D+01*OneOvT                         
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            G3=G2*0.2500000000000000D+01*OneOvT                         
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            G4=G3*0.3500000000000000D+01*OneOvT                         
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            G5=G4*0.4500000000000000D+01*OneOvT                         
            AuxR(5,I)=o1*G5                                             
            o1=o2*o1                                                    
            G6=G5*0.5500000000000000D+01*OneOvT                         
            AuxR(6,I)=o1*G6                                             
            o1=o2*o1                                                    
            G7=G6*0.6500000000000000D+01*OneOvT                         
            AuxR(7,I)=o1*G7                                             
         ENDIF                                                          
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen8(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:8,*)                                
      REAL*8 ET,TwoT,T,T2,T3,T4,RTE,RPE,Upq,TwoPi5x2                    
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(SqrtPi=1.7724538509055160273D0)                         
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)                       
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      INCLUDE "Gamma_Asymptotics.Inc"                                   
      INCLUDE "Mesh.Inc"                                                
      REAL*8 F8_0(0:Mesh)                                               
      REAL*8 F8_1(0:Mesh)                                               
      REAL*8 F8_2(0:Mesh)                                               
      REAL*8 F8_3(0:Mesh)                                               
      INCLUDE "Gamma_8.Inc"                                             
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=TwoPi5x2/(RTE*SQRT(RPE))                                   
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         IF(T.LT.Switch)THEN                                            
            T2=T*T                                                      
            T3=T*T2                                                     
            J=AINT(T*Grid)                                              
            G8=(F8_0(J)+T*F8_1(J)+T2*F8_2(J)+T3*F8_3(J))                
            ET=DEXP(-T)                                                 
            TwoT=2.0D0*T                                                
            G7=0.6666666666666667D-01*(TwoT*G8+ET)                      
            G6=0.7692307692307693D-01*(TwoT*G7+ET)                      
            G5=0.9090909090909091D-01*(TwoT*G6+ET)                      
            G4=0.1111111111111111D+00*(TwoT*G5+ET)                      
            G3=0.1428571428571428D+00*(TwoT*G4+ET)                      
            G2=0.2000000000000000D+00*(TwoT*G3+ET)                      
            G1=0.3333333333333333D+00*(TwoT*G2+ET)                      
            G0=TwoT*G1+ET                                               
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            AuxR(5,I)=o1*G5                                             
            o1=o2*o1                                                    
            AuxR(6,I)=o1*G6                                             
            o1=o2*o1                                                    
            AuxR(7,I)=o1*G7                                             
            o1=o2*o1                                                    
            AuxR(8,I)=o1*G8                                             
         ELSE                                                           
            SqrtT=DSQRT(T)                                              
            OneOvT=One/T                                                
            G0=SqrtPi/(Two*SqrtT)                                       
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            G1=G0*0.5000000000000000D+00*OneOvT                         
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            G2=G1*0.1500000000000000D+01*OneOvT                         
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            G3=G2*0.2500000000000000D+01*OneOvT                         
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            G4=G3*0.3500000000000000D+01*OneOvT                         
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            G5=G4*0.4500000000000000D+01*OneOvT                         
            AuxR(5,I)=o1*G5                                             
            o1=o2*o1                                                    
            G6=G5*0.5500000000000000D+01*OneOvT                         
            AuxR(6,I)=o1*G6                                             
            o1=o2*o1                                                    
            G7=G6*0.6500000000000000D+01*OneOvT                         
            AuxR(7,I)=o1*G7                                             
            o1=o2*o1                                                    
            G8=G7*0.7500000000000000D+01*OneOvT                         
            AuxR(8,I)=o1*G8                                             
         ENDIF                                                          
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen9(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:9,*)                                
      REAL*8 ET,TwoT,T,T2,T3,T4,RTE,RPE,Upq,TwoPi5x2                    
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(SqrtPi=1.7724538509055160273D0)                         
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)                       
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      INCLUDE "Gamma_Asymptotics.Inc"                                   
      INCLUDE "Mesh.Inc"                                                
      REAL*8 F9_0(0:Mesh)                                               
      REAL*8 F9_1(0:Mesh)                                               
      REAL*8 F9_2(0:Mesh)                                               
      REAL*8 F9_3(0:Mesh)                                               
      INCLUDE "Gamma_9.Inc"                                             
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=TwoPi5x2/(RTE*SQRT(RPE))                                   
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         IF(T.LT.Switch)THEN                                            
            T2=T*T                                                      
            T3=T*T2                                                     
            J=AINT(T*Grid)                                              
            G9=(F9_0(J)+T*F9_1(J)+T2*F9_2(J)+T3*F9_3(J))                
            ET=DEXP(-T)                                                 
            TwoT=2.0D0*T                                                
            G8=0.5882352941176471D-01*(TwoT*G9+ET)                      
            G7=0.6666666666666667D-01*(TwoT*G8+ET)                      
            G6=0.7692307692307693D-01*(TwoT*G7+ET)                      
            G5=0.9090909090909091D-01*(TwoT*G6+ET)                      
            G4=0.1111111111111111D+00*(TwoT*G5+ET)                      
            G3=0.1428571428571428D+00*(TwoT*G4+ET)                      
            G2=0.2000000000000000D+00*(TwoT*G3+ET)                      
            G1=0.3333333333333333D+00*(TwoT*G2+ET)                      
            G0=TwoT*G1+ET                                               
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            AuxR(5,I)=o1*G5                                             
            o1=o2*o1                                                    
            AuxR(6,I)=o1*G6                                             
            o1=o2*o1                                                    
            AuxR(7,I)=o1*G7                                             
            o1=o2*o1                                                    
            AuxR(8,I)=o1*G8                                             
            o1=o2*o1                                                    
            AuxR(9,I)=o1*G9                                             
         ELSE                                                           
            SqrtT=DSQRT(T)                                              
            OneOvT=One/T                                                
            G0=SqrtPi/(Two*SqrtT)                                       
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            G1=G0*0.5000000000000000D+00*OneOvT                         
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            G2=G1*0.1500000000000000D+01*OneOvT                         
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            G3=G2*0.2500000000000000D+01*OneOvT                         
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            G4=G3*0.3500000000000000D+01*OneOvT                         
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            G5=G4*0.4500000000000000D+01*OneOvT                         
            AuxR(5,I)=o1*G5                                             
            o1=o2*o1                                                    
            G6=G5*0.5500000000000000D+01*OneOvT                         
            AuxR(6,I)=o1*G6                                             
            o1=o2*o1                                                    
            G7=G6*0.6500000000000000D+01*OneOvT                         
            AuxR(7,I)=o1*G7                                             
            o1=o2*o1                                                    
            G8=G7*0.7500000000000000D+01*OneOvT                         
            AuxR(8,I)=o1*G8                                             
            o1=o2*o1                                                    
            G9=G8*0.8500000000000000D+01*OneOvT                         
            AuxR(9,I)=o1*G9                                             
         ENDIF                                                          
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen10(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)         
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:10,*)                               
      REAL*8 ET,TwoT,T,T2,T3,T4,RTE,RPE,Upq,TwoPi5x2                    
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(SqrtPi=1.7724538509055160273D0)                         
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)                       
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      INCLUDE "Gamma_Asymptotics.Inc"                                   
      INCLUDE "Mesh.Inc"                                                
      REAL*8 F10_0(0:Mesh)                                              
      REAL*8 F10_1(0:Mesh)                                              
      REAL*8 F10_2(0:Mesh)                                              
      REAL*8 F10_3(0:Mesh)                                              
      INCLUDE "Gamma_10.Inc"                                            
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=TwoPi5x2/(RTE*SQRT(RPE))                                   
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         IF(T.LT.Switch)THEN                                            
            T2=T*T                                                      
            T3=T*T2                                                     
            J=AINT(T*Grid)                                              
            G10=(F10_0(J)+T*F10_1(J)+T2*F10_2(J)+T3*F10_3(J))           
            ET=DEXP(-T)                                                 
            TwoT=2.0D0*T                                                
            G9=0.5263157894736842D-01*(TwoT*G10+ET)                     
            G8=0.5882352941176471D-01*(TwoT*G9+ET)                      
            G7=0.6666666666666667D-01*(TwoT*G8+ET)                      
            G6=0.7692307692307693D-01*(TwoT*G7+ET)                      
            G5=0.9090909090909091D-01*(TwoT*G6+ET)                      
            G4=0.1111111111111111D+00*(TwoT*G5+ET)                      
            G3=0.1428571428571428D+00*(TwoT*G4+ET)                      
            G2=0.2000000000000000D+00*(TwoT*G3+ET)                      
            G1=0.3333333333333333D+00*(TwoT*G2+ET)                      
            G0=TwoT*G1+ET                                               
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            AuxR(5,I)=o1*G5                                             
            o1=o2*o1                                                    
            AuxR(6,I)=o1*G6                                             
            o1=o2*o1                                                    
            AuxR(7,I)=o1*G7                                             
            o1=o2*o1                                                    
            AuxR(8,I)=o1*G8                                             
            o1=o2*o1                                                    
            AuxR(9,I)=o1*G9                                             
            o1=o2*o1                                                    
            AuxR(10,I)=o1*G10                                           
         ELSE                                                           
            SqrtT=DSQRT(T)                                              
            OneOvT=One/T                                                
            G0=SqrtPi/(Two*SqrtT)                                       
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            G1=G0*0.5000000000000000D+00*OneOvT                         
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            G2=G1*0.1500000000000000D+01*OneOvT                         
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            G3=G2*0.2500000000000000D+01*OneOvT                         
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            G4=G3*0.3500000000000000D+01*OneOvT                         
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            G5=G4*0.4500000000000000D+01*OneOvT                         
            AuxR(5,I)=o1*G5                                             
            o1=o2*o1                                                    
            G6=G5*0.5500000000000000D+01*OneOvT                         
            AuxR(6,I)=o1*G6                                             
            o1=o2*o1                                                    
            G7=G6*0.6500000000000000D+01*OneOvT                         
            AuxR(7,I)=o1*G7                                             
            o1=o2*o1                                                    
            G8=G7*0.7500000000000000D+01*OneOvT                         
            AuxR(8,I)=o1*G8                                             
            o1=o2*o1                                                    
            G9=G8*0.8500000000000000D+01*OneOvT                         
            AuxR(9,I)=o1*G9                                             
            o1=o2*o1                                                    
            G10=G9*0.9500000000000000D+01*OneOvT                        
            AuxR(10,I)=o1*G10                                           
         ENDIF                                                          
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen11(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)         
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:11,*)                               
      REAL*8 ET,TwoT,T,T2,T3,T4,RTE,RPE,Upq,TwoPi5x2                    
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(SqrtPi=1.7724538509055160273D0)                         
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)                       
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      INCLUDE "Gamma_Asymptotics.Inc"                                   
      INCLUDE "Mesh.Inc"                                                
      REAL*8 F11_0(0:Mesh)                                              
      REAL*8 F11_1(0:Mesh)                                              
      REAL*8 F11_2(0:Mesh)                                              
      REAL*8 F11_3(0:Mesh)                                              
      INCLUDE "Gamma_11.Inc"                                            
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=TwoPi5x2/(RTE*SQRT(RPE))                                   
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         IF(T.LT.Switch)THEN                                            
            T2=T*T                                                      
            T3=T*T2                                                     
            J=AINT(T*Grid)                                              
            G11=(F11_0(J)+T*F11_1(J)+T2*F11_2(J)+T3*F11_3(J))           
            ET=DEXP(-T)                                                 
            TwoT=2.0D0*T                                                
            G10=0.4761904761904762D-01*(TwoT*G11+ET)                    
            G9=0.5263157894736842D-01*(TwoT*G10+ET)                     
            G8=0.5882352941176471D-01*(TwoT*G9+ET)                      
            G7=0.6666666666666667D-01*(TwoT*G8+ET)                      
            G6=0.7692307692307693D-01*(TwoT*G7+ET)                      
            G5=0.9090909090909091D-01*(TwoT*G6+ET)                      
            G4=0.1111111111111111D+00*(TwoT*G5+ET)                      
            G3=0.1428571428571428D+00*(TwoT*G4+ET)                      
            G2=0.2000000000000000D+00*(TwoT*G3+ET)                      
            G1=0.3333333333333333D+00*(TwoT*G2+ET)                      
            G0=TwoT*G1+ET                                               
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            AuxR(5,I)=o1*G5                                             
            o1=o2*o1                                                    
            AuxR(6,I)=o1*G6                                             
            o1=o2*o1                                                    
            AuxR(7,I)=o1*G7                                             
            o1=o2*o1                                                    
            AuxR(8,I)=o1*G8                                             
            o1=o2*o1                                                    
            AuxR(9,I)=o1*G9                                             
            o1=o2*o1                                                    
            AuxR(10,I)=o1*G10                                           
            o1=o2*o1                                                    
            AuxR(11,I)=o1*G11                                           
         ELSE                                                           
            SqrtT=DSQRT(T)                                              
            OneOvT=One/T                                                
            G0=SqrtPi/(Two*SqrtT)                                       
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            G1=G0*0.5000000000000000D+00*OneOvT                         
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            G2=G1*0.1500000000000000D+01*OneOvT                         
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            G3=G2*0.2500000000000000D+01*OneOvT                         
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            G4=G3*0.3500000000000000D+01*OneOvT                         
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            G5=G4*0.4500000000000000D+01*OneOvT                         
            AuxR(5,I)=o1*G5                                             
            o1=o2*o1                                                    
            G6=G5*0.5500000000000000D+01*OneOvT                         
            AuxR(6,I)=o1*G6                                             
            o1=o2*o1                                                    
            G7=G6*0.6500000000000000D+01*OneOvT                         
            AuxR(7,I)=o1*G7                                             
            o1=o2*o1                                                    
            G8=G7*0.7500000000000000D+01*OneOvT                         
            AuxR(8,I)=o1*G8                                             
            o1=o2*o1                                                    
            G9=G8*0.8500000000000000D+01*OneOvT                         
            AuxR(9,I)=o1*G9                                             
            o1=o2*o1                                                    
            G10=G9*0.9500000000000000D+01*OneOvT                        
            AuxR(10,I)=o1*G10                                           
            o1=o2*o1                                                    
            G11=G10*0.1050000000000000D+02*OneOvT                       
            AuxR(11,I)=o1*G11                                           
         ENDIF                                                          
 100  CONTINUE
      RETURN
      END
      SUBROUTINE AuxIGen12(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)         
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:12,*)                               
      REAL*8 ET,TwoT,T,T2,T3,T4,RTE,RPE,Upq,TwoPi5x2                    
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(SqrtPi=1.7724538509055160273D0)                         
      PARAMETER(TwoPi5x2=3.4986836655249725693D1)                       
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      INCLUDE "Gamma_Asymptotics.Inc"                                   
      INCLUDE "Mesh.Inc"                                                
      REAL*8 F12_0(0:Mesh)                                              
      REAL*8 F12_1(0:Mesh)                                              
      REAL*8 F12_2(0:Mesh)                                              
      REAL*8 F12_3(0:Mesh)                                              
      INCLUDE "Gamma_12.Inc"                                            
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=TwoPi5x2/(RTE*SQRT(RPE))                                   
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         IF(T.LT.Switch)THEN                                            
            T2=T*T                                                      
            T3=T*T2                                                     
            J=AINT(T*Grid)                                              
            G12=(F12_0(J)+T*F12_1(J)+T2*F12_2(J)+T3*F12_3(J))           
            ET=DEXP(-T)                                                 
            TwoT=2.0D0*T                                                
            G11=0.4347826086956522D-01*(TwoT*G12+ET)                    
            G10=0.4761904761904762D-01*(TwoT*G11+ET)                    
            G9=0.5263157894736842D-01*(TwoT*G10+ET)                     
            G8=0.5882352941176471D-01*(TwoT*G9+ET)                      
            G7=0.6666666666666667D-01*(TwoT*G8+ET)                      
            G6=0.7692307692307693D-01*(TwoT*G7+ET)                      
            G5=0.9090909090909091D-01*(TwoT*G6+ET)                      
            G4=0.1111111111111111D+00*(TwoT*G5+ET)                      
            G3=0.1428571428571428D+00*(TwoT*G4+ET)                      
            G2=0.2000000000000000D+00*(TwoT*G3+ET)                      
            G1=0.3333333333333333D+00*(TwoT*G2+ET)                      
            G0=TwoT*G1+ET                                               
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            AuxR(5,I)=o1*G5                                             
            o1=o2*o1                                                    
            AuxR(6,I)=o1*G6                                             
            o1=o2*o1                                                    
            AuxR(7,I)=o1*G7                                             
            o1=o2*o1                                                    
            AuxR(8,I)=o1*G8                                             
            o1=o2*o1                                                    
            AuxR(9,I)=o1*G9                                             
            o1=o2*o1                                                    
            AuxR(10,I)=o1*G10                                           
            o1=o2*o1                                                    
            AuxR(11,I)=o1*G11                                           
            o1=o2*o1                                                    
            AuxR(12,I)=o1*G12                                           
         ELSE                                                           
            SqrtT=DSQRT(T)                                              
            OneOvT=One/T                                                
            G0=SqrtPi/(Two*SqrtT)                                       
            o1=Upq                                                      
            o2=-2D0*Omega                                               
            AuxR(0,I)=o1*G0                                             
            o1=o2*o1                                                    
            G1=G0*0.5000000000000000D+00*OneOvT                         
            AuxR(1,I)=o1*G1                                             
            o1=o2*o1                                                    
            G2=G1*0.1500000000000000D+01*OneOvT                         
            AuxR(2,I)=o1*G2                                             
            o1=o2*o1                                                    
            G3=G2*0.2500000000000000D+01*OneOvT                         
            AuxR(3,I)=o1*G3                                             
            o1=o2*o1                                                    
            G4=G3*0.3500000000000000D+01*OneOvT                         
            AuxR(4,I)=o1*G4                                             
            o1=o2*o1                                                    
            G5=G4*0.4500000000000000D+01*OneOvT                         
            AuxR(5,I)=o1*G5                                             
            o1=o2*o1                                                    
            G6=G5*0.5500000000000000D+01*OneOvT                         
            AuxR(6,I)=o1*G6                                             
            o1=o2*o1                                                    
            G7=G6*0.6500000000000000D+01*OneOvT                         
            AuxR(7,I)=o1*G7                                             
            o1=o2*o1                                                    
            G8=G7*0.7500000000000000D+01*OneOvT                         
            AuxR(8,I)=o1*G8                                             
            o1=o2*o1                                                    
            G9=G8*0.8500000000000000D+01*OneOvT                         
            AuxR(9,I)=o1*G9                                             
            o1=o2*o1                                                    
            G10=G9*0.9500000000000000D+01*OneOvT                        
            AuxR(10,I)=o1*G10                                           
            o1=o2*o1                                                    
            G11=G10*0.1050000000000000D+02*OneOvT                       
            AuxR(11,I)=o1*G11                                           
            o1=o2*o1                                                    
            G12=G11*0.1150000000000000D+02*OneOvT                       
            AuxR(12,I)=o1*G12                                           
         ENDIF                                                          
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen0(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:0,*)                                
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi                                   
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(Pi=3.1415926535897932385D0)                             
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=(Pi/RPE)**15D-1                                            
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         ET=EXP(-T)                                                     
         TwoT=2.0D0*T                                                   
         o1=Upq                                                         
         AuxR(0,I)=o1*ET                                                
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen1(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:1,*)                                
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi                                   
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(Pi=3.1415926535897932385D0)                             
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=(Pi/RPE)**15D-1                                            
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         ET=EXP(-T)                                                     
         TwoT=2.0D0*T                                                   
         o1=Upq                                                         
         o2=-2D0*Omega                                                  
         AuxR(0,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(1,I)=o1*ET                                                
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen2(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:2,*)                                
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi                                   
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(Pi=3.1415926535897932385D0)                             
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=(Pi/RPE)**15D-1                                            
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         ET=EXP(-T)                                                     
         TwoT=2.0D0*T                                                   
         o1=Upq                                                         
         o2=-2D0*Omega                                                  
         AuxR(0,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(1,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(2,I)=o1*ET                                                
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen3(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:3,*)                                
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi                                   
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(Pi=3.1415926535897932385D0)                             
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=(Pi/RPE)**15D-1                                            
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         ET=EXP(-T)                                                     
         TwoT=2.0D0*T                                                   
         o1=Upq                                                         
         o2=-2D0*Omega                                                  
         AuxR(0,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(1,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(2,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(3,I)=o1*ET                                                
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen4(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:4,*)                                
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi                                   
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(Pi=3.1415926535897932385D0)                             
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=(Pi/RPE)**15D-1                                            
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         ET=EXP(-T)                                                     
         TwoT=2.0D0*T                                                   
         o1=Upq                                                         
         o2=-2D0*Omega                                                  
         AuxR(0,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(1,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(2,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(3,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(4,I)=o1*ET                                                
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen5(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:5,*)                                
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi                                   
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(Pi=3.1415926535897932385D0)                             
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=(Pi/RPE)**15D-1                                            
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         ET=EXP(-T)                                                     
         TwoT=2.0D0*T                                                   
         o1=Upq                                                         
         o2=-2D0*Omega                                                  
         AuxR(0,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(1,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(2,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(3,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(4,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(5,I)=o1*ET                                                
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen6(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:6,*)                                
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi                                   
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(Pi=3.1415926535897932385D0)                             
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=(Pi/RPE)**15D-1                                            
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         ET=EXP(-T)                                                     
         TwoT=2.0D0*T                                                   
         o1=Upq                                                         
         o2=-2D0*Omega                                                  
         AuxR(0,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(1,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(2,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(3,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(4,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(5,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(6,I)=o1*ET                                                
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen7(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:7,*)                                
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi                                   
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(Pi=3.1415926535897932385D0)                             
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=(Pi/RPE)**15D-1                                            
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         ET=EXP(-T)                                                     
         TwoT=2.0D0*T                                                   
         o1=Upq                                                         
         o2=-2D0*Omega                                                  
         AuxR(0,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(1,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(2,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(3,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(4,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(5,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(6,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(7,I)=o1*ET                                                
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen8(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:8,*)                                
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi                                   
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(Pi=3.1415926535897932385D0)                             
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=(Pi/RPE)**15D-1                                            
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         ET=EXP(-T)                                                     
         TwoT=2.0D0*T                                                   
         o1=Upq                                                         
         o2=-2D0*Omega                                                  
         AuxR(0,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(1,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(2,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(3,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(4,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(5,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(6,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(7,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(8,I)=o1*ET                                                
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen9(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)          
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:9,*)                                
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi                                   
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(Pi=3.1415926535897932385D0)                             
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=(Pi/RPE)**15D-1                                            
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         ET=EXP(-T)                                                     
         TwoT=2.0D0*T                                                   
         o1=Upq                                                         
         o2=-2D0*Omega                                                  
         AuxR(0,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(1,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(2,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(3,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(4,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(5,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(6,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(7,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(8,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(9,I)=o1*ET                                                
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen10(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)         
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:10,*)                               
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi                                   
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(Pi=3.1415926535897932385D0)                             
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=(Pi/RPE)**15D-1                                            
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         ET=EXP(-T)                                                     
         TwoT=2.0D0*T                                                   
         o1=Upq                                                         
         o2=-2D0*Omega                                                  
         AuxR(0,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(1,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(2,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(3,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(4,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(5,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(6,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(7,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(8,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(9,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(10,I)=o1*ET                                               
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen11(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)         
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:11,*)                               
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi                                   
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(Pi=3.1415926535897932385D0)                             
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=(Pi/RPE)**15D-1                                            
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         ET=EXP(-T)                                                     
         TwoT=2.0D0*T                                                   
         o1=Upq                                                         
         o2=-2D0*Omega                                                  
         AuxR(0,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(1,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(2,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(3,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(4,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(5,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(6,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(7,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(8,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(9,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(10,I)=o1*ET                                               
         o1=o2*o1                                                       
         AuxR(11,I)=o1*ET                                               
 100  CONTINUE
      RETURN
      END
      SUBROUTINE OvrIGen12(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)         
      IMPLICIT DOUBLE PRECISION(G,O,S,F)                                
      INTEGER Nc,Ns,Mesh,I,J,K                                          
      REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)                          
      REAL*8  ZetaP,ZetaQ(*),AuxR(0:12,*)                               
      REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi                                   
      REAL*8 SqrtT,SqrtPi,One,Two                                       
      PARAMETER(Pi=3.1415926535897932385D0)                             
      PARAMETER(One=1D0)                                                
      PARAMETER(Two=2D0)                                                
      DO 100 I=1,Nc                                                     
         PQx(I)=-(P(1)-Q(1,I))                                          
         PQy(I)=-(P(2)-Q(2,I))                                          
         PQz(I)=-(P(3)-Q(3,I))                                          
         RTE=ZetaP*ZetaQ(I)                                             
         RPE=ZetaP+ZetaQ(I)                                             
         Omega=RTE/RPE                                                  
         Upq=(Pi/RPE)**15D-1                                            
         T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))            
         ET=EXP(-T)                                                     
         TwoT=2.0D0*T                                                   
         o1=Upq                                                         
         o2=-2D0*Omega                                                  
         AuxR(0,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(1,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(2,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(3,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(4,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(5,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(6,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(7,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(8,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(9,I)=o1*ET                                                
         o1=o2*o1                                                       
         AuxR(10,I)=o1*ET                                               
         o1=o2*o1                                                       
         AuxR(11,I)=o1*ET                                               
         o1=o2*o1                                                       
         AuxR(12,I)=o1*ET                                               
 100  CONTINUE
      RETURN
      END
