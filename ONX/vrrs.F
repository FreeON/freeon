      SUBROUTINE VRRs(N,LBra,LKet,W,R,WR,WZ)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 W(N,*),R(N,*),WR(12,*),WZ(5,*)

      LCode=10*LBra+LKet
      IF(LCode.EQ.01)THEN                                          
         DO I=1,N                                                         
            W(I,1)=R(I,1)                                                    
            W(I,4)=R(I,2)                                                    
            W(I,2)=WR(2,I)*W(I,1)+WR(8,I)*W(I,4)                             
            W(I,3)=WR(4,I)*W(I,1)+WR(10,I)*W(I,4)                            
            W(I,4)=WR(6,I)*W(I,1)+WR(12,I)*W(I,4)                            
         ENDDO                                                            
         RETURN                                                           
      ELSEIF(LCode.EQ.02)THEN                                          
         DO I=1,N                                                         
            W(I,1)=R(I,1)                                                    
            W(I,6)=R(I,2)                                                    
            W(I,5)=R(I,3)                                                    
            W(I,2)=WR(2,I)*W(I,1)+WR(8,I)*W(I,6)                             
            W(I,3)=WR(4,I)*W(I,1)+WR(10,I)*W(I,6)                            
            W(I,4)=WR(6,I)*W(I,1)+WR(12,I)*W(I,6)                            
            W(I,7)=WR(2,I)*W(I,6)+WR(8,I)*W(I,5)                             
            W(I,8)=WR(4,I)*W(I,6)+WR(10,I)*W(I,5)                            
            W(I,9)=WR(6,I)*W(I,6)+WR(12,I)*W(I,5)                            
            Ts0=WZ(1,I)*(W(I,1)-WZ(3,I)*W(I,6))
            W(I,5)=WR(2,I)*W(I,2)+WR(8,I)*W(I,7)+ Ts0
            W(I,7)=WR(4,I)*W(I,3)+WR(10,I)*W(I,8)+ Ts0
            W(I,10)=WR(6,I)*W(I,4)+WR(12,I)*W(I,9)+ Ts0
            W(I,6)=WR(2,I)*W(I,3)+WR(8,I)*W(I,8)                             
            W(I,8)=WR(2,I)*W(I,4)+WR(8,I)*W(I,9)                             
            W(I,9)=WR(4,I)*W(I,4)+WR(10,I)*W(I,9)                            
         ENDDO                                                            
         RETURN                                                           
      ELSEIF(LCode.EQ.03)THEN                                          
         DO I=1,N                                                         
            W(I,1)=R(I,1)                                                    
            W(I,11)=R(I,2)                                                   
            W(I,12)=R(I,3)                                                   
            W(I,5)=R(I,4)                                                    
            W(I,2)=WR(2,I)*W(I,1)+WR(8,I)*W(I,11)                            
            W(I,3)=WR(4,I)*W(I,1)+WR(10,I)*W(I,11)                           
            W(I,4)=WR(6,I)*W(I,1)+WR(12,I)*W(I,11)                           
            W(I,13)=WR(2,I)*W(I,11)+WR(8,I)*W(I,12)                          
            W(I,14)=WR(4,I)*W(I,11)+WR(10,I)*W(I,12)                         
            W(I,16)=WR(6,I)*W(I,11)+WR(12,I)*W(I,12)                         
            W(I,15)=WR(2,I)*W(I,12)+WR(8,I)*W(I,5)                           
            W(I,17)=WR(4,I)*W(I,12)+WR(10,I)*W(I,5)                          
            W(I,18)=WR(6,I)*W(I,12)+WR(12,I)*W(I,5)                          
            Ts0=WZ(1,I)*(W(I,1)-WZ(3,I)*W(I,11))
            W(I,5)=WR(2,I)*W(I,2)+WR(8,I)*W(I,13)+ Ts0
            W(I,6)=WR(2,I)*W(I,3)+WR(8,I)*W(I,14)                            
            W(I,7)=WR(4,I)*W(I,3)+WR(10,I)*W(I,14)+ Ts0
            W(I,8)=WR(2,I)*W(I,4)+WR(8,I)*W(I,16)                            
            W(I,9)=WR(4,I)*W(I,4)+WR(10,I)*W(I,16)                           
            W(I,10)=WR(6,I)*W(I,4)+WR(12,I)*W(I,16)+ Ts0
            Ts1=WZ(1,I)*(W(I,11)-WZ(3,I)*W(I,12))
            W(I,19)=WR(2,I)*W(I,13)+WR(8,I)*W(I,15)+ Ts1
            W(I,20)=WR(2,I)*W(I,14)+WR(8,I)*W(I,17)                          
            W(I,21)=WR(4,I)*W(I,14)+WR(10,I)*W(I,17)+ Ts1
            W(I,22)=WR(2,I)*W(I,16)+WR(8,I)*W(I,18)                          
            W(I,23)=WR(4,I)*W(I,16)+WR(10,I)*W(I,18)                         
            W(I,24)=WR(6,I)*W(I,16)+WR(12,I)*W(I,18)+ Ts1
            W(I,11)=WR(2,I)*W(I,5)+WR(8,I)*W(I,19)+                          
     &        2.0D0*WZ(1,I)*(W(I,2)-WZ(3,I)*W(I,13))                           
            W(I,12)=WR(2,I)*W(I,6)+WR(8,I)*W(I,20)+                          
     &        WZ(1,I)*(W(I,3)-WZ(3,I)*W(I,14))                                 
            W(I,14)=WR(4,I)*W(I,7)+WR(10,I)*W(I,21)+                         
     &        2.0D0*WZ(1,I)*(W(I,3)-WZ(3,I)*W(I,14))                           
            W(I,15)=WR(2,I)*W(I,8)+WR(8,I)*W(I,22)+                          
     &        WZ(1,I)*(W(I,4)-WZ(3,I)*W(I,16))                                 
            W(I,17)=WR(4,I)*W(I,9)+WR(10,I)*W(I,23)+                         
     &        WZ(1,I)*(W(I,4)-WZ(3,I)*W(I,16))                                 
            W(I,20)=WR(6,I)*W(I,10)+WR(12,I)*W(I,24)+                        
     &        2.0D0*WZ(1,I)*(W(I,4)-WZ(3,I)*W(I,16))                           
            W(I,13)=WR(2,I)*W(I,7)+WR(8,I)*W(I,21)                           
            W(I,16)=WR(2,I)*W(I,9)+WR(8,I)*W(I,23)                           
            W(I,18)=WR(2,I)*W(I,10)+WR(8,I)*W(I,24)                          
            W(I,19)=WR(4,I)*W(I,10)+WR(10,I)*W(I,24)                         
         ENDDO                                                            
         RETURN                                                           
      ELSEIF(LCode.EQ.04)THEN                                          
         DO I=1,N                                                         
            W(I,1)=R(I,1)                                                    
            W(I,11)=R(I,2)                                                   
            W(I,12)=R(I,3)                                                   
            W(I,13)=R(I,4)                                                   
            W(I,5)=R(I,5)                                                    
            W(I,2)=WR(2,I)*W(I,1)+WR(8,I)*W(I,11)                            
            W(I,3)=WR(4,I)*W(I,1)+WR(10,I)*W(I,11)                           
            W(I,4)=WR(6,I)*W(I,1)+WR(12,I)*W(I,11)                           
            W(I,21)=WR(2,I)*W(I,11)+WR(8,I)*W(I,12)                          
            W(I,22)=WR(4,I)*W(I,11)+WR(10,I)*W(I,12)                         
            W(I,23)=WR(6,I)*W(I,11)+WR(12,I)*W(I,12)                         
            W(I,24)=WR(2,I)*W(I,12)+WR(8,I)*W(I,13)                          
            W(I,25)=WR(4,I)*W(I,12)+WR(10,I)*W(I,13)                         
            W(I,26)=WR(6,I)*W(I,12)+WR(12,I)*W(I,13)                         
            W(I,14)=WR(2,I)*W(I,13)+WR(8,I)*W(I,5)                           
            W(I,15)=WR(4,I)*W(I,13)+WR(10,I)*W(I,5)                          
            W(I,16)=WR(6,I)*W(I,13)+WR(12,I)*W(I,5)                          
            Ts0=WZ(1,I)*(W(I,1)-WZ(3,I)*W(I,11))
            W(I,5)=WR(2,I)*W(I,2)+WR(8,I)*W(I,21)+ Ts0
            W(I,6)=WR(2,I)*W(I,3)+WR(8,I)*W(I,22)                            
            W(I,7)=WR(4,I)*W(I,3)+WR(10,I)*W(I,22)+ Ts0
            W(I,8)=WR(2,I)*W(I,4)+WR(8,I)*W(I,23)                            
            W(I,9)=WR(4,I)*W(I,4)+WR(10,I)*W(I,23)                           
            W(I,10)=WR(6,I)*W(I,4)+WR(12,I)*W(I,23)+ Ts0
            Ts1=WZ(1,I)*(W(I,11)-WZ(3,I)*W(I,12))
            W(I,27)=WR(2,I)*W(I,21)+WR(8,I)*W(I,24)+ Ts1
            W(I,28)=WR(2,I)*W(I,22)+WR(8,I)*W(I,25)                          
            W(I,29)=WR(4,I)*W(I,22)+WR(10,I)*W(I,25)+ Ts1
            W(I,30)=WR(2,I)*W(I,23)+WR(8,I)*W(I,26)                          
            W(I,31)=WR(4,I)*W(I,23)+WR(10,I)*W(I,26)                         
            W(I,33)=WR(6,I)*W(I,23)+WR(12,I)*W(I,26)+ Ts1
            Ts2=WZ(1,I)*(W(I,12)-WZ(3,I)*W(I,13))
            W(I,32)=WR(2,I)*W(I,24)+WR(8,I)*W(I,14)+ Ts2
            W(I,34)=WR(2,I)*W(I,25)+WR(8,I)*W(I,15)                          
            W(I,35)=WR(4,I)*W(I,25)+WR(10,I)*W(I,15)+ Ts2
            W(I,36)=WR(2,I)*W(I,26)+WR(8,I)*W(I,16)                          
            W(I,37)=WR(4,I)*W(I,26)+WR(10,I)*W(I,16)                         
            W(I,38)=WR(6,I)*W(I,26)+WR(12,I)*W(I,16)+ Ts2
            W(I,11)=WR(2,I)*W(I,5)+WR(8,I)*W(I,27)+                          
     &        2.0D0*WZ(1,I)*(W(I,2)-WZ(3,I)*W(I,21))                           
            W(I,12)=WR(2,I)*W(I,6)+WR(8,I)*W(I,28)+                          
     &        WZ(1,I)*(W(I,3)-WZ(3,I)*W(I,22))                                 
            W(I,13)=WR(2,I)*W(I,7)+WR(8,I)*W(I,29)                           
            W(I,14)=WR(4,I)*W(I,7)+WR(10,I)*W(I,29)+                         
     &        2.0D0*WZ(1,I)*(W(I,3)-WZ(3,I)*W(I,22))                           
            W(I,15)=WR(2,I)*W(I,8)+WR(8,I)*W(I,30)+                          
     &        WZ(1,I)*(W(I,4)-WZ(3,I)*W(I,23))                                 
            W(I,16)=WR(2,I)*W(I,9)+WR(8,I)*W(I,31)                           
            W(I,17)=WR(4,I)*W(I,9)+WR(10,I)*W(I,31)+                         
     &        WZ(1,I)*(W(I,4)-WZ(3,I)*W(I,23))                                 
            W(I,18)=WR(2,I)*W(I,10)+WR(8,I)*W(I,33)                          
            W(I,19)=WR(4,I)*W(I,10)+WR(10,I)*W(I,33)                         
            W(I,20)=WR(6,I)*W(I,10)+WR(12,I)*W(I,33)+                        
     &        2.0D0*WZ(1,I)*(W(I,4)-WZ(3,I)*W(I,23))                           
            W(I,39)=WR(2,I)*W(I,27)+WR(8,I)*W(I,32)+                         
     &        2.0D0*WZ(1,I)*(W(I,21)-WZ(3,I)*W(I,24))                          
            W(I,40)=WR(2,I)*W(I,28)+WR(8,I)*W(I,34)+                         
     &        WZ(1,I)*(W(I,22)-WZ(3,I)*W(I,25))                                
            W(I,41)=WR(2,I)*W(I,29)+WR(8,I)*W(I,35)                          
            W(I,42)=WR(4,I)*W(I,29)+WR(10,I)*W(I,35)+                        
     &        2.0D0*WZ(1,I)*(W(I,22)-WZ(3,I)*W(I,25))                          
            W(I,43)=WR(2,I)*W(I,30)+WR(8,I)*W(I,36)+                         
     &        WZ(1,I)*(W(I,23)-WZ(3,I)*W(I,26))                                
            W(I,44)=WR(2,I)*W(I,31)+WR(8,I)*W(I,37)                          
            W(I,45)=WR(4,I)*W(I,31)+WR(10,I)*W(I,37)+                        
     &        WZ(1,I)*(W(I,23)-WZ(3,I)*W(I,26))                                
            W(I,46)=WR(2,I)*W(I,33)+WR(8,I)*W(I,38)                          
            W(I,47)=WR(4,I)*W(I,33)+WR(10,I)*W(I,38)                         
            W(I,48)=WR(6,I)*W(I,33)+WR(12,I)*W(I,38)+                        
     &        2.0D0*WZ(1,I)*(W(I,23)-WZ(3,I)*W(I,26))                          
            W(I,21)=WR(2,I)*W(I,11)+WR(8,I)*W(I,39)+                         
     &        3.0D0*WZ(1,I)*(W(I,5)-WZ(3,I)*W(I,27))                           
            W(I,22)=WR(2,I)*W(I,12)+WR(8,I)*W(I,40)+                         
     &        2.0D0*WZ(1,I)*(W(I,6)-WZ(3,I)*W(I,28))                           
            W(I,23)=WR(2,I)*W(I,13)+WR(8,I)*W(I,41)+                         
     &        WZ(1,I)*(W(I,7)-WZ(3,I)*W(I,29))                                 
            W(I,25)=WR(4,I)*W(I,14)+WR(10,I)*W(I,42)+                        
     &        3.0D0*WZ(1,I)*(W(I,7)-WZ(3,I)*W(I,29))                           
            W(I,26)=WR(2,I)*W(I,15)+WR(8,I)*W(I,43)+                         
     &        2.0D0*WZ(1,I)*(W(I,8)-WZ(3,I)*W(I,30))                           
            W(I,27)=WR(2,I)*W(I,16)+WR(8,I)*W(I,44)+                         
     &        WZ(1,I)*(W(I,9)-WZ(3,I)*W(I,31))                                 
            W(I,29)=WR(4,I)*W(I,17)+WR(10,I)*W(I,45)+                        
     &        2.0D0*WZ(1,I)*(W(I,9)-WZ(3,I)*W(I,31))                           
            W(I,30)=WR(2,I)*W(I,18)+WR(8,I)*W(I,46)+                         
     &        WZ(1,I)*(W(I,10)-WZ(3,I)*W(I,33))                                
            W(I,32)=WR(4,I)*W(I,19)+WR(10,I)*W(I,47)+                        
     &        WZ(1,I)*(W(I,10)-WZ(3,I)*W(I,33))                                
            W(I,35)=WR(6,I)*W(I,20)+WR(12,I)*W(I,48)+                        
     &        3.0D0*WZ(1,I)*(W(I,10)-WZ(3,I)*W(I,33))                          
            W(I,24)=WR(2,I)*W(I,14)+WR(8,I)*W(I,42)                          
            W(I,28)=WR(2,I)*W(I,17)+WR(8,I)*W(I,45)                          
            W(I,31)=WR(2,I)*W(I,19)+WR(8,I)*W(I,47)                          
            W(I,33)=WR(2,I)*W(I,20)+WR(8,I)*W(I,48)                          
            W(I,34)=WR(4,I)*W(I,20)+WR(10,I)*W(I,48)                         
         ENDDO                                                            
         RETURN                                                           
      ELSEIF(LCode.EQ.10)THEN                                          
         DO I=1,N                                                         
            W(I,1)=R(I,1)                                                    
            W(I,4)=R(I,2)                                                    
            W(I,2)=WR(1,I)*W(I,1)+WR(7,I)*W(I,4)                             
            W(I,3)=WR(3,I)*W(I,1)+WR(9,I)*W(I,4)                             
            W(I,4)=WR(5,I)*W(I,1)+WR(11,I)*W(I,4)                            
         ENDDO                                                            
         RETURN                                                           
      ELSEIF(LCode.EQ.11)THEN                                          
         DO I=1,N                                                         
            W(I,1)=R(I,1)                                                    
            W(I,4)=R(I,2)                                                    
            W(I,2)=R(I,3)                                                    
            W(I,5)=WR(2,I)*W(I,1)+WR(8,I)*W(I,4)                             
            W(I,9)=WR(4,I)*W(I,1)+WR(10,I)*W(I,4)                            
            W(I,13)=WR(6,I)*W(I,1)+WR(12,I)*W(I,4)                           
            W(I,8)=WR(2,I)*W(I,4)+WR(8,I)*W(I,2)                             
            W(I,12)=WR(4,I)*W(I,4)+WR(10,I)*W(I,2)                           
            W(I,15)=WR(6,I)*W(I,4)+WR(12,I)*W(I,2)                           
            Ts0=WZ(5,I)*W(I,4)
            W(I,6)=WR(1,I)*W(I,5)+WR(7,I)*W(I,8)+ Ts0
            W(I,11)=WR(3,I)*W(I,9)+WR(9,I)*W(I,12)+ Ts0
            W(I,16)=WR(5,I)*W(I,13)+WR(11,I)*W(I,15)+ Ts0
            W(I,2)=WR(1,I)*W(I,1)+WR(7,I)*W(I,4)                             
            W(I,3)=WR(3,I)*W(I,1)+WR(9,I)*W(I,4)                             
            W(I,4)=WR(5,I)*W(I,1)+WR(11,I)*W(I,4)                            
            W(I,7)=WR(3,I)*W(I,5)+WR(9,I)*W(I,8)                             
            W(I,8)=WR(5,I)*W(I,5)+WR(11,I)*W(I,8)                            
            W(I,10)=WR(1,I)*W(I,9)+WR(7,I)*W(I,12)                           
            W(I,12)=WR(5,I)*W(I,9)+WR(11,I)*W(I,12)                          
            W(I,14)=WR(1,I)*W(I,13)+WR(7,I)*W(I,15)                          
            W(I,15)=WR(3,I)*W(I,13)+WR(9,I)*W(I,15)                          
         ENDDO                                                            
         RETURN                                                           
      ELSEIF(LCode.EQ.12)THEN                                          
         DO I=1,N                                                         
            W(I,1)=R(I,1)                                                    
            W(I,4)=R(I,2)                                                    
            W(I,2)=R(I,3)                                                    
            W(I,3)=R(I,4)                                                    
            W(I,5)=WR(2,I)*W(I,1)+WR(8,I)*W(I,4)                             
            W(I,9)=WR(4,I)*W(I,1)+WR(10,I)*W(I,4)                            
            W(I,13)=WR(6,I)*W(I,1)+WR(12,I)*W(I,4)                           
            W(I,8)=WR(2,I)*W(I,4)+WR(8,I)*W(I,2)                             
            W(I,12)=WR(4,I)*W(I,4)+WR(10,I)*W(I,2)                           
            W(I,15)=WR(6,I)*W(I,4)+WR(12,I)*W(I,2)                           
            W(I,6)=WR(2,I)*W(I,2)+WR(8,I)*W(I,3)                             
            W(I,7)=WR(4,I)*W(I,2)+WR(10,I)*W(I,3)                            
            W(I,10)=WR(6,I)*W(I,2)+WR(12,I)*W(I,3)                           
            Ts0=WZ(1,I)*(W(I,1)-WZ(3,I)*W(I,4))
            W(I,17)=WR(2,I)*W(I,5)+WR(8,I)*W(I,8)+ Ts0
            W(I,21)=WR(2,I)*W(I,9)+WR(8,I)*W(I,12)                           
            W(I,25)=WR(4,I)*W(I,9)+WR(10,I)*W(I,12)+ Ts0
            W(I,29)=WR(2,I)*W(I,13)+WR(8,I)*W(I,15)                          
            W(I,33)=WR(4,I)*W(I,13)+WR(10,I)*W(I,15)                         
            W(I,37)=WR(6,I)*W(I,13)+WR(12,I)*W(I,15)+ Ts0
            Ts1=WZ(1,I)*(W(I,4)-WZ(3,I)*W(I,2))
            W(I,20)=WR(2,I)*W(I,8)+WR(8,I)*W(I,6)+ Ts1
            W(I,24)=WR(2,I)*W(I,12)+WR(8,I)*W(I,7)                           
            W(I,28)=WR(4,I)*W(I,12)+WR(10,I)*W(I,7)+ Ts1
            W(I,31)=WR(2,I)*W(I,15)+WR(8,I)*W(I,10)                          
            W(I,34)=WR(4,I)*W(I,15)+WR(10,I)*W(I,10)                         
            W(I,39)=WR(6,I)*W(I,15)+WR(12,I)*W(I,10)+ Ts1
            Ts2=WZ(5,I)*W(I,4)
            W(I,6)=WR(1,I)*W(I,5)+WR(7,I)*W(I,8)+ Ts2
            W(I,11)=WR(3,I)*W(I,9)+WR(9,I)*W(I,12)+ Ts2
            W(I,16)=WR(5,I)*W(I,13)+WR(11,I)*W(I,15)+ Ts2
            W(I,18)=WR(1,I)*W(I,17)+WR(7,I)*W(I,20)+                         
     &        2.0D0*WZ(5,I)*W(I,8)                                             
            W(I,22)=WR(1,I)*W(I,21)+WR(7,I)*W(I,24)+                         
     &        WZ(5,I)*W(I,12)                                                  
            W(I,23)=WR(3,I)*W(I,21)+WR(9,I)*W(I,24)+                         
     &        WZ(5,I)*W(I,8)                                                   
            W(I,27)=WR(3,I)*W(I,25)+WR(9,I)*W(I,28)+                         
     &        2.0D0*WZ(5,I)*W(I,12)                                            
            W(I,30)=WR(1,I)*W(I,29)+WR(7,I)*W(I,31)+                         
     &        WZ(5,I)*W(I,15)                                                  
            W(I,32)=WR(5,I)*W(I,29)+WR(11,I)*W(I,31)+                        
     &        WZ(5,I)*W(I,8)                                                   
            W(I,35)=WR(3,I)*W(I,33)+WR(9,I)*W(I,34)+                         
     &        WZ(5,I)*W(I,15)                                                  
            W(I,36)=WR(5,I)*W(I,33)+WR(11,I)*W(I,34)+                        
     &        WZ(5,I)*W(I,12)                                                  
            W(I,40)=WR(5,I)*W(I,37)+WR(11,I)*W(I,39)+                        
     &        2.0D0*WZ(5,I)*W(I,15)                                            
            W(I,2)=WR(1,I)*W(I,1)+WR(7,I)*W(I,4)                             
            W(I,3)=WR(3,I)*W(I,1)+WR(9,I)*W(I,4)                             
            W(I,4)=WR(5,I)*W(I,1)+WR(11,I)*W(I,4)                            
            W(I,7)=WR(3,I)*W(I,5)+WR(9,I)*W(I,8)                             
            W(I,8)=WR(5,I)*W(I,5)+WR(11,I)*W(I,8)                            
            W(I,10)=WR(1,I)*W(I,9)+WR(7,I)*W(I,12)                           
            W(I,12)=WR(5,I)*W(I,9)+WR(11,I)*W(I,12)                          
            W(I,14)=WR(1,I)*W(I,13)+WR(7,I)*W(I,15)                          
            W(I,15)=WR(3,I)*W(I,13)+WR(9,I)*W(I,15)                          
            W(I,19)=WR(3,I)*W(I,17)+WR(9,I)*W(I,20)                          
            W(I,20)=WR(5,I)*W(I,17)+WR(11,I)*W(I,20)                         
            W(I,24)=WR(5,I)*W(I,21)+WR(11,I)*W(I,24)                         
            W(I,26)=WR(1,I)*W(I,25)+WR(7,I)*W(I,28)                          
            W(I,28)=WR(5,I)*W(I,25)+WR(11,I)*W(I,28)                         
            W(I,31)=WR(3,I)*W(I,29)+WR(9,I)*W(I,31)                          
            W(I,34)=WR(1,I)*W(I,33)+WR(7,I)*W(I,34)                          
            W(I,38)=WR(1,I)*W(I,37)+WR(7,I)*W(I,39)                          
            W(I,39)=WR(3,I)*W(I,37)+WR(9,I)*W(I,39)                          
         ENDDO                                                            
         RETURN                                                           
      ELSEIF(LCode.EQ.13)THEN                                          
         DO I=1,N                                                         
            W(I,1)=R(I,1)                                                    
            W(I,4)=R(I,2)                                                    
            W(I,2)=R(I,3)                                                    
            W(I,3)=R(I,4)                                                    
            W(I,6)=R(I,5)                                                    
            W(I,5)=WR(2,I)*W(I,1)+WR(8,I)*W(I,4)                             
            W(I,9)=WR(4,I)*W(I,1)+WR(10,I)*W(I,4)                            
            W(I,13)=WR(6,I)*W(I,1)+WR(12,I)*W(I,4)                           
            W(I,8)=WR(2,I)*W(I,4)+WR(8,I)*W(I,2)                             
            W(I,12)=WR(4,I)*W(I,4)+WR(10,I)*W(I,2)                           
            W(I,15)=WR(6,I)*W(I,4)+WR(12,I)*W(I,2)                           
            W(I,7)=WR(2,I)*W(I,2)+WR(8,I)*W(I,3)                             
            W(I,10)=WR(4,I)*W(I,2)+WR(10,I)*W(I,3)                           
            W(I,11)=WR(6,I)*W(I,2)+WR(12,I)*W(I,3)                           
            W(I,14)=WR(2,I)*W(I,3)+WR(8,I)*W(I,6)                            
            W(I,16)=WR(4,I)*W(I,3)+WR(10,I)*W(I,6)                           
            W(I,18)=WR(6,I)*W(I,3)+WR(12,I)*W(I,6)                           
            W(I,17)=WR(2,I)*W(I,5)+WR(8,I)*W(I,8)+                           
     &        WZ(1,I)*(W(I,1)-WZ(3,I)*W(I,4))                                  
            W(I,21)=WR(2,I)*W(I,9)+WR(8,I)*W(I,12)                           
            W(I,25)=WR(4,I)*W(I,9)+WR(10,I)*W(I,12)+                         
     &        WZ(1,I)*(W(I,1)-WZ(3,I)*W(I,4))                                  
            W(I,29)=WR(2,I)*W(I,13)+WR(8,I)*W(I,15)                          
            W(I,33)=WR(4,I)*W(I,13)+WR(10,I)*W(I,15)                         
            W(I,37)=WR(6,I)*W(I,13)+WR(12,I)*W(I,15)+                        
     &        WZ(1,I)*(W(I,1)-WZ(3,I)*W(I,4))                                  
            W(I,20)=WR(2,I)*W(I,8)+WR(8,I)*W(I,7)+                           
     &        WZ(1,I)*(W(I,4)-WZ(3,I)*W(I,2))                                  
            W(I,24)=WR(2,I)*W(I,12)+WR(8,I)*W(I,10)                          
            W(I,28)=WR(4,I)*W(I,12)+WR(10,I)*W(I,10)+                        
     &        WZ(1,I)*(W(I,4)-WZ(3,I)*W(I,2))                                  
            W(I,31)=WR(2,I)*W(I,15)+WR(8,I)*W(I,11)                          
            W(I,34)=WR(4,I)*W(I,15)+WR(10,I)*W(I,11)                         
            W(I,39)=WR(6,I)*W(I,15)+WR(12,I)*W(I,11)+                        
     &        WZ(1,I)*(W(I,4)-WZ(3,I)*W(I,2))                                  
            W(I,19)=WR(2,I)*W(I,7)+WR(8,I)*W(I,14)+                          
     &        WZ(1,I)*(W(I,2)-WZ(3,I)*W(I,3))                                  
            W(I,22)=WR(2,I)*W(I,10)+WR(8,I)*W(I,16)                          
            W(I,23)=WR(4,I)*W(I,10)+WR(10,I)*W(I,16)+                        
     &        WZ(1,I)*(W(I,2)-WZ(3,I)*W(I,3))                                  
            W(I,26)=WR(2,I)*W(I,11)+WR(8,I)*W(I,18)                          
            W(I,27)=WR(4,I)*W(I,11)+WR(10,I)*W(I,18)                         
            W(I,30)=WR(6,I)*W(I,11)+WR(12,I)*W(I,18)+                        
     &        WZ(1,I)*(W(I,2)-WZ(3,I)*W(I,3))                                  
            W(I,41)=WR(2,I)*W(I,17)+WR(8,I)*W(I,20)+                         
     &        2.0D0*WZ(1,I)*(W(I,5)-WZ(3,I)*W(I,8))                            
            W(I,45)=WR(2,I)*W(I,21)+WR(8,I)*W(I,24)+                         
     &        WZ(1,I)*(W(I,9)-WZ(3,I)*W(I,12))                                 
            W(I,49)=WR(2,I)*W(I,25)+WR(8,I)*W(I,28)                          
            W(I,53)=WR(4,I)*W(I,25)+WR(10,I)*W(I,28)+                        
     &        2.0D0*WZ(1,I)*(W(I,9)-WZ(3,I)*W(I,12))                           
            W(I,57)=WR(2,I)*W(I,29)+WR(8,I)*W(I,31)+                         
     &        WZ(1,I)*(W(I,13)-WZ(3,I)*W(I,15))                                
            W(I,61)=WR(2,I)*W(I,33)+WR(8,I)*W(I,34)                          
            W(I,65)=WR(4,I)*W(I,33)+WR(10,I)*W(I,34)+                        
     &        WZ(1,I)*(W(I,13)-WZ(3,I)*W(I,15))                                
            W(I,69)=WR(2,I)*W(I,37)+WR(8,I)*W(I,39)                          
            W(I,73)=WR(4,I)*W(I,37)+WR(10,I)*W(I,39)                         
            W(I,77)=WR(6,I)*W(I,37)+WR(12,I)*W(I,39)+                        
     &        2.0D0*WZ(1,I)*(W(I,13)-WZ(3,I)*W(I,15))                          
            W(I,44)=WR(2,I)*W(I,20)+WR(8,I)*W(I,19)+                         
     &        2.0D0*WZ(1,I)*(W(I,8)-WZ(3,I)*W(I,7))                            
            W(I,48)=WR(2,I)*W(I,24)+WR(8,I)*W(I,22)+                         
     &        WZ(1,I)*(W(I,12)-WZ(3,I)*W(I,10))                                
            W(I,52)=WR(2,I)*W(I,28)+WR(8,I)*W(I,23)                          
            W(I,56)=WR(4,I)*W(I,28)+WR(10,I)*W(I,23)+                        
     &        2.0D0*WZ(1,I)*(W(I,12)-WZ(3,I)*W(I,10))                          
            W(I,59)=WR(2,I)*W(I,31)+WR(8,I)*W(I,26)+                         
     &        WZ(1,I)*(W(I,15)-WZ(3,I)*W(I,11))                                
            W(I,38)=WR(2,I)*W(I,34)+WR(8,I)*W(I,27)                          
            W(I,66)=WR(4,I)*W(I,34)+WR(10,I)*W(I,27)+                        
     &        WZ(1,I)*(W(I,15)-WZ(3,I)*W(I,11))                                
            W(I,71)=WR(2,I)*W(I,39)+WR(8,I)*W(I,30)                          
            W(I,74)=WR(4,I)*W(I,39)+WR(10,I)*W(I,30)                         
            W(I,79)=WR(6,I)*W(I,39)+WR(12,I)*W(I,30)+                        
     &        2.0D0*WZ(1,I)*(W(I,15)-WZ(3,I)*W(I,11))                          
            W(I,6)=WR(1,I)*W(I,5)+WR(7,I)*W(I,8)+                            
     &        WZ(5,I)*W(I,4)                                                   
            W(I,11)=WR(3,I)*W(I,9)+WR(9,I)*W(I,12)+                          
     &        WZ(5,I)*W(I,4)                                                   
            W(I,16)=WR(5,I)*W(I,13)+WR(11,I)*W(I,15)+                        
     &        WZ(5,I)*W(I,4)                                                   
            W(I,18)=WR(1,I)*W(I,17)+WR(7,I)*W(I,20)+                         
     &        2.0D0*WZ(5,I)*W(I,8)                                             
            W(I,22)=WR(1,I)*W(I,21)+WR(7,I)*W(I,24)+                         
     &        WZ(5,I)*W(I,12)                                                  
            W(I,23)=WR(3,I)*W(I,21)+WR(9,I)*W(I,24)+                         
     &        WZ(5,I)*W(I,8)                                                   
            W(I,27)=WR(3,I)*W(I,25)+WR(9,I)*W(I,28)+                         
     &        2.0D0*WZ(5,I)*W(I,12)                                            
            W(I,30)=WR(1,I)*W(I,29)+WR(7,I)*W(I,31)+                         
     &        WZ(5,I)*W(I,15)                                                  
            W(I,32)=WR(5,I)*W(I,29)+WR(11,I)*W(I,31)+                        
     &        WZ(5,I)*W(I,8)                                                   
            W(I,35)=WR(3,I)*W(I,33)+WR(9,I)*W(I,34)+                         
     &        WZ(5,I)*W(I,15)                                                  
            W(I,36)=WR(5,I)*W(I,33)+WR(11,I)*W(I,34)+                        
     &        WZ(5,I)*W(I,12)                                                  
            W(I,40)=WR(5,I)*W(I,37)+WR(11,I)*W(I,39)+                        
     &        2.0D0*WZ(5,I)*W(I,15)                                            
            W(I,42)=WR(1,I)*W(I,41)+WR(7,I)*W(I,44)+                         
     &        3.0D0*WZ(5,I)*W(I,20)                                            
            W(I,46)=WR(1,I)*W(I,45)+WR(7,I)*W(I,48)+                         
     &        2.0D0*WZ(5,I)*W(I,24)                                            
            W(I,47)=WR(3,I)*W(I,45)+WR(9,I)*W(I,48)+                         
     &        WZ(5,I)*W(I,20)                                                  
            W(I,50)=WR(1,I)*W(I,49)+WR(7,I)*W(I,52)+                         
     &        WZ(5,I)*W(I,28)                                                  
            W(I,51)=WR(3,I)*W(I,49)+WR(9,I)*W(I,52)+                         
     &        2.0D0*WZ(5,I)*W(I,24)                                            
            W(I,55)=WR(3,I)*W(I,53)+WR(9,I)*W(I,56)+                         
     &        3.0D0*WZ(5,I)*W(I,28)                                            
            W(I,58)=WR(1,I)*W(I,57)+WR(7,I)*W(I,59)+                         
     &        2.0D0*WZ(5,I)*W(I,31)                                            
            W(I,60)=WR(5,I)*W(I,57)+WR(11,I)*W(I,59)+                        
     &        WZ(5,I)*W(I,20)                                                  
            W(I,62)=WR(1,I)*W(I,61)+WR(7,I)*W(I,38)+                         
     &        WZ(5,I)*W(I,34)                                                  
            W(I,63)=WR(3,I)*W(I,61)+WR(9,I)*W(I,38)+                         
     &        WZ(5,I)*W(I,31)                                                  
            W(I,64)=WR(5,I)*W(I,61)+WR(11,I)*W(I,38)+                        
     &        WZ(5,I)*W(I,24)                                                  
            W(I,67)=WR(3,I)*W(I,65)+WR(9,I)*W(I,66)+                         
     &        2.0D0*WZ(5,I)*W(I,34)                                            
            W(I,68)=WR(5,I)*W(I,65)+WR(11,I)*W(I,66)+                        
     &        WZ(5,I)*W(I,28)                                                  
            W(I,70)=WR(1,I)*W(I,69)+WR(7,I)*W(I,71)+                         
     &        WZ(5,I)*W(I,39)                                                  
            W(I,72)=WR(5,I)*W(I,69)+WR(11,I)*W(I,71)+                        
     &        2.0D0*WZ(5,I)*W(I,31)                                            
            W(I,75)=WR(3,I)*W(I,73)+WR(9,I)*W(I,74)+                         
     &        WZ(5,I)*W(I,39)                                                  
            W(I,76)=WR(5,I)*W(I,73)+WR(11,I)*W(I,74)+                        
     &        2.0D0*WZ(5,I)*W(I,34)                                            
            W(I,80)=WR(5,I)*W(I,77)+WR(11,I)*W(I,79)+                        
     &        3.0D0*WZ(5,I)*W(I,39)                                            
            W(I,2)=WR(1,I)*W(I,1)+WR(7,I)*W(I,4)                             
            W(I,3)=WR(3,I)*W(I,1)+WR(9,I)*W(I,4)                             
            W(I,4)=WR(5,I)*W(I,1)+WR(11,I)*W(I,4)                            
            W(I,7)=WR(3,I)*W(I,5)+WR(9,I)*W(I,8)                             
            W(I,8)=WR(5,I)*W(I,5)+WR(11,I)*W(I,8)                            
            W(I,10)=WR(1,I)*W(I,9)+WR(7,I)*W(I,12)                           
            W(I,12)=WR(5,I)*W(I,9)+WR(11,I)*W(I,12)                          
            W(I,14)=WR(1,I)*W(I,13)+WR(7,I)*W(I,15)                          
            W(I,15)=WR(3,I)*W(I,13)+WR(9,I)*W(I,15)                          
            W(I,19)=WR(3,I)*W(I,17)+WR(9,I)*W(I,20)                          
            W(I,20)=WR(5,I)*W(I,17)+WR(11,I)*W(I,20)                         
            W(I,24)=WR(5,I)*W(I,21)+WR(11,I)*W(I,24)                         
            W(I,26)=WR(1,I)*W(I,25)+WR(7,I)*W(I,28)                          
            W(I,28)=WR(5,I)*W(I,25)+WR(11,I)*W(I,28)                         
            W(I,31)=WR(3,I)*W(I,29)+WR(9,I)*W(I,31)                          
            W(I,34)=WR(1,I)*W(I,33)+WR(7,I)*W(I,34)                          
            W(I,38)=WR(1,I)*W(I,37)+WR(7,I)*W(I,39)                          
            W(I,39)=WR(3,I)*W(I,37)+WR(9,I)*W(I,39)                          
            W(I,43)=WR(3,I)*W(I,41)+WR(9,I)*W(I,44)                          
            W(I,44)=WR(5,I)*W(I,41)+WR(11,I)*W(I,44)                         
            W(I,48)=WR(5,I)*W(I,45)+WR(11,I)*W(I,48)                         
            W(I,52)=WR(5,I)*W(I,49)+WR(11,I)*W(I,52)                         
            W(I,54)=WR(1,I)*W(I,53)+WR(7,I)*W(I,56)                          
            W(I,56)=WR(5,I)*W(I,53)+WR(11,I)*W(I,56)                         
            W(I,59)=WR(3,I)*W(I,57)+WR(9,I)*W(I,59)                          
            W(I,66)=WR(1,I)*W(I,65)+WR(7,I)*W(I,66)                          
            W(I,71)=WR(3,I)*W(I,69)+WR(9,I)*W(I,71)                          
            W(I,74)=WR(1,I)*W(I,73)+WR(7,I)*W(I,74)                          
            W(I,78)=WR(1,I)*W(I,77)+WR(7,I)*W(I,79)                          
            W(I,79)=WR(3,I)*W(I,77)+WR(9,I)*W(I,79)                          
         ENDDO                                                            
         RETURN                                                           
      ELSEIF(LCode.EQ.20)THEN                                          
         DO I=1,N                                                         
            W(I,1)=R(I,1)                                                    
            W(I,6)=R(I,2)                                                    
            W(I,5)=R(I,3)                                                    
            W(I,2)=WR(1,I)*W(I,1)+WR(7,I)*W(I,6)                             
            W(I,3)=WR(3,I)*W(I,1)+WR(9,I)*W(I,6)                             
            W(I,4)=WR(5,I)*W(I,1)+WR(11,I)*W(I,6)                            
            W(I,7)=WR(1,I)*W(I,6)+WR(7,I)*W(I,5)                             
            W(I,8)=WR(3,I)*W(I,6)+WR(9,I)*W(I,5)                             
            W(I,9)=WR(5,I)*W(I,6)+WR(11,I)*W(I,5)                            
            Ts0=WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,6))
            W(I,5)=WR(1,I)*W(I,2)+WR(7,I)*W(I,7)+ Ts0
            W(I,7)=WR(3,I)*W(I,3)+WR(9,I)*W(I,8)+ Ts0
            W(I,10)=WR(5,I)*W(I,4)+WR(11,I)*W(I,9)+ Ts0
            W(I,6)=WR(1,I)*W(I,3)+WR(7,I)*W(I,8)                             
            W(I,8)=WR(1,I)*W(I,4)+WR(7,I)*W(I,9)                             
            W(I,9)=WR(3,I)*W(I,4)+WR(9,I)*W(I,9)                             
         ENDDO                                                            
         RETURN                                                           
      ELSEIF(LCode.EQ.21)THEN                                          
         DO I=1,N                                                         
            W(I,1)=R(I,1)                                                    
            W(I,6)=R(I,2)                                                    
            W(I,5)=R(I,3)                                                    
            W(I,2)=R(I,4)                                                    
            W(I,18)=WR(2,I)*W(I,1)+WR(8,I)*W(I,6)                            
            W(I,35)=WR(4,I)*W(I,1)+WR(10,I)*W(I,6)                           
            W(I,52)=WR(6,I)*W(I,1)+WR(12,I)*W(I,6)                           
            W(I,8)=WR(2,I)*W(I,6)+WR(8,I)*W(I,5)                             
            W(I,9)=WR(4,I)*W(I,6)+WR(10,I)*W(I,5)                            
            W(I,11)=WR(6,I)*W(I,6)+WR(12,I)*W(I,5)                           
            W(I,7)=WR(2,I)*W(I,5)+WR(8,I)*W(I,2)                             
            W(I,10)=WR(4,I)*W(I,5)+WR(10,I)*W(I,2)                           
            W(I,12)=WR(6,I)*W(I,5)+WR(12,I)*W(I,2)                           
            W(I,2)=WR(1,I)*W(I,1)+WR(7,I)*W(I,6)                             
            W(I,3)=WR(3,I)*W(I,1)+WR(9,I)*W(I,6)                             
            W(I,4)=WR(5,I)*W(I,1)+WR(11,I)*W(I,6)                            
            W(I,13)=WR(1,I)*W(I,6)+WR(7,I)*W(I,5)                            
            W(I,14)=WR(3,I)*W(I,6)+WR(9,I)*W(I,5)                            
            W(I,15)=WR(5,I)*W(I,6)+WR(11,I)*W(I,5)                           
            W(I,19)=WR(1,I)*W(I,18)+WR(7,I)*W(I,8)+                          
     &        WZ(5,I)*W(I,6)                                                   
            W(I,20)=WR(3,I)*W(I,18)+WR(9,I)*W(I,8)                           
            W(I,21)=WR(5,I)*W(I,18)+WR(11,I)*W(I,8)                          
            W(I,36)=WR(1,I)*W(I,35)+WR(7,I)*W(I,9)                           
            W(I,37)=WR(3,I)*W(I,35)+WR(9,I)*W(I,9)+                          
     &        WZ(5,I)*W(I,6)                                                   
            W(I,38)=WR(5,I)*W(I,35)+WR(11,I)*W(I,9)                          
            W(I,53)=WR(1,I)*W(I,52)+WR(7,I)*W(I,11)                          
            W(I,54)=WR(3,I)*W(I,52)+WR(9,I)*W(I,11)                          
            W(I,55)=WR(5,I)*W(I,52)+WR(11,I)*W(I,11)+                        
     &        WZ(5,I)*W(I,6)                                                   
            W(I,16)=WR(1,I)*W(I,8)+WR(7,I)*W(I,7)+                           
     &        WZ(5,I)*W(I,5)                                                   
            W(I,17)=WR(3,I)*W(I,8)+WR(9,I)*W(I,7)                            
            W(I,26)=WR(5,I)*W(I,8)+WR(11,I)*W(I,7)                           
            W(I,28)=WR(1,I)*W(I,9)+WR(7,I)*W(I,10)                           
            W(I,29)=WR(3,I)*W(I,9)+WR(9,I)*W(I,10)+                          
     &        WZ(5,I)*W(I,5)                                                   
            W(I,30)=WR(5,I)*W(I,9)+WR(11,I)*W(I,10)                          
            W(I,31)=WR(1,I)*W(I,11)+WR(7,I)*W(I,12)                          
            W(I,32)=WR(3,I)*W(I,11)+WR(9,I)*W(I,12)                          
            W(I,33)=WR(5,I)*W(I,11)+WR(11,I)*W(I,12)+                        
     &        WZ(5,I)*W(I,5)                                                   
            W(I,22)=WR(1,I)*W(I,19)+WR(7,I)*W(I,16)+                         
     &        WZ(2,I)*(W(I,18)-WZ(4,I)*W(I,8))+                                
     &        WZ(5,I)*W(I,13)                                                  
            W(I,23)=WR(1,I)*W(I,20)+WR(7,I)*W(I,17)+                         
     &        WZ(5,I)*W(I,14)                                                  
            W(I,25)=WR(1,I)*W(I,21)+WR(7,I)*W(I,26)+                         
     &        WZ(5,I)*W(I,15)                                                  
            W(I,41)=WR(3,I)*W(I,37)+WR(9,I)*W(I,29)+                         
     &        WZ(2,I)*(W(I,35)-WZ(4,I)*W(I,9))+                                
     &        WZ(5,I)*W(I,14)                                                  
            W(I,43)=WR(3,I)*W(I,38)+WR(9,I)*W(I,30)+                         
     &        WZ(5,I)*W(I,15)                                                  
            W(I,61)=WR(5,I)*W(I,55)+WR(11,I)*W(I,33)+                        
     &        WZ(2,I)*(W(I,52)-WZ(4,I)*W(I,11))+                               
     &        WZ(5,I)*W(I,15)                                                  
            W(I,5)=WR(1,I)*W(I,2)+WR(7,I)*W(I,13)+                           
     &        WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,6))                                  
            W(I,7)=WR(3,I)*W(I,3)+WR(9,I)*W(I,14)+                           
     &        WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,6))                                  
            W(I,10)=WR(5,I)*W(I,4)+WR(11,I)*W(I,15)+                         
     &        WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,6))                                  
            W(I,24)=WR(3,I)*W(I,20)+WR(9,I)*W(I,17)+                         
     &        WZ(2,I)*(W(I,18)-WZ(4,I)*W(I,8))                                 
            W(I,27)=WR(5,I)*W(I,21)+WR(11,I)*W(I,26)+                        
     &        WZ(2,I)*(W(I,18)-WZ(4,I)*W(I,8))                                 
            W(I,39)=WR(1,I)*W(I,36)+WR(7,I)*W(I,28)+                         
     &        WZ(2,I)*(W(I,35)-WZ(4,I)*W(I,9))                                 
            W(I,44)=WR(5,I)*W(I,38)+WR(11,I)*W(I,30)+                        
     &        WZ(2,I)*(W(I,35)-WZ(4,I)*W(I,9))                                 
            W(I,56)=WR(1,I)*W(I,53)+WR(7,I)*W(I,31)+                         
     &        WZ(2,I)*(W(I,52)-WZ(4,I)*W(I,11))                                
            W(I,58)=WR(3,I)*W(I,54)+WR(9,I)*W(I,32)+                         
     &        WZ(2,I)*(W(I,52)-WZ(4,I)*W(I,11))                                
            W(I,6)=WR(1,I)*W(I,3)+WR(7,I)*W(I,14)                            
            W(I,8)=WR(1,I)*W(I,4)+WR(7,I)*W(I,15)                            
            W(I,9)=WR(3,I)*W(I,4)+WR(9,I)*W(I,15)                            
            W(I,26)=WR(3,I)*W(I,21)+WR(9,I)*W(I,26)                          
            W(I,40)=WR(1,I)*W(I,37)+WR(7,I)*W(I,29)                          
            W(I,42)=WR(1,I)*W(I,38)+WR(7,I)*W(I,30)                          
            W(I,57)=WR(1,I)*W(I,54)+WR(7,I)*W(I,32)                          
            W(I,59)=WR(1,I)*W(I,55)+WR(7,I)*W(I,33)                          
            W(I,60)=WR(3,I)*W(I,55)+WR(9,I)*W(I,33)                          
         ENDDO                                                            
         RETURN                                                           
      ELSEIF(LCode.EQ.22)THEN                                          
         DO I=1,N                                                         
            W(I,1)=R(I,1)                                                    
            W(I,6)=R(I,2)                                                    
            W(I,5)=R(I,3)                                                    
            W(I,2)=R(I,4)                                                    
            W(I,3)=R(I,5)                                                    
            W(I,18)=WR(2,I)*W(I,1)+WR(8,I)*W(I,6)                            
            W(I,35)=WR(4,I)*W(I,1)+WR(10,I)*W(I,6)                           
            W(I,52)=WR(6,I)*W(I,1)+WR(12,I)*W(I,6)                           
            W(I,8)=WR(2,I)*W(I,6)+WR(8,I)*W(I,5)                             
            W(I,9)=WR(4,I)*W(I,6)+WR(10,I)*W(I,5)                            
            W(I,11)=WR(6,I)*W(I,6)+WR(12,I)*W(I,5)                           
            W(I,7)=WR(2,I)*W(I,5)+WR(8,I)*W(I,2)                             
            W(I,10)=WR(4,I)*W(I,5)+WR(10,I)*W(I,2)                           
            W(I,12)=WR(6,I)*W(I,5)+WR(12,I)*W(I,2)                           
            W(I,4)=WR(2,I)*W(I,2)+WR(8,I)*W(I,3)                             
            W(I,13)=WR(4,I)*W(I,2)+WR(10,I)*W(I,3)                           
            W(I,14)=WR(6,I)*W(I,2)+WR(12,I)*W(I,3)                           
            W(I,69)=WR(2,I)*W(I,18)+WR(8,I)*W(I,8)+                          
     &        WZ(1,I)*(W(I,1)-WZ(3,I)*W(I,6))                                  
            W(I,86)=WR(2,I)*W(I,35)+WR(8,I)*W(I,9)                           
            W(I,103)=WR(4,I)*W(I,35)+WR(10,I)*W(I,9)+                        
     &        WZ(1,I)*(W(I,1)-WZ(3,I)*W(I,6))                                  
            W(I,120)=WR(2,I)*W(I,52)+WR(8,I)*W(I,11)                         
            W(I,137)=WR(4,I)*W(I,52)+WR(10,I)*W(I,11)                        
            W(I,154)=WR(6,I)*W(I,52)+WR(12,I)*W(I,11)+                       
     &        WZ(1,I)*(W(I,1)-WZ(3,I)*W(I,6))                                  
            W(I,15)=WR(2,I)*W(I,8)+WR(8,I)*W(I,7)+                           
     &        WZ(1,I)*(W(I,6)-WZ(3,I)*W(I,5))                                  
            W(I,16)=WR(2,I)*W(I,9)+WR(8,I)*W(I,10)                           
            W(I,17)=WR(4,I)*W(I,9)+WR(10,I)*W(I,10)+                         
     &        WZ(1,I)*(W(I,6)-WZ(3,I)*W(I,5))                                  
            W(I,26)=WR(2,I)*W(I,11)+WR(8,I)*W(I,12)                          
            W(I,28)=WR(4,I)*W(I,11)+WR(10,I)*W(I,12)                         
            W(I,29)=WR(6,I)*W(I,11)+WR(12,I)*W(I,12)+                        
     &        WZ(1,I)*(W(I,6)-WZ(3,I)*W(I,5))                                  
            W(I,22)=WR(2,I)*W(I,7)+WR(8,I)*W(I,4)+                           
     &        WZ(1,I)*(W(I,5)-WZ(3,I)*W(I,2))                                  
            W(I,23)=WR(2,I)*W(I,10)+WR(8,I)*W(I,13)                          
            W(I,24)=WR(4,I)*W(I,10)+WR(10,I)*W(I,13)+                        
     &        WZ(1,I)*(W(I,5)-WZ(3,I)*W(I,2))                                  
            W(I,25)=WR(2,I)*W(I,12)+WR(8,I)*W(I,14)                          
            W(I,27)=WR(4,I)*W(I,12)+WR(10,I)*W(I,14)                         
            W(I,30)=WR(6,I)*W(I,12)+WR(12,I)*W(I,14)+                        
     &        WZ(1,I)*(W(I,5)-WZ(3,I)*W(I,2))                                  
            W(I,2)=WR(1,I)*W(I,1)+WR(7,I)*W(I,6)                             
            W(I,3)=WR(3,I)*W(I,1)+WR(9,I)*W(I,6)                             
            W(I,4)=WR(5,I)*W(I,1)+WR(11,I)*W(I,6)                            
            W(I,31)=WR(1,I)*W(I,6)+WR(7,I)*W(I,5)                            
            W(I,32)=WR(3,I)*W(I,6)+WR(9,I)*W(I,5)                            
            W(I,33)=WR(5,I)*W(I,6)+WR(11,I)*W(I,5)                           
            W(I,19)=WR(1,I)*W(I,18)+WR(7,I)*W(I,8)+                          
     &        WZ(5,I)*W(I,6)                                                   
            W(I,20)=WR(3,I)*W(I,18)+WR(9,I)*W(I,8)                           
            W(I,21)=WR(5,I)*W(I,18)+WR(11,I)*W(I,8)                          
            W(I,36)=WR(1,I)*W(I,35)+WR(7,I)*W(I,9)                           
            W(I,37)=WR(3,I)*W(I,35)+WR(9,I)*W(I,9)+                          
     &        WZ(5,I)*W(I,6)                                                   
            W(I,38)=WR(5,I)*W(I,35)+WR(11,I)*W(I,9)                          
            W(I,53)=WR(1,I)*W(I,52)+WR(7,I)*W(I,11)                          
            W(I,54)=WR(3,I)*W(I,52)+WR(9,I)*W(I,11)                          
            W(I,55)=WR(5,I)*W(I,52)+WR(11,I)*W(I,11)+                        
     &        WZ(5,I)*W(I,6)                                                   
            W(I,34)=WR(1,I)*W(I,8)+WR(7,I)*W(I,7)+                           
     &        WZ(5,I)*W(I,5)                                                   
            W(I,39)=WR(3,I)*W(I,8)+WR(9,I)*W(I,7)                            
            W(I,40)=WR(5,I)*W(I,8)+WR(11,I)*W(I,7)                           
            W(I,42)=WR(1,I)*W(I,9)+WR(7,I)*W(I,10)                           
            W(I,45)=WR(3,I)*W(I,9)+WR(9,I)*W(I,10)+                          
     &        WZ(5,I)*W(I,5)                                                   
            W(I,46)=WR(5,I)*W(I,9)+WR(11,I)*W(I,10)                          
            W(I,47)=WR(1,I)*W(I,11)+WR(7,I)*W(I,12)                          
            W(I,48)=WR(3,I)*W(I,11)+WR(9,I)*W(I,12)                          
            W(I,49)=WR(5,I)*W(I,11)+WR(11,I)*W(I,12)+                        
     &        WZ(5,I)*W(I,5)                                                   
            W(I,70)=WR(1,I)*W(I,69)+WR(7,I)*W(I,15)+                         
     &        2.0D0*WZ(5,I)*W(I,8)                                             
            W(I,71)=WR(3,I)*W(I,69)+WR(9,I)*W(I,15)                          
            W(I,72)=WR(5,I)*W(I,69)+WR(11,I)*W(I,15)                         
            W(I,87)=WR(1,I)*W(I,86)+WR(7,I)*W(I,16)+                         
     &        WZ(5,I)*W(I,9)                                                   
            W(I,88)=WR(3,I)*W(I,86)+WR(9,I)*W(I,16)+                         
     &        WZ(5,I)*W(I,8)                                                   
            W(I,89)=WR(5,I)*W(I,86)+WR(11,I)*W(I,16)                         
            W(I,104)=WR(1,I)*W(I,103)+WR(7,I)*W(I,17)                        
            W(I,105)=WR(3,I)*W(I,103)+WR(9,I)*W(I,17)+                       
     &        2.0D0*WZ(5,I)*W(I,9)                                             
            W(I,106)=WR(5,I)*W(I,103)+WR(11,I)*W(I,17)                       
            W(I,121)=WR(1,I)*W(I,120)+WR(7,I)*W(I,26)+                       
     &        WZ(5,I)*W(I,11)                                                  
            W(I,122)=WR(3,I)*W(I,120)+WR(9,I)*W(I,26)                        
            W(I,123)=WR(5,I)*W(I,120)+WR(11,I)*W(I,26)+                      
     &        WZ(5,I)*W(I,8)                                                   
            W(I,138)=WR(1,I)*W(I,137)+WR(7,I)*W(I,28)                        
            W(I,139)=WR(3,I)*W(I,137)+WR(9,I)*W(I,28)+                       
     &        WZ(5,I)*W(I,11)                                                  
            W(I,140)=WR(5,I)*W(I,137)+WR(11,I)*W(I,28)+                      
     &        WZ(5,I)*W(I,9)                                                   
            W(I,155)=WR(1,I)*W(I,154)+WR(7,I)*W(I,29)                        
            W(I,156)=WR(3,I)*W(I,154)+WR(9,I)*W(I,29)                        
            W(I,157)=WR(5,I)*W(I,154)+WR(11,I)*W(I,29)+                      
     &        2.0D0*WZ(5,I)*W(I,11)                                            
            W(I,44)=WR(1,I)*W(I,15)+WR(7,I)*W(I,22)+                         
     &        2.0D0*WZ(5,I)*W(I,7)                                             
            W(I,50)=WR(3,I)*W(I,15)+WR(9,I)*W(I,22)                          
            W(I,51)=WR(5,I)*W(I,15)+WR(11,I)*W(I,22)                         
            W(I,56)=WR(1,I)*W(I,16)+WR(7,I)*W(I,23)+                         
     &        WZ(5,I)*W(I,10)                                                  
            W(I,57)=WR(3,I)*W(I,16)+WR(9,I)*W(I,23)+                         
     &        WZ(5,I)*W(I,7)                                                   
            W(I,59)=WR(5,I)*W(I,16)+WR(11,I)*W(I,23)                         
            W(I,60)=WR(1,I)*W(I,17)+WR(7,I)*W(I,24)                          
            W(I,62)=WR(3,I)*W(I,17)+WR(9,I)*W(I,24)+                         
     &        2.0D0*WZ(5,I)*W(I,10)                                            
            W(I,63)=WR(5,I)*W(I,17)+WR(11,I)*W(I,24)                         
            W(I,58)=WR(1,I)*W(I,26)+WR(7,I)*W(I,25)+                         
     &        WZ(5,I)*W(I,12)                                                  
            W(I,64)=WR(3,I)*W(I,26)+WR(9,I)*W(I,25)                          
            W(I,65)=WR(5,I)*W(I,26)+WR(11,I)*W(I,25)+                        
     &        WZ(5,I)*W(I,7)                                                   
            W(I,66)=WR(1,I)*W(I,28)+WR(7,I)*W(I,27)                          
            W(I,67)=WR(3,I)*W(I,28)+WR(9,I)*W(I,27)+                         
     &        WZ(5,I)*W(I,12)                                                  
            W(I,68)=WR(5,I)*W(I,28)+WR(11,I)*W(I,27)+                        
     &        WZ(5,I)*W(I,10)                                                  
            W(I,77)=WR(1,I)*W(I,29)+WR(7,I)*W(I,30)                          
            W(I,79)=WR(3,I)*W(I,29)+WR(9,I)*W(I,30)                          
            W(I,80)=WR(5,I)*W(I,29)+WR(11,I)*W(I,30)+                        
     &        2.0D0*WZ(5,I)*W(I,12)                                            
            W(I,22)=WR(1,I)*W(I,19)+WR(7,I)*W(I,34)+                         
     &        WZ(2,I)*(W(I,18)-WZ(4,I)*W(I,8))+                                
     &        WZ(5,I)*W(I,31)                                                  
            W(I,23)=WR(1,I)*W(I,20)+WR(7,I)*W(I,39)+                         
     &        WZ(5,I)*W(I,32)                                                  
            W(I,25)=WR(1,I)*W(I,21)+WR(7,I)*W(I,40)+                         
     &        WZ(5,I)*W(I,33)                                                  
            W(I,41)=WR(3,I)*W(I,37)+WR(9,I)*W(I,45)+                         
     &        WZ(2,I)*(W(I,35)-WZ(4,I)*W(I,9))+                                
     &        WZ(5,I)*W(I,32)                                                  
            W(I,43)=WR(3,I)*W(I,38)+WR(9,I)*W(I,46)+                         
     &        WZ(5,I)*W(I,33)                                                  
            W(I,61)=WR(5,I)*W(I,55)+WR(11,I)*W(I,49)+                        
     &        WZ(2,I)*(W(I,52)-WZ(4,I)*W(I,11))+                               
     &        WZ(5,I)*W(I,33)                                                  
            W(I,73)=WR(1,I)*W(I,70)+WR(7,I)*W(I,44)+                         
     &        WZ(2,I)*(W(I,69)-WZ(4,I)*W(I,15))+                               
     &        2.0D0*WZ(5,I)*W(I,34)                                            
            W(I,74)=WR(1,I)*W(I,71)+WR(7,I)*W(I,50)+                         
     &        2.0D0*WZ(5,I)*W(I,39)                                            
            W(I,76)=WR(1,I)*W(I,72)+WR(7,I)*W(I,51)+                         
     &        2.0D0*WZ(5,I)*W(I,40)                                            
            W(I,90)=WR(1,I)*W(I,87)+WR(7,I)*W(I,56)+                         
     &        WZ(2,I)*(W(I,86)-WZ(4,I)*W(I,16))+                               
     &        WZ(5,I)*W(I,42)                                                  
            W(I,91)=WR(1,I)*W(I,88)+WR(7,I)*W(I,57)+                         
     &        WZ(5,I)*W(I,45)                                                  
            W(I,92)=WR(3,I)*W(I,88)+WR(9,I)*W(I,57)+                         
     &        WZ(2,I)*(W(I,86)-WZ(4,I)*W(I,16))+                               
     &        WZ(5,I)*W(I,39)                                                  
            W(I,93)=WR(1,I)*W(I,89)+WR(7,I)*W(I,59)+                         
     &        WZ(5,I)*W(I,46)                                                  
            W(I,94)=WR(3,I)*W(I,89)+WR(9,I)*W(I,59)+                         
     &        WZ(5,I)*W(I,40)                                                  
            W(I,109)=WR(3,I)*W(I,105)+WR(9,I)*W(I,62)+                       
     &        WZ(2,I)*(W(I,103)-WZ(4,I)*W(I,17))+                              
     &        2.0D0*WZ(5,I)*W(I,45)                                            
            W(I,111)=WR(3,I)*W(I,106)+WR(9,I)*W(I,63)+                       
     &        2.0D0*WZ(5,I)*W(I,46)                                            
            W(I,124)=WR(1,I)*W(I,121)+WR(7,I)*W(I,58)+                       
     &        WZ(2,I)*(W(I,120)-WZ(4,I)*W(I,26))+                              
     &        WZ(5,I)*W(I,47)                                                  
            W(I,125)=WR(1,I)*W(I,122)+WR(7,I)*W(I,64)+                       
     &        WZ(5,I)*W(I,48)                                                  
            W(I,127)=WR(1,I)*W(I,123)+WR(7,I)*W(I,65)+                       
     &        WZ(5,I)*W(I,49)                                                  
            W(I,129)=WR(5,I)*W(I,123)+WR(11,I)*W(I,65)+                      
     &        WZ(2,I)*(W(I,120)-WZ(4,I)*W(I,26))+                              
     &        WZ(5,I)*W(I,40)                                                  
            W(I,143)=WR(3,I)*W(I,139)+WR(9,I)*W(I,67)+                       
     &        WZ(2,I)*(W(I,137)-WZ(4,I)*W(I,28))+                              
     &        WZ(5,I)*W(I,48)                                                  
            W(I,145)=WR(3,I)*W(I,140)+WR(9,I)*W(I,68)+                       
     &        WZ(5,I)*W(I,49)                                                  
            W(I,146)=WR(5,I)*W(I,140)+WR(11,I)*W(I,68)+                      
     &        WZ(2,I)*(W(I,137)-WZ(4,I)*W(I,28))+                              
     &        WZ(5,I)*W(I,46)                                                  
            W(I,163)=WR(5,I)*W(I,157)+WR(11,I)*W(I,80)+                      
     &        WZ(2,I)*(W(I,154)-WZ(4,I)*W(I,29))+                              
     &        2.0D0*WZ(5,I)*W(I,49)                                            
            W(I,5)=WR(1,I)*W(I,2)+WR(7,I)*W(I,31)+                           
     &        WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,6))                                  
            W(I,7)=WR(3,I)*W(I,3)+WR(9,I)*W(I,32)+                           
     &        WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,6))                                  
            W(I,10)=WR(5,I)*W(I,4)+WR(11,I)*W(I,33)+                         
     &        WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,6))                                  
            W(I,24)=WR(3,I)*W(I,20)+WR(9,I)*W(I,39)+                         
     &        WZ(2,I)*(W(I,18)-WZ(4,I)*W(I,8))                                 
            W(I,27)=WR(5,I)*W(I,21)+WR(11,I)*W(I,40)+                        
     &        WZ(2,I)*(W(I,18)-WZ(4,I)*W(I,8))                                 
            W(I,39)=WR(1,I)*W(I,36)+WR(7,I)*W(I,42)+                         
     &        WZ(2,I)*(W(I,35)-WZ(4,I)*W(I,9))                                 
            W(I,44)=WR(5,I)*W(I,38)+WR(11,I)*W(I,46)+                        
     &        WZ(2,I)*(W(I,35)-WZ(4,I)*W(I,9))                                 
            W(I,56)=WR(1,I)*W(I,53)+WR(7,I)*W(I,47)+                         
     &        WZ(2,I)*(W(I,52)-WZ(4,I)*W(I,11))                                
            W(I,58)=WR(3,I)*W(I,54)+WR(9,I)*W(I,48)+                         
     &        WZ(2,I)*(W(I,52)-WZ(4,I)*W(I,11))                                
            W(I,75)=WR(3,I)*W(I,71)+WR(9,I)*W(I,50)+                         
     &        WZ(2,I)*(W(I,69)-WZ(4,I)*W(I,15))                                
            W(I,78)=WR(5,I)*W(I,72)+WR(11,I)*W(I,51)+                        
     &        WZ(2,I)*(W(I,69)-WZ(4,I)*W(I,15))                                
            W(I,95)=WR(5,I)*W(I,89)+WR(11,I)*W(I,59)+                        
     &        WZ(2,I)*(W(I,86)-WZ(4,I)*W(I,16))                                
            W(I,107)=WR(1,I)*W(I,104)+WR(7,I)*W(I,60)+                       
     &        WZ(2,I)*(W(I,103)-WZ(4,I)*W(I,17))                               
            W(I,112)=WR(5,I)*W(I,106)+WR(11,I)*W(I,63)+                      
     &        WZ(2,I)*(W(I,103)-WZ(4,I)*W(I,17))                               
            W(I,126)=WR(3,I)*W(I,122)+WR(9,I)*W(I,64)+                       
     &        WZ(2,I)*(W(I,120)-WZ(4,I)*W(I,26))                               
            W(I,141)=WR(1,I)*W(I,138)+WR(7,I)*W(I,66)+                       
     &        WZ(2,I)*(W(I,137)-WZ(4,I)*W(I,28))                               
            W(I,158)=WR(1,I)*W(I,155)+WR(7,I)*W(I,77)+                       
     &        WZ(2,I)*(W(I,154)-WZ(4,I)*W(I,29))                               
            W(I,160)=WR(3,I)*W(I,156)+WR(9,I)*W(I,79)+                       
     &        WZ(2,I)*(W(I,154)-WZ(4,I)*W(I,29))                               
            W(I,6)=WR(1,I)*W(I,3)+WR(7,I)*W(I,32)                            
            W(I,8)=WR(1,I)*W(I,4)+WR(7,I)*W(I,33)                            
            W(I,9)=WR(3,I)*W(I,4)+WR(9,I)*W(I,33)                            
            W(I,26)=WR(3,I)*W(I,21)+WR(9,I)*W(I,40)                          
            W(I,40)=WR(1,I)*W(I,37)+WR(7,I)*W(I,45)                          
            W(I,42)=WR(1,I)*W(I,38)+WR(7,I)*W(I,46)                          
            W(I,57)=WR(1,I)*W(I,54)+WR(7,I)*W(I,48)                          
            W(I,59)=WR(1,I)*W(I,55)+WR(7,I)*W(I,49)                          
            W(I,60)=WR(3,I)*W(I,55)+WR(9,I)*W(I,49)                          
            W(I,77)=WR(3,I)*W(I,72)+WR(9,I)*W(I,51)                          
            W(I,108)=WR(1,I)*W(I,105)+WR(7,I)*W(I,62)                        
            W(I,110)=WR(1,I)*W(I,106)+WR(7,I)*W(I,63)                        
            W(I,128)=WR(3,I)*W(I,123)+WR(9,I)*W(I,65)                        
            W(I,142)=WR(1,I)*W(I,139)+WR(7,I)*W(I,67)                        
            W(I,144)=WR(1,I)*W(I,140)+WR(7,I)*W(I,68)                        
            W(I,159)=WR(1,I)*W(I,156)+WR(7,I)*W(I,79)                        
            W(I,161)=WR(1,I)*W(I,157)+WR(7,I)*W(I,80)                        
            W(I,162)=WR(3,I)*W(I,157)+WR(9,I)*W(I,80)                        
         ENDDO                                                            
         RETURN                                                           
      ELSEIF(LCode.EQ.30)THEN                                          
         DO I=1,N                                                         
            W(I,1)=R(I,1)                                                    
            W(I,11)=R(I,2)                                                   
            W(I,12)=R(I,3)                                                   
            W(I,5)=R(I,4)                                                    
            W(I,2)=WR(1,I)*W(I,1)+WR(7,I)*W(I,11)                            
            W(I,3)=WR(3,I)*W(I,1)+WR(9,I)*W(I,11)                            
            W(I,4)=WR(5,I)*W(I,1)+WR(11,I)*W(I,11)                           
            W(I,13)=WR(1,I)*W(I,11)+WR(7,I)*W(I,12)                          
            W(I,14)=WR(3,I)*W(I,11)+WR(9,I)*W(I,12)                          
            W(I,16)=WR(5,I)*W(I,11)+WR(11,I)*W(I,12)                         
            W(I,15)=WR(1,I)*W(I,12)+WR(7,I)*W(I,5)                           
            W(I,17)=WR(3,I)*W(I,12)+WR(9,I)*W(I,5)                           
            W(I,18)=WR(5,I)*W(I,12)+WR(11,I)*W(I,5)                          
            W(I,5)=WR(1,I)*W(I,2)+WR(7,I)*W(I,13)+                           
     &        WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,11))                                 
            W(I,6)=WR(1,I)*W(I,3)+WR(7,I)*W(I,14)                            
            W(I,7)=WR(3,I)*W(I,3)+WR(9,I)*W(I,14)+                           
     &        WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,11))                                 
            W(I,8)=WR(1,I)*W(I,4)+WR(7,I)*W(I,16)                            
            W(I,9)=WR(3,I)*W(I,4)+WR(9,I)*W(I,16)                            
            W(I,10)=WR(5,I)*W(I,4)+WR(11,I)*W(I,16)+                         
     &        WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,11))                                 
            W(I,19)=WR(1,I)*W(I,13)+WR(7,I)*W(I,15)+                         
     &        WZ(2,I)*(W(I,11)-WZ(4,I)*W(I,12))                                
            W(I,20)=WR(1,I)*W(I,14)+WR(7,I)*W(I,17)                          
            W(I,21)=WR(3,I)*W(I,14)+WR(9,I)*W(I,17)+                         
     &        WZ(2,I)*(W(I,11)-WZ(4,I)*W(I,12))                                
            W(I,22)=WR(1,I)*W(I,16)+WR(7,I)*W(I,18)                          
            W(I,23)=WR(3,I)*W(I,16)+WR(9,I)*W(I,18)                          
            W(I,24)=WR(5,I)*W(I,16)+WR(11,I)*W(I,18)+                        
     &        WZ(2,I)*(W(I,11)-WZ(4,I)*W(I,12))                                
            W(I,11)=WR(1,I)*W(I,5)+WR(7,I)*W(I,19)+                          
     &        2.0D0*WZ(2,I)*(W(I,2)-WZ(4,I)*W(I,13))                           
            W(I,12)=WR(1,I)*W(I,6)+WR(7,I)*W(I,20)+                          
     &        WZ(2,I)*(W(I,3)-WZ(4,I)*W(I,14))                                 
            W(I,14)=WR(3,I)*W(I,7)+WR(9,I)*W(I,21)+                          
     &        2.0D0*WZ(2,I)*(W(I,3)-WZ(4,I)*W(I,14))                           
            W(I,15)=WR(1,I)*W(I,8)+WR(7,I)*W(I,22)+                          
     &        WZ(2,I)*(W(I,4)-WZ(4,I)*W(I,16))                                 
            W(I,17)=WR(3,I)*W(I,9)+WR(9,I)*W(I,23)+                          
     &        WZ(2,I)*(W(I,4)-WZ(4,I)*W(I,16))                                 
            W(I,20)=WR(5,I)*W(I,10)+WR(11,I)*W(I,24)+                        
     &        2.0D0*WZ(2,I)*(W(I,4)-WZ(4,I)*W(I,16))                           
            W(I,13)=WR(1,I)*W(I,7)+WR(7,I)*W(I,21)                           
            W(I,16)=WR(1,I)*W(I,9)+WR(7,I)*W(I,23)                           
            W(I,18)=WR(1,I)*W(I,10)+WR(7,I)*W(I,24)                          
            W(I,19)=WR(3,I)*W(I,10)+WR(9,I)*W(I,24)                          
         ENDDO                                                            
         RETURN                                                           
      ELSEIF(LCode.EQ.31)THEN                                          
         DO I=1,N                                                         
            W(I,1)=R(I,1)                                                    
            W(I,11)=R(I,2)                                                   
            W(I,12)=R(I,3)                                                   
            W(I,5)=R(I,4)                                                    
            W(I,2)=R(I,5)                                                    
            W(I,25)=WR(2,I)*W(I,1)+WR(8,I)*W(I,11)                           
            W(I,49)=WR(4,I)*W(I,1)+WR(10,I)*W(I,11)                          
            W(I,73)=WR(6,I)*W(I,1)+WR(12,I)*W(I,11)                          
            W(I,13)=WR(2,I)*W(I,11)+WR(8,I)*W(I,12)                          
            W(I,14)=WR(4,I)*W(I,11)+WR(10,I)*W(I,12)                         
            W(I,15)=WR(6,I)*W(I,11)+WR(12,I)*W(I,12)                         
            W(I,16)=WR(2,I)*W(I,12)+WR(8,I)*W(I,5)                           
            W(I,17)=WR(4,I)*W(I,12)+WR(10,I)*W(I,5)                          
            W(I,18)=WR(6,I)*W(I,12)+WR(12,I)*W(I,5)                          
            W(I,6)=WR(2,I)*W(I,5)+WR(8,I)*W(I,2)                             
            W(I,7)=WR(4,I)*W(I,5)+WR(10,I)*W(I,2)                            
            W(I,8)=WR(6,I)*W(I,5)+WR(12,I)*W(I,2)                            
            W(I,2)=WR(1,I)*W(I,1)+WR(7,I)*W(I,11)                            
            W(I,3)=WR(3,I)*W(I,1)+WR(9,I)*W(I,11)                            
            W(I,4)=WR(5,I)*W(I,1)+WR(11,I)*W(I,11)                           
            W(I,19)=WR(1,I)*W(I,11)+WR(7,I)*W(I,12)                          
            W(I,20)=WR(3,I)*W(I,11)+WR(9,I)*W(I,12)                          
            W(I,21)=WR(5,I)*W(I,11)+WR(11,I)*W(I,12)                         
            W(I,22)=WR(1,I)*W(I,12)+WR(7,I)*W(I,5)                           
            W(I,23)=WR(3,I)*W(I,12)+WR(9,I)*W(I,5)                           
            W(I,24)=WR(5,I)*W(I,12)+WR(11,I)*W(I,5)                          
            W(I,26)=WR(1,I)*W(I,25)+WR(7,I)*W(I,13)+                         
     &        WZ(5,I)*W(I,11)                                                  
            W(I,27)=WR(3,I)*W(I,25)+WR(9,I)*W(I,13)                          
            W(I,28)=WR(5,I)*W(I,25)+WR(11,I)*W(I,13)                         
            W(I,50)=WR(1,I)*W(I,49)+WR(7,I)*W(I,14)                          
            W(I,51)=WR(3,I)*W(I,49)+WR(9,I)*W(I,14)+                         
     &        WZ(5,I)*W(I,11)                                                  
            W(I,52)=WR(5,I)*W(I,49)+WR(11,I)*W(I,14)                         
            W(I,74)=WR(1,I)*W(I,73)+WR(7,I)*W(I,15)                          
            W(I,75)=WR(3,I)*W(I,73)+WR(9,I)*W(I,15)                          
            W(I,76)=WR(5,I)*W(I,73)+WR(11,I)*W(I,15)+                        
     &        WZ(5,I)*W(I,11)                                                  
            W(I,35)=WR(1,I)*W(I,13)+WR(7,I)*W(I,16)+                         
     &        WZ(5,I)*W(I,12)                                                  
            W(I,38)=WR(3,I)*W(I,13)+WR(9,I)*W(I,16)                          
            W(I,43)=WR(5,I)*W(I,13)+WR(11,I)*W(I,16)                         
            W(I,45)=WR(1,I)*W(I,14)+WR(7,I)*W(I,17)                          
            W(I,46)=WR(3,I)*W(I,14)+WR(9,I)*W(I,17)+                         
     &        WZ(5,I)*W(I,12)                                                  
            W(I,47)=WR(5,I)*W(I,14)+WR(11,I)*W(I,17)                         
            W(I,48)=WR(1,I)*W(I,15)+WR(7,I)*W(I,18)                          
            W(I,61)=WR(3,I)*W(I,15)+WR(9,I)*W(I,18)                          
            W(I,64)=WR(5,I)*W(I,15)+WR(11,I)*W(I,18)+                        
     &        WZ(5,I)*W(I,12)                                                  
            W(I,36)=WR(1,I)*W(I,16)+WR(7,I)*W(I,6)+                          
     &        WZ(5,I)*W(I,5)                                                   
            W(I,37)=WR(3,I)*W(I,16)+WR(9,I)*W(I,6)                           
            W(I,39)=WR(5,I)*W(I,16)+WR(11,I)*W(I,6)                          
            W(I,40)=WR(1,I)*W(I,17)+WR(7,I)*W(I,7)                           
            W(I,41)=WR(3,I)*W(I,17)+WR(9,I)*W(I,7)+                          
     &        WZ(5,I)*W(I,5)                                                   
            W(I,42)=WR(5,I)*W(I,17)+WR(11,I)*W(I,7)                          
            W(I,44)=WR(1,I)*W(I,18)+WR(7,I)*W(I,8)                           
            W(I,59)=WR(3,I)*W(I,18)+WR(9,I)*W(I,8)                           
            W(I,60)=WR(5,I)*W(I,18)+WR(11,I)*W(I,8)+                         
     &        WZ(5,I)*W(I,5)                                                   
            W(I,5)=WR(1,I)*W(I,2)+WR(7,I)*W(I,19)+                           
     &        WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,11))                                 
            W(I,6)=WR(1,I)*W(I,3)+WR(7,I)*W(I,20)                            
            W(I,7)=WR(3,I)*W(I,3)+WR(9,I)*W(I,20)+                           
     &        WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,11))                                 
            W(I,8)=WR(1,I)*W(I,4)+WR(7,I)*W(I,21)                            
            W(I,9)=WR(3,I)*W(I,4)+WR(9,I)*W(I,21)                            
            W(I,10)=WR(5,I)*W(I,4)+WR(11,I)*W(I,21)+                         
     &        WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,11))                                 
            W(I,63)=WR(1,I)*W(I,19)+WR(7,I)*W(I,22)+                         
     &        WZ(2,I)*(W(I,11)-WZ(4,I)*W(I,12))                                
            W(I,66)=WR(1,I)*W(I,20)+WR(7,I)*W(I,23)                          
            W(I,69)=WR(3,I)*W(I,20)+WR(9,I)*W(I,23)+                         
     &        WZ(2,I)*(W(I,11)-WZ(4,I)*W(I,12))                                
            W(I,68)=WR(1,I)*W(I,21)+WR(7,I)*W(I,24)                          
            W(I,70)=WR(3,I)*W(I,21)+WR(9,I)*W(I,24)                          
            W(I,71)=WR(5,I)*W(I,21)+WR(11,I)*W(I,24)+                        
     &        WZ(2,I)*(W(I,11)-WZ(4,I)*W(I,12))                                
            W(I,29)=WR(1,I)*W(I,26)+WR(7,I)*W(I,35)+                         
     &        WZ(2,I)*(W(I,25)-WZ(4,I)*W(I,13))+                               
     &        WZ(5,I)*W(I,19)                                                  
            W(I,30)=WR(1,I)*W(I,27)+WR(7,I)*W(I,38)+                         
     &        WZ(5,I)*W(I,20)                                                  
            W(I,31)=WR(3,I)*W(I,27)+WR(9,I)*W(I,38)+                         
     &        WZ(2,I)*(W(I,25)-WZ(4,I)*W(I,13))                                
            W(I,32)=WR(1,I)*W(I,28)+WR(7,I)*W(I,43)+                         
     &        WZ(5,I)*W(I,21)                                                  
            W(I,33)=WR(3,I)*W(I,28)+WR(9,I)*W(I,43)                          
            W(I,34)=WR(5,I)*W(I,28)+WR(11,I)*W(I,43)+                        
     &        WZ(2,I)*(W(I,25)-WZ(4,I)*W(I,13))                                
            W(I,53)=WR(1,I)*W(I,50)+WR(7,I)*W(I,45)+                         
     &        WZ(2,I)*(W(I,49)-WZ(4,I)*W(I,14))                                
            W(I,54)=WR(1,I)*W(I,51)+WR(7,I)*W(I,46)                          
            W(I,55)=WR(3,I)*W(I,51)+WR(9,I)*W(I,46)+                         
     &        WZ(2,I)*(W(I,49)-WZ(4,I)*W(I,14))+                               
     &        WZ(5,I)*W(I,20)                                                  
            W(I,56)=WR(1,I)*W(I,52)+WR(7,I)*W(I,47)                          
            W(I,57)=WR(3,I)*W(I,52)+WR(9,I)*W(I,47)+                         
     &        WZ(5,I)*W(I,21)                                                  
            W(I,58)=WR(5,I)*W(I,52)+WR(11,I)*W(I,47)+                        
     &        WZ(2,I)*(W(I,49)-WZ(4,I)*W(I,14))                                
            W(I,77)=WR(1,I)*W(I,74)+WR(7,I)*W(I,48)+                         
     &        WZ(2,I)*(W(I,73)-WZ(4,I)*W(I,15))                                
            W(I,78)=WR(1,I)*W(I,75)+WR(7,I)*W(I,61)                          
            W(I,79)=WR(3,I)*W(I,75)+WR(9,I)*W(I,61)+                         
     &        WZ(2,I)*(W(I,73)-WZ(4,I)*W(I,15))                                
            W(I,80)=WR(1,I)*W(I,76)+WR(7,I)*W(I,64)                          
            W(I,81)=WR(3,I)*W(I,76)+WR(9,I)*W(I,64)                          
            W(I,82)=WR(5,I)*W(I,76)+WR(11,I)*W(I,64)+                        
     &        WZ(2,I)*(W(I,73)-WZ(4,I)*W(I,15))+                               
     &        WZ(5,I)*W(I,21)                                                  
            W(I,62)=WR(1,I)*W(I,35)+WR(7,I)*W(I,36)+                         
     &        WZ(2,I)*(W(I,13)-WZ(4,I)*W(I,16))+                               
     &        WZ(5,I)*W(I,22)                                                  
            W(I,65)=WR(1,I)*W(I,38)+WR(7,I)*W(I,37)+                         
     &        WZ(5,I)*W(I,23)                                                  
            W(I,72)=WR(3,I)*W(I,38)+WR(9,I)*W(I,37)+                         
     &        WZ(2,I)*(W(I,13)-WZ(4,I)*W(I,16))                                
            W(I,67)=WR(1,I)*W(I,43)+WR(7,I)*W(I,39)+                         
     &        WZ(5,I)*W(I,24)                                                  
            W(I,83)=WR(3,I)*W(I,43)+WR(9,I)*W(I,39)                          
            W(I,85)=WR(5,I)*W(I,43)+WR(11,I)*W(I,39)+                        
     &        WZ(2,I)*(W(I,13)-WZ(4,I)*W(I,16))                                
            W(I,84)=WR(1,I)*W(I,45)+WR(7,I)*W(I,40)+                         
     &        WZ(2,I)*(W(I,14)-WZ(4,I)*W(I,17))                                
            W(I,86)=WR(1,I)*W(I,46)+WR(7,I)*W(I,41)                          
            W(I,88)=WR(3,I)*W(I,46)+WR(9,I)*W(I,41)+                         
     &        WZ(2,I)*(W(I,14)-WZ(4,I)*W(I,17))+                               
     &        WZ(5,I)*W(I,23)                                                  
            W(I,87)=WR(1,I)*W(I,47)+WR(7,I)*W(I,42)                          
            W(I,90)=WR(3,I)*W(I,47)+WR(9,I)*W(I,42)+                         
     &        WZ(5,I)*W(I,24)                                                  
            W(I,91)=WR(5,I)*W(I,47)+WR(11,I)*W(I,42)+                        
     &        WZ(2,I)*(W(I,14)-WZ(4,I)*W(I,17))                                
            W(I,89)=WR(1,I)*W(I,48)+WR(7,I)*W(I,44)+                         
     &        WZ(2,I)*(W(I,15)-WZ(4,I)*W(I,18))                                
            W(I,93)=WR(1,I)*W(I,61)+WR(7,I)*W(I,59)                          
            W(I,94)=WR(3,I)*W(I,61)+WR(9,I)*W(I,59)+                         
     &        WZ(2,I)*(W(I,15)-WZ(4,I)*W(I,18))                                
            W(I,95)=WR(1,I)*W(I,64)+WR(7,I)*W(I,60)                          
            W(I,96)=WR(3,I)*W(I,64)+WR(9,I)*W(I,60)                          
            W(I,97)=WR(5,I)*W(I,64)+WR(11,I)*W(I,60)+                        
     &        WZ(2,I)*(W(I,15)-WZ(4,I)*W(I,18))+                               
     &        WZ(5,I)*W(I,24)                                                  
            W(I,35)=WR(1,I)*W(I,29)+WR(7,I)*W(I,62)+                         
     &        2.0D0*WZ(2,I)*(W(I,26)-WZ(4,I)*W(I,35))+                         
     &        WZ(5,I)*W(I,63)                                                  
            W(I,36)=WR(1,I)*W(I,30)+WR(7,I)*W(I,65)+                         
     &        WZ(2,I)*(W(I,27)-WZ(4,I)*W(I,38))+                               
     &        WZ(5,I)*W(I,66)                                                  
            W(I,37)=WR(1,I)*W(I,31)+WR(7,I)*W(I,72)+                         
     &        WZ(5,I)*W(I,69)                                                  
            W(I,39)=WR(1,I)*W(I,32)+WR(7,I)*W(I,67)+                         
     &        WZ(2,I)*(W(I,28)-WZ(4,I)*W(I,43))+                               
     &        WZ(5,I)*W(I,68)                                                  
            W(I,40)=WR(1,I)*W(I,33)+WR(7,I)*W(I,83)+                         
     &        WZ(5,I)*W(I,70)                                                  
            W(I,42)=WR(1,I)*W(I,34)+WR(7,I)*W(I,85)+                         
     &        WZ(5,I)*W(I,71)                                                  
            W(I,62)=WR(3,I)*W(I,55)+WR(9,I)*W(I,88)+                         
     &        2.0D0*WZ(2,I)*(W(I,51)-WZ(4,I)*W(I,46))+                         
     &        WZ(5,I)*W(I,69)                                                  
            W(I,65)=WR(3,I)*W(I,57)+WR(9,I)*W(I,90)+                         
     &        WZ(2,I)*(W(I,52)-WZ(4,I)*W(I,47))+                               
     &        WZ(5,I)*W(I,70)                                                  
            W(I,67)=WR(3,I)*W(I,58)+WR(9,I)*W(I,91)+                         
     &        WZ(5,I)*W(I,71)                                                  
            W(I,92)=WR(5,I)*W(I,82)+WR(11,I)*W(I,97)+                        
     &        2.0D0*WZ(2,I)*(W(I,76)-WZ(4,I)*W(I,64))+                         
     &        WZ(5,I)*W(I,71)                                                  
            W(I,11)=WR(1,I)*W(I,5)+WR(7,I)*W(I,63)+                          
     &        2.0D0*WZ(2,I)*(W(I,2)-WZ(4,I)*W(I,19))                           
            W(I,12)=WR(1,I)*W(I,6)+WR(7,I)*W(I,66)+                          
     &        WZ(2,I)*(W(I,3)-WZ(4,I)*W(I,20))                                 
            W(I,14)=WR(3,I)*W(I,7)+WR(9,I)*W(I,69)+                          
     &        2.0D0*WZ(2,I)*(W(I,3)-WZ(4,I)*W(I,20))                           
            W(I,15)=WR(1,I)*W(I,8)+WR(7,I)*W(I,68)+                          
     &        WZ(2,I)*(W(I,4)-WZ(4,I)*W(I,21))                                 
            W(I,17)=WR(3,I)*W(I,9)+WR(9,I)*W(I,70)+                          
     &        WZ(2,I)*(W(I,4)-WZ(4,I)*W(I,21))                                 
            W(I,20)=WR(5,I)*W(I,10)+WR(11,I)*W(I,71)+                        
     &        2.0D0*WZ(2,I)*(W(I,4)-WZ(4,I)*W(I,21))                           
            W(I,38)=WR(3,I)*W(I,31)+WR(9,I)*W(I,72)+                         
     &        2.0D0*WZ(2,I)*(W(I,27)-WZ(4,I)*W(I,38))                          
            W(I,41)=WR(3,I)*W(I,33)+WR(9,I)*W(I,83)+                         
     &        WZ(2,I)*(W(I,28)-WZ(4,I)*W(I,43))                                
            W(I,44)=WR(5,I)*W(I,34)+WR(11,I)*W(I,85)+                        
     &        2.0D0*WZ(2,I)*(W(I,28)-WZ(4,I)*W(I,43))                          
            W(I,59)=WR(1,I)*W(I,53)+WR(7,I)*W(I,84)+                         
     &        2.0D0*WZ(2,I)*(W(I,50)-WZ(4,I)*W(I,45))                          
            W(I,60)=WR(1,I)*W(I,54)+WR(7,I)*W(I,86)+                         
     &        WZ(2,I)*(W(I,51)-WZ(4,I)*W(I,46))                                
            W(I,63)=WR(1,I)*W(I,56)+WR(7,I)*W(I,87)+                         
     &        WZ(2,I)*(W(I,52)-WZ(4,I)*W(I,47))                                
            W(I,68)=WR(5,I)*W(I,58)+WR(11,I)*W(I,91)+                        
     &        2.0D0*WZ(2,I)*(W(I,52)-WZ(4,I)*W(I,47))                          
            W(I,83)=WR(1,I)*W(I,77)+WR(7,I)*W(I,89)+                         
     &        2.0D0*WZ(2,I)*(W(I,74)-WZ(4,I)*W(I,48))                          
            W(I,84)=WR(1,I)*W(I,78)+WR(7,I)*W(I,93)+                         
     &        WZ(2,I)*(W(I,75)-WZ(4,I)*W(I,61))                                
            W(I,86)=WR(3,I)*W(I,79)+WR(9,I)*W(I,94)+                         
     &        2.0D0*WZ(2,I)*(W(I,75)-WZ(4,I)*W(I,61))                          
            W(I,87)=WR(1,I)*W(I,80)+WR(7,I)*W(I,95)+                         
     &        WZ(2,I)*(W(I,76)-WZ(4,I)*W(I,64))                                
            W(I,89)=WR(3,I)*W(I,81)+WR(9,I)*W(I,96)+                         
     &        WZ(2,I)*(W(I,76)-WZ(4,I)*W(I,64))                                
            W(I,13)=WR(1,I)*W(I,7)+WR(7,I)*W(I,69)                           
            W(I,16)=WR(1,I)*W(I,9)+WR(7,I)*W(I,70)                           
            W(I,18)=WR(1,I)*W(I,10)+WR(7,I)*W(I,71)                          
            W(I,19)=WR(3,I)*W(I,10)+WR(9,I)*W(I,71)                          
            W(I,43)=WR(3,I)*W(I,34)+WR(9,I)*W(I,85)                          
            W(I,61)=WR(1,I)*W(I,55)+WR(7,I)*W(I,88)                          
            W(I,64)=WR(1,I)*W(I,57)+WR(7,I)*W(I,90)                          
            W(I,66)=WR(1,I)*W(I,58)+WR(7,I)*W(I,91)                          
            W(I,85)=WR(1,I)*W(I,79)+WR(7,I)*W(I,94)                          
            W(I,88)=WR(1,I)*W(I,81)+WR(7,I)*W(I,96)                          
            W(I,90)=WR(1,I)*W(I,82)+WR(7,I)*W(I,97)                          
            W(I,91)=WR(3,I)*W(I,82)+WR(9,I)*W(I,97)                          
         ENDDO                                                            
         RETURN                                                           
      ELSEIF(LCode.EQ.40)THEN                                          
         DO I=1,N                                                         
            W(I,1)=R(I,1)                                                    
            W(I,11)=R(I,2)                                                   
            W(I,12)=R(I,3)                                                   
            W(I,13)=R(I,4)                                                   
            W(I,5)=R(I,5)                                                    
            W(I,2)=WR(1,I)*W(I,1)+WR(7,I)*W(I,11)                            
            W(I,3)=WR(3,I)*W(I,1)+WR(9,I)*W(I,11)                            
            W(I,4)=WR(5,I)*W(I,1)+WR(11,I)*W(I,11)                           
            W(I,21)=WR(1,I)*W(I,11)+WR(7,I)*W(I,12)                          
            W(I,22)=WR(3,I)*W(I,11)+WR(9,I)*W(I,12)                          
            W(I,23)=WR(5,I)*W(I,11)+WR(11,I)*W(I,12)                         
            W(I,24)=WR(1,I)*W(I,12)+WR(7,I)*W(I,13)                          
            W(I,25)=WR(3,I)*W(I,12)+WR(9,I)*W(I,13)                          
            W(I,26)=WR(5,I)*W(I,12)+WR(11,I)*W(I,13)                         
            W(I,14)=WR(1,I)*W(I,13)+WR(7,I)*W(I,5)                           
            W(I,15)=WR(3,I)*W(I,13)+WR(9,I)*W(I,5)                           
            W(I,16)=WR(5,I)*W(I,13)+WR(11,I)*W(I,5)                          
            W(I,5)=WR(1,I)*W(I,2)+WR(7,I)*W(I,21)+                           
     &        WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,11))                                 
            W(I,6)=WR(1,I)*W(I,3)+WR(7,I)*W(I,22)                            
            W(I,7)=WR(3,I)*W(I,3)+WR(9,I)*W(I,22)+                           
     &        WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,11))                                 
            W(I,8)=WR(1,I)*W(I,4)+WR(7,I)*W(I,23)                            
            W(I,9)=WR(3,I)*W(I,4)+WR(9,I)*W(I,23)                            
            W(I,10)=WR(5,I)*W(I,4)+WR(11,I)*W(I,23)+                         
     &        WZ(2,I)*(W(I,1)-WZ(4,I)*W(I,11))                                 
            W(I,27)=WR(1,I)*W(I,21)+WR(7,I)*W(I,24)+                         
     &        WZ(2,I)*(W(I,11)-WZ(4,I)*W(I,12))                                
            W(I,28)=WR(1,I)*W(I,22)+WR(7,I)*W(I,25)                          
            W(I,29)=WR(3,I)*W(I,22)+WR(9,I)*W(I,25)+                         
     &        WZ(2,I)*(W(I,11)-WZ(4,I)*W(I,12))                                
            W(I,30)=WR(1,I)*W(I,23)+WR(7,I)*W(I,26)                          
            W(I,31)=WR(3,I)*W(I,23)+WR(9,I)*W(I,26)                          
            W(I,33)=WR(5,I)*W(I,23)+WR(11,I)*W(I,26)+                        
     &        WZ(2,I)*(W(I,11)-WZ(4,I)*W(I,12))                                
            W(I,32)=WR(1,I)*W(I,24)+WR(7,I)*W(I,14)+                         
     &        WZ(2,I)*(W(I,12)-WZ(4,I)*W(I,13))                                
            W(I,34)=WR(1,I)*W(I,25)+WR(7,I)*W(I,15)                          
            W(I,35)=WR(3,I)*W(I,25)+WR(9,I)*W(I,15)+                         
     &        WZ(2,I)*(W(I,12)-WZ(4,I)*W(I,13))                                
            W(I,36)=WR(1,I)*W(I,26)+WR(7,I)*W(I,16)                          
            W(I,37)=WR(3,I)*W(I,26)+WR(9,I)*W(I,16)                          
            W(I,38)=WR(5,I)*W(I,26)+WR(11,I)*W(I,16)+                        
     &        WZ(2,I)*(W(I,12)-WZ(4,I)*W(I,13))                                
            W(I,11)=WR(1,I)*W(I,5)+WR(7,I)*W(I,27)+                          
     &        2.0D0*WZ(2,I)*(W(I,2)-WZ(4,I)*W(I,21))                           
            W(I,12)=WR(1,I)*W(I,6)+WR(7,I)*W(I,28)+                          
     &        WZ(2,I)*(W(I,3)-WZ(4,I)*W(I,22))                                 
            W(I,13)=WR(1,I)*W(I,7)+WR(7,I)*W(I,29)                           
            W(I,14)=WR(3,I)*W(I,7)+WR(9,I)*W(I,29)+                          
     &        2.0D0*WZ(2,I)*(W(I,3)-WZ(4,I)*W(I,22))                           
            W(I,15)=WR(1,I)*W(I,8)+WR(7,I)*W(I,30)+                          
     &        WZ(2,I)*(W(I,4)-WZ(4,I)*W(I,23))                                 
            W(I,16)=WR(1,I)*W(I,9)+WR(7,I)*W(I,31)                           
            W(I,17)=WR(3,I)*W(I,9)+WR(9,I)*W(I,31)+                          
     &        WZ(2,I)*(W(I,4)-WZ(4,I)*W(I,23))                                 
            W(I,18)=WR(1,I)*W(I,10)+WR(7,I)*W(I,33)                          
            W(I,19)=WR(3,I)*W(I,10)+WR(9,I)*W(I,33)                          
            W(I,20)=WR(5,I)*W(I,10)+WR(11,I)*W(I,33)+                        
     &        2.0D0*WZ(2,I)*(W(I,4)-WZ(4,I)*W(I,23))                           
            W(I,39)=WR(1,I)*W(I,27)+WR(7,I)*W(I,32)+                         
     &        2.0D0*WZ(2,I)*(W(I,21)-WZ(4,I)*W(I,24))                          
            W(I,40)=WR(1,I)*W(I,28)+WR(7,I)*W(I,34)+                         
     &        WZ(2,I)*(W(I,22)-WZ(4,I)*W(I,25))                                
            W(I,41)=WR(1,I)*W(I,29)+WR(7,I)*W(I,35)                          
            W(I,42)=WR(3,I)*W(I,29)+WR(9,I)*W(I,35)+                         
     &        2.0D0*WZ(2,I)*(W(I,22)-WZ(4,I)*W(I,25))                          
            W(I,43)=WR(1,I)*W(I,30)+WR(7,I)*W(I,36)+                         
     &        WZ(2,I)*(W(I,23)-WZ(4,I)*W(I,26))                                
            W(I,44)=WR(1,I)*W(I,31)+WR(7,I)*W(I,37)                          
            W(I,45)=WR(3,I)*W(I,31)+WR(9,I)*W(I,37)+                         
     &        WZ(2,I)*(W(I,23)-WZ(4,I)*W(I,26))                                
            W(I,46)=WR(1,I)*W(I,33)+WR(7,I)*W(I,38)                          
            W(I,47)=WR(3,I)*W(I,33)+WR(9,I)*W(I,38)                          
            W(I,48)=WR(5,I)*W(I,33)+WR(11,I)*W(I,38)+                        
     &        2.0D0*WZ(2,I)*(W(I,23)-WZ(4,I)*W(I,26))                          
            W(I,21)=WR(1,I)*W(I,11)+WR(7,I)*W(I,39)+                         
     &        3.0D0*WZ(2,I)*(W(I,5)-WZ(4,I)*W(I,27))                           
            W(I,22)=WR(1,I)*W(I,12)+WR(7,I)*W(I,40)+                         
     &        2.0D0*WZ(2,I)*(W(I,6)-WZ(4,I)*W(I,28))                           
            W(I,23)=WR(1,I)*W(I,13)+WR(7,I)*W(I,41)+                         
     &        WZ(2,I)*(W(I,7)-WZ(4,I)*W(I,29))                                 
            W(I,25)=WR(3,I)*W(I,14)+WR(9,I)*W(I,42)+                         
     &        3.0D0*WZ(2,I)*(W(I,7)-WZ(4,I)*W(I,29))                           
            W(I,26)=WR(1,I)*W(I,15)+WR(7,I)*W(I,43)+                         
     &        2.0D0*WZ(2,I)*(W(I,8)-WZ(4,I)*W(I,30))                           
            W(I,27)=WR(1,I)*W(I,16)+WR(7,I)*W(I,44)+                         
     &        WZ(2,I)*(W(I,9)-WZ(4,I)*W(I,31))                                 
            W(I,29)=WR(3,I)*W(I,17)+WR(9,I)*W(I,45)+                         
     &        2.0D0*WZ(2,I)*(W(I,9)-WZ(4,I)*W(I,31))                           
            W(I,30)=WR(1,I)*W(I,18)+WR(7,I)*W(I,46)+                         
     &        WZ(2,I)*(W(I,10)-WZ(4,I)*W(I,33))                                
            W(I,32)=WR(3,I)*W(I,19)+WR(9,I)*W(I,47)+                         
     &        WZ(2,I)*(W(I,10)-WZ(4,I)*W(I,33))                                
            W(I,35)=WR(5,I)*W(I,20)+WR(11,I)*W(I,48)+                        
     &        3.0D0*WZ(2,I)*(W(I,10)-WZ(4,I)*W(I,33))                          
            W(I,24)=WR(1,I)*W(I,14)+WR(7,I)*W(I,42)                          
            W(I,28)=WR(1,I)*W(I,17)+WR(7,I)*W(I,45)                          
            W(I,31)=WR(1,I)*W(I,19)+WR(7,I)*W(I,47)                          
            W(I,33)=WR(1,I)*W(I,20)+WR(7,I)*W(I,48)                          
            W(I,34)=WR(3,I)*W(I,20)+WR(9,I)*W(I,48)                          
         ENDDO                                                            
      ENDIF
      RETURN                                                           
      END
