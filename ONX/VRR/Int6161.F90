SUBROUTINE Int6161(N,IntCode,CBra,CKet,DisBufB,PrmBufB,DB,IB,SB,C,U)
 USE DerivedTypes
 USE GlobalScalars
 USE PrettyPrint
 USE InvExp
 INCLUDE "Int6161.Inc"
 REAL(DOUBLE)  :: Cx,Cy,Cz
 REAL(DOUBLE)  :: Zeta,Up,Px,Py,Pz,cp1,cp2,cp3
 REAL(DOUBLE)  :: Eta, Uq,Qx,Qy,Qz,cq1,cq2,cq3
 REAL(DOUBLE)  :: EtaQx,EtaQy,EtaQz
 REAL(DOUBLE)  :: r1xZpE,ZpE,Rkk,r1xZ,r1x2Z,rpxZ
 REAL(DOUBLE)  :: PAx,PAy,PAz,PQx,PQy,PQz
 REAL(DOUBLE)  :: QBx,QBy,QBz,WQx,WQy,WQz
 REAL(DOUBLE)  :: Wx,Wy,Wz,WPx,WPy,WPz
 REAL(DOUBLE)  :: T1,T2,T3,T4,TwoT
 REAL(DOUBLE)  :: R1,R2,R3,R4,R5,R6,R7,R8,R9
 REAL(DOUBLE)  :: G0,G1,G2,G3,G4,G5,G6,G7,G8,ET
 REAL(DOUBLE)  :: Ts0,Ts1,r1xE,r1x2E,rpxE,r1x2ZE
 INTEGER       :: I,J,K,ME,MG,l
 INTEGER       :: I0,I1,I2
 IF(N.EQ.0) RETURN
 C=0.0D0
 Cx=DisBufB( 8)
 Cy=DisBufB( 9)
 Cz=DisBufB(10)
 IF(IntCode.EQ. 6010601) THEN
 DO I=1,N
 I0 = SB%SLPrm%I(I)
 I2 = SB%SLDis%I(I)-4
 DO J=1,CKet
 I1=I0+(J-1)*DB%MAXP
 Eta = DB%PrmBuf%D(I1)
 Qx  = DB%PrmBuf%D(I1+1)
 Qy  = DB%PrmBuf%D(I1+2)
 Qz  = DB%PrmBuf%D(I1+3)
 Uq  = DB%PrmBuf%D(I1+4)
 cq1 = DB%PrmBuf%D(I1+5)
 cq2 = DB%PrmBuf%D(I1+6)
 cq3 = DB%PrmBuf%D(I1+7)
 r1xE  = 1.0D0/Eta
 r1x2E = 0.5D0*r1xE
 EtaQx = Eta*Qx
 EtaQy = Eta*Qy
 EtaQz = Eta*Qz
 QBx = Qx-DB%DisBuf%D(I2+7)
 QBy = Qy-DB%DisBuf%D(I2+8)
 QBz = Qz-DB%DisBuf%D(I2+9)
 DO K=1,CBra
 Zeta = PrmBufB(1,K)
 Px   = PrmBufB(2,K)
 Py   = PrmBufB(3,K)
 Pz   = PrmBufB(4,K)
 Up   = PrmBufB(5,K)
 cp1  = PrmBufB(6,K)
 cp2  = PrmBufB(7,K)
 cp3  = PrmBufB(8,K)
 ZpE  = Zeta+Eta
 r1xZpE = 1.0D0/ZpE
 r1x2ZE = 0.5D0*r1xZpE
 r1xZ   = 1.0D0/Zeta
 r1x2Z  = 0.5D0*r1xZ
 rpxZ   = Eta*r1xZpE
 rpxE   = Zeta*r1xZpE
 Rkk  = Up*Uq*DSQRT(r1xZpE)
 PAx = Px-Cx
 PAy = Py-Cy
 PAz = Pz-Cz
 Wx  = (Zeta*Px+EtaQx)*r1xZpE
 Wy  = (Zeta*Py+EtaQy)*r1xZpE
 Wz  = (Zeta*Pz+EtaQz)*r1xZpE
 WPx = Wx-Px
 WPy = Wy-Py
 WPz = Wz-Pz
 WQx = Wx - Qx
 WQy = Wy - Qy
 WQz = Wz - Qz
 PQx = Px - Qx
 PQy = Py - Qy
 PQz = Pz - Qz
 T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
 INCLUDE "RGen4.Inc"
 U(1)=R1                                                          
 U(6)=R2                                                          
 U(5)=R3                                                          
 U(2)=R4                                                          
 U(3)=R5                                                          
 U(18)=QBx*U(1)+WQx*U(6)                                          
 U(35)=QBy*U(1)+WQy*U(6)                                          
 U(52)=QBz*U(1)+WQz*U(6)                                          
 U(8)=QBx*U(6)+WQx*U(5)                                           
 U(9)=QBy*U(6)+WQy*U(5)                                           
 U(11)=QBz*U(6)+WQz*U(5)                                          
 U(7)=QBx*U(5)+WQx*U(2)                                           
 U(10)=QBy*U(5)+WQy*U(2)                                          
 U(12)=QBz*U(5)+WQz*U(2)                                          
 U(4)=QBx*U(2)+WQx*U(3)                                           
 U(13)=QBy*U(2)+WQy*U(3)                                          
 U(14)=QBz*U(2)+WQz*U(3)                                          
 U(69)=QBx*U(18)+WQx*U(8)+r1x2E*(U(1)-rpxE*U(6))                                           
 U(86)=QBx*U(35)+WQx*U(9)                                         
 U(103)=QBy*U(35)+WQy*U(9)+r1x2E*(U(1)-rpxE*U(6))                                          
 U(120)=QBx*U(52)+WQx*U(11)                                       
 U(137)=QBy*U(52)+WQy*U(11)                                       
 U(154)=QBz*U(52)+WQz*U(11)+r1x2E*(U(1)-rpxE*U(6))                                         
 U(15)=QBx*U(8)+WQx*U(7)+r1x2E*(U(6)-rpxE*U(5))                                            
 U(16)=QBx*U(9)+WQx*U(10)                                         
 U(17)=QBy*U(9)+WQy*U(10)+r1x2E*(U(6)-rpxE*U(5))                                           
 U(26)=QBx*U(11)+WQx*U(12)                                        
 U(28)=QBy*U(11)+WQy*U(12)                                        
 U(29)=QBz*U(11)+WQz*U(12)+r1x2E*(U(6)-rpxE*U(5))                                          
 U(22)=QBx*U(7)+WQx*U(4)+r1x2E*(U(5)-rpxE*U(2))                                            
 U(23)=QBx*U(10)+WQx*U(13)                                        
 U(24)=QBy*U(10)+WQy*U(13)+r1x2E*(U(5)-rpxE*U(2))                                          
 U(25)=QBx*U(12)+WQx*U(14)                                        
 U(27)=QBy*U(12)+WQy*U(14)                                        
 U(30)=QBz*U(12)+WQz*U(14)+r1x2E*(U(5)-rpxE*U(2))                                          
 U(2)=PAx*U(1)+WPx*U(6)                                           
 U(3)=PAy*U(1)+WPy*U(6)                                           
 U(4)=PAz*U(1)+WPz*U(6)                                           
 U(31)=PAx*U(6)+WPx*U(5)                                          
 U(32)=PAy*U(6)+WPy*U(5)                                          
 U(33)=PAz*U(6)+WPz*U(5)                                          
 U(19)=PAx*U(18)+WPx*U(8)+r1x2ZE*U(6)                                                      
 U(20)=PAy*U(18)+WPy*U(8)                                         
 U(21)=PAz*U(18)+WPz*U(8)                                         
 U(36)=PAx*U(35)+WPx*U(9)                                         
 U(37)=PAy*U(35)+WPy*U(9)+r1x2ZE*U(6)                                                      
 U(38)=PAz*U(35)+WPz*U(9)                                         
 U(53)=PAx*U(52)+WPx*U(11)                                        
 U(54)=PAy*U(52)+WPy*U(11)                                        
 U(55)=PAz*U(52)+WPz*U(11)+r1x2ZE*U(6)                                                     
 U(34)=PAx*U(8)+WPx*U(7)+r1x2ZE*U(5)                                                       
 U(39)=PAy*U(8)+WPy*U(7)                                          
 U(40)=PAz*U(8)+WPz*U(7)                                          
 U(42)=PAx*U(9)+WPx*U(10)                                         
 U(45)=PAy*U(9)+WPy*U(10)+r1x2ZE*U(5)                                                      
 U(46)=PAz*U(9)+WPz*U(10)                                         
 U(47)=PAx*U(11)+WPx*U(12)                                        
 U(48)=PAy*U(11)+WPy*U(12)                                        
 U(49)=PAz*U(11)+WPz*U(12)+r1x2ZE*U(5)                                                     
 U(70)=PAx*U(69)+WPx*U(15)+2.0D0*r1x2ZE*U(8)                                               
 U(71)=PAy*U(69)+WPy*U(15)                                        
 U(72)=PAz*U(69)+WPz*U(15)                                        
 U(87)=PAx*U(86)+WPx*U(16)+r1x2ZE*U(9)                                                     
 U(88)=PAy*U(86)+WPy*U(16)+r1x2ZE*U(8)                                                     
 U(89)=PAz*U(86)+WPz*U(16)                                        
 U(104)=PAx*U(103)+WPx*U(17)                                      
 U(105)=PAy*U(103)+WPy*U(17)+2.0D0*r1x2ZE*U(9)                                             
 U(106)=PAz*U(103)+WPz*U(17)                                      
 U(121)=PAx*U(120)+WPx*U(26)+r1x2ZE*U(11)                                                  
 U(122)=PAy*U(120)+WPy*U(26)                                      
 U(123)=PAz*U(120)+WPz*U(26)+r1x2ZE*U(8)                                                   
 U(138)=PAx*U(137)+WPx*U(28)                                      
 U(139)=PAy*U(137)+WPy*U(28)+r1x2ZE*U(11)                                                  
 U(140)=PAz*U(137)+WPz*U(28)+r1x2ZE*U(9)                                                   
 U(155)=PAx*U(154)+WPx*U(29)                                      
 U(156)=PAy*U(154)+WPy*U(29)                                      
 U(157)=PAz*U(154)+WPz*U(29)+2.0D0*r1x2ZE*U(11)                                            
 U(44)=PAx*U(15)+WPx*U(22)+2.0D0*r1x2ZE*U(7)                                               
 U(50)=PAy*U(15)+WPy*U(22)                                        
 U(51)=PAz*U(15)+WPz*U(22)                                        
 U(56)=PAx*U(16)+WPx*U(23)+r1x2ZE*U(10)                                                    
 U(57)=PAy*U(16)+WPy*U(23)+r1x2ZE*U(7)                                                     
 U(59)=PAz*U(16)+WPz*U(23)                                        
 U(60)=PAx*U(17)+WPx*U(24)                                        
 U(62)=PAy*U(17)+WPy*U(24)+2.0D0*r1x2ZE*U(10)                                              
 U(63)=PAz*U(17)+WPz*U(24)                                        
 U(58)=PAx*U(26)+WPx*U(25)+r1x2ZE*U(12)                                                    
 U(64)=PAy*U(26)+WPy*U(25)                                        
 U(65)=PAz*U(26)+WPz*U(25)+r1x2ZE*U(7)                                                     
 U(66)=PAx*U(28)+WPx*U(27)                                        
 U(67)=PAy*U(28)+WPy*U(27)+r1x2ZE*U(12)                                                    
 U(68)=PAz*U(28)+WPz*U(27)+r1x2ZE*U(10)                                                    
 U(77)=PAx*U(29)+WPx*U(30)                                        
 U(79)=PAy*U(29)+WPy*U(30)                                        
 U(80)=PAz*U(29)+WPz*U(30)+2.0D0*r1x2ZE*U(12)                                              
 U(22)=PAx*U(19)+WPx*U(34)+r1x2Z*(U(18)-rpxZ*U(8))+r1x2ZE*U(31)                            
 U(23)=PAx*U(20)+WPx*U(39)+r1x2ZE*U(32)                                                    
 U(25)=PAx*U(21)+WPx*U(40)+r1x2ZE*U(33)                                                    
 U(41)=PAy*U(37)+WPy*U(45)+r1x2Z*(U(35)-rpxZ*U(9))+r1x2ZE*U(32)                            
 U(43)=PAy*U(38)+WPy*U(46)+r1x2ZE*U(33)                                                    
 U(61)=PAz*U(55)+WPz*U(49)+r1x2Z*(U(52)-rpxZ*U(11))+r1x2ZE*U(33)                           
 U(73)=PAx*U(70)+WPx*U(44)+r1x2Z*(U(69)-rpxZ*U(15))+2.0D0*r1x2ZE*U(34)                     
 U(74)=PAx*U(71)+WPx*U(50)+2.0D0*r1x2ZE*U(39)                                              
 U(76)=PAx*U(72)+WPx*U(51)+2.0D0*r1x2ZE*U(40)                                              
 U(90)=PAx*U(87)+WPx*U(56)+r1x2Z*(U(86)-rpxZ*U(16))+r1x2ZE*U(42)                           
 U(91)=PAx*U(88)+WPx*U(57)+r1x2ZE*U(45)                                                    
 U(92)=PAy*U(88)+WPy*U(57)+r1x2Z*(U(86)-rpxZ*U(16))+r1x2ZE*U(39)                           
 U(93)=PAx*U(89)+WPx*U(59)+r1x2ZE*U(46)                                                    
 U(94)=PAy*U(89)+WPy*U(59)+r1x2ZE*U(40)                                                    
 U(109)=PAy*U(105)+WPy*U(62)+r1x2Z*(U(103)-rpxZ*U(17))+2.0D0*r1x2ZE*U(45)                  
 U(111)=PAy*U(106)+WPy*U(63)+2.0D0*r1x2ZE*U(46)                                            
 U(124)=PAx*U(121)+WPx*U(58)+r1x2Z*(U(120)-rpxZ*U(26))+r1x2ZE*U(47)                        
 U(125)=PAx*U(122)+WPx*U(64)+r1x2ZE*U(48)                                                  
 U(127)=PAx*U(123)+WPx*U(65)+r1x2ZE*U(49)                                                  
 U(129)=PAz*U(123)+WPz*U(65)+r1x2Z*(U(120)-rpxZ*U(26))+r1x2ZE*U(40)                        
 U(143)=PAy*U(139)+WPy*U(67)+r1x2Z*(U(137)-rpxZ*U(28))+r1x2ZE*U(48)                        
 U(145)=PAy*U(140)+WPy*U(68)+r1x2ZE*U(49)                                                  
 U(146)=PAz*U(140)+WPz*U(68)+r1x2Z*(U(137)-rpxZ*U(28))+r1x2ZE*U(46)                        
 U(163)=PAz*U(157)+WPz*U(80)+r1x2Z*(U(154)-rpxZ*U(29))+2.0D0*r1x2ZE*U(49)                  
 U(5)=PAx*U(2)+WPx*U(31)+r1x2Z*(U(1)-rpxZ*U(6))                                            
 U(7)=PAy*U(3)+WPy*U(32)+r1x2Z*(U(1)-rpxZ*U(6))                                            
 U(10)=PAz*U(4)+WPz*U(33)+r1x2Z*(U(1)-rpxZ*U(6))                                           
 U(24)=PAy*U(20)+WPy*U(39)+r1x2Z*(U(18)-rpxZ*U(8))                                         
 U(27)=PAz*U(21)+WPz*U(40)+r1x2Z*(U(18)-rpxZ*U(8))                                         
 U(39)=PAx*U(36)+WPx*U(42)+r1x2Z*(U(35)-rpxZ*U(9))                                         
 U(44)=PAz*U(38)+WPz*U(46)+r1x2Z*(U(35)-rpxZ*U(9))                                         
 U(56)=PAx*U(53)+WPx*U(47)+r1x2Z*(U(52)-rpxZ*U(11))                                        
 U(58)=PAy*U(54)+WPy*U(48)+r1x2Z*(U(52)-rpxZ*U(11))                                        
 U(75)=PAy*U(71)+WPy*U(50)+r1x2Z*(U(69)-rpxZ*U(15))                                        
 U(78)=PAz*U(72)+WPz*U(51)+r1x2Z*(U(69)-rpxZ*U(15))                                        
 U(95)=PAz*U(89)+WPz*U(59)+r1x2Z*(U(86)-rpxZ*U(16))                                        
 U(107)=PAx*U(104)+WPx*U(60)+r1x2Z*(U(103)-rpxZ*U(17))                                     
 U(112)=PAz*U(106)+WPz*U(63)+r1x2Z*(U(103)-rpxZ*U(17))                                     
 U(126)=PAy*U(122)+WPy*U(64)+r1x2Z*(U(120)-rpxZ*U(26))                                     
 U(141)=PAx*U(138)+WPx*U(66)+r1x2Z*(U(137)-rpxZ*U(28))                                     
 U(158)=PAx*U(155)+WPx*U(77)+r1x2Z*(U(154)-rpxZ*U(29))                                     
 U(160)=PAy*U(156)+WPy*U(79)+r1x2Z*(U(154)-rpxZ*U(29))                                     
 U(6)=PAx*U(3)+WPx*U(32)                                          
 U(8)=PAx*U(4)+WPx*U(33)                                          
 U(9)=PAy*U(4)+WPy*U(33)                                          
 U(26)=PAy*U(21)+WPy*U(40)                                        
 U(40)=PAx*U(37)+WPx*U(45)                                        
 U(42)=PAx*U(38)+WPx*U(46)                                        
 U(57)=PAx*U(54)+WPx*U(48)                                        
 U(59)=PAx*U(55)+WPx*U(49)                                        
 U(60)=PAy*U(55)+WPy*U(49)                                        
 U(77)=PAy*U(72)+WPy*U(51)                                        
 U(108)=PAx*U(105)+WPx*U(62)                                      
 U(110)=PAx*U(106)+WPx*U(63)                                      
 U(128)=PAy*U(123)+WPy*U(65)                                      
 U(142)=PAx*U(139)+WPx*U(67)                                      
 U(144)=PAx*U(140)+WPx*U(68)                                      
 U(159)=PAx*U(156)+WPx*U(79)                                      
 U(161)=PAx*U(157)+WPx*U(80)                                      
 U(162)=PAy*U(157)+WPy*U(80)                                      
 C(I,1)=C(I,1)+U(73)                                              
 C(I,2)=C(I,2)+U(74)                                              
 C(I,3)=C(I,3)+U(75)                                              
 C(I,4)=C(I,4)+U(76)                                              
 C(I,5)=C(I,5)+U(77)                                              
 C(I,6)=C(I,6)+U(78)                                              
 C(I,18)=C(I,18)+U(90)                                            
 C(I,19)=C(I,19)+U(91)                                            
 C(I,20)=C(I,20)+U(92)                                            
 C(I,21)=C(I,21)+U(93)                                            
 C(I,22)=C(I,22)+U(94)                                            
 C(I,23)=C(I,23)+U(95)                                            
 C(I,35)=C(I,35)+U(107)                                           
 C(I,36)=C(I,36)+U(108)                                           
 C(I,37)=C(I,37)+U(109)                                           
 C(I,38)=C(I,38)+U(110)                                           
 C(I,39)=C(I,39)+U(111)                                           
 C(I,40)=C(I,40)+U(112)                                           
 C(I,52)=C(I,52)+U(124)                                           
 C(I,53)=C(I,53)+U(125)                                           
 C(I,54)=C(I,54)+U(126)                                           
 C(I,55)=C(I,55)+U(127)                                           
 C(I,56)=C(I,56)+U(128)                                           
 C(I,57)=C(I,57)+U(129)                                           
 C(I,69)=C(I,69)+U(141)                                           
 C(I,70)=C(I,70)+U(142)                                           
 C(I,71)=C(I,71)+U(143)                                           
 C(I,72)=C(I,72)+U(144)                                           
 C(I,73)=C(I,73)+U(145)                                           
 C(I,74)=C(I,74)+U(146)                                           
 C(I,86)=C(I,86)+U(158)                                           
 C(I,87)=C(I,87)+U(159)                                           
 C(I,88)=C(I,88)+U(160)                                           
 C(I,89)=C(I,89)+U(161)                                           
 C(I,90)=C(I,90)+U(162)                                           
 C(I,91)=C(I,91)+U(163)                                           
 ENDDO
 ENDDO
 ENDDO
 ELSE
  CALL Halt("Illegal IntCode")
 ENDIF
END SUBROUTINE Int6161
