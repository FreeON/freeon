            IF(LCode==0) THEN
               CALL VMD0(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                         Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                         Rho%Est%D(oq),Rho%Co%D(or),HGKet%D) 
            ELSEIF(LCode== 1) THEN               
               CALL VMD1(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                         Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                         Rho%Est%D(oq),Rho%Co%D(or),HGKet%D) 
            ELSEIF(LCode== 2) THEN  
               CALL VMD2(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                         Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                         Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode== 3) THEN  
               CALL VMD3(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                         Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                         Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode== 4) THEN  
               CALL VMD4(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                         Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                         Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode==10) THEN  
               CALL VMD10(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode==11) THEN  
               CALL VMD11(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode==12) THEN 
               CALL VMD12(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)   
            ELSEIF(LCode==13) THEN  
               CALL VMD13(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode==14) THEN  
               CALL VMD14(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode==20) THEN  
               CALL VMD20(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode==21) THEN  
               CALL VMD21(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode==22) THEN  
               CALL VMD22(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode==23) THEN  
               CALL VMD23(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode==24) THEN 
               CALL VMD24(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D) 
            ELSEIF(LCode==30) THEN  
               CALL VMD30(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                         Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                         Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)      
            ELSEIF(LCode==31) THEN  
               CALL VMD31(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                         Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                         Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode==32) THEN  
               CALL VMD32(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode==33) THEN  
               CALL VMD33(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode==34) THEN  
               CALL VMD34(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)    
            ELSEIF(LCode==40) THEN  
               CALL VMD40(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode==41) THEN  
               CALL VMD41(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode==42) THEN  
               CALL VMD42(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSEIF(LCode==43) THEN    
               CALL VMD43(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D) 
            ELSEIF(LCode==44) THEN    
               CALL VMD44(NQ,Px,Py,Pz,Omega,Upq,Tol, &
                          Rho%Qx%D(oq),Rho%Qy%D(oq),Rho%Qz%D(oq), &
                          Rho%Est%D(oq),Rho%Co%D(or),HGKet%D)  
            ELSE
                CALL Halt(' No code generated for LCode = '//TRIM(IntToChar(LCode)))
!               CALL VMD(zq,Px,Py,Pz,Omega,Upq,Tol,Rho)  
            ENDIF
