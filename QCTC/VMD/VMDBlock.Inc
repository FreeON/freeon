         SELECT CASE(LCode)
            ELSEIF(LBra.EQ.0.AND.LKet.EQ.0)THEN                          
               CALL VMD0(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                  
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.0.AND.LKet.EQ.1)THEN                          
               CALL VMD1(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                  
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.0.AND.LKet.EQ.2)THEN                          
               CALL VMD2(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                  
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.0.AND.LKet.EQ.3)THEN                          
               CALL VMD3(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                  
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.0.AND.LKet.EQ.4)THEN                          
               CALL VMD4(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                  
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.1.AND.LKet.EQ.0)THEN                          
               CALL VMD10(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.1.AND.LKet.EQ.1)THEN                          
               CALL VMD11(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.1.AND.LKet.EQ.2)THEN                          
               CALL VMD12(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.1.AND.LKet.EQ.3)THEN                          
               CALL VMD13(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.1.AND.LKet.EQ.4)THEN                          
               CALL VMD14(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.2.AND.LKet.EQ.0)THEN                          
               CALL VMD20(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.2.AND.LKet.EQ.1)THEN                          
               CALL VMD21(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.2.AND.LKet.EQ.2)THEN                          
               CALL VMD22(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.2.AND.LKet.EQ.3)THEN                          
               CALL VMD23(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.2.AND.LKet.EQ.4)THEN                          
               CALL VMD24(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.3.AND.LKet.EQ.0)THEN                          
               CALL VMD30(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.3.AND.LKet.EQ.1)THEN                          
               CALL VMD31(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.3.AND.LKet.EQ.2)THEN                          
               CALL VMD32(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.3.AND.LKet.EQ.3)THEN                          
               CALL VMD33(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.3.AND.LKet.EQ.4)THEN                          
               CALL VMD34(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.4.AND.LKet.EQ.0)THEN                          
               CALL VMD40(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.4.AND.LKet.EQ.1)THEN                          
               CALL VMD41(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.4.AND.LKet.EQ.2)THEN                          
               CALL VMD42(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.4.AND.LKet.EQ.3)THEN                          
               CALL VMD43(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.4.AND.LKet.EQ.4)THEN                          
               CALL VMD44(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                 
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
           CASE(45:)                                                     
              CALL QCrash(
           END SELECT
