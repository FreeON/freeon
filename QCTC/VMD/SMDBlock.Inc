         SELECT CASE(LCode)
            ELSEIF(LBra.EQ.0.AND.LKet.EQ.0)THEN                          
               CALL VMD0(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                  
     $              rQ(oq+1,1),rQ(oq+1,2),rQ(oq+1,3),rQ(oq+1,4),         
     $              RhoCo(or+1),HGKet)                                   
            ELSEIF(LBra.EQ.0.AND.LKet.EQ.1)THEN                          
               CALL VMD1(NQ(zq),Px,Py,Pz,Omega,Upq,Tol,                  
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
           CASE(12:)                                                     
              CALL QCrash(
           END SELECT
