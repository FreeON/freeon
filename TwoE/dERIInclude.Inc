   SELECT CASE(IntType)
   CASE(1010101)
   CALL dInt1111(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(3010101)
   CALL dInt3111(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(6010101)
   CALL dInt6111(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(3030101)
   CALL dInt3311(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(6030101)
   CALL dInt6311(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(6060101)
   CALL dInt6611(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(3010301)
   CALL dInt3131(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(6010301)
   CALL dInt6131(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(3030301)
   CALL dInt3331(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(6030301)
   CALL dInt6331(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(6060301)
   CALL dInt6631(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(6010601)
   CALL dInt6161(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(6030601)
   CALL dInt6361(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(6060601)
   CALL dInt6661(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(6010303)
   CALL dInt6133(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(3030303)
   CALL dInt3333(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(6030303)
   CALL dInt6333(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(6060303)
   CALL dInt6633(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(6030603)
   CALL dInt6363(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(6060603)
   CALL dInt6663(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(6060606)
   CALL dInt6666(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE DEFAULT
     STOP 'MISSED AN INTEGRAL' 
   END SELECT 
