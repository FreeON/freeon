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

   CASE(1030101)
   CALL dInt1311(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
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

   CASE(1010301)
   CALL dInt1131(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
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

   CASE(1030301)
   CALL dInt1331(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
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

   CASE(1010103)
   CALL dInt1113(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(3010103)
   CALL dInt3113(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(1030103)
   CALL dInt1313(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(3030103)
   CALL dInt3313(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(1010303)
   CALL dInt1133(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(3010303)
   CALL dInt3133(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,1,4,7,10,NIntBlk,GM%PBC,dC(1,1)) 

   CASE(1030303)
   CALL dInt1333(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
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

   CASE DEFAULT
     STOP 'MISSED AN INTEGRAL' 
   END SELECT 
