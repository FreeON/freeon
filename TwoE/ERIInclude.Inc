   SELECT CASE(IntType)
   CASE(1010101)
   CALL Int1111(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3010101)
   CALL Int3111(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6010101)
   CALL Int6111(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3030101)
   CALL Int3311(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6030101)
   CALL Int6311(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060101)
   CALL Int6611(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3010301)
   CALL Int3131(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6010301)
   CALL Int6131(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3030301)
   CALL Int3331(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6030301)
   CALL Int6331(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060301)
   CALL Int6631(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6010601)
   CALL Int6161(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6030601)
   CALL Int6361(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060601)
   CALL Int6661(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6010303)
   CALL Int6133(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3030303)
   CALL Int3333(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6030303)
   CALL Int6333(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060303)
   CALL Int6633(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6030603)
   CALL Int6363(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060603)
   CALL Int6663(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060606)
   CALL Int6666(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE DEFAULT
     STOP 'MISSED AN INTEGRAL' 
   END SELECT 
