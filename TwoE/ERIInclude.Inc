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

   CASE(1030101)
   CALL Int1311(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
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

   CASE(1010301)
   CALL Int1131(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
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

   CASE(1030301)
   CALL Int1331(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
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

   CASE(1010103)
   CALL Int1113(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3010103)
   CALL Int3113(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(1030103)
   CALL Int1313(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3030103)
   CALL Int3313(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(1010303)
   CALL Int1133(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3010303)
   CALL Int3133(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                       OffSet%A  ,1              , & 
                       OffSet%C-1,NBFA           , & 
                       OffSet%B-1,NBFA*NBFC      , & 
                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(1030303)
   CALL Int1333(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
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

   CASE DEFAULT
     STOP 'MISSED AN INTEGRAL' 
   END SELECT 
