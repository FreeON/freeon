   SELECT CASE(IntType)
   CASE(1010101)
   CALL IntB1010101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(2010101)
   CALL IntB2010101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3010101)
   CALL IntB3010101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6010101)
   CALL IntB6010101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10010101)
   CALL IntB10010101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(2020101)
   CALL IntB2020101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3020101)
   CALL IntB3020101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6020101)
   CALL IntB6020101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10020101)
   CALL IntB10020101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3030101)
   CALL IntB3030101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6030101)
   CALL IntB6030101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10030101)
   CALL IntB10030101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060101)
   CALL IntB6060101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10060101)
   CALL IntB10060101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10100101)
   CALL IntB10100101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(2010201)
   CALL IntB2010201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3010201)
   CALL IntB3010201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6010201)
   CALL IntB6010201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10010201)
   CALL IntB10010201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(2020201)
   CALL IntB2020201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3020201)
   CALL IntB3020201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6020201)
   CALL IntB6020201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10020201)
   CALL IntB10020201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3030201)
   CALL IntB3030201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6030201)
   CALL IntB6030201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10030201)
   CALL IntB10030201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060201)
   CALL IntB6060201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10060201)
   CALL IntB10060201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10100201)
   CALL IntB10100201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3010301)
   CALL IntB3010301(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6010301)
   CALL IntB6010301(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10010301)
   CALL IntB10010301(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3020301)
   CALL IntB3020301(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6020301)
   CALL IntB6020301(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10020301)
   CALL IntB10020301(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3030301)
   CALL IntB3030301(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6030301)
   CALL IntB6030301(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10030301)
   CALL IntB10030301(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060301)
   CALL IntB6060301(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10060301)
   CALL IntB10060301(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10100301)
   CALL IntB10100301(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6010601)
   CALL IntB6010601(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10010601)
   CALL IntB10010601(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6020601)
   CALL IntB6020601(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10020601)
   CALL IntB10020601(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6030601)
   CALL IntB6030601(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10030601)
   CALL IntB10030601(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060601)
   CALL IntB6060601(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10060601)
   CALL IntB10060601(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10100601)
   CALL IntB10100601(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10011001)
   CALL IntB10011001(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10021001)
   CALL IntB10021001(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10031001)
   CALL IntB10031001(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10061001)
   CALL IntB10061001(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10101001)
   CALL IntB10101001(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3010202)
   CALL IntB3010202(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6010202)
   CALL IntB6010202(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10010202)
   CALL IntB10010202(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(2020202)
   CALL IntB2020202(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3020202)
   CALL IntB3020202(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6020202)
   CALL IntB6020202(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10020202)
   CALL IntB10020202(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3030202)
   CALL IntB3030202(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6030202)
   CALL IntB6030202(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10030202)
   CALL IntB10030202(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060202)
   CALL IntB6060202(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10060202)
   CALL IntB10060202(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10100202)
   CALL IntB10100202(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6010302)
   CALL IntB6010302(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10010302)
   CALL IntB10010302(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3020302)
   CALL IntB3020302(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6020302)
   CALL IntB6020302(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10020302)
   CALL IntB10020302(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3030302)
   CALL IntB3030302(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6030302)
   CALL IntB6030302(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10030302)
   CALL IntB10030302(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060302)
   CALL IntB6060302(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10060302)
   CALL IntB10060302(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10100302)
   CALL IntB10100302(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10010602)
   CALL IntB10010602(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6020602)
   CALL IntB6020602(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10020602)
   CALL IntB10020602(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6030602)
   CALL IntB6030602(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10030602)
   CALL IntB10030602(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060602)
   CALL IntB6060602(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10060602)
   CALL IntB10060602(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10100602)
   CALL IntB10100602(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10021002)
   CALL IntB10021002(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10031002)
   CALL IntB10031002(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10061002)
   CALL IntB10061002(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10101002)
   CALL IntB10101002(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6010303)
   CALL IntB6010303(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10010303)
   CALL IntB10010303(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6020303)
   CALL IntB6020303(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10020303)
   CALL IntB10020303(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(3030303)
   CALL IntB3030303(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6030303)
   CALL IntB6030303(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10030303)
   CALL IntB10030303(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060303)
   CALL IntB6060303(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10060303)
   CALL IntB10060303(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10100303)
   CALL IntB10100303(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10010603)
   CALL IntB10010603(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10020603)
   CALL IntB10020603(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6030603)
   CALL IntB6030603(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10030603)
   CALL IntB10030603(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060603)
   CALL IntB6060603(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10060603)
   CALL IntB10060603(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10100603)
   CALL IntB10100603(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10031003)
   CALL IntB10031003(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10061003)
   CALL IntB10061003(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10101003)
   CALL IntB10101003(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10010606)
   CALL IntB10010606(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10020606)
   CALL IntB10020606(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10030606)
   CALL IntB10030606(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(6060606)
   CALL IntB6060606(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10060606)
   CALL IntB10060606(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10100606)
   CALL IntB10100606(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10061006)
   CALL IntB10061006(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10101006)
   CALL IntB10101006(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE(10101010)
   CALL IntB10101010(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                        OffSet%A  ,1              , & 
                        OffSet%C-1,NBFA           , & 
                        OffSet%B-1,NBFA*NBFC      , & 
                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) 

   CASE DEFAULT
    ! STOP 'MISSED AN INTEGRAL' 
   END SELECT 
