   SELECT CASE(IntType)
   CASE(1010101)
     CALL dInt1111(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  1,4,7,10,NIntBlk,GM%PBC,C(1))
     LocNInt=12
   CASE(3010101)
     CALL dInt3111(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  1,4,7,10,NIntBlk,GM%PBC,C(1))
     LocNInt=36
   CASE(1030101)
     CALL dInt3111(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  4,1,7,10,NIntBlk,GM%PBC,C(1))
     LocNInt=36
   CASE(3030101)
     CALL dInt3311(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  1,4,7,10,NIntBlk,GM%PBC,C(1))
     LocNInt=108
   CASE(1010301)
     CALL dInt3111(BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  7,10,1,4,NIntBlk,GM%PBC,C(1))
     LocNInt=36
   CASE(3010301)
     CALL dInt3131(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  1,4,7,10,NIntBlk,GM%PBC,C(1))
     LocNInt=108
   CASE(1030301)
     CALL dInt3131(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  4,1,7,10,NIntBlk,GM%PBC,C(1))
     LocNInt=108
   CASE(3030301)
     CALL dInt3331(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  1,4,7,10,NIntBlk,GM%PBC,C(1))
     LocNInt=324
   CASE(1010103)
     CALL dInt3111(BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  10,7,1,4,NIntBlk,GM%PBC,C(1))
     LocNInt=36
   CASE(3010103)
     CALL dInt3131(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  1,4,10,7,NIntBlk,GM%PBC,C(1))
     LocNInt=108
   CASE(1030103)
     CALL dInt3131(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  4,1,10,7,NIntBlk,GM%PBC,C(1))
     LocNInt=108
   CASE(3030103)
     CALL dInt3331(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  1,4,10,7,NIntBlk,GM%PBC,C(1))
     LocNInt=324
   CASE(1010303)
     CALL dInt3311(BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  7,10,1,4,NIntBlk,GM%PBC,C(1))
     LocNInt=108
   CASE(3010303)
     CALL dInt3331(BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  7,10,1,4,NIntBlk,GM%PBC,C(1))
     LocNInt=324
   CASE(1030303)
     CALL dInt3331(BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  7,10,4,1,NIntBlk,GM%PBC,C(1))
     LocNInt=324
   CASE(3030303)
     CALL dInt3333(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  1,4,7,10,NIntBlk,GM%PBC,C(1))
     LocNInt=972
   CASE DEFAULT
      STOP 'MISS AN INTEGRAL'
   END SELECT
