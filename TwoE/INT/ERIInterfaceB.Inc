   SELECT CASE(IntType)
   CASE(1010101)
     CALL IntB1010101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1
   CASE(2010101)
     CALL IntB2010101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=4
   CASE(3010101)
     CALL IntB3010101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=3
   CASE(6010101)
     CALL IntB6010101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=6
   CASE(1020101)
     CALL IntB2010101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=4
   CASE(2020101)
     CALL IntB2020101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=16
   CASE(3020101)
     CALL IntB3020101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=12
   CASE(6020101)
     CALL IntB6020101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=24
   CASE(1030101)
     CALL IntB3010101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=3
   CASE(2030101)
     CALL IntB3020101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=12
   CASE(3030101)
     CALL IntB3030101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=9
   CASE(6030101)
     CALL IntB6030101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=18
   CASE(1060101)
     CALL IntB6010101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=6
   CASE(2060101)
     CALL IntB6020101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=24
   CASE(3060101)
     CALL IntB6030101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=18
   CASE(6060101)
     CALL IntB6060101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=36
   CASE(1010201)
     CALL IntB2010101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=4
   CASE(2010201)
     CALL IntB2010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=16
   CASE(3010201)
     CALL IntB3010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=12
   CASE(6010201)
     CALL IntB6010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=24
   CASE(1020201)
     CALL IntB2010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=16
   CASE(2020201)
     CALL IntB2020201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=64
   CASE(3020201)
     CALL IntB3020201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=48
   CASE(6020201)
     CALL IntB6020201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=96
   CASE(1030201)
     CALL IntB3010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=12
   CASE(2030201)
     CALL IntB3020201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=48
   CASE(3030201)
     CALL IntB3030201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=36
   CASE(6030201)
     CALL IntB6030201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=72
   CASE(1060201)
     CALL IntB6010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=24
   CASE(2060201)
     CALL IntB6020201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=96
   CASE(3060201)
     CALL IntB6030201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=72
   CASE(6060201)
     CALL IntB6060201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=144
   CASE(1010301)
     CALL IntB3010101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=3
   CASE(2010301)
     CALL IntB3010201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=12
   CASE(3010301)
     CALL IntB3010301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=9
   CASE(6010301)
     CALL IntB6010301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=18
   CASE(1020301)
     CALL IntB3010201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=12
   CASE(2020301)
     CALL IntB3010202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=48
   CASE(3020301)
     CALL IntB3020301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=36
   CASE(6020301)
     CALL IntB6020301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=72
   CASE(1030301)
     CALL IntB3010301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=9
   CASE(2030301)
     CALL IntB3020301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=36
   CASE(3030301)
     CALL IntB3030301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=27
   CASE(6030301)
     CALL IntB6030301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=54
   CASE(1060301)
     CALL IntB6010301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=18
   CASE(2060301)
     CALL IntB6020301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=72
   CASE(3060301)
     CALL IntB6030301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=54
   CASE(6060301)
     CALL IntB6060301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=108
   CASE(1010601)
     CALL IntB6010101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=6
   CASE(2010601)
     CALL IntB6010201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=24
   CASE(3010601)
     CALL IntB6010301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=18
   CASE(6010601)
     CALL IntB6010601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=36
   CASE(1020601)
     CALL IntB6010201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=24
   CASE(2020601)
     CALL IntB6010202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=96
   CASE(3020601)
     CALL IntB6010302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=72
   CASE(6020601)
     CALL IntB6020601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=144
   CASE(1030601)
     CALL IntB6010301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=18
   CASE(2030601)
     CALL IntB6010302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=72
   CASE(3030601)
     CALL IntB6010303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=54
   CASE(6030601)
     CALL IntB6030601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=108
   CASE(1060601)
     CALL IntB6010601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=36
   CASE(2060601)
     CALL IntB6020601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=144
   CASE(3060601)
     CALL IntB6030601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=108
   CASE(6060601)
     CALL IntB6060601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=216
   CASE(1010102)
     CALL IntB2010101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=4
   CASE(2010102)
     CALL IntB2010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=16
   CASE(3010102)
     CALL IntB3010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=12
   CASE(6010102)
     CALL IntB6010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=24
   CASE(1020102)
     CALL IntB2010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=16
   CASE(2020102)
     CALL IntB2020201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=64
   CASE(3020102)
     CALL IntB3020201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=48
   CASE(6020102)
     CALL IntB6020201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=96
   CASE(1030102)
     CALL IntB3010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=12
   CASE(2030102)
     CALL IntB3020201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=48
   CASE(3030102)
     CALL IntB3030201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=36
   CASE(6030102)
     CALL IntB6030201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=72
   CASE(1060102)
     CALL IntB6010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=24
   CASE(2060102)
     CALL IntB6020201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=96
   CASE(3060102)
     CALL IntB6030201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=72
   CASE(6060102)
     CALL IntB6060201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=144
   CASE(1010202)
     CALL IntB2020101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=16
   CASE(2010202)
     CALL IntB2020201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=64
   CASE(3010202)
     CALL IntB3010202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=48
   CASE(6010202)
     CALL IntB6010202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=96
   CASE(1020202)
     CALL IntB2020201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=64
   CASE(2020202)
     CALL IntB2020202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=256
   CASE(3020202)
     CALL IntB3020202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=192
   CASE(6020202)
     CALL IntB6020202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=384
   CASE(1030202)
     CALL IntB3010202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=48
   CASE(2030202)
     CALL IntB3020202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=192
   CASE(3030202)
     CALL IntB3030202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=144
   CASE(6030202)
     CALL IntB6030202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=288
   CASE(1060202)
     CALL IntB6010202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=96
   CASE(2060202)
     CALL IntB6020202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=384
   CASE(3060202)
     CALL IntB6030202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=288
   CASE(6060202)
     CALL IntB6060202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=576
   CASE(1010302)
     CALL IntB3020101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=12
   CASE(2010302)
     CALL IntB3020201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=48
   CASE(3010302)
     CALL IntB3020301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=36
   CASE(6010302)
     CALL IntB6010302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=72
   CASE(1020302)
     CALL IntB3020201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=48
   CASE(2020302)
     CALL IntB3020202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=192
   CASE(3020302)
     CALL IntB3020302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=144
   CASE(6020302)
     CALL IntB6020302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=288
   CASE(1030302)
     CALL IntB3020301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=36
   CASE(2030302)
     CALL IntB3020302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=144
   CASE(3030302)
     CALL IntB3030302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=108
   CASE(6030302)
     CALL IntB6030302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=216
   CASE(1060302)
     CALL IntB6010302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=72
   CASE(2060302)
     CALL IntB6020302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=288
   CASE(3060302)
     CALL IntB6030302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=216
   CASE(6060302)
     CALL IntB6060302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=432
   CASE(1010602)
     CALL IntB6020101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=24
   CASE(2010602)
     CALL IntB6020201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=96
   CASE(3010602)
     CALL IntB6020301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=72
   CASE(6010602)
     CALL IntB6020601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=144
   CASE(1020602)
     CALL IntB6020201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=96
   CASE(2020602)
     CALL IntB6020202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=384
   CASE(3020602)
     CALL IntB6020302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=288
   CASE(6020602)
     CALL IntB6020602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=576
   CASE(1030602)
     CALL IntB6020301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=72
   CASE(2030602)
     CALL IntB6020302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=288
   CASE(3030602)
     CALL IntB6020303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=216
   CASE(6030602)
     CALL IntB6030602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=432
   CASE(1060602)
     CALL IntB6020601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=144
   CASE(2060602)
     CALL IntB6020602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=576
   CASE(3060602)
     CALL IntB6030602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=432
   CASE(6060602)
     CALL IntB6060602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=864
   CASE(1010103)
     CALL IntB3010101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=3
   CASE(2010103)
     CALL IntB3010201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=12
   CASE(3010103)
     CALL IntB3010301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=9
   CASE(6010103)
     CALL IntB6010301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=18
   CASE(1020103)
     CALL IntB3010201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=12
   CASE(2020103)
     CALL IntB3010202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=48
   CASE(3020103)
     CALL IntB3020301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=36
   CASE(6020103)
     CALL IntB6020301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=72
   CASE(1030103)
     CALL IntB3010301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=9
   CASE(2030103)
     CALL IntB3020301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=36
   CASE(3030103)
     CALL IntB3030301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=27
   CASE(6030103)
     CALL IntB6030301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=54
   CASE(1060103)
     CALL IntB6010301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=18
   CASE(2060103)
     CALL IntB6020301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=72
   CASE(3060103)
     CALL IntB6030301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=54
   CASE(6060103)
     CALL IntB6060301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=108
   CASE(1010203)
     CALL IntB3020101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=12
   CASE(2010203)
     CALL IntB3020201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=48
   CASE(3010203)
     CALL IntB3020301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=36
   CASE(6010203)
     CALL IntB6010302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=72
   CASE(1020203)
     CALL IntB3020201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=48
   CASE(2020203)
     CALL IntB3020202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=192
   CASE(3020203)
     CALL IntB3020302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=144
   CASE(6020203)
     CALL IntB6020302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=288
   CASE(1030203)
     CALL IntB3020301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=36
   CASE(2030203)
     CALL IntB3020302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=144
   CASE(3030203)
     CALL IntB3030302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=108
   CASE(6030203)
     CALL IntB6030302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=216
   CASE(1060203)
     CALL IntB6010302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=72
   CASE(2060203)
     CALL IntB6020302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=288
   CASE(3060203)
     CALL IntB6030302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=216
   CASE(6060203)
     CALL IntB6060302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=432
   CASE(1010303)
     CALL IntB3030101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=9
   CASE(2010303)
     CALL IntB3030201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=36
   CASE(3010303)
     CALL IntB3030301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=27
   CASE(6010303)
     CALL IntB6010303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=54
   CASE(1020303)
     CALL IntB3030201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=36
   CASE(2020303)
     CALL IntB3030202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=144
   CASE(3020303)
     CALL IntB3030302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=108
   CASE(6020303)
     CALL IntB6020303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=216
   CASE(1030303)
     CALL IntB3030301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=27
   CASE(2030303)
     CALL IntB3030302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=108
   CASE(3030303)
     CALL IntB3030303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=81
   CASE(6030303)
     CALL IntB6030303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=162
   CASE(1060303)
     CALL IntB6010303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=54
   CASE(2060303)
     CALL IntB6020303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=216
   CASE(3060303)
     CALL IntB6030303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=162
   CASE(6060303)
     CALL IntB6060303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=324
   CASE(1010603)
     CALL IntB6030101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=18
   CASE(2010603)
     CALL IntB6030201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=72
   CASE(3010603)
     CALL IntB6030301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=54
   CASE(6010603)
     CALL IntB6030601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=108
   CASE(1020603)
     CALL IntB6030201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=72
   CASE(2020603)
     CALL IntB6030202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=288
   CASE(3020603)
     CALL IntB6030302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=216
   CASE(6020603)
     CALL IntB6030602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=432
   CASE(1030603)
     CALL IntB6030301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=54
   CASE(2030603)
     CALL IntB6030302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=216
   CASE(3030603)
     CALL IntB6030303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=162
   CASE(6030603)
     CALL IntB6030603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=324
   CASE(1060603)
     CALL IntB6030601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=108
   CASE(2060603)
     CALL IntB6030602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=432
   CASE(3060603)
     CALL IntB6030603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=324
   CASE(6060603)
     CALL IntB6060603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=648
   CASE(1010106)
     CALL IntB6010101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=6
   CASE(2010106)
     CALL IntB6010201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=24
   CASE(3010106)
     CALL IntB6010301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=18
   CASE(6010106)
     CALL IntB6010601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=36
   CASE(1020106)
     CALL IntB6010201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=24
   CASE(2020106)
     CALL IntB6010202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=96
   CASE(3020106)
     CALL IntB6010302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=72
   CASE(6020106)
     CALL IntB6020601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=144
   CASE(1030106)
     CALL IntB6010301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=18
   CASE(2030106)
     CALL IntB6010302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=72
   CASE(3030106)
     CALL IntB6010303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=54
   CASE(6030106)
     CALL IntB6030601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=108
   CASE(1060106)
     CALL IntB6010601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=36
   CASE(2060106)
     CALL IntB6020601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=144
   CASE(3060106)
     CALL IntB6030601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=108
   CASE(6060106)
     CALL IntB6060601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=216
   CASE(1010206)
     CALL IntB6020101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=24
   CASE(2010206)
     CALL IntB6020201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=96
   CASE(3010206)
     CALL IntB6020301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=72
   CASE(6010206)
     CALL IntB6020601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=144
   CASE(1020206)
     CALL IntB6020201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=96
   CASE(2020206)
     CALL IntB6020202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=384
   CASE(3020206)
     CALL IntB6020302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=288
   CASE(6020206)
     CALL IntB6020602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=576
   CASE(1030206)
     CALL IntB6020301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=72
   CASE(2030206)
     CALL IntB6020302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=288
   CASE(3030206)
     CALL IntB6020303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=216
   CASE(6030206)
     CALL IntB6030602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=432
   CASE(1060206)
     CALL IntB6020601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=144
   CASE(2060206)
     CALL IntB6020602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=576
   CASE(3060206)
     CALL IntB6030602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=432
   CASE(6060206)
     CALL IntB6060602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=864
   CASE(1010306)
     CALL IntB6030101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=18
   CASE(2010306)
     CALL IntB6030201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=72
   CASE(3010306)
     CALL IntB6030301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=54
   CASE(6010306)
     CALL IntB6030601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=108
   CASE(1020306)
     CALL IntB6030201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=72
   CASE(2020306)
     CALL IntB6030202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=288
   CASE(3020306)
     CALL IntB6030302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=216
   CASE(6020306)
     CALL IntB6030602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=432
   CASE(1030306)
     CALL IntB6030301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=54
   CASE(2030306)
     CALL IntB6030302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=216
   CASE(3030306)
     CALL IntB6030303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=162
   CASE(6030306)
     CALL IntB6030603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=324
   CASE(1060306)
     CALL IntB6030601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=108
   CASE(2060306)
     CALL IntB6030602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=432
   CASE(3060306)
     CALL IntB6030603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=324
   CASE(6060306)
     CALL IntB6060603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=648
   CASE(1010606)
     CALL IntB6060101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=36
   CASE(2010606)
     CALL IntB6060201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=144
   CASE(3010606)
     CALL IntB6060301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=108
   CASE(6010606)
     CALL IntB6060601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=216
   CASE(1020606)
     CALL IntB6060201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=144
   CASE(2020606)
     CALL IntB6060202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=576
   CASE(3020606)
     CALL IntB6060302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=432
   CASE(6020606)
     CALL IntB6060602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=864
   CASE(1030606)
     CALL IntB6060301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=108
   CASE(2030606)
     CALL IntB6060302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=432
   CASE(3030606)
     CALL IntB6060303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=324
   CASE(6030606)
     CALL IntB6060603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=648
   CASE(1060606)
     CALL IntB6060601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=216
   CASE(2060606)
     CALL IntB6060602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=864
   CASE(3060606)
     CALL IntB6060603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=648
   CASE(6060606)
     CALL IntB6060606(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1296
   CASE DEFAULT
      WRITE(*,*) 'We are in ERIInterface.Inc'
      WRITE(*,*) 'IntType=',IntType
      STOP 'MISS AN INTEGRAL'
   END SELECT
