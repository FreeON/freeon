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
   CASE(10010101)
     CALL IntB10010101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=10
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
   CASE(10020101)
     CALL IntB10020101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=40
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
   CASE(10030101)
     CALL IntB10030101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=30
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
   CASE(10060101)
     CALL IntB10060101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=60
   CASE(1100101)
     CALL IntB10010101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=10
   CASE(2100101)
     CALL IntB10020101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=40
   CASE(3100101)
     CALL IntB10030101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=30
   CASE(6100101)
     CALL IntB10060101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=60
   CASE(10100101)
     CALL IntB10100101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=100
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
   CASE(10010201)
     CALL IntB10010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=40
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
   CASE(10020201)
     CALL IntB10020201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=160
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
   CASE(10030201)
     CALL IntB10030201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=120
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
   CASE(10060201)
     CALL IntB10060201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=240
   CASE(1100201)
     CALL IntB10010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=40
   CASE(2100201)
     CALL IntB10020201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=160
   CASE(3100201)
     CALL IntB10030201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=120
   CASE(6100201)
     CALL IntB10060201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=240
   CASE(10100201)
     CALL IntB10100201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=400
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
   CASE(10010301)
     CALL IntB10010301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=30
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
   CASE(10020301)
     CALL IntB10020301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=120
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
   CASE(10030301)
     CALL IntB10030301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=90
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
   CASE(10060301)
     CALL IntB10060301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=180
   CASE(1100301)
     CALL IntB10010301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=30
   CASE(2100301)
     CALL IntB10020301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=120
   CASE(3100301)
     CALL IntB10030301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=90
   CASE(6100301)
     CALL IntB10060301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=180
   CASE(10100301)
     CALL IntB10100301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=300
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
   CASE(10010601)
     CALL IntB10010601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=60
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
   CASE(10020601)
     CALL IntB10020601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=240
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
   CASE(10030601)
     CALL IntB10030601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=180
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
   CASE(10060601)
     CALL IntB10060601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=360
   CASE(1100601)
     CALL IntB10010601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=60
   CASE(2100601)
     CALL IntB10020601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=240
   CASE(3100601)
     CALL IntB10030601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=180
   CASE(6100601)
     CALL IntB10060601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=360
   CASE(10100601)
     CALL IntB10100601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=600
   CASE(1011001)
     CALL IntB10010101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=10
   CASE(2011001)
     CALL IntB10010201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=40
   CASE(3011001)
     CALL IntB10010301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=30
   CASE(6011001)
     CALL IntB10010601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=60
   CASE(10011001)
     CALL IntB10011001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=100
   CASE(1021001)
     CALL IntB10010201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=40
   CASE(2021001)
     CALL IntB10010202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=160
   CASE(3021001)
     CALL IntB10010302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=120
   CASE(6021001)
     CALL IntB10010602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=240
   CASE(10021001)
     CALL IntB10021001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=400
   CASE(1031001)
     CALL IntB10010301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=30
   CASE(2031001)
     CALL IntB10010302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=120
   CASE(3031001)
     CALL IntB10010303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=90
   CASE(6031001)
     CALL IntB10010603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=180
   CASE(10031001)
     CALL IntB10031001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=300
   CASE(1061001)
     CALL IntB10010601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=60
   CASE(2061001)
     CALL IntB10010602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=240
   CASE(3061001)
     CALL IntB10010603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=180
   CASE(6061001)
     CALL IntB10010606(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=360
   CASE(10061001)
     CALL IntB10061001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=600
   CASE(1101001)
     CALL IntB10011001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=100
   CASE(2101001)
     CALL IntB10021001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=400
   CASE(3101001)
     CALL IntB10031001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=300
   CASE(6101001)
     CALL IntB10061001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=600
   CASE(10101001)
     CALL IntB10101001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1000
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
   CASE(10010102)
     CALL IntB10010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=40
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
   CASE(10020102)
     CALL IntB10020201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=160
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
   CASE(10030102)
     CALL IntB10030201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=120
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
   CASE(10060102)
     CALL IntB10060201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=240
   CASE(1100102)
     CALL IntB10010201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=40
   CASE(2100102)
     CALL IntB10020201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=160
   CASE(3100102)
     CALL IntB10030201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=120
   CASE(6100102)
     CALL IntB10060201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=240
   CASE(10100102)
     CALL IntB10100201(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=400
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
   CASE(10010202)
     CALL IntB10010202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=160
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
   CASE(10020202)
     CALL IntB10020202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=640
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
   CASE(10030202)
     CALL IntB10030202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=480
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
   CASE(10060202)
     CALL IntB10060202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=960
   CASE(1100202)
     CALL IntB10010202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=160
   CASE(2100202)
     CALL IntB10020202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=640
   CASE(3100202)
     CALL IntB10030202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=480
   CASE(6100202)
     CALL IntB10060202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=960
   CASE(10100202)
     CALL IntB10100202(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1600
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
   CASE(10010302)
     CALL IntB10010302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=120
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
   CASE(10020302)
     CALL IntB10020302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=480
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
   CASE(10030302)
     CALL IntB10030302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=360
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
   CASE(10060302)
     CALL IntB10060302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=720
   CASE(1100302)
     CALL IntB10010302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=120
   CASE(2100302)
     CALL IntB10020302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=480
   CASE(3100302)
     CALL IntB10030302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=360
   CASE(6100302)
     CALL IntB10060302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=720
   CASE(10100302)
     CALL IntB10100302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1200
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
   CASE(10010602)
     CALL IntB10010602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=240
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
   CASE(10020602)
     CALL IntB10020602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=960
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
   CASE(10030602)
     CALL IntB10030602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=720
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
   CASE(10060602)
     CALL IntB10060602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1440
   CASE(1100602)
     CALL IntB10010602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=240
   CASE(2100602)
     CALL IntB10020602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=960
   CASE(3100602)
     CALL IntB10030602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=720
   CASE(6100602)
     CALL IntB10060602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1440
   CASE(10100602)
     CALL IntB10100602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=2400
   CASE(1011002)
     CALL IntB10020101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=40
   CASE(2011002)
     CALL IntB10020201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=160
   CASE(3011002)
     CALL IntB10020301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=120
   CASE(6011002)
     CALL IntB10020601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=240
   CASE(10011002)
     CALL IntB10021001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=400
   CASE(1021002)
     CALL IntB10020201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=160
   CASE(2021002)
     CALL IntB10020202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=640
   CASE(3021002)
     CALL IntB10020302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=480
   CASE(6021002)
     CALL IntB10020602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=960
   CASE(10021002)
     CALL IntB10021002(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1600
   CASE(1031002)
     CALL IntB10020301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=120
   CASE(2031002)
     CALL IntB10020302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=480
   CASE(3031002)
     CALL IntB10020303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=360
   CASE(6031002)
     CALL IntB10020603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=720
   CASE(10031002)
     CALL IntB10031002(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1200
   CASE(1061002)
     CALL IntB10020601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=240
   CASE(2061002)
     CALL IntB10020602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=960
   CASE(3061002)
     CALL IntB10020603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=720
   CASE(6061002)
     CALL IntB10020606(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1440
   CASE(10061002)
     CALL IntB10061002(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=2400
   CASE(1101002)
     CALL IntB10021001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=400
   CASE(2101002)
     CALL IntB10021002(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1600
   CASE(3101002)
     CALL IntB10031002(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1200
   CASE(6101002)
     CALL IntB10061002(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=2400
   CASE(10101002)
     CALL IntB10101002(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=4000
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
   CASE(10010103)
     CALL IntB10010301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=30
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
   CASE(10020103)
     CALL IntB10020301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=120
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
   CASE(10030103)
     CALL IntB10030301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=90
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
   CASE(10060103)
     CALL IntB10060301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=180
   CASE(1100103)
     CALL IntB10010301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=30
   CASE(2100103)
     CALL IntB10020301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=120
   CASE(3100103)
     CALL IntB10030301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=90
   CASE(6100103)
     CALL IntB10060301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=180
   CASE(10100103)
     CALL IntB10100301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=300
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
   CASE(10010203)
     CALL IntB10010302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=120
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
   CASE(10020203)
     CALL IntB10020302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=480
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
   CASE(10030203)
     CALL IntB10030302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=360
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
   CASE(10060203)
     CALL IntB10060302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=720
   CASE(1100203)
     CALL IntB10010302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=120
   CASE(2100203)
     CALL IntB10020302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=480
   CASE(3100203)
     CALL IntB10030302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=360
   CASE(6100203)
     CALL IntB10060302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=720
   CASE(10100203)
     CALL IntB10100302(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=1200
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
   CASE(10010303)
     CALL IntB10010303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=90
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
   CASE(10020303)
     CALL IntB10020303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=360
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
   CASE(10030303)
     CALL IntB10030303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=270
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
   CASE(10060303)
     CALL IntB10060303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=540
   CASE(1100303)
     CALL IntB10010303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=90
   CASE(2100303)
     CALL IntB10020303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=360
   CASE(3100303)
     CALL IntB10030303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=270
   CASE(6100303)
     CALL IntB10060303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=540
   CASE(10100303)
     CALL IntB10100303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=900
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
   CASE(10010603)
     CALL IntB10010603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=180
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
   CASE(10020603)
     CALL IntB10020603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=720
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
   CASE(10030603)
     CALL IntB10030603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=540
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
   CASE(10060603)
     CALL IntB10060603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1080
   CASE(1100603)
     CALL IntB10010603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=180
   CASE(2100603)
     CALL IntB10020603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=720
   CASE(3100603)
     CALL IntB10030603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=540
   CASE(6100603)
     CALL IntB10060603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1080
   CASE(10100603)
     CALL IntB10100603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1800
   CASE(1011003)
     CALL IntB10030101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=30
   CASE(2011003)
     CALL IntB10030201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=120
   CASE(3011003)
     CALL IntB10030301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=90
   CASE(6011003)
     CALL IntB10030601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=180
   CASE(10011003)
     CALL IntB10031001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=300
   CASE(1021003)
     CALL IntB10030201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=120
   CASE(2021003)
     CALL IntB10030202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=480
   CASE(3021003)
     CALL IntB10030302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=360
   CASE(6021003)
     CALL IntB10030602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=720
   CASE(10021003)
     CALL IntB10031002(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1200
   CASE(1031003)
     CALL IntB10030301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=90
   CASE(2031003)
     CALL IntB10030302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=360
   CASE(3031003)
     CALL IntB10030303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=270
   CASE(6031003)
     CALL IntB10030603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=540
   CASE(10031003)
     CALL IntB10031003(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=900
   CASE(1061003)
     CALL IntB10030601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=180
   CASE(2061003)
     CALL IntB10030602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=720
   CASE(3061003)
     CALL IntB10030603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=540
   CASE(6061003)
     CALL IntB10030606(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1080
   CASE(10061003)
     CALL IntB10061003(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1800
   CASE(1101003)
     CALL IntB10031001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=300
   CASE(2101003)
     CALL IntB10031002(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=1200
   CASE(3101003)
     CALL IntB10031003(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=900
   CASE(6101003)
     CALL IntB10061003(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1800
   CASE(10101003)
     CALL IntB10101003(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=3000
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
   CASE(10010106)
     CALL IntB10010601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=60
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
   CASE(10020106)
     CALL IntB10020601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=240
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
   CASE(10030106)
     CALL IntB10030601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=180
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
   CASE(10060106)
     CALL IntB10060601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=360
   CASE(1100106)
     CALL IntB10010601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=60
   CASE(2100106)
     CALL IntB10020601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=240
   CASE(3100106)
     CALL IntB10030601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=180
   CASE(6100106)
     CALL IntB10060601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=360
   CASE(10100106)
     CALL IntB10100601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=600
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
   CASE(10010206)
     CALL IntB10010602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=240
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
   CASE(10020206)
     CALL IntB10020602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=960
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
   CASE(10030206)
     CALL IntB10030602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=720
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
   CASE(10060206)
     CALL IntB10060602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=1440
   CASE(1100206)
     CALL IntB10010602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=240
   CASE(2100206)
     CALL IntB10020602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=960
   CASE(3100206)
     CALL IntB10030602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=720
   CASE(6100206)
     CALL IntB10060602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=1440
   CASE(10100206)
     CALL IntB10100602(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=2400
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
   CASE(10010306)
     CALL IntB10010603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=180
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
   CASE(10020306)
     CALL IntB10020603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=720
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
   CASE(10030306)
     CALL IntB10030603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=540
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
   CASE(10060306)
     CALL IntB10060603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=1080
   CASE(1100306)
     CALL IntB10010603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=180
   CASE(2100306)
     CALL IntB10020603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=720
   CASE(3100306)
     CALL IntB10030603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=540
   CASE(6100306)
     CALL IntB10060603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=1080
   CASE(10100306)
     CALL IntB10100603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=1800
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
   CASE(10010606)
     CALL IntB10010606(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=360
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
   CASE(10020606)
     CALL IntB10020606(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1440
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
   CASE(10030606)
     CALL IntB10030606(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1080
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
   CASE(10060606)
     CALL IntB10060606(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=2160
   CASE(1100606)
     CALL IntB10010606(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=360
   CASE(2100606)
     CALL IntB10020606(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1440
   CASE(3100606)
     CALL IntB10030606(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=1080
   CASE(6100606)
     CALL IntB10060606(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=2160
   CASE(10100606)
     CALL IntB10100606(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=3600
   CASE(1011006)
     CALL IntB10060101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=60
   CASE(2011006)
     CALL IntB10060201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=240
   CASE(3011006)
     CALL IntB10060301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=180
   CASE(6011006)
     CALL IntB10060601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=360
   CASE(10011006)
     CALL IntB10061001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=600
   CASE(1021006)
     CALL IntB10060201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=240
   CASE(2021006)
     CALL IntB10060202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=960
   CASE(3021006)
     CALL IntB10060302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=720
   CASE(6021006)
     CALL IntB10060602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1440
   CASE(10021006)
     CALL IntB10061002(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=2400
   CASE(1031006)
     CALL IntB10060301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=180
   CASE(2031006)
     CALL IntB10060302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=720
   CASE(3031006)
     CALL IntB10060303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=540
   CASE(6031006)
     CALL IntB10060603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1080
   CASE(10031006)
     CALL IntB10061003(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1800
   CASE(1061006)
     CALL IntB10060601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=360
   CASE(2061006)
     CALL IntB10060602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=1440
   CASE(3061006)
     CALL IntB10060603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=1080
   CASE(6061006)
     CALL IntB10060606(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=2160
   CASE(10061006)
     CALL IntB10061006(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=3600
   CASE(1101006)
     CALL IntB10061001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=600
   CASE(2101006)
     CALL IntB10061002(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=2400
   CASE(3101006)
     CALL IntB10061003(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=1800
   CASE(6101006)
     CALL IntB10061006(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=3600
   CASE(10101006)
     CALL IntB10101006(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=6000
   CASE(1010110)
     CALL IntB10010101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=10
   CASE(2010110)
     CALL IntB10010201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=40
   CASE(3010110)
     CALL IntB10010301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=30
   CASE(6010110)
     CALL IntB10010601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=60
   CASE(10010110)
     CALL IntB10011001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=100
   CASE(1020110)
     CALL IntB10010201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=40
   CASE(2020110)
     CALL IntB10010202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=160
   CASE(3020110)
     CALL IntB10010302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=120
   CASE(6020110)
     CALL IntB10010602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=240
   CASE(10020110)
     CALL IntB10021001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=400
   CASE(1030110)
     CALL IntB10010301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=30
   CASE(2030110)
     CALL IntB10010302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=120
   CASE(3030110)
     CALL IntB10010303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=90
   CASE(6030110)
     CALL IntB10010603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=180
   CASE(10030110)
     CALL IntB10031001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=300
   CASE(1060110)
     CALL IntB10010601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=60
   CASE(2060110)
     CALL IntB10010602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=240
   CASE(3060110)
     CALL IntB10010603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=180
   CASE(6060110)
     CALL IntB10010606(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=360
   CASE(10060110)
     CALL IntB10061001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=600
   CASE(1100110)
     CALL IntB10011001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=100
   CASE(2100110)
     CALL IntB10021001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=400
   CASE(3100110)
     CALL IntB10031001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=300
   CASE(6100110)
     CALL IntB10061001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=600
   CASE(10100110)
     CALL IntB10101001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=1000
   CASE(1010210)
     CALL IntB10020101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=40
   CASE(2010210)
     CALL IntB10020201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=160
   CASE(3010210)
     CALL IntB10020301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=120
   CASE(6010210)
     CALL IntB10020601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=240
   CASE(10010210)
     CALL IntB10021001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=400
   CASE(1020210)
     CALL IntB10020201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=160
   CASE(2020210)
     CALL IntB10020202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=640
   CASE(3020210)
     CALL IntB10020302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=480
   CASE(6020210)
     CALL IntB10020602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=960
   CASE(10020210)
     CALL IntB10021002(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=1600
   CASE(1030210)
     CALL IntB10020301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=120
   CASE(2030210)
     CALL IntB10020302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=480
   CASE(3030210)
     CALL IntB10020303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=360
   CASE(6030210)
     CALL IntB10020603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=720
   CASE(10030210)
     CALL IntB10031002(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=1200
   CASE(1060210)
     CALL IntB10020601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=240
   CASE(2060210)
     CALL IntB10020602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=960
   CASE(3060210)
     CALL IntB10020603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=720
   CASE(6060210)
     CALL IntB10020606(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1440
   CASE(10060210)
     CALL IntB10061002(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=2400
   CASE(1100210)
     CALL IntB10021001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=400
   CASE(2100210)
     CALL IntB10021002(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=1600
   CASE(3100210)
     CALL IntB10031002(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=1200
   CASE(6100210)
     CALL IntB10061002(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=2400
   CASE(10100210)
     CALL IntB10101002(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=4000
   CASE(1010310)
     CALL IntB10030101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=30
   CASE(2010310)
     CALL IntB10030201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=120
   CASE(3010310)
     CALL IntB10030301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=90
   CASE(6010310)
     CALL IntB10030601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=180
   CASE(10010310)
     CALL IntB10031001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=300
   CASE(1020310)
     CALL IntB10030201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=120
   CASE(2020310)
     CALL IntB10030202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=480
   CASE(3020310)
     CALL IntB10030302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=360
   CASE(6020310)
     CALL IntB10030602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=720
   CASE(10020310)
     CALL IntB10031002(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1200
   CASE(1030310)
     CALL IntB10030301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=90
   CASE(2030310)
     CALL IntB10030302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=360
   CASE(3030310)
     CALL IntB10030303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=270
   CASE(6030310)
     CALL IntB10030603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=540
   CASE(10030310)
     CALL IntB10031003(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=900
   CASE(1060310)
     CALL IntB10030601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=180
   CASE(2060310)
     CALL IntB10030602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=720
   CASE(3060310)
     CALL IntB10030603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=540
   CASE(6060310)
     CALL IntB10030606(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1080
   CASE(10060310)
     CALL IntB10061003(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=1800
   CASE(1100310)
     CALL IntB10031001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=300
   CASE(2100310)
     CALL IntB10031002(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=1200
   CASE(3100310)
     CALL IntB10031003(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=900
   CASE(6100310)
     CALL IntB10061003(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=1800
   CASE(10100310)
     CALL IntB10101003(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=3000
   CASE(1010610)
     CALL IntB10060101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=60
   CASE(2010610)
     CALL IntB10060201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=240
   CASE(3010610)
     CALL IntB10060301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=180
   CASE(6010610)
     CALL IntB10060601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=360
   CASE(10010610)
     CALL IntB10061001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=600
   CASE(1020610)
     CALL IntB10060201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=240
   CASE(2020610)
     CALL IntB10060202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=960
   CASE(3020610)
     CALL IntB10060302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=720
   CASE(6020610)
     CALL IntB10060602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1440
   CASE(10020610)
     CALL IntB10061002(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=2400
   CASE(1030610)
     CALL IntB10060301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=180
   CASE(2030610)
     CALL IntB10060302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=720
   CASE(3030610)
     CALL IntB10060303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=540
   CASE(6030610)
     CALL IntB10060603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1080
   CASE(10030610)
     CALL IntB10061003(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1800
   CASE(1060610)
     CALL IntB10060601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=360
   CASE(2060610)
     CALL IntB10060602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=1440
   CASE(3060610)
     CALL IntB10060603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=1080
   CASE(6060610)
     CALL IntB10060606(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=2160
   CASE(10060610)
     CALL IntB10061006(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=3600
   CASE(1100610)
     CALL IntB10061001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=600
   CASE(2100610)
     CALL IntB10061002(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=2400
   CASE(3100610)
     CALL IntB10061003(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=1800
   CASE(6100610)
     CALL IntB10061006(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=3600
   CASE(10100610)
     CALL IntB10101006(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=6000
   CASE(1011010)
     CALL IntB10100101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=100
   CASE(2011010)
     CALL IntB10100201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=400
   CASE(3011010)
     CALL IntB10100301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=300
   CASE(6011010)
     CALL IntB10100601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=600
   CASE(10011010)
     CALL IntB10101001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1000
   CASE(1021010)
     CALL IntB10100201(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=400
   CASE(2021010)
     CALL IntB10100202(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1600
   CASE(3021010)
     CALL IntB10100302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1200
   CASE(6021010)
     CALL IntB10100602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=2400
   CASE(10021010)
     CALL IntB10101002(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=4000
   CASE(1031010)
     CALL IntB10100301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=300
   CASE(2031010)
     CALL IntB10100302(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=1200
   CASE(3031010)
     CALL IntB10100303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=900
   CASE(6031010)
     CALL IntB10100603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=1800
   CASE(10031010)
     CALL IntB10101003(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=3000
   CASE(1061010)
     CALL IntB10100601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=600
   CASE(2061010)
     CALL IntB10100602(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=2400
   CASE(3061010)
     CALL IntB10100603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=1800
   CASE(6061010)
     CALL IntB10100606(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=3600
   CASE(10061010)
     CALL IntB10101006(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=6000
   CASE(1101010)
     CALL IntB10101001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=1000
   CASE(2101010)
     CALL IntB10101002(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=4000
   CASE(3101010)
     CALL IntB10101003(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=3000
   CASE(6101010)
     CALL IntB10101006(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=6000
   CASE(10101010)
     CALL IntB10101010(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=10000
   CASE DEFAULT
      WRITE(*,*) 'We are in ERIInterface.Inc'
      WRITE(*,*) 'IntType=',IntType
      STOP 'MISS AN INTEGRAL'
   END SELECT
