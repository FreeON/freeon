   SELECT CASE(IntType)
   CASE(1010101)
     CALL dIntB1010101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=12
   CASE(3010101)
     CALL dIntB3010101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=36
   CASE(6010101)
     CALL dIntB6010101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=72
   CASE(10010101)
     CALL dIntB10010101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=120
   CASE(3030101)
     CALL dIntB3030101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=108
   CASE(6030101)
     CALL dIntB6030101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=216
   CASE(10030101)
     CALL dIntB10030101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=360
   CASE(6060101)
     CALL dIntB6060101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=432
   CASE(10060101)
     CALL dIntB10060101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=720
   CASE(10100101)
     CALL dIntB10100101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=1200
   CASE(1010301)
     CALL dIntB3010101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=36
   CASE(3010301)
     CALL dIntB3010301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=108
   CASE(6010301)
     CALL dIntB6010301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=216
   CASE(10010301)
     CALL dIntB10010301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=360
   CASE(3030301)
     CALL dIntB3030301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=324
   CASE(6030301)
     CALL dIntB6030301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=648
   CASE(10030301)
     CALL dIntB10030301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=1080
   CASE(6060301)
     CALL dIntB6060301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=1296
   CASE(10060301)
     CALL dIntB10060301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=2160
   CASE(10100301)
     CALL dIntB10100301(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=3600
   CASE(1010601)
     CALL dIntB6010101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=72
   CASE(3010601)
     CALL dIntB6010301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=216
   CASE(6010601)
     CALL dIntB6010601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=432
   CASE(10010601)
     CALL dIntB10010601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=720
   CASE(3030601)
     CALL dIntB6010303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=648
   CASE(6030601)
     CALL dIntB6030601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=1296
   CASE(10030601)
     CALL dIntB10030601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=2160
   CASE(6060601)
     CALL dIntB6060601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=2592
   CASE(10060601)
     CALL dIntB10060601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=4320
   CASE(10100601)
     CALL dIntB10100601(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=7200
   CASE(1011001)
     CALL dIntB10010101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=120
   CASE(3011001)
     CALL dIntB10010301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=360
   CASE(6011001)
     CALL dIntB10010601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=720
   CASE(10011001)
     CALL dIntB10011001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=1200
   CASE(3031001)
     CALL dIntB10010303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=1080
   CASE(6031001)
     CALL dIntB10010603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=2160
   CASE(10031001)
     CALL dIntB10031001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=3600
   CASE(6061001)
     CALL dIntB10010606(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=4320
   CASE(10061001)
     CALL dIntB10061001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=7200
   CASE(10101001)
     CALL dIntB10101001(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=12000
   CASE(1010303)
     CALL dIntB3030101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=108
   CASE(3010303)
     CALL dIntB3030301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=324
   CASE(6010303)
     CALL dIntB6010303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=648
   CASE(10010303)
     CALL dIntB10010303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=1080
   CASE(3030303)
     CALL dIntB3030303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=972
   CASE(6030303)
     CALL dIntB6030303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=1944
   CASE(10030303)
     CALL dIntB10030303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=3240
   CASE(6060303)
     CALL dIntB6060303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=3888
   CASE(10060303)
     CALL dIntB10060303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=6480
   CASE(10100303)
     CALL dIntB10100303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=10800
   CASE(1010603)
     CALL dIntB6030101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=216
   CASE(3010603)
     CALL dIntB6030301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=648
   CASE(6010603)
     CALL dIntB6030601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=1296
   CASE(10010603)
     CALL dIntB10010603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=2160
   CASE(3030603)
     CALL dIntB6030303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=1944
   CASE(6030603)
     CALL dIntB6030603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=3888
   CASE(10030603)
     CALL dIntB10030603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=6480
   CASE(6060603)
     CALL dIntB6060603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=7776
   CASE(10060603)
     CALL dIntB10060603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=12960
   CASE(10100603)
     CALL dIntB10100603(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=21600
   CASE(1011003)
     CALL dIntB10030101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=360
   CASE(3011003)
     CALL dIntB10030301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=1080
   CASE(6011003)
     CALL dIntB10030601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=2160
   CASE(10011003)
     CALL dIntB10031001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=3600
   CASE(3031003)
     CALL dIntB10030303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=3240
   CASE(6031003)
     CALL dIntB10030603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=6480
   CASE(10031003)
     CALL dIntB10031003(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=10800
   CASE(6061003)
     CALL dIntB10030606(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=12960
   CASE(10061003)
     CALL dIntB10061003(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=21600
   CASE(10101003)
     CALL dIntB10101003(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=36000
   CASE(1010606)
     CALL dIntB6060101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=432
   CASE(3010606)
     CALL dIntB6060301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=1296
   CASE(6010606)
     CALL dIntB6060601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=2592
   CASE(10010606)
     CALL dIntB10010606(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=4320
   CASE(3030606)
     CALL dIntB6060303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=3888
   CASE(6030606)
     CALL dIntB6060603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=7776
   CASE(10030606)
     CALL dIntB10030606(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=12960
   CASE(6060606)
     CALL dIntB6060606(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=15552
   CASE(10060606)
     CALL dIntB10060606(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=25920
   CASE(10100606)
     CALL dIntB10100606(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=43200
   CASE(1011006)
     CALL dIntB10060101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=720
   CASE(3011006)
     CALL dIntB10060301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=2160
   CASE(6011006)
     CALL dIntB10060601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=4320
   CASE(10011006)
     CALL dIntB10061001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=7200
   CASE(3031006)
     CALL dIntB10060303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=6480
   CASE(6031006)
     CALL dIntB10060603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=12960
   CASE(10031006)
     CALL dIntB10061003(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=21600
   CASE(6061006)
     CALL dIntB10060606(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=25920
   CASE(10061006)
     CALL dIntB10061006(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=43200
   CASE(10101006)
     CALL dIntB10101006(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=72000
   CASE(1011010)
     CALL dIntB10100101(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=1200
   CASE(3011010)
     CALL dIntB10100301(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=3600
   CASE(6011010)
     CALL dIntB10100601(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=7200
   CASE(10011010)
     CALL dIntB10101001(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=12000
   CASE(3031010)
     CALL dIntB10100303(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=10800
   CASE(6031010)
     CALL dIntB10100603(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=21600
   CASE(10031010)
     CALL dIntB10101003(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=36000
   CASE(6061010)
     CALL dIntB10100606(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=43200
   CASE(10061010)
     CALL dIntB10101006(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=72000
   CASE(10101010)
     CALL dIntB10101010(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk,GMc%PBC,C(1))
     LocNInt=120000
   CASE DEFAULT
      WRITE(*,*) 'We are in dERIInterfaceB.Inc'
      WRITE(*,*) 'IntType=',IntType
      STOP 'MISS AN INTEGRAL'
   END SELECT
