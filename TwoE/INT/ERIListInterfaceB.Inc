   SELECT CASE(IntType)
   CASE(101)
     OA=1;LDA=1;OB=0;LDB=1;OC=0;LDC=1;OD=0;LDD=1;
     CALL IntB1010101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GMc%PBC,C(1))
     LocNInt=1
   CASE(201)
     OA=1;LDA=1;OB=0;LDB=4;OC=0;LDC=4;OD=0;LDD=16;
     CALL IntB2010201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GMc%PBC,C(1))
     LocNInt=16
   CASE(301)
     OA=1;LDA=1;OB=0;LDB=3;OC=0;LDC=3;OD=0;LDD=9;
     CALL IntB3010301(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GMc%PBC,C(1))
     LocNInt=9
   CASE(601)
     OA=1;LDA=1;OB=0;LDB=6;OC=0;LDC=6;OD=0;LDD=36;
     CALL IntB6010601(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GMc%PBC,C(1))
     LocNInt=36
   CASE(1001)
     OA=1;LDA=1;OB=0;LDB=10;OC=0;LDC=10;OD=0;LDD=100;
     CALL IntB10011001(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GMc%PBC,C(1))
     LocNInt=100
   CASE(202)
     OA=1;LDA=1;OB=0;LDB=4;OC=0;LDC=16;OD=0;LDD=64;
     CALL IntB2020202(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GMc%PBC,C(1))
     LocNInt=256
   CASE(302)
     OA=1;LDA=1;OB=0;LDB=3;OC=0;LDC=12;OD=0;LDD=36;
     CALL IntB3020302(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GMc%PBC,C(1))
     LocNInt=144
   CASE(602)
     OA=1;LDA=1;OB=0;LDB=6;OC=0;LDC=24;OD=0;LDD=144;
     CALL IntB6020602(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GMc%PBC,C(1))
     LocNInt=576
   CASE(1002)
     OA=1;LDA=1;OB=0;LDB=10;OC=0;LDC=40;OD=0;LDD=400;
     CALL IntB10021002(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GMc%PBC,C(1))
     LocNInt=1600
   CASE(303)
     OA=1;LDA=1;OB=0;LDB=3;OC=0;LDC=9;OD=0;LDD=27;
     CALL IntB3030303(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GMc%PBC,C(1))
     LocNInt=81
   CASE(603)
     OA=1;LDA=1;OB=0;LDB=6;OC=0;LDC=18;OD=0;LDD=108;
     CALL IntB6030603(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GMc%PBC,C(1))
     LocNInt=324
   CASE(1003)
     OA=1;LDA=1;OB=0;LDB=10;OC=0;LDC=30;OD=0;LDD=300;
     CALL IntB10031003(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GMc%PBC,C(1))
     LocNInt=900
   CASE(606)
     OA=1;LDA=1;OB=0;LDB=6;OC=0;LDC=36;OD=0;LDD=216;
     CALL IntB6060606(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GMc%PBC,C(1))
     LocNInt=1296
   CASE(1006)
     OA=1;LDA=1;OB=0;LDB=10;OC=0;LDC=60;OD=0;LDD=600;
     CALL IntB10061006(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GMc%PBC,C(1))
     LocNInt=3600
   CASE(1010)
     OA=1;LDA=1;OB=0;LDB=10;OC=0;LDC=100;OD=0;LDD=1000;
     CALL IntB10101010(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  OA,LDA,OB,LDB,OC,LDC,OD,LDD,GMc%PBC,C(1))
     LocNInt=10000
   CASE DEFAULT
      WRITE(*,*) 'We are in ERIListInterface.Inc'
      WRITE(*,*) 'IntType=',IntType
      STOP 'MISS AN INTEGRAL'
   END SELECT
