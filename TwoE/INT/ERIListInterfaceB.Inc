   SELECT CASE(IntType)
   CASE(101)
     CALL IntB1010101(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,1,0,1,0,1,GMc%PBC,C(1))
     LocNInt=1
   CASE(201)
     CALL IntB2010201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,4,0,4,0,16,GMc%PBC,C(1))
     LocNInt=16
   CASE(301)
     CALL IntB3010301(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,3,0,3,0,9,GMc%PBC,C(1))
     LocNInt=9
   CASE(601)
     CALL IntB6010601(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,6,0,6,0,36,GMc%PBC,C(1))
     LocNInt=36
   CASE(102)
     CALL IntB2010201(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,4,0,4,0,16,GMc%PBC,C(1))
     LocNInt=16
   CASE(202)
     CALL IntB2020202(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,4,0,16,0,64,GMc%PBC,C(1))
     LocNInt=256
   CASE(302)
     CALL IntB3020302(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,3,0,12,0,36,GMc%PBC,C(1))
     LocNInt=144
   CASE(602)
     CALL IntB6020602(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,6,0,24,0,144,GMc%PBC,C(1))
     LocNInt=576
   CASE(103)
     CALL IntB3010301(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,3,0,3,0,9,GMc%PBC,C(1))
     LocNInt=9
   CASE(203)
     CALL IntB3020302(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,3,0,12,0,36,GMc%PBC,C(1))
     LocNInt=144
   CASE(303)
     CALL IntB3030303(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,3,0,9,0,27,GMc%PBC,C(1))
     LocNInt=81
   CASE(603)
     CALL IntB6030603(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,6,0,18,0,108,GMc%PBC,C(1))
     LocNInt=324
   CASE(106)
     CALL IntB6010601(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,6,0,6,0,36,GMc%PBC,C(1))
     LocNInt=36
   CASE(206)
     CALL IntB6020602(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,6,0,24,0,144,GMc%PBC,C(1))
     LocNInt=576
   CASE(306)
     CALL IntB6030603(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,6,0,18,0,108,GMc%PBC,C(1))
     LocNInt=324
   CASE(606)
     CALL IntB6060606(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,6,0,36,0,216,GMc%PBC,C(1))
     LocNInt=1296
   CASE DEFAULT
      WRITE(*,*) 'We are in ERIListInterface.Inc'
      WRITE(*,*) 'IntType=',IntType
      STOP 'MISS AN INTEGRAL'
   END SELECT
