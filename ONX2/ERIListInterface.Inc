   SELECT CASE(IntType)
   CASE(101)
     CALL IntSSSS5(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, &
           &        ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L,C(1),GMc%PBC)
     LocNInt=1
   CASE(301)
     CALL Int3131(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,3,0,3,0,9,GMc%PBC,C(1))
     LocNInt=9
   CASE(103)
     CALL Int3131(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,3,0,3,0,9,GMc%PBC,C(1))
     LocNInt=9
   CASE(303)
     CALL Int3333(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & 
                  ACAtmPair(CFAC)%SP%AtmInfo,ACAtmPair(CFAC)%SP%AtmInfo, & 
                  1,1,0,3,0,9,0,27,GMc%PBC,C(1))
     LocNInt=81
   CASE DEFAULT
      STOP 'MISS AN INTEGRAL'
   END SELECT
