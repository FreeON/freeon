   SELECT CASE(IntType)
   CASE(1010101)
      Off=    OffSet%A+ &
           & (OffSet%B-1)*NBFA+ &
           & (OffSet%C-1)*NBFA*NBFB+ &
           & (OffSet%D-1)*NBFA*NBFB*NBFC
      CALL IntSSSS5(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, &
           &        BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L,C(Off),GMc%PBC)
     LocNInt=1
   CASE(3010101)
     CALL Int3111(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=3
   CASE(1030101)
     CALL Int3111(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=3
   CASE(3030101)
     CALL Int3311(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=9
   CASE(1010301)
     CALL Int3111(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=3
   CASE(3010301)
     CALL Int3131(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=9
   CASE(1030301)
     CALL Int3131(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=9
   CASE(3030301)
     CALL Int3331(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=27
   CASE(1010103)
     CALL Int3111(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=3
   CASE(3010103)
     CALL Int3131(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=9
   CASE(1030103)
     CALL Int3131(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=9
   CASE(3030103)
     CALL Int3331(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%B-1,NBFA          ,GMc%PBC,C(1))
     LocNInt=27
   CASE(1010303)
     CALL Int3311(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=9
   CASE(3010303)
     CALL Int3331(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     ,GMc%PBC,C(1))
     LocNInt=27
   CASE(1030303)
     CALL Int3331(BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC, & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%A  ,1             ,GMc%PBC,C(1))
     LocNInt=27
   CASE(3030303)
     CALL Int3333(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=81
   CASE DEFAULT
      STOP 'MISS AN INTEGRAL'
   END SELECT
