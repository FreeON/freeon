   INTEGER,     PARAMETER                    :: NGrid = 35
   REAL(DOUBLE),PARAMETER                    :: Grd0=Zero
   REAL(DOUBLE),PARAMETER                    :: Grd1=0.8864052604279183D0
   REAL(DOUBLE),PARAMETER                    :: Grd2=One
   REAL(DOUBLE),PARAMETER                    :: Grd3=0.6158817620514397D0
   REAL(DOUBLE),PARAMETER                    :: Grd4=One
   REAL(DOUBLE),PARAMETER                    :: Wht0=0.9246403428221610D0
   REAL(DOUBLE),PARAMETER                    :: Wht1=0.3839608203244567D0
   REAL(DOUBLE),PARAMETER                    :: Wht2=0.5925925925925926D-1
   REAL(DOUBLE),PARAMETER                    :: Wht3=0.4976839098051219D0
   REAL(DOUBLE),PARAMETER                    :: Wht4=0.9876543209876543D-2
   REAL(DOUBLE),PARAMETER,DIMENSION(NGrid)   :: CubeRuleWght=(/             Wht0, &
                                                   Wht1,Wht1,Wht1,Wht1,Wht1,Wht1, &
                     Wht2,Wht2,Wht2,Wht2,Wht2,Wht2,Wht2,Wht2,Wht2,Wht2,Wht2,Wht2, &
                                         Wht3,Wht3,Wht3,Wht3,Wht3,Wht3,Wht3,Wht3, &
                                         Wht4,Wht4,Wht4,Wht4,Wht4,Wht4,Wht4,Wht4/)
   REAL(DOUBLE),PARAMETER,DIMENSION(3,NGrid) :: CubeRuleGrid=RESHAPE( (/       & 
                                                      Zero, Zero, Zero, & ! 0
                                                      Grd1, Zero, Zero, & ! 1
                                                     -Grd1, Zero, Zero, & ! 2
                                                      Zero, Grd1, Zero, & ! 3
                                                      Zero,-Grd1, Zero, & ! 4
                                                      Zero, Zero, Grd1, & ! 5
                                                      Zero, Zero,-Grd1, & ! 6
                                                      Zero, Grd2, Grd2, & ! 7
                                                      Grd2, Zero, Grd2, & ! 8
                                                      Grd2, Grd2, Zero, & ! 9
                                                      Zero,-Grd2, Grd2, & ! 10
                                                      Zero, Grd2,-Grd2, & ! 11
                                                      Zero,-Grd2,-Grd2, & ! 12
                                                     -Grd2, Zero, Grd2, & ! 13
                                                      Grd2, Zero,-Grd2, & ! 14
                                                     -Grd2, Zero,-Grd2, & ! 15
                                                     -Grd2, Grd2, Zero, & ! 16
                                                      Grd2,-Grd2, Zero, & ! 17
                                                     -Grd2,-Grd2, Zero, & ! 18
                                                      Grd3, Grd3, Grd3, & ! 19
                                                     -Grd3, Grd3, Grd3, & ! 20
                                                      Grd3,-Grd3, Grd3, & ! 21
                                                      Grd3, Grd3,-Grd3, & ! 22
                                                     -Grd3, Grd3,-Grd3, & ! 23
                                                      Grd3,-Grd3,-Grd3, & ! 24
                                                     -Grd3,-Grd3, Grd3, & ! 25
                                                     -Grd3,-Grd3,-Grd3, & ! 26
                                                      Grd4, Grd4, Grd4, & ! 27
                                                     -Grd4, Grd4, Grd4, & ! 28
                                                      Grd4,-Grd4, Grd4, & ! 29
                                                      Grd4, Grd4,-Grd4, & ! 30
                                                     -Grd4, Grd4,-Grd4, & ! 31
                                                      Grd4,-Grd4,-Grd4, & ! 32
                                                     -Grd4,-Grd4, Grd4, & ! 33
                                                     -Grd4,-Grd4,-Grd4  & ! 34
                                                     /), (/3,NGrid/) )
