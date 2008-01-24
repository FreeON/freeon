   INTEGER,     PARAMETER                    :: NGrid = 14
   REAL(DOUBLE),PARAMETER                    :: Grd1=0.79582242575422146D0
   REAL(DOUBLE),PARAMETER                    :: Grd2=0.75878691063932815D0
   REAL(DOUBLE),PARAMETER                    :: Wht1=0.88642659279778393D0
   REAL(DOUBLE),PARAMETER                    :: Wht2=0.33518005540166205D0
   REAL(DOUBLE),PARAMETER,DIMENSION(NGrid)   :: CubeRuleWght=(/Wht1,Wht1,Wht1,Wht1,Wht1,Wht1, &
                                                     Wht2,Wht2,Wht2,Wht2,Wht2,Wht2,Wht2,Wht2/)
   REAL(DOUBLE),PARAMETER,DIMENSION(3,NGrid) :: CubeRuleGrid=RESHAPE( (/       & 
                                                      Grd1, Zero, Zero, & ! 1
                                                     -Grd1, Zero, Zero, & ! 2
                                                      Zero, Grd1, Zero, & ! 3
                                                      Zero,-Grd1, Zero, & ! 4
                                                      Zero, Zero, Grd1, & ! 5
                                                      Zero, Zero,-Grd1, & ! 6
                                                      Grd2, Grd2, Grd2, & ! 7
                                                     -Grd2, Grd2, Grd2, & ! 8
                                                      Grd2,-Grd2, Grd2, & ! 9
                                                      Grd2, Grd2,-Grd2, & ! 10
                                                     -Grd2, Grd2,-Grd2, & ! 11
                                                      Grd2,-Grd2,-Grd2, & ! 12
                                                     -Grd2,-Grd2, Grd2, & ! 13
                                                     -Grd2,-Grd2,-Grd2  & ! 14
                                                     /), (/3,NGrid/) )
