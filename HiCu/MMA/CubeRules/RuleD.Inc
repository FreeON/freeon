   INTEGER,    PARAMETER                     :: NGrid = 6
   REAL(DOUBLE),PARAMETER                    :: Grd1=1.0D0
   REAL(DOUBLE),PARAMETER                    :: Wht1=4.0D0/3.0D0
   REAL(DOUBLE),PARAMETER,DIMENSION(NGrid)   :: CubeRuleWght=(/Wht1,Wht1,Wht1,Wht1,Wht1,Wht1/)
   REAL(DOUBLE),PARAMETER,DIMENSION(3,NGrid) :: CubeRuleGrid=RESHAPE( (/       & 
                                                      Grd1, Zero, Zero, & ! 1
                                                     -Grd1, Zero, Zero, & ! 2
                                                      Zero, Grd1, Zero, & ! 3
                                                      Zero,-Grd1, Zero, & ! 4
                                                      Zero, Zero, Grd1, & ! 5
                                                      Zero, Zero,-Grd1  & ! 6
                                                     /), (/3,NGrid/) )

