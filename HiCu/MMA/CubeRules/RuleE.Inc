   INTEGER,    PARAMETER                     :: NGrid = 8
   REAL(DOUBLE),PARAMETER                    :: Grd2=0.5773502691896257645D0
   REAL(DOUBLE),PARAMETER                    :: Wht2=1.0D0
   REAL(DOUBLE),PARAMETER,DIMENSION(NGrid)   :: CubeRuleWght=(/Wht2,Wht2,Wht2,Wht2,Wht2,Wht2,Wht2,Wht2/)

  REAL(DOUBLE),PARAMETER,DIMENSION(3,NGrid) :: CubeRuleGrid=RESHAPE( (/       & 
                                                      Grd2, Grd2, Grd2, & ! 7
                                                     -Grd2, Grd2, Grd2, & ! 8
                                                      Grd2,-Grd2, Grd2, & ! 9
                                                      Grd2, Grd2,-Grd2, & ! 10
                                                     -Grd2, Grd2,-Grd2, & ! 11
                                                      Grd2,-Grd2,-Grd2, & ! 12
                                                     -Grd2,-Grd2, Grd2, & ! 13
                                                     -Grd2,-Grd2,-Grd2  & ! 14
                                                     /), (/3,NGrid/) )

