      FUNCTION ChebErf(X)
      REAL*8 ChebErf,X,Switch,Grid
      REAL*8 Erf_0(0:1000),Erf_1(0:1000),
     >       Erf_2(0:1000),Erf_3(0:1000)
      INCLUDE 'Erf.Inc'
      PARAMETER(Switch=5.5D0,Grid=1.82D+02)
      IF(X.GE.Switch)THEN
         ChebErf=1.0D0
      ELSE
         J=AINT(X*Grid)
         ChebErf=Erf_0(J)+X*(Erf_1(J)+X*(Erf_2(J)+X*Erf_3(J))) 
      ENDIF                                                      
      RETURN
      END


