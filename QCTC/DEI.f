      FUNCTION DEI ( X1 )
C
C  AN EXPONENTIAL INTEGRAL ROUTINE.
C  FOR X GREATER THAN 0, THE EXPONENTIAL INTEGRAL, EI, IS DEFINED BY
C  EI(X) = INTEGRAL ( EXP ( T ) / T DT ), FROM T = -INFINITY TO T = X
C  WHERE THE INTEGRAL IS TO BE INTERPRETED AS THE CAUCHY PRINCIPAL
C  VALUE.  FOR X LESS THAN 0, EI(X) = -E1(-X), WHERE
C  E1(Z) = INTEGRAL ( EXP ( -T ) / T DT ) FROM T = Z TO T = INFINITY.
C
C
C  Modified:
C
C    04 October 2006
C
C  Reference:
C
C    Kathleen Paciorek,
C    Algorithm 385:
C    Exponential Integral Ei(x),
C    Communications of the ACM,
C    Volume 13, Number 7, July 1970, pages 446-447.
C

      IMPLICIT NONE

      DOUBLE PRECISION A(6)
      DOUBLE PRECISION B(6)
      DOUBLE PRECISION C(8)
      DOUBLE PRECISION D(8)
      DOUBLE PRECISION DEI
      DOUBLE PRECISION DENM
      DOUBLE PRECISION E(8)
      DOUBLE PRECISION F(8)
      DOUBLE PRECISION FRAC
      INTEGER I
      INTEGER J
      INTEGER L
      DOUBLE PRECISION P0(6)
      DOUBLE PRECISION P1(9)
      DOUBLE PRECISION P2(9)
      DOUBLE PRECISION P3(10)
      DOUBLE PRECISION P4(10)
      DOUBLE PRECISION PX(10)
      DOUBLE PRECISION Q0(6)
      DOUBLE PRECISION Q1(9)
      DOUBLE PRECISION Q2(8)
      DOUBLE PRECISION Q3(9)
      DOUBLE PRECISION Q4(9)
      DOUBLE PRECISION QX(10)
      DOUBLE PRECISION R
      DOUBLE PRECISION SUMP
      DOUBLE PRECISION SUMQ
      DOUBLE PRECISION T
      DOUBLE PRECISION W
      DOUBLE PRECISION X
      DOUBLE PRECISION X0
      DOUBLE PRECISION X1
      DOUBLE PRECISION XMX0
      DOUBLE PRECISION Y, MAXEXP

      SAVE A
      SAVE B
      SAVE C
      SAVE D
      SAVE E
      SAVE F
      SAVE P0
      SAVE P1
      SAVE P2
      SAVE P3
      SAVE P4
      SAVE Q0
      SAVE Q1
      SAVE Q2
      SAVE Q3
      SAVE Q4
      SAVE X0

      DATA A /
     &  -5.77215664901532863D-01,
     &   7.54164313663016620D-01,
     &   1.29849232927373234D-01,
     &   2.40681355683977413D-02,
     &   1.32084309209609371D-03,
     &   6.57739399753264501D-05 /
      DATA B /
     &  1.0D+00,
     &  4.25899193811589822D-01,
     &  7.9779471841022822D-02,
     &  8.30208476098771677D-03,
     &  4.86427138393016416D-04,
     &  1.30655195822848878D-05 /
      DATA C /
     &  8.67745954838443744D-08,
     &  9.99995519301390302D-01,
     &  1.18483105554945844D+01,
     &  4.55930644253389823D+01,
     &  6.99279451291003023D+01,
     &  4.25202034768840779D+01,
     &  8.83671808803843939D+00,
     &  4.01377664940664720D-01 /
      DATA D /
     &  1.0D+00,
     &  1.28481935379156650D+01,
     &  5.64433569561803199D+01,
     &  1.06645183769913883D+02,
     &  8.97311097125289802D+01,
     &  3.14971849170440750D+01,
     &  3.79559003762122243D+00,
     &  9.08804569188869219D-02 /
      DATA E /
     &  -9.99999999999973414D-01,
     &  -3.44061995006684895D+01,
     &  -4.27532671201988539D+02,
     &  -2.39601943247490540D+03,
     &  -6.16885210055476351D+03,
     &  -6.57609698748021179D+03,
     &  -2.10607737142633289D+03,
     &  -1.48990849972948169D+01 /
      DATA F /
     &  1.0D+00,
     &  3.64061995006459804D+01,
     &  4.94345070209903645D+02,
     &  3.19027237489543304D+03,
     &  1.03370753085840977D+04,
     &  1.63241453557783503D+04,
     &  1.11497752871096620D+04,
     &  2.37813899102160221D+03 /
      DATA P0 /
     &  1.0D+00,
     &  2.23069937666899751D+00,
     &  1.70277059606809295D+00,
     &  5.10499279623219400D-01,
     &  4.89089253789279154D-02,
     &  3.65462224132368429D-04 /
      DATA P1 /
     &  5.99569946892370010D+09,
     & -2.50389994886351362D+08,
     &  7.05921609590056747D+08,
     & -3.36899564201591901D+06,
     &  8.98683291643758313D+06,
     &  7.37147790184657443D+04,
     &  2.85446881813647015D+04,
     &  4.12626667248911939D+02,
     &  1.10639547241639580D+01 /
      DATA P2 /
     &  9.98957666516551704D-01,
     &  5.73116705744508018D+00,
     &  4.18102422562856622D+00,
     &  5.88658240753281111D+00,
     & -1.94132967514430702D+01,
     &  7.89472209294457221D+00,
     &  2.32730233839039141D+01,
     & -3.67783113478311458D+01,
     & -2.46940983448361265D+00 /
      DATA P3 /
     &  9.99993310616056874D-01,
     & -1.84508623239127867D+00,
     &  2.65257581845279982D+01,
     &  2.49548773040205944D+01,
     & -3.32361257934396228D+01,
     & -9.13483569999874255D-01,
     & -2.10574079954804045D+01,
     & -1.00064191398928483D+01,
     & -1.86009212172643758D+01,
     & -1.64772117246346314D+00 /
      DATA P4 /
     &  1.00000000000000486D+00,
     & -3.00000000320981266D+00,
     & -5.00006640413131002D+00,
     & -7.06810977895029359D+00,
     & -1.52856623636929637D+01,
     & -7.63147701620253631D+00,
     & -2.79798528624305389D+01,
     & -1.81949664929868906D+01,
     & -2.23127670777632410D+02,
     &  1.75338801265465972D+02 /
      DATA Q0 /
     &  1.0D+00,
     &  2.73069937666899751D+00,
     &  2.73478695106925836D+00,
     &  1.21765962960151532D+00,
     &  2.28817933990526412D-01,
     &  1.31114151194977706D-02 /   
      DATA Q1 /
     &  2.55926497607616350D+09,
     & -2.79673351122984591D+09,
     &  8.02827782946956507D+08,
     & -1.44980714393023883D+08,
     &  1.77158308010799884D+07,
     & -1.49575457202559218D+06,
     &  8.53771000180749097D+04,
     & -3.02523682238227410D+03,
     &  5.12578125D+01 /
      DATA Q2 /
     &  1.14625253249016191D+00,
     & -1.99149600231235164D+02,
     &  3.41365212524375539D+02,
     &  5.23165568734558614D+01,
     &  3.17279489254369328D+02,
     & -8.38767084189640707D+00,
     &  9.65405217429280303D+02,
     &  2.63983007318024593D+00 /
      DATA Q3 /
     &  1.00153385204534270D+00,
     & -1.09355619539109124D+01,
     &  1.99100447081774247D+02,
     &  1.19283242396860101D+03,
     &  4.42941317833792840D+01,
     &  2.53881931563070803D+02,
     &  5.99493232566740736D+01,
     &  6.40380040535241555D+01,
     &  9.79240359921729030D+01 /
      DATA Q4 /
     &  1.99999999999048104D+00,
     & -2.99999894040324960D+00,
     & -7.99243595776339741D+00,
     & -1.20187763547154743D+01,
     &  7.04831847180424676D+01,
     &  1.17179220502086455D+02,
     &  1.37790390235747999D+02,
     &  3.97277109100414518D+00,
     &  3.97845977167414721D+04 /
      DATA X0 / 0.372507410781366634D+00 /

C MAXEXP needs to be set to the largest argument of exp
C that will not cause an overflow. This is computed here
c but could be embedded as a constant for efficiency reasons.
      MAXEXP = (INT(LOG(HUGE(0.0D0))*100))/100.0D0

      X = X1
1     IF ( X .LE. 0.0D+00 ) GO TO 100
      IF ( X .GE. 12.0D+00 ) GO TO 60
      IF ( X .GE. 6.0D+00 ) GO TO 40
C
C  X IN (0,6).
C
      T = X + X
      T = T / 3.0D+00 - 2.0D+00
      PX(10) = 0.0D+00
      QX(10) = 0.0D+00
      PX(9) = P1(9)
      QX(9) = Q1(9)
C
C  THE RATIONAL FUNCTION IS EXPRESSED AS A RATIO OF FINITE SUMS OF
C  SHIFTED CHEBYSHEV POLYNOMIALS, AND IS EVALUATED BY NOTING THAT
C  T*(X) = T(2*X-1) AND USING THE CLENSHAW-RICE ALGORITHM FOUND IN
C  REFERENCE (4).
C
      DO 10 L = 2, 8
        I = 10 - L
        PX(I) = T * PX(I+1) - PX(I+2) + P1(I)
        QX(I) = T * QX(I+1) - QX(I+2) + Q1(I)
10    CONTINUE

      R = ( 0.5D+00 * T * PX(2) - PX(3) + P1(1) ) 
     &  / ( 0.5D+00 * T * QX(2) - QX(3) + Q1(1) )
C
C  ( X - X0 ) = ( X - X1 ) - X2, WHERE X1 = 409576229586. / 2**40 AND
C  X2 = -.7671772501993940D-12.
C
      XMX0 = ( X - 409576229586.0D+00 / 1099511627776.0D+00 ) 
     & - 0.7671772501993940D-12
      IF ( DABS ( XMX0 ) .LT. 0.037D+00 ) GO TO 15
      DEI = LOG ( X / X0 ) + XMX0 * R
      RETURN
15    Y = XMX0 / X0
C
C  A RATIONAL APPROXIMATION TO LOG ( X / X0 ) * LOG ( 1 + Y ), 
C  WHERE Y = ( X - X0 ) / X0, AND DABS ( Y ) IS LESS THAN 0.1,
C  THAT IS FOR DABS ( X - X0 ) LESS THAN 0.037.
C
      SUMP = (((( P0(6) 
     &      * Y + P0(5) ) 
     &      * Y + P0(4) ) 
     &      * Y + P0(3) ) 
     &      * Y + P0(2) ) 
     &      * Y + P0(1)

      SUMQ = (((( Q0(6) 
     &      * Y + Q0(5) ) 
     &      * Y + Q0(4) ) 
     &      * Y + Q0(3) ) 
     &      * Y + Q0(2) ) 
     &      * Y + Q0(1)

      DEI = ( SUMP / ( SUMQ * X0 ) + R ) * XMX0
      RETURN
C
C  X IN (6,12).
C
40    DENM = P2(9) + X
      FRAC = Q2(8) + X
C
C  THE RATIONAL FUNCTION IS EXPRESSED AS A J-FRACTION.
C
      DO 25 J = 2, 8
        I = 9 - J
        DENM = P2(I+1) + X + FRAC
        FRAC = Q2(I) / DENM
25    CONTINUE

      DEI = EXP ( X ) * ( ( P2(1) + FRAC ) / X )
      RETURN

60    IF ( X .GE. 24.0D+00 ) GO TO 80
C
C  X IN (12,24).
C
      DENM = P3(10) + X
      FRAC = Q3(9) / DENM
C
C  THE RATIONAL FUNCTION IS EXPRESSED AS A J-FRACTION.
C
      DO 26 J = 2, 9
        I = 10 - J
        DENM = P3(I+1) + X + FRAC
        FRAC = Q3(I) / DENM
26    CONTINUE

      DEI = EXP ( X ) * ( ( P3(1) + FRAC ) / X )
      RETURN
C
C  X GREATER THAN 24.
C
80    IF ( X .LE. MAXEXP ) GO TO 90
C
C  X IS GREATER THAN MAXEXP AND DEI IS SET TO INFINITY.
C
      DEI = HUGE(0.0D0)
      RETURN
90    Y = 1.0D+00 / X
      DENM = P4(10) + X
      FRAC = Q4(9) / DENM
C
C  THE RATIONAL FUNCTION IS EXPRESSED AS A J-FRACTION.
C
      DO 28 J = 2, 9
        I = 10 - J
        DENM = P4(I+1) + X + FRAC
        FRAC = Q4(I) / DENM
28    CONTINUE

      DEI = EXP ( X ) * ( Y + Y * Y * ( P4(1) + FRAC ) )
      RETURN

100   IF ( X .NE. 0.0D+00 ) GO TO 101
C
C  X = 0 AND DEI IS SET TO -INFINITY.
C
      DEI = -HUGE(0.0D0)
      WRITE(*,500)
500   FORMAT ( 
     & ' DEI CALLED WITH A ZERO ARGUMENT, RESULT SET TO -INFINITY')
      RETURN
101   Y = -X
110   W = 1.0D+00 / Y
      IF ( Y .GT. 4.0D+00 ) GO TO 300
      IF ( Y .GT. 1.0D+00 ) GO TO 200
C
C  X IN (-1,0).
C
      DEI = LOG ( Y ) - (((((
     &        A(6)
     &  * Y + A(5) )
     &  * Y + A(4) )
     &  * Y + A(4) )
     &  * Y + A(2) )
     &  * Y + A(1) ) / (((((
     &        B(6)
     &  * Y + B(5) )
     &  * Y + B(4) )
     &  * Y + B(4) )
     &  * Y + B(2) )
     &  * Y + B(1) )
      RETURN
C
C  X IN (-4,-1).
C
200   DEI = -EXP ( -Y ) * ((((((((
     &        C(8)
     &  * W + C(7) )
     &  * W + C(6) )
     &  * W + C(5) )
     &  * W + C(4) )
     &  * W + C(3) )
     &  * W + C(2) )
     &  * W + C(1) ) / (((((((
     &        D(8)
     &  * W + D(7) )
     &  * W + D(6) )
     &  * W + D(5) )
     &  * W + D(4) )
     &  * W + D(3) )
     &  * W + D(2) )
     &  * W + D(1) ) )
      RETURN
C
C  X LESS THAN -4.
C
300   DEI = -EXP ( -Y ) * ( W * ( 1.0D+00 + W * (((((((
     &        E(8)
     &  * W + E(7) )
     &  * W + E(6) )
     &  * W + E(5) )
     &  * W + E(4) )
     &  * W + E(3) )
     &  * W + E(2) )
     &  * W + E(1) ) / (((((((
     &        F(8)
     &  * W + F(7) )
     &  * W + F(6) )
     &  * W + F(5) )
     &  * W + F(4) )
     &  * W + F(3) )
     &  * W + F(2) )
     &  * W + F(1) ) ) )

      RETURN
      END
