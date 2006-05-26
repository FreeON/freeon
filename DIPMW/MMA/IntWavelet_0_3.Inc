    IntD0L0 = One/(2.D0*SqZeta*XZeta1)
    IntD0L1 = IntD0L0*SqZeta
    IntD0L2 = IntD0L1*SqZeta
    IntD0L3 = IntD0L2*SqZeta
    SELECT CASE(SType)
    CASE(-1)
        RZeta2 = RZeta1*RZeta1
        t1 = -1.D0*Exp0 + ExpRmX0 - 1.772453850905516D0*Erf0*RZeta1
        t2 = 1.772453850905516D0*ErfRmX0*RZeta1 + 1.772453850905516D0*Erf0*XZeta1 - 1.772453850905516D0*ErfRmX0*XZeta1
        WV(0,0) = IntD0L0*(t1 + t2)
        WV(0,1) = IntD0L1*(1.772453850905516D0*Erf0 - 1.772453850905516D0*ErfRmX0 - 2.D0*Exp0*XZeta1)
        WV(0,2) = IntD0L2*(-2.D0*Exp0 + 2.D0*ExpRmX0 - 4.D0*Exp0*RZeta1*XZeta1)
        t1 = -4.D0*Exp0*RZeta1 + 4.D0*ExpRmX0*RZeta1
        t2 = 4.D0*Exp0*XZeta1 - 4.D0*ExpRmX0*XZeta1 - 8.D0*Exp0*RZeta2*XZeta1
        WV(0,3) = IntD0L3*(t1 + t2)
    CASE(0)
        t1 = -2.D0*Exp0 + ExpRmX0 + ExpRpX0 - 3.544907701811032D0*Erf0*RZeta1
        t2 = 1.772453850905516D0*ErfRmX0*RZeta1 + 1.772453850905516D0*ErfRpX0*RZeta1
        t3 = -1.772453850905516D0*ErfRmX0*XZeta1 + 1.772453850905516D0*ErfRpX0*XZeta1
        WV(0,0) = IntD0L0*(t1 + t2 + t3)
        WV(0,1) = (3.544907701811032D0*Erf0 - 1.772453850905516D0*ErfRmX0 - 1.772453850905516D0*ErfRpX0)*IntD0L1
        WV(0,2) = (-4.D0*Exp0 + 2.D0*ExpRmX0 + 2.D0*ExpRpX0)*IntD0L2
        t1 = -8.D0*Exp0*RZeta1 + 4.D0*ExpRmX0*RZeta1
        t2 = 4.D0*ExpRpX0*RZeta1 - 4.D0*ExpRmX0*XZeta1 + 4.D0*ExpRpX0*XZeta1
        WV(0,3) = IntD0L3*(t1 + t2)
    CASE(1)
        RZeta2 = RZeta1*RZeta1
        t1 = -1.D0*Exp0 + ExpRpX0 - 1.772453850905516D0*Erf0*RZeta1
        t2 = 1.772453850905516D0*ErfRpX0*RZeta1 - 1.772453850905516D0*Erf0*XZeta1 + 1.772453850905516D0*ErfRpX0*XZeta1
        WV(0,0) = IntD0L0*(t1 + t2)
        WV(0,1) = IntD0L1*(1.772453850905516D0*Erf0 - 1.772453850905516D0*ErfRpX0 + 2.D0*Exp0*XZeta1)
        WV(0,2) = IntD0L2*(-2.D0*Exp0 + 2.D0*ExpRpX0 + 4.D0*Exp0*RZeta1*XZeta1)
        t1 = -4.D0*Exp0*RZeta1 + 4.D0*ExpRpX0*RZeta1
        t2 = -4.D0*Exp0*XZeta1 + 4.D0*ExpRpX0*XZeta1 + 8.D0*Exp0*RZeta2*XZeta1
        WV(0,3) = IntD0L3*(t1 + t2)
    END SELECT
