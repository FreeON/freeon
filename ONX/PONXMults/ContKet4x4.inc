c       ------------------------------ 
        t1=Cpq2*CBufK(2,KK,J)+Cpq4*CBufK(4,KK,J)
        t2=Cpq5*CBufK(5,KK,J)+Cpq7*CBufK(7,KK,J)+Cpq8*CBufK(8,KK,J)
        t3=Cpq9*CBufK(9,KK,J)+Cpq11*CBufK(11,KK,J)
        t4=Cpq12*CBufK(12,KK,J)+Cpq13*CBufK(13,KK,J)+Cpq14*CBufK(14,K
     $   K,J)
        t5=Spq2*SBufK(2,KK,J)+Spq4*SBufK(4,KK,J)
        t6=Spq5*SBufK(5,KK,J)+Spq7*SBufK(7,KK,J)+Spq8*SBufK(8,KK,J)
        t7=Spq9*SBufK(9,KK,J)+Spq11*SBufK(11,KK,J)
        t8=Spq12*SBufK(12,KK,J)+Spq13*SBufK(13,KK,J)+Spq14*SBufK(14,K
     $   K,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        CTwo0=t9+t10
        t1=-(Cpq4*CBufK(2,KK,J))-Cpq7*CBufK(4,KK,J)
        t2=-(Cpq8*CBufK(5,KK,J))
        t3=-(Cpq11*CBufK(7,KK,J))-Cpq12*CBufK(8,KK,J)
        t4=-(Cpq13*CBufK(9,KK,J))-Cpq16*CBufK(11,KK,J)
        t5=-(Cpq17*CBufK(12,KK,J))
        t6=-(Cpq18*CBufK(13,KK,J))-Cpq19*CBufK(14,KK,J)
        t7=-(Spq4*SBufK(2,KK,J))-Spq7*SBufK(4,KK,J)
        t8=-(Spq8*SBufK(5,KK,J))
        t9=-(Spq11*SBufK(7,KK,J))-Spq12*SBufK(8,KK,J)
        t10=-(Spq13*SBufK(9,KK,J))-Spq16*SBufK(11,KK,J)
        t11=-(Spq17*SBufK(12,KK,J))
        t12=-(Spq18*SBufK(13,KK,J))-Spq19*SBufK(14,KK,J)
        t13=t1+t2+t3+t4+t5+t6
        t14=t7+t8+t9+t10+t11+t12
        CTwo1=t13+t14
        t1=-(Cpq2*CBufK(0,KK,J))
        t2=-(Cpq4*CBufK(1,KK,J))+(Cpq3-Cpq5)*CBufK(2,KK,J)
        t3=-(Cpq7*CBufK(3,KK,J))
        t4=(Cpq6-Cpq8)*CBufK(4,KK,J)+(Cpq7-Cpq9)*CBufK(5,KK,J)
        t5=-(Cpq11*CBufK(6,KK,J))
        t6=(Cpq10-Cpq12)*CBufK(7,KK,J)+(Cpq11-Cpq13)*CBufK(8,KK,J)
        t7=(Cpq12-Cpq14)*CBufK(9,KK,J)
        t8=-(Cpq16*CBufK(10,KK,J))+(Cpq15-Cpq17)*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=(Cpq16-Cpq18)*CBufK(12,KK,J)
        t12=(Cpq17-Cpq19)*CBufK(13,KK,J)+(Cpq18-Cpq20)*CBufK(14,KK,J)
        t13=-(Spq5*SBufK(2,KK,J))
        t14=-(Spq8*SBufK(4,KK,J))+(Spq7-Spq9)*SBufK(5,KK,J)
        t15=-(Spq12*SBufK(7,KK,J))
        t16=(Spq11-Spq13)*SBufK(8,KK,J)+(Spq12-Spq14)*SBufK(9,KK,J)
        t17=-(Spq17*SBufK(11,KK,J))+(Spq16-Spq18)*SBufK(12,KK,J)
        t18=(Spq17-Spq19)*SBufK(13,KK,J)+(Spq18-Spq20)*SBufK(14,KK,J)
        t19=t11+t12+t13+t14
        t20=t15+t16+t17+t18
        CTwo2=t9+t10+t19+t20
        t1=-(Spq2*CBufK(0,KK,J))
        t2=-(Spq4*CBufK(1,KK,J))-Spq5*CBufK(2,KK,J)
        t3=-(Spq7*CBufK(3,KK,J))
        t4=-(Spq8*CBufK(4,KK,J))-(Spq7+Spq9)*CBufK(5,KK,J)
        t5=-(Spq11*CBufK(6,KK,J))
        t6=-(Spq12*CBufK(7,KK,J))-(Spq11+Spq13)*CBufK(8,KK,J)
        t7=-((Spq12+Spq14)*CBufK(9,KK,J))
        t8=-(Spq16*CBufK(10,KK,J))-Spq17*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=-((Spq16+Spq18)*CBufK(12,KK,J))
        t12=-((Spq17+Spq19)*CBufK(13,KK,J))-(Spq18+Spq20)*CBufK(14,KK
     $   ,J)
        t13=(Cpq3+Cpq5)*SBufK(2,KK,J)
        t14=(Cpq6+Cpq8)*SBufK(4,KK,J)+(Cpq7+Cpq9)*SBufK(5,KK,J)
        t15=(Cpq10+Cpq12)*SBufK(7,KK,J)
        t16=(Cpq11+Cpq13)*SBufK(8,KK,J)+(Cpq12+Cpq14)*SBufK(9,KK,J)
        t17=(Cpq15+Cpq17)*SBufK(11,KK,J)+(Cpq16+Cpq18)*SBufK(12,KK,J)
        t18=(Cpq17+Cpq19)*SBufK(13,KK,J)+(Cpq18+Cpq20)*SBufK(14,KK,J)
        t19=t11+t12+t13+t14
        t20=t15+t16+t17+t18
        STwo2=t9+t10+t19+t20
        t1=Cpq7*CBufK(2,KK,J)+Cpq11*CBufK(4,KK,J)
        t2=Cpq12*CBufK(5,KK,J)+Cpq16*CBufK(7,KK,J)+Cpq17*CBufK(8,KK,J
     $   )
        t3=Cpq18*CBufK(9,KK,J)+Cpq22*CBufK(11,KK,J)
        t4=Cpq23*CBufK(12,KK,J)+Cpq24*CBufK(13,KK,J)+Cpq25*CBufK(14,K
     $   K,J)
        t5=Spq7*SBufK(2,KK,J)+Spq11*SBufK(4,KK,J)
        t6=Spq12*SBufK(5,KK,J)+Spq16*SBufK(7,KK,J)+Spq17*SBufK(8,KK,J
     $   )
        t7=Spq18*SBufK(9,KK,J)+Spq22*SBufK(11,KK,J)
        t8=Spq23*SBufK(12,KK,J)+Spq24*SBufK(13,KK,J)+Spq25*SBufK(14,K
     $   K,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        CTwo3=t9+t10
        t1=Cpq4*CBufK(0,KK,J)
        t2=Cpq7*CBufK(1,KK,J)+(-Cpq6+Cpq8)*CBufK(2,KK,J)
        t3=Cpq11*CBufK(3,KK,J)
        t4=(-Cpq10+Cpq12)*CBufK(4,KK,J)+(-Cpq11+Cpq13)*CBufK(5,KK,J)
        t5=Cpq16*CBufK(6,KK,J)
        t6=(-Cpq15+Cpq17)*CBufK(7,KK,J)+(-Cpq16+Cpq18)*CBufK(8,KK,J)
        t7=(-Cpq17+Cpq19)*CBufK(9,KK,J)
        t8=Cpq22*CBufK(10,KK,J)+(-Cpq21+Cpq23)*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=(-Cpq22+Cpq24)*CBufK(12,KK,J)
        t12=(-Cpq23+Cpq25)*CBufK(13,KK,J)+(-Cpq24+Cpq26)*CBufK(14,KK,
     $   J)
        t13=Spq8*SBufK(2,KK,J)
        t14=Spq12*SBufK(4,KK,J)+(-Spq11+Spq13)*SBufK(5,KK,J)
        t15=Spq17*SBufK(7,KK,J)
        t16=(-Spq16+Spq18)*SBufK(8,KK,J)+(-Spq17+Spq19)*SBufK(9,KK,J)
        t17=Spq23*SBufK(11,KK,J)+(-Spq22+Spq24)*SBufK(12,KK,J)
        t18=(-Spq23+Spq25)*SBufK(13,KK,J)+(-Spq24+Spq26)*SBufK(14,KK,
     $   J)
        t19=t11+t12+t13+t14
        t20=t15+t16+t17+t18
        CTwo4=t9+t10+t19+t20
        t1=Spq4*CBufK(0,KK,J)+Spq7*CBufK(1,KK,J)+Spq8*CBufK(2,KK,J)
        t2=Spq11*CBufK(3,KK,J)+Spq12*CBufK(4,KK,J)+(Spq11+Spq13)*CBuf
     $   K(5,KK,J)
        t3=Spq16*CBufK(6,KK,J)+Spq17*CBufK(7,KK,J)+(Spq16+Spq18)*CBuf
     $   K(8,KK,J)
        t4=(Spq17+Spq19)*CBufK(9,KK,J)+Spq22*CBufK(10,KK,J)+Spq23*CBu
     $   fK(11,KK,J)
        t5=(Spq22+Spq24)*CBufK(12,KK,J)
        t6=(Spq23+Spq25)*CBufK(13,KK,J)+(Spq24+Spq26)*CBufK(14,KK,J)
        t7=-((Cpq6+Cpq8)*SBufK(2,KK,J))
        t8=-((Cpq10+Cpq12)*SBufK(4,KK,J))-(Cpq11+Cpq13)*SBufK(5,KK,J)
        t9=-((Cpq15+Cpq17)*SBufK(7,KK,J))
        t10=-((Cpq16+Cpq18)*SBufK(8,KK,J))-(Cpq17+Cpq19)*SBufK(9,KK,J
     $   )
        t11=-((Cpq21+Cpq23)*SBufK(11,KK,J))-(Cpq22+Cpq24)*SBufK(12,KK
     $   ,J)
        t12=-((Cpq23+Cpq25)*SBufK(13,KK,J))-(Cpq24+Cpq26)*SBufK(14,KK
     $   ,J)
        t13=t5+t6+t7+t8
        t14=t9+t10+t11+t12
        STwo4=t1+t2+t3+t4+t13+t14
        t1=Cpq5*CBufK(0,KK,J)
        t2=Cpq8*CBufK(1,KK,J)+(-Cpq7+Cpq9)*CBufK(2,KK,J)
        t3=Cpq12*CBufK(3,KK,J)
        t4=(-Cpq11+Cpq13)*CBufK(4,KK,J)+(Cpq10+Cpq14)*CBufK(5,KK,J)
        t5=Cpq17*CBufK(6,KK,J)
        t6=(-Cpq16+Cpq18)*CBufK(7,KK,J)+(Cpq15+Cpq19)*CBufK(8,KK,J)
        t7=(Cpq16+Cpq20)*CBufK(9,KK,J)
        t8=Cpq23*CBufK(10,KK,J)+(-Cpq22+Cpq24)*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=(Cpq21+Cpq25)*CBufK(12,KK,J)
        t12=(Cpq22+Cpq26)*CBufK(13,KK,J)+(Cpq23+Cpq27)*CBufK(14,KK,J)
        t13=(Spq7+Spq9)*SBufK(2,KK,J)+(Spq11+Spq13)*SBufK(4,KK,J)+Spq
     $   14*SBufK(5,KK,J)
        t14=(Spq16+Spq18)*SBufK(7,KK,J)+Spq19*SBufK(8,KK,J)+(Spq16+Sp
     $   q20)*SBufK(9,KK,J)
        t15=(Spq22+Spq24)*SBufK(11,KK,J)+Spq25*SBufK(12,KK,J)
        t16=(Spq22+Spq26)*SBufK(13,KK,J)+(Spq23+Spq27)*SBufK(14,KK,J)
        t17=t9+t10+t11+t12
        t18=t13+t14+t15+t16
        CTwo5=t17+t18
        t1=Spq5*CBufK(0,KK,J)
        t2=Spq8*CBufK(1,KK,J)+(-Spq7+Spq9)*CBufK(2,KK,J)
        t3=Spq12*CBufK(3,KK,J)
        t4=(-Spq11+Spq13)*CBufK(4,KK,J)+Spq14*CBufK(5,KK,J)
        t5=Spq17*CBufK(6,KK,J)
        t6=(-Spq16+Spq18)*CBufK(7,KK,J)+Spq19*CBufK(8,KK,J)
        t7=(-Spq16+Spq20)*CBufK(9,KK,J)
        t8=Spq23*CBufK(10,KK,J)+(-Spq22+Spq24)*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=Spq25*CBufK(12,KK,J)
        t12=(-Spq22+Spq26)*CBufK(13,KK,J)+(-Spq23+Spq27)*CBufK(14,KK,
     $   J)
        t13=-((Cpq7+Cpq9)*SBufK(2,KK,J))
        t14=-((Cpq11+Cpq13)*SBufK(4,KK,J))+(Cpq10-Cpq14)*SBufK(5,KK,J
     $   )
        t15=-((Cpq16+Cpq18)*SBufK(7,KK,J))
        t16=(Cpq15-Cpq19)*SBufK(8,KK,J)+(Cpq16-Cpq20)*SBufK(9,KK,J)
        t17=-((Cpq22+Cpq24)*SBufK(11,KK,J))+(Cpq21-Cpq25)*SBufK(12,KK
     $   ,J)
        t18=(Cpq22-Cpq26)*SBufK(13,KK,J)+(Cpq23-Cpq27)*SBufK(14,KK,J)
        t19=t11+t12+t13+t14
        t20=t15+t16+t17+t18
        STwo5=t9+t10+t19+t20
        t1=-(Cpq11*CBufK(2,KK,J))-Cpq16*CBufK(4,KK,J)
        t2=-(Cpq17*CBufK(5,KK,J))
        t3=-(Cpq22*CBufK(7,KK,J))-Cpq23*CBufK(8,KK,J)
        t4=-(Cpq24*CBufK(9,KK,J))-Cpq29*CBufK(11,KK,J)
        t5=-(Cpq30*CBufK(12,KK,J))
        t6=-(Cpq31*CBufK(13,KK,J))-Cpq32*CBufK(14,KK,J)
        t7=-(Spq11*SBufK(2,KK,J))-Spq16*SBufK(4,KK,J)
        t8=-(Spq17*SBufK(5,KK,J))
        t9=-(Spq22*SBufK(7,KK,J))-Spq23*SBufK(8,KK,J)
        t10=-(Spq24*SBufK(9,KK,J))-Spq29*SBufK(11,KK,J)
        t11=-(Spq30*SBufK(12,KK,J))
        t12=-(Spq31*SBufK(13,KK,J))-Spq32*SBufK(14,KK,J)
        t13=t1+t2+t3+t4+t5+t6
        t14=t7+t8+t9+t10+t11+t12
        CTwo6=t13+t14
        t1=-(Cpq7*CBufK(0,KK,J))
        t2=-(Cpq11*CBufK(1,KK,J))+(Cpq10-Cpq12)*CBufK(2,KK,J)
        t3=-(Cpq16*CBufK(3,KK,J))
        t4=(Cpq15-Cpq17)*CBufK(4,KK,J)+(Cpq16-Cpq18)*CBufK(5,KK,J)
        t5=-(Cpq22*CBufK(6,KK,J))
        t6=(Cpq21-Cpq23)*CBufK(7,KK,J)+(Cpq22-Cpq24)*CBufK(8,KK,J)
        t7=(Cpq23-Cpq25)*CBufK(9,KK,J)
        t8=-(Cpq29*CBufK(10,KK,J))+(Cpq28-Cpq30)*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=(Cpq29-Cpq31)*CBufK(12,KK,J)
        t12=(Cpq30-Cpq32)*CBufK(13,KK,J)+(Cpq31-Cpq33)*CBufK(14,KK,J)
        t13=-(Spq12*SBufK(2,KK,J))
        t14=-(Spq17*SBufK(4,KK,J))+(Spq16-Spq18)*SBufK(5,KK,J)
        t15=-(Spq23*SBufK(7,KK,J))
        t16=(Spq22-Spq24)*SBufK(8,KK,J)+(Spq23-Spq25)*SBufK(9,KK,J)
        t17=-(Spq30*SBufK(11,KK,J))+(Spq29-Spq31)*SBufK(12,KK,J)
        t18=(Spq30-Spq32)*SBufK(13,KK,J)+(Spq31-Spq33)*SBufK(14,KK,J)
        t19=t11+t12+t13+t14
        t20=t15+t16+t17+t18
        CTwo7=t9+t10+t19+t20
        t1=-(Spq7*CBufK(0,KK,J))
        t2=-(Spq11*CBufK(1,KK,J))-Spq12*CBufK(2,KK,J)
        t3=-(Spq16*CBufK(3,KK,J))
        t4=-(Spq17*CBufK(4,KK,J))-(Spq16+Spq18)*CBufK(5,KK,J)
        t5=-(Spq22*CBufK(6,KK,J))
        t6=-(Spq23*CBufK(7,KK,J))-(Spq22+Spq24)*CBufK(8,KK,J)
        t7=-((Spq23+Spq25)*CBufK(9,KK,J))
        t8=-(Spq29*CBufK(10,KK,J))-Spq30*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=-((Spq29+Spq31)*CBufK(12,KK,J))
        t12=-((Spq30+Spq32)*CBufK(13,KK,J))-(Spq31+Spq33)*CBufK(14,KK
     $   ,J)
        t13=(Cpq10+Cpq12)*SBufK(2,KK,J)
        t14=(Cpq15+Cpq17)*SBufK(4,KK,J)+(Cpq16+Cpq18)*SBufK(5,KK,J)
        t15=(Cpq21+Cpq23)*SBufK(7,KK,J)
        t16=(Cpq22+Cpq24)*SBufK(8,KK,J)+(Cpq23+Cpq25)*SBufK(9,KK,J)
        t17=(Cpq28+Cpq30)*SBufK(11,KK,J)+(Cpq29+Cpq31)*SBufK(12,KK,J)
        t18=(Cpq30+Cpq32)*SBufK(13,KK,J)+(Cpq31+Cpq33)*SBufK(14,KK,J)
        t19=t11+t12+t13+t14
        t20=t15+t16+t17+t18
        STwo7=t9+t10+t19+t20
        t1=-(Cpq8*CBufK(0,KK,J))-Cpq12*CBufK(1,KK,J)
        t2=Cpq11*CBufK(2,KK,J)-Cpq13*CBufK(2,KK,J)-Cpq17*CBufK(3,KK,J
     $   )
        t3=Cpq16*CBufK(4,KK,J)-Cpq18*CBufK(4,KK,J)
        t4=-(Cpq15*CBufK(5,KK,J))
        t5=-(Cpq19*CBufK(5,KK,J))-Cpq23*CBufK(6,KK,J)
        t6=Cpq22*CBufK(7,KK,J)-Cpq24*CBufK(7,KK,J)
        t7=-(Cpq21*CBufK(8,KK,J))
        t8=-(Cpq25*CBufK(8,KK,J))-Cpq22*CBufK(9,KK,J)
        t9=-(Cpq26*CBufK(9,KK,J))-Cpq30*CBufK(10,KK,J)+Cpq29*CBufK(11
     $   ,KK,J)
        t10=-(Cpq31*CBufK(11,KK,J))
        t11=-(Cpq28*CBufK(12,KK,J))-Cpq32*CBufK(12,KK,J)
        t12=t1+t2+t3+t4+t5
        t13=t6+t7+t8+t9+t10+t11
        t14=-(Cpq29*CBufK(13,KK,J))-Cpq33*CBufK(13,KK,J)
        t15=-(Cpq30*CBufK(14,KK,J))
        t16=-(Cpq34*CBufK(14,KK,J))-Spq11*SBufK(2,KK,J)
        t17=-(Spq13*SBufK(2,KK,J))-Spq16*SBufK(4,KK,J)
        t18=-(Spq18*SBufK(4,KK,J))
        t19=-(Spq19*SBufK(5,KK,J))-Spq22*SBufK(7,KK,J)
        t20=-(Spq24*SBufK(7,KK,J))-Spq25*SBufK(8,KK,J)
        t21=-(Spq22*SBufK(9,KK,J))
        t22=-(Spq26*SBufK(9,KK,J))-Spq29*SBufK(11,KK,J)
        t23=-(Spq31*SBufK(11,KK,J))
        t24=-(Spq32*SBufK(12,KK,J))-Spq29*SBufK(13,KK,J)
        t25=-(Spq33*SBufK(13,KK,J))
        t26=-(Spq30*SBufK(14,KK,J))-Spq34*SBufK(14,KK,J)
        t27=t20+t21+t22
        t28=t23+t24+t25+t26
        t29=t14+t15+t16+t17
        t30=t18+t19+t27+t28
        CTwo8=t12+t13+t29+t30
        t1=-(Spq8*CBufK(0,KK,J))
        t2=-(Spq12*CBufK(1,KK,J))+(Spq11-Spq13)*CBufK(2,KK,J)
        t3=-(Spq17*CBufK(3,KK,J))
        t4=(Spq16-Spq18)*CBufK(4,KK,J)-Spq19*CBufK(5,KK,J)
        t5=-(Spq23*CBufK(6,KK,J))
        t6=(Spq22-Spq24)*CBufK(7,KK,J)-Spq25*CBufK(8,KK,J)
        t7=(Spq22-Spq26)*CBufK(9,KK,J)
        t8=-(Spq30*CBufK(10,KK,J))+(Spq29-Spq31)*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=-(Spq32*CBufK(12,KK,J))
        t12=(Spq29-Spq33)*CBufK(13,KK,J)+(Spq30-Spq34)*CBufK(14,KK,J)
        t13=(Cpq11+Cpq13)*SBufK(2,KK,J)
        t14=(Cpq16+Cpq18)*SBufK(4,KK,J)+(-Cpq15+Cpq19)*SBufK(5,KK,J)
        t15=(Cpq22+Cpq24)*SBufK(7,KK,J)
        t16=(-Cpq21+Cpq25)*SBufK(8,KK,J)+(-Cpq22+Cpq26)*SBufK(9,KK,J)
        t17=(Cpq29+Cpq31)*SBufK(11,KK,J)+(-Cpq28+Cpq32)*SBufK(12,KK,J
     $   )
        t18=(-Cpq29+Cpq33)*SBufK(13,KK,J)+(-Cpq30+Cpq34)*SBufK(14,KK,
     $   J)
        t19=t11+t12+t13+t14
        t20=t15+t16+t17+t18
        STwo8=t9+t10+t19+t20
        t1=-(Cpq9*CBufK(0,KK,J))
        t2=-(Cpq13*CBufK(1,KK,J))+(Cpq12-Cpq14)*CBufK(2,KK,J)
        t3=-(Cpq18*CBufK(3,KK,J))
        t4=(Cpq17-Cpq19)*CBufK(4,KK,J)-(Cpq16+Cpq20)*CBufK(5,KK,J)
        t5=-(Cpq24*CBufK(6,KK,J))
        t6=(Cpq23-Cpq25)*CBufK(7,KK,J)-(Cpq22+Cpq26)*CBufK(8,KK,J)
        t7=(Cpq21-Cpq27)*CBufK(9,KK,J)
        t8=-(Cpq31*CBufK(10,KK,J))+(Cpq30-Cpq32)*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=-((Cpq29+Cpq33)*CBufK(12,KK,J))
        t12=(Cpq28-Cpq34)*CBufK(13,KK,J)+(Cpq29-Cpq35)*CBufK(14,KK,J)
        t13=-((Spq12+Spq14)*SBufK(2,KK,J))
        t14=-((Spq17+Spq19)*SBufK(4,KK,J))+(Spq16-Spq20)*SBufK(5,KK,J
     $   )
        t15=-((Spq23+Spq25)*SBufK(7,KK,J))
        t16=(Spq22-Spq26)*SBufK(8,KK,J)-Spq27*SBufK(9,KK,J)
        t17=-((Spq30+Spq32)*SBufK(11,KK,J))+(Spq29-Spq33)*SBufK(12,KK
     $   ,J)
        t18=-(Spq34*SBufK(13,KK,J))+(Spq29-Spq35)*SBufK(14,KK,J)
        t19=t11+t12+t13+t14
        t20=t15+t16+t17+t18
        CTwo9=t9+t10+t19+t20
        t1=-(Spq9*CBufK(0,KK,J))
        t2=-(Spq13*CBufK(1,KK,J))+(Spq12-Spq14)*CBufK(2,KK,J)
        t3=-(Spq18*CBufK(3,KK,J))
        t4=(Spq17-Spq19)*CBufK(4,KK,J)-(Spq16+Spq20)*CBufK(5,KK,J)
        t5=-(Spq24*CBufK(6,KK,J))
        t6=(Spq23-Spq25)*CBufK(7,KK,J)-(Spq22+Spq26)*CBufK(8,KK,J)
        t7=-(Spq27*CBufK(9,KK,J))
        t8=-(Spq31*CBufK(10,KK,J))+(Spq30-Spq32)*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=-((Spq29+Spq33)*CBufK(12,KK,J))
        t12=-(Spq34*CBufK(13,KK,J))-(Spq29+Spq35)*CBufK(14,KK,J)
        t13=(Cpq12+Cpq14)*SBufK(2,KK,J)
        t14=(Cpq17+Cpq19)*SBufK(4,KK,J)+(-Cpq16+Cpq20)*SBufK(5,KK,J)
        t15=(Cpq23+Cpq25)*SBufK(7,KK,J)
        t16=(-Cpq22+Cpq26)*SBufK(8,KK,J)+(Cpq21+Cpq27)*SBufK(9,KK,J)
        t17=(Cpq30+Cpq32)*SBufK(11,KK,J)+(-Cpq29+Cpq33)*SBufK(12,KK,J
     $   )
        t18=(Cpq28+Cpq34)*SBufK(13,KK,J)+(Cpq29+Cpq35)*SBufK(14,KK,J)
        t19=t11+t12+t13+t14
        t20=t15+t16+t17+t18
        STwo9=t9+t10+t19+t20
        t1=Cpq16*CBufK(2,KK,J)+Cpq22*CBufK(4,KK,J)
        t2=Cpq23*CBufK(5,KK,J)+Cpq29*CBufK(7,KK,J)+Cpq30*CBufK(8,KK,J
     $   )
        t3=Cpq31*CBufK(9,KK,J)+Cpq37*CBufK(11,KK,J)
        t4=Cpq38*CBufK(12,KK,J)+Cpq39*CBufK(13,KK,J)+Cpq40*CBufK(14,K
     $   K,J)
        t5=Spq16*SBufK(2,KK,J)+Spq22*SBufK(4,KK,J)
        t6=Spq23*SBufK(5,KK,J)+Spq29*SBufK(7,KK,J)+Spq30*SBufK(8,KK,J
     $   )
        t7=Spq31*SBufK(9,KK,J)+Spq37*SBufK(11,KK,J)
        t8=Spq38*SBufK(12,KK,J)+Spq39*SBufK(13,KK,J)+Spq40*SBufK(14,K
     $   K,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        CTwo10=t9+t10
        t1=Cpq11*CBufK(0,KK,J)
        t2=Cpq16*CBufK(1,KK,J)+(-Cpq15+Cpq17)*CBufK(2,KK,J)
        t3=Cpq22*CBufK(3,KK,J)
        t4=(-Cpq21+Cpq23)*CBufK(4,KK,J)+(-Cpq22+Cpq24)*CBufK(5,KK,J)
        t5=Cpq29*CBufK(6,KK,J)
        t6=(-Cpq28+Cpq30)*CBufK(7,KK,J)+(-Cpq29+Cpq31)*CBufK(8,KK,J)
        t7=(-Cpq30+Cpq32)*CBufK(9,KK,J)
        t8=Cpq37*CBufK(10,KK,J)+(-Cpq36+Cpq38)*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=(-Cpq37+Cpq39)*CBufK(12,KK,J)
        t12=(-Cpq38+Cpq40)*CBufK(13,KK,J)+(-Cpq39+Cpq41)*CBufK(14,KK,
     $   J)
        t13=Spq17*SBufK(2,KK,J)
        t14=Spq23*SBufK(4,KK,J)+(-Spq22+Spq24)*SBufK(5,KK,J)
        t15=Spq30*SBufK(7,KK,J)
        t16=(-Spq29+Spq31)*SBufK(8,KK,J)+(-Spq30+Spq32)*SBufK(9,KK,J)
        t17=Spq38*SBufK(11,KK,J)+(-Spq37+Spq39)*SBufK(12,KK,J)
        t18=(-Spq38+Spq40)*SBufK(13,KK,J)+(-Spq39+Spq41)*SBufK(14,KK,
     $   J)
        t19=t11+t12+t13+t14
        t20=t15+t16+t17+t18
        CTwo11=t9+t10+t19+t20
        t1=Spq11*CBufK(0,KK,J)+Spq16*CBufK(1,KK,J)+Spq17*CBufK(2,KK,J
     $   )
        t2=Spq22*CBufK(3,KK,J)+Spq23*CBufK(4,KK,J)+(Spq22+Spq24)*CBuf
     $   K(5,KK,J)
        t3=Spq29*CBufK(6,KK,J)+Spq30*CBufK(7,KK,J)+(Spq29+Spq31)*CBuf
     $   K(8,KK,J)
        t4=(Spq30+Spq32)*CBufK(9,KK,J)+Spq37*CBufK(10,KK,J)+Spq38*CBu
     $   fK(11,KK,J)
        t5=(Spq37+Spq39)*CBufK(12,KK,J)
        t6=(Spq38+Spq40)*CBufK(13,KK,J)+(Spq39+Spq41)*CBufK(14,KK,J)
        t7=-((Cpq15+Cpq17)*SBufK(2,KK,J))
        t8=-((Cpq21+Cpq23)*SBufK(4,KK,J))-(Cpq22+Cpq24)*SBufK(5,KK,J)
        t9=-((Cpq28+Cpq30)*SBufK(7,KK,J))
        t10=-((Cpq29+Cpq31)*SBufK(8,KK,J))-(Cpq30+Cpq32)*SBufK(9,KK,J
     $   )
        t11=-((Cpq36+Cpq38)*SBufK(11,KK,J))-(Cpq37+Cpq39)*SBufK(12,KK
     $   ,J)
        t12=-((Cpq38+Cpq40)*SBufK(13,KK,J))-(Cpq39+Cpq41)*SBufK(14,KK
     $   ,J)
        t13=t5+t6+t7+t8
        t14=t9+t10+t11+t12
        STwo11=t1+t2+t3+t4+t13+t14
        t1=Cpq12*CBufK(0,KK,J)
        t2=Cpq17*CBufK(1,KK,J)+(-Cpq16+Cpq18)*CBufK(2,KK,J)
        t3=Cpq23*CBufK(3,KK,J)
        t4=(-Cpq22+Cpq24)*CBufK(4,KK,J)+(Cpq21+Cpq25)*CBufK(5,KK,J)
        t5=Cpq30*CBufK(6,KK,J)
        t6=(-Cpq29+Cpq31)*CBufK(7,KK,J)+(Cpq28+Cpq32)*CBufK(8,KK,J)
        t7=(Cpq29+Cpq33)*CBufK(9,KK,J)
        t8=Cpq38*CBufK(10,KK,J)+(-Cpq37+Cpq39)*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=(Cpq36+Cpq40)*CBufK(12,KK,J)
        t12=(Cpq37+Cpq41)*CBufK(13,KK,J)+(Cpq38+Cpq42)*CBufK(14,KK,J)
        t13=(Spq16+Spq18)*SBufK(2,KK,J)+(Spq22+Spq24)*SBufK(4,KK,J)+S
     $   pq25*SBufK(5,KK,J)
        t14=(Spq29+Spq31)*SBufK(7,KK,J)+Spq32*SBufK(8,KK,J)+(Spq29+Sp
     $   q33)*SBufK(9,KK,J)
        t15=(Spq37+Spq39)*SBufK(11,KK,J)+Spq40*SBufK(12,KK,J)
        t16=(Spq37+Spq41)*SBufK(13,KK,J)+(Spq38+Spq42)*SBufK(14,KK,J)
        t17=t9+t10+t11+t12
        t18=t13+t14+t15+t16
        CTwo12=t17+t18
        t1=Spq12*CBufK(0,KK,J)
        t2=Spq17*CBufK(1,KK,J)+(-Spq16+Spq18)*CBufK(2,KK,J)
        t3=Spq23*CBufK(3,KK,J)
        t4=(-Spq22+Spq24)*CBufK(4,KK,J)+Spq25*CBufK(5,KK,J)
        t5=Spq30*CBufK(6,KK,J)
        t6=(-Spq29+Spq31)*CBufK(7,KK,J)+Spq32*CBufK(8,KK,J)
        t7=(-Spq29+Spq33)*CBufK(9,KK,J)
        t8=Spq38*CBufK(10,KK,J)+(-Spq37+Spq39)*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=Spq40*CBufK(12,KK,J)
        t12=(-Spq37+Spq41)*CBufK(13,KK,J)+(-Spq38+Spq42)*CBufK(14,KK,
     $   J)
        t13=-((Cpq16+Cpq18)*SBufK(2,KK,J))
        t14=-((Cpq22+Cpq24)*SBufK(4,KK,J))+(Cpq21-Cpq25)*SBufK(5,KK,J
     $   )
        t15=-((Cpq29+Cpq31)*SBufK(7,KK,J))
        t16=(Cpq28-Cpq32)*SBufK(8,KK,J)+(Cpq29-Cpq33)*SBufK(9,KK,J)
        t17=-((Cpq37+Cpq39)*SBufK(11,KK,J))+(Cpq36-Cpq40)*SBufK(12,KK
     $   ,J)
        t18=(Cpq37-Cpq41)*SBufK(13,KK,J)+(Cpq38-Cpq42)*SBufK(14,KK,J)
        t19=t11+t12+t13+t14
        t20=t15+t16+t17+t18
        STwo12=t9+t10+t19+t20
        t1=Cpq13*CBufK(0,KK,J)
        t2=Cpq18*CBufK(1,KK,J)+(-Cpq17+Cpq19)*CBufK(2,KK,J)
        t3=Cpq24*CBufK(3,KK,J)
        t4=(-Cpq23+Cpq25)*CBufK(4,KK,J)+(Cpq22+Cpq26)*CBufK(5,KK,J)
        t5=Cpq31*CBufK(6,KK,J)
        t6=(-Cpq30+Cpq32)*CBufK(7,KK,J)+(Cpq29+Cpq33)*CBufK(8,KK,J)
        t7=(-Cpq28+Cpq34)*CBufK(9,KK,J)
        t8=Cpq39*CBufK(10,KK,J)+(-Cpq38+Cpq40)*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=(Cpq37+Cpq41)*CBufK(12,KK,J)
        t12=(-Cpq36+Cpq42)*CBufK(13,KK,J)+(-Cpq37+Cpq43)*CBufK(14,KK,
     $   J)
        t13=(Spq17+Spq19)*SBufK(2,KK,J)
        t14=(Spq23+Spq25)*SBufK(4,KK,J)+(-Spq22+Spq26)*SBufK(5,KK,J)
        t15=(Spq30+Spq32)*SBufK(7,KK,J)
        t16=(-Spq29+Spq33)*SBufK(8,KK,J)+Spq34*SBufK(9,KK,J)
        t17=(Spq38+Spq40)*SBufK(11,KK,J)+(-Spq37+Spq41)*SBufK(12,KK,J
     $   )
        t18=Spq42*SBufK(13,KK,J)+(-Spq37+Spq43)*SBufK(14,KK,J)
        t19=t11+t12+t13+t14
        t20=t15+t16+t17+t18
        CTwo13=t9+t10+t19+t20
        t1=Spq13*CBufK(0,KK,J)
        t2=Spq18*CBufK(1,KK,J)+(-Spq17+Spq19)*CBufK(2,KK,J)
        t3=Spq24*CBufK(3,KK,J)
        t4=(-Spq23+Spq25)*CBufK(4,KK,J)+(Spq22+Spq26)*CBufK(5,KK,J)
        t5=Spq31*CBufK(6,KK,J)
        t6=(-Spq30+Spq32)*CBufK(7,KK,J)+(Spq29+Spq33)*CBufK(8,KK,J)
        t7=Spq34*CBufK(9,KK,J)
        t8=Spq39*CBufK(10,KK,J)+(-Spq38+Spq40)*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=(Spq37+Spq41)*CBufK(12,KK,J)+Spq42*CBufK(13,KK,J)+(Spq37+
     $   Spq43)*CBufK(14,KK,J)
        t12=-((Cpq17+Cpq19)*SBufK(2,KK,J))
        t13=-((Cpq23+Cpq25)*SBufK(4,KK,J))+(Cpq22-Cpq26)*SBufK(5,KK,J
     $   )
        t14=-((Cpq30+Cpq32)*SBufK(7,KK,J))
        t15=(Cpq29-Cpq33)*SBufK(8,KK,J)-(Cpq28+Cpq34)*SBufK(9,KK,J)
        t16=-((Cpq38+Cpq40)*SBufK(11,KK,J))+(Cpq37-Cpq41)*SBufK(12,KK
     $   ,J)
        t17=-((Cpq36+Cpq42)*SBufK(13,KK,J))-(Cpq37+Cpq43)*SBufK(14,KK
     $   ,J)
        t18=t11+t12+t13
        t19=t14+t15+t16+t17
        STwo13=t9+t10+t18+t19
        t1=Cpq14*CBufK(0,KK,J)
        t2=Cpq19*CBufK(1,KK,J)+(-Cpq18+Cpq20)*CBufK(2,KK,J)
        t3=Cpq25*CBufK(3,KK,J)
        t4=(-Cpq24+Cpq26)*CBufK(4,KK,J)+(Cpq23+Cpq27)*CBufK(5,KK,J)
        t5=Cpq32*CBufK(6,KK,J)
        t6=(-Cpq31+Cpq33)*CBufK(7,KK,J)+(Cpq30+Cpq34)*CBufK(8,KK,J)
        t7=(-Cpq29+Cpq35)*CBufK(9,KK,J)
        t8=Cpq40*CBufK(10,KK,J)+(-Cpq39+Cpq41)*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=(Cpq38+Cpq42)*CBufK(12,KK,J)
        t12=(-Cpq37+Cpq43)*CBufK(13,KK,J)+(Cpq36+Cpq44)*CBufK(14,KK,J
     $   )
        t13=(Spq18+Spq20)*SBufK(2,KK,J)
        t14=(Spq24+Spq26)*SBufK(4,KK,J)+(-Spq23+Spq27)*SBufK(5,KK,J)
        t15=(Spq31+Spq33)*SBufK(7,KK,J)
        t16=(-Spq30+Spq34)*SBufK(8,KK,J)+(Spq29+Spq35)*SBufK(9,KK,J)
        t17=(Spq39+Spq41)*SBufK(11,KK,J)+(-Spq38+Spq42)*SBufK(12,KK,J
     $   )
        t18=(Spq37+Spq43)*SBufK(13,KK,J)+Spq44*SBufK(14,KK,J)
        t19=t11+t12+t13+t14
        t20=t15+t16+t17+t18
        CTwo14=t9+t10+t19+t20
        t1=Spq14*CBufK(0,KK,J)
        t2=Spq19*CBufK(1,KK,J)+(-Spq18+Spq20)*CBufK(2,KK,J)
        t3=Spq25*CBufK(3,KK,J)
        t4=(-Spq24+Spq26)*CBufK(4,KK,J)+(Spq23+Spq27)*CBufK(5,KK,J)
        t5=Spq32*CBufK(6,KK,J)
        t6=(-Spq31+Spq33)*CBufK(7,KK,J)+(Spq30+Spq34)*CBufK(8,KK,J)
        t7=(-Spq29+Spq35)*CBufK(9,KK,J)
        t8=Spq40*CBufK(10,KK,J)+(-Spq39+Spq41)*CBufK(11,KK,J)
        t9=t1+t2+t3+t4
        t10=t5+t6+t7+t8
        t11=(Spq38+Spq42)*CBufK(12,KK,J)
        t12=(-Spq37+Spq43)*CBufK(13,KK,J)+Spq44*CBufK(14,KK,J)
        t13=-((Cpq18+Cpq20)*SBufK(2,KK,J))
        t14=-((Cpq24+Cpq26)*SBufK(4,KK,J))+(Cpq23-Cpq27)*SBufK(5,KK,J
     $   )
        t15=-((Cpq31+Cpq33)*SBufK(7,KK,J))
        t16=(Cpq30-Cpq34)*SBufK(8,KK,J)-(Cpq29+Cpq35)*SBufK(9,KK,J)
        t17=-((Cpq39+Cpq41)*SBufK(11,KK,J))+(Cpq38-Cpq42)*SBufK(12,KK
     $   ,J)
        t18=-((Cpq37+Cpq43)*SBufK(13,KK,J))+(Cpq36-Cpq44)*SBufK(14,KK
     $   ,J)
        t19=t11+t12+t13+t14
        t20=t15+t16+t17+t18
        STwo14=t9+t10+t19+t20
c       ------------------------------ 
        t1=Cpq0*CBufK(0,KK,J)+Cpq1*CBufK(1,KK,J)
        t2=Cpq3*CBufK(3,KK,J)+Cpq6*CBufK(6,KK,J)+Cpq10*CBufK(10,KK,J)
        Cone0=t1+t2
        t1=-(Cpq1*CBufK(0,KK,J))-Cpq3*CBufK(1,KK,J)
        t2=-(Cpq6*CBufK(3,KK,J))
        t3=-(Cpq10*CBufK(6,KK,J))-Cpq15*CBufK(10,KK,J)
        Cone1=t1+t2+t3
        t1=Cpq3*CBufK(0,KK,J)+Cpq6*CBufK(1,KK,J)
        t2=Cpq10*CBufK(3,KK,J)+Cpq15*CBufK(6,KK,J)+Cpq21*CBufK(10,KK,
     >   J)
        Cone3=t1+t2
        t1=-(Cpq6*CBufK(0,KK,J))-Cpq10*CBufK(1,KK,J)
        t2=-(Cpq15*CBufK(3,KK,J))
        t3=-(Cpq21*CBufK(6,KK,J))-Cpq28*CBufK(10,KK,J)
        Cone6=t1+t2+t3
        t1=Cpq10*CBufK(0,KK,J)+Cpq15*CBufK(1,KK,J)
        t2=Cpq21*CBufK(3,KK,J)+Cpq28*CBufK(6,KK,J)+Cpq36*CBufK(10,KK,
     >   J)
        Cone10=t1+t2
c       ------------------------------ 
        CKet(0)=COne0+2.D0*CTwo0
        SKet(0)=Zero
        CKet(1)=COne1+2.D0*CTwo1
        SKet(1)=Zero
        CKet(2)=2.D0*CTwo2
        SKet(2)=2.D0*STwo2
        CKet(3)=COne3+2.D0*CTwo3
        SKet(3)=Zero
        CKet(4)=2.D0*CTwo4
        SKet(4)=2.D0*STwo4
        CKet(5)=2.D0*CTwo5
        SKet(5)=2.D0*STwo5
        CKet(6)=COne6+2.D0*CTwo6
        SKet(6)=Zero
        CKet(7)=2.D0*CTwo7
        SKet(7)=2.D0*STwo7
        CKet(8)=2.D0*CTwo8
        SKet(8)=2.D0*STwo8
        CKet(9)=2.D0*CTwo9
        SKet(9)=2.D0*STwo9
        CKet(10)=COne10+2.D0*CTwo10
        SKet(10)=Zero
        CKet(11)=2.D0*CTwo11
        SKet(11)=2.D0*STwo11
        CKet(12)=2.D0*CTwo12
        SKet(12)=2.D0*STwo12
        CKet(13)=2.D0*CTwo13
        SKet(13)=2.D0*STwo13
        CKet(14)=2.D0*CTwo14
        SKet(14)=2.D0*STwo14
c       ------------------------------ 
         
