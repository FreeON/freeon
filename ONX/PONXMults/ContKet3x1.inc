c       ------------------------------ 
        CTwo0=Cpq2*CBufK(2,KK,J)+Spq2*SBufK(2,KK,J)
        CTwo1=-(Cpq4*CBufK(2,KK,J))-Spq4*SBufK(2,KK,J)
        t1=-(Cpq2*CBufK(0,KK,J))-Cpq4*CBufK(1,KK,J)
        t2=Cpq3*CBufK(2,KK,J)-Cpq5*CBufK(2,KK,J)-Spq5*SBufK(2,KK,J)
        CTwo2=t1+t2
        t1=-(Spq2*CBufK(0,KK,J))-Spq4*CBufK(1,KK,J)
        t2=-(Spq5*CBufK(2,KK,J))+(Cpq3+Cpq5)*SBufK(2,KK,J)
        STwo2=t1+t2
        CTwo3=Cpq7*CBufK(2,KK,J)+Spq7*SBufK(2,KK,J)
        t1=Cpq4*CBufK(0,KK,J)+Cpq7*CBufK(1,KK,J)
        t2=(-Cpq6+Cpq8)*CBufK(2,KK,J)+Spq8*SBufK(2,KK,J)
        CTwo4=t1+t2
        t1=Spq4*CBufK(0,KK,J)+Spq7*CBufK(1,KK,J)
        t2=Spq8*CBufK(2,KK,J)-(Cpq6+Cpq8)*SBufK(2,KK,J)
        STwo4=t1+t2
        t1=Cpq5*CBufK(0,KK,J)+Cpq8*CBufK(1,KK,J)
        t2=(-Cpq7+Cpq9)*CBufK(2,KK,J)+(Spq7+Spq9)*SBufK(2,KK,J)
        CTwo5=t1+t2
        t1=Spq5*CBufK(0,KK,J)+Spq8*CBufK(1,KK,J)
        t2=(-Spq7+Spq9)*CBufK(2,KK,J)-(Cpq7+Cpq9)*SBufK(2,KK,J)
        STwo5=t1+t2
        CTwo6=-(Cpq11*CBufK(2,KK,J))-Spq11*SBufK(2,KK,J)
        t1=-(Cpq7*CBufK(0,KK,J))-Cpq11*CBufK(1,KK,J)
        t2=Cpq10*CBufK(2,KK,J)-Cpq12*CBufK(2,KK,J)-Spq12*SBufK(2,KK,J
     $   )
        CTwo7=t1+t2
        t1=-(Spq7*CBufK(0,KK,J))-Spq11*CBufK(1,KK,J)
        t2=-(Spq12*CBufK(2,KK,J))+(Cpq10+Cpq12)*SBufK(2,KK,J)
        STwo7=t1+t2
        t1=-(Cpq8*CBufK(0,KK,J))-Cpq12*CBufK(1,KK,J)+Cpq11*CBufK(2,KK
     $   ,J)
        t2=-(Cpq13*CBufK(2,KK,J))
        t3=-(Spq11*SBufK(2,KK,J))-Spq13*SBufK(2,KK,J)
        CTwo8=t1+t2+t3
        t1=-(Spq8*CBufK(0,KK,J))-Spq12*CBufK(1,KK,J)
        t2=(Spq11-Spq13)*CBufK(2,KK,J)+(Cpq11+Cpq13)*SBufK(2,KK,J)
        STwo8=t1+t2
        t1=-(Cpq9*CBufK(0,KK,J))-Cpq13*CBufK(1,KK,J)+Cpq12*CBufK(2,KK
     $   ,J)
        t2=-(Cpq14*CBufK(2,KK,J))
        t3=-(Spq12*SBufK(2,KK,J))-Spq14*SBufK(2,KK,J)
        CTwo9=t1+t2+t3
        t1=-(Spq9*CBufK(0,KK,J))-Spq13*CBufK(1,KK,J)
        t2=(Spq12-Spq14)*CBufK(2,KK,J)+(Cpq12+Cpq14)*SBufK(2,KK,J)
        STwo9=t1+t2
c       ------------------------------ 
        Cone0=Cpq0*CBufK(0,KK,J)+Cpq1*CBufK(1,KK,J)
        Cone1=-(Cpq1*CBufK(0,KK,J))-Cpq3*CBufK(1,KK,J)
        Cone3=Cpq3*CBufK(0,KK,J)+Cpq6*CBufK(1,KK,J)
        Cone6=-(Cpq6*CBufK(0,KK,J))-Cpq10*CBufK(1,KK,J)
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
c       ------------------------------ 
         
