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
c       ------------------------------ 
        Cone0=Cpq0*CBufK(0,KK,J)+Cpq1*CBufK(1,KK,J)
        Cone1=-(Cpq1*CBufK(0,KK,J))-Cpq3*CBufK(1,KK,J)
        Cone3=Cpq3*CBufK(0,KK,J)+Cpq6*CBufK(1,KK,J)
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
c       ------------------------------ 
         
