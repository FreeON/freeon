c       ------------------------------ 
        CTwo0=Cpq2*CBufK(2,KK,J)+Spq2*SBufK(2,KK,J)
        CTwo1=-(Cpq4*CBufK(2,KK,J))-Spq4*SBufK(2,KK,J)
        t1=-(Cpq2*CBufK(0,KK,J))-Cpq4*CBufK(1,KK,J)
        t2=Cpq3*CBufK(2,KK,J)-Cpq5*CBufK(2,KK,J)-Spq5*SBufK(2,KK,J)
        CTwo2=t1+t2
        t1=-(Spq2*CBufK(0,KK,J))-Spq4*CBufK(1,KK,J)
        t2=-(Spq5*CBufK(2,KK,J))+(Cpq3+Cpq5)*SBufK(2,KK,J)
        STwo2=t1+t2
c       ------------------------------ 
        Cone0=Cpq0*CBufK(0,KK,J)+Cpq1*CBufK(1,KK,J)
        Cone1=-(Cpq1*CBufK(0,KK,J))-Cpq3*CBufK(1,KK,J)
c       ------------------------------ 
        CKet(0)=COne0+2.D0*CTwo0
        SKet(0)=Zero
        CKet(1)=COne1+2.D0*CTwo1
        SKet(1)=Zero
        CKet(2)=2.D0*CTwo2
        SKet(2)=2.D0*STwo2
c       ------------------------------ 
         
