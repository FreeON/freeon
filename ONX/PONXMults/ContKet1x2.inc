c       ------------------------------ 
        t1=Cpq2*CBufK(2,KK,J)+Cpq4*CBufK(4,KK,J)+Cpq5*CBufK(5,KK,J)
        t2=Spq2*SBufK(2,KK,J)+Spq4*SBufK(4,KK,J)+Spq5*SBufK(5,KK,J)
        CTwo0=t1+t2
        t1=-(Cpq4*CBufK(2,KK,J))
        t2=-(Cpq7*CBufK(4,KK,J))-Cpq8*CBufK(5,KK,J)
        t3=-(Spq4*SBufK(2,KK,J))
        t4=-(Spq7*SBufK(4,KK,J))-Spq8*SBufK(5,KK,J)
        CTwo1=t1+t2+t3+t4
        t1=-(Cpq2*CBufK(0,KK,J))-Cpq4*CBufK(1,KK,J)
        t2=(Cpq3-Cpq5)*CBufK(2,KK,J)-Cpq7*CBufK(3,KK,J)
        t3=(Cpq6-Cpq8)*CBufK(4,KK,J)+(Cpq7-Cpq9)*CBufK(5,KK,J)
        t4=-(Spq5*SBufK(2,KK,J))
        t5=-(Spq8*SBufK(4,KK,J))+(Spq7-Spq9)*SBufK(5,KK,J)
        CTwo2=t1+t2+t3+t4+t5
        t1=-(Spq2*CBufK(0,KK,J))-Spq4*CBufK(1,KK,J)
        t2=-(Spq5*CBufK(2,KK,J))-Spq7*CBufK(3,KK,J)
        t3=-(Spq8*CBufK(4,KK,J))-(Spq7+Spq9)*CBufK(5,KK,J)
        t4=(Cpq3+Cpq5)*SBufK(2,KK,J)
        t5=(Cpq6+Cpq8)*SBufK(4,KK,J)+(Cpq7+Cpq9)*SBufK(5,KK,J)
        STwo2=t1+t2+t3+t4+t5
c       ------------------------------ 
        Cone0=Cpq0*CBufK(0,KK,J)+Cpq1*CBufK(1,KK,J)+Cpq3*CBufK(3,KK,J
     >   )
        t1=-(Cpq1*CBufK(0,KK,J))
        t2=-(Cpq3*CBufK(1,KK,J))-Cpq6*CBufK(3,KK,J)
        Cone1=t1+t2
c       ------------------------------ 
        CKet(0)=COne0+2.D0*CTwo0
        SKet(0)=Zero
        CKet(1)=COne1+2.D0*CTwo1
        SKet(1)=Zero
        CKet(2)=2.D0*CTwo2
        SKet(2)=2.D0*STwo2
c       ------------------------------ 
         
