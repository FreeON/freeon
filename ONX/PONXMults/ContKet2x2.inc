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
        t1=Cpq7*CBufK(2,KK,J)+Cpq11*CBufK(4,KK,J)+Cpq12*CBufK(5,KK,J)
        t2=Spq7*SBufK(2,KK,J)+Spq11*SBufK(4,KK,J)+Spq12*SBufK(5,KK,J)
        CTwo3=t1+t2
        t1=Cpq4*CBufK(0,KK,J)+Cpq7*CBufK(1,KK,J)
        t2=(-Cpq6+Cpq8)*CBufK(2,KK,J)+Cpq11*CBufK(3,KK,J)
        t3=(-Cpq10+Cpq12)*CBufK(4,KK,J)+(-Cpq11+Cpq13)*CBufK(5,KK,J)
        t4=Spq8*SBufK(2,KK,J)
        t5=Spq12*SBufK(4,KK,J)+(-Spq11+Spq13)*SBufK(5,KK,J)
        CTwo4=t1+t2+t3+t4+t5
        t1=Spq4*CBufK(0,KK,J)+Spq7*CBufK(1,KK,J)
        t2=Spq8*CBufK(2,KK,J)+Spq11*CBufK(3,KK,J)
        t3=Spq12*CBufK(4,KK,J)+(Spq11+Spq13)*CBufK(5,KK,J)
        t4=-((Cpq6+Cpq8)*SBufK(2,KK,J))
        t5=-((Cpq10+Cpq12)*SBufK(4,KK,J))-(Cpq11+Cpq13)*SBufK(5,KK,J)
        STwo4=t1+t2+t3+t4+t5
        t1=Cpq5*CBufK(0,KK,J)+Cpq8*CBufK(1,KK,J)
        t2=(-Cpq7+Cpq9)*CBufK(2,KK,J)+Cpq12*CBufK(3,KK,J)
        t3=(-Cpq11+Cpq13)*CBufK(4,KK,J)+(Cpq10+Cpq14)*CBufK(5,KK,J)
        t4=(Spq7+Spq9)*SBufK(2,KK,J)+(Spq11+Spq13)*SBufK(4,KK,J)+Spq1
     $   4*SBufK(5,KK,J)
        CTwo5=t1+t2+t3+t4
        t1=Spq5*CBufK(0,KK,J)+Spq8*CBufK(1,KK,J)
        t2=(-Spq7+Spq9)*CBufK(2,KK,J)+Spq12*CBufK(3,KK,J)
        t3=(-Spq11+Spq13)*CBufK(4,KK,J)+Spq14*CBufK(5,KK,J)
        t4=-((Cpq7+Cpq9)*SBufK(2,KK,J))
        t5=-((Cpq11+Cpq13)*SBufK(4,KK,J))+(Cpq10-Cpq14)*SBufK(5,KK,J)
        STwo5=t1+t2+t3+t4+t5
c       ------------------------------ 
        Cone0=Cpq0*CBufK(0,KK,J)+Cpq1*CBufK(1,KK,J)+Cpq3*CBufK(3,KK,J
     >   )
        t1=-(Cpq1*CBufK(0,KK,J))
        t2=-(Cpq3*CBufK(1,KK,J))-Cpq6*CBufK(3,KK,J)
        Cone1=t1+t2
        Cone3=Cpq3*CBufK(0,KK,J)+Cpq6*CBufK(1,KK,J)+Cpq10*CBufK(3,KK,
     >   J)
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
         
