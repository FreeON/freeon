c       ------------------------------ 
        t1=Cpq2*CBufK(2,KK,J)+Cpq4*CBufK(4,KK,J)+Cpq5*CBufK(5,KK,J)
        t2=Cpq7*CBufK(7,KK,J)+Cpq8*CBufK(8,KK,J)+Cpq9*CBufK(9,KK,J)
        t3=Spq2*SBufK(2,KK,J)+Spq4*SBufK(4,KK,J)+Spq5*SBufK(5,KK,J)
        t4=Spq7*SBufK(7,KK,J)+Spq8*SBufK(8,KK,J)+Spq9*SBufK(9,KK,J)
        CTwo0=t1+t2+t3+t4
c       ------------------------------ 
        t1=Cpq0*CBufK(0,KK,J)+Cpq1*CBufK(1,KK,J)
        t2=Cpq3*CBufK(3,KK,J)+Cpq6*CBufK(6,KK,J)
        Cone0=t1+t2
c       ------------------------------ 
        CKet(0)=COne0+2.D0*CTwo0
        SKet(0)=Zero
c       ------------------------------ 
         
