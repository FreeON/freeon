c       ------------------------------ 
        t1=Cpq2*CBufK(2,KK,J)+Cpq4*CBufK(4,KK,J)+Cpq5*CBufK(5,KK,J)
        t2=Spq2*SBufK(2,KK,J)+Spq4*SBufK(4,KK,J)+Spq5*SBufK(5,KK,J)
        CTwo0=t1+t2
c       ------------------------------ 
        Cone0=Cpq0*CBufK(0,KK,J)+Cpq1*CBufK(1,KK,J)+Cpq3*CBufK(3,KK,J
     >   )
c       ------------------------------ 
        CKet(0)=COne0+2.D0*CTwo0
        SKet(0)=Zero
c       ------------------------------ 
         
