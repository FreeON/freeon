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
c       ------------------------------ 
        t1=Cpq0*CBufK(0,KK,J)+Cpq1*CBufK(1,KK,J)
        t2=Cpq3*CBufK(3,KK,J)+Cpq6*CBufK(6,KK,J)+Cpq10*CBufK(10,KK,J)
        Cone0=t1+t2
c       ------------------------------ 
        CKet(0)=COne0+2.D0*CTwo0
        SKet(0)=Zero
c       ------------------------------ 
         
