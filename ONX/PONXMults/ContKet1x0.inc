c       ------------------------------ 
        CTwo2=-(Cpq2*CBufK(0,KK,J))
        STwo2=-(Spq2*CBufK(0,KK,J))
c       ------------------------------ 
        Cone0=Cpq0*CBufK(0,KK,J)
        Cone1=-(Cpq1*CBufK(0,KK,J))
c       ------------------------------ 
        CKet(0)=COne0
        SKet(0)=Zero
        CKet(1)=COne1
        SKet(1)=Zero
        CKet(2)=2.D0*CTwo2
        SKet(2)=2.D0*STwo2
c       ------------------------------ 
         
