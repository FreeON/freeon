c       ------------------------------ 
        CTwo2=-(Cpq2*CBufK(0,KK,J))
        STwo2=-(Spq2*CBufK(0,KK,J))
        CTwo4=Cpq4*CBufK(0,KK,J)
        STwo4=Spq4*CBufK(0,KK,J)
        CTwo5=Cpq5*CBufK(0,KK,J)
        STwo5=Spq5*CBufK(0,KK,J)
c       ------------------------------ 
        Cone0=Cpq0*CBufK(0,KK,J)
        Cone1=-(Cpq1*CBufK(0,KK,J))
        Cone3=Cpq3*CBufK(0,KK,J)
c       ------------------------------ 
        CKet(0)=COne0
        SKet(0)=Zero
        CKet(1)=COne1
        SKet(1)=Zero
        CKet(2)=2.D0*CTwo2
        SKet(2)=2.D0*STwo2
        CKet(3)=COne3
        SKet(3)=Zero
        CKet(4)=2.D0*CTwo4
        SKet(4)=2.D0*STwo4
        CKet(5)=2.D0*CTwo5
        SKet(5)=2.D0*STwo5
c       ------------------------------ 
         
