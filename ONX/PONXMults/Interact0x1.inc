c       ------------------------------- 
        PQx2=PQx*PQx 
        PQy2=PQy*PQy 
        PQ2=PQx2+PQy2+PQz*PQz 
        PQ=DSQRT(PQ2) 
        OneOverPQ=1.D0/PQ
        Zpq=OneOverPQ*PQz
        Cp0=1.0D0 
        Sp0=0.0D0 
        Rxy=DSQRT(PQx2+PQy2) 
        IF(Rxy.NE.0.0D0)THEN  
           OneOvRxy=1.0D0/Rxy 
           Cp1=PQy*OneOvRxy 
           Sp1=PQx*OneOvRxy  
        ELSE 
           Cp1=0.70710678118654752D0 
           Sp1=0.70710678118654752D0 
        ENDIF 
        TwoC=2.0D0*Cp1 
c       ------------------------------- 
        SqZpq=DSQRT(1.D0-Zpq*Zpq) 
        rss=1.D0 
        ALP0=1.D0*rss
        rss=rss*sqZpq
        ALP2=-rss
c       ------------------------------- 
        ALP1=ALP0*Zpq
c       ------------------------------- 
c       ------------------------------- 
        rr=1.D0
        rr=OneOverPQ*rr
        tt=ALP0*rr
        Cpq0=Cp0*tt
        Spq0=Sp0*tt
        rr=OneOverPQ*rr
        tt=ALP1*rr
        Cpq1=Cp0*tt
        Spq1=Sp0*tt
        tt=ALP2*rr
        Cpq2=Cp1*tt
        Spq2=Sp1*tt
c       ------------------------------ 
         
