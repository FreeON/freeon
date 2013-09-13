              LP=SP%Ell
              LQ=Ell
              LPQ=MAX(LP,LQ)
              QP=HG%Cent(Ell)%D(:,I)-SP%Center
              CALL Regular(LPQ,QP(1),QP(2),QP(3))
              CALL XLate77(LP,LQ,SP%C(0),SP%S(0),Cpq(0),S(0),C(0),S(0))
