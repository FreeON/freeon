trace[m_List]:=Sum[ m[[i,i]],{i,Length[m]}]
mdot[x_List,y_List]:=Flatten[x].Flatten[y];

evec=Eigenvectors[S];
eval=Eigenvalues[S];
ev=Table[Table[0,{j,NBasF}],{i,NBasF}];
Do[ev[[i,i]]=1/Sqrt[eval[[i]]],{i,NBasF}];
X=Transpose[evec].ev.evec;

it=3; MaxIts=8;

n=NBasF; xi=Inverse[X]; f=X.(F[it]).X; t=X.T.X; v=X.V.X;  

p=xi.P[it].xi;
p=IdentityMatrix[ n ]*(NEl/(2*NBasF)); 

GradO=3*p.f+3*f.p-2*p.p.f-2*p.f.p-2*f.p.p;

mu=-trace[GradO] / n ;
gg=(GradO-mu*IdentityMatrix[ n ]);
g[0]=-gg;

h=g[0];

Do[
   Print[" = = = = = = = = = = = = ",k," = = = = = = = = = = = = \n "];

   a=trace[(3 p.p - 2 p.p.p).f];
   b=trace[h.gg];    
   c=trace[ 3 h.h.f-2 p.h.h.f-4 p.h.h.f ];
   d=trace[-2*h.h.h.f];

   Print[" c = ",c," c2 = ",c2];
   disc=c^2 - 3*b*d;
   If[disc>=0,
      x1=(-c - Sqrt[c^2 - 3*b*d])/(3*d);
      x2=(-c + Sqrt[c^2 - 3*b*d])/(3*d);
      f1=(a + b x + c x^2 + d x^3)/.{x->x1};
      f2=(a + b x + c x^2 + d x^3)/.{x->x2};
      f3=Min[f1,f2];
      If[f2==f3,lambda=x2,lambda=x1];
         Print[" \n "];
         Print[" lambda = ",N[lambda,20]," Min = ",N[f3,20]];
         Print[" \n "];
        ,      
         Abort[];
        ];

   p=p+lambda*h ;

   GradO=3*p.f+3*f.p-2*p.p.f-2*p.f.p-2*f.p.p;

   mu=-trace[GradO] / n ;

   gg=(GradO-mu*IdentityMatrix[ n ]);
   g[k+1]=-gg;

   beta  = mdot[g[k+1]-g[k],g[k+1]]/mdot[g[k],g[k]]; 

   h = g[k+1] + beta * h ;

,{k,0,MaxIts}];

Do[ p=3 p.p-2 p.p.p ,{i,10}];

p=X.p.X;
Do[ p=3 p.S.p-2 p.S.p.S.p ,{i,10}];

Print["Achieved ETot  = ",N[trace[p.(F[it]+T+V)]+ENucNuc,20]];
Print["Objective ETot = ",N[trace[P[it+1].(F[it]+T+V)]+ENucNuc,20]];
