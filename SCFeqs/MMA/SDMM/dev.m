trace[m_List]:=Sum[ m[[i,i]],{i,Length[m]}]
mdot[x_List,y_List]:=Flatten[x].Flatten[y];
SpectBounds[m_List]:=Chop[{Min[Eigenvalues[m]],Max[Eigenvalues[m]]},10^(-5)];

<<scf.in

WP=14;
MaxIts=40;

Print["NEl = ",NEl," NBasF = ",NBasF," \n "];

Print[Sort[Eigenvalues[f]]];

p =IdentityMatrix[NBasF]*(NEl/(2*NBasF)); 

Print["NEl = ",NEl," NBasF = ",NBasF," Tr P = ",trace[p]," \n "];

(*
p =SetPrecision[p, WP];
f =SetPrecision[f, WP];
pp=SetPrecision[pp,WP];
*)

GradO=3*p.f+3*f.p-2*p.p.f-2*p.f.p-2*f.p.p;

mu=-trace[GradO]/NBasF;

gg=(GradO+mu*IdentityMatrix[NBasF]);

g[0]=-gg;

Print["Trace GradO = ",trace[GradO]," \n"];
Print["Trace     I = ",trace[IdentityMatrix[NBasF]]," \n"];
Print["mu          = ",mu," \n"];
Print["Trace gg = ",trace[gg]," \n"];

h=g[0];

EOld=trace[f.p];

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
         Print[" lambda = ",N[lambda,6]," Min = ",N[f3,6],"\n"];
         Print[" \n "];
        , 
         Print["ABORTING"];     
         Abort[];
        ];
   p=p+lambda*h;

(*   p=SetPrecision[p,WP] ;*)

   ENew=trace[f.p];
   Pcnt=Abs[(ENew-EOld)/EOld];
   EOld=ENew;

   GradO=3*p.f+3*f.p-2*p.p.f-2*p.f.p-2*f.p.p;

   mu=-trace[GradO]/NBasF;

   gg=(GradO+mu*IdentityMatrix[NBasF]);
   g[k+1]=-gg;


   Print[" Percent         = ",Pcnt," \n "];
   Print[" Spectral Bounds = ",N[SpectBounds[p],4],"\n"];
   Print[" Tr{P}           = ",N[trace[p],WP]];
   Print[" Tr{F.P}         = ",N[trace[p.f],WP]];
   Print[" Tr{G}           = ",N[trace[gg],WP],"\n"];


   beta  = mdot[g[k+1]-g[k],g[k+1]]/mdot[g[k],g[k]]; 

  
   h = g[k+1] + beta * h ;

,{k,0,MaxIts}];

Do[ 
   p=3 p.p-2 p.p.p;
(*   p=SetPrecision[p,WP];*)
   Print[" Tr P = ",trace[p]," Tr F.P = ",trace[p.f]];
,{i,10}];

Print["Achieved ETot  = ",N[trace[p.f],20]];
Print["Objective ETot = ",N[trace[pp.f],20]];
