<<scf.in;

trace[m_List]:=Sum[ m[[i,i]],{i,Length[m]}]

mdot[x_List,y_List]:=Flatten[x].Flatten[y];

SpectBounds[m_List]:=Chop[{Min[Eigenvalues[m]],Max[Eigenvalues[m]]},10^(-5)];

fxpt[p_,sw_]:=Module[{p2,p3}, p2=p.p; p3=p.p2; (1/2)sw+(1-sw)trace[p2-p3]/trace[p-p2]];

coef[p_,sw_]:=Module[{c,u,v,w},
                    c=fxpt[p,sw]; 
                    If[c<=1/2, u=(1-2c)/(1-c); v=(1+c)/(1-c); w=-1/(1-c);,
                               u=0;            v=(1+c)/c;     w=-1/c;];
                {u,v,w}];
 
grad[p_,f_,uvw_]:=Module[{u,v,w,G0,mu},
                          u=uvw[[1]]; v=uvw[[2]]; w=uvw[[3]];
                          G0=(u*Id+(2*v*Id +3 w p).p).f;
                          mu=-trace[G0]/NBasF;
                          G0+mu*IdentityMatrix[NBasF]];

line[p_,f_,h_,uvw_]:=Module[{a,b,c,d,u,v,w,disc,x1,x2,f1,f2,f3,lambda},
                             u=uvw[[1]]; v=uvw[[2]]; w=uvw[[3]]; 
                             a=trace[(u p  + v p.p + w p.p.p).f];
                             b=trace[(u Id + 2 v p + 3 w p.p).h.f];
                             c=trace[(v Id + 3 w p).h.h.f];
                             d=trace[ w h.h.h.f ];
                             disc=c^2 - 3*b*d;
                             If[disc>=0,
                                x1=(-c - Sqrt[c^2 - 3*b*d])/(3*d);
                                x2=(-c + Sqrt[c^2 - 3*b*d])/(3*d);
                                f1=(a + b x + c x^2 + d x^3)/.{x->x1};
                                f2=(a + b x + c x^2 + d x^3)/.{x->x2};
                                f3=Min[f1,f2];
                                If[f2==f3,lambda=x2,lambda=x1];
                                (*
                                Print[" \n "];
                                Print[" lambda = ",N[lambda,6]," Min = ",N[f3,6],"\n"];
                                Print[" \n "];
                                *)
                               , 
                                Print["ABORTING"];     
                                Abort[];
                               ];
                            lambda]


WP=14;
MaxCG=11;
MaxPur=3;

Print["NEl = ",NEl," NBasF = ",NBasF," \n "];

Id=IdentityMatrix[NBasF];

p=Id*(NEl/(2*NBasF)); 

uvw=coef[p,1];

g[0]=-grad[p,f,uvw];

h=g[0];

Do[
   Print[" = = = = = = = = = = = = ",k," = = = = = = = = = = = = \n "];

   lambda=line[p,f,h,uvw];

   p=p+lambda*h;

   uvw=coef[p,1];

   g[k+1]=-grad[p,f,uvw];

   beta=mdot[g[k+1]-g[k],g[k+1]]/mdot[g[k],g[k]];

   h=g[k+1]+beta*h;

   Print[" lambda = ",lambda," beta = ",beta];
   Print[" Spectral Bounds   = ",N[SpectBounds[p],4]];
   Print[" FixedPointC       = ",N[ fxpt[p,0],WP]];
   Print[" Tr{p.f}           = ",N[trace[p.f],WP]];

,{k,0,MaxCG}];


Do[ 
   c=fxpt[p,0]; uvw=coef[p,0]; 
   Print[" c = ",c];
   Print[" uvw = ",uvw];
   p=uvw[[1]]*p+uvw[[2]]*p.p+uvw[[3]]*p.p.p;
   Print[" dotp = ",mdot[p,p]];
   Print[i," c = ",c," Tr P = ",trace[p]," Tr F.P = ",trace[p.f]];
,{i,MaxPur}];


Print["Achieved ETot  = ",N[trace[p.f],20]];
Print["Objective ETot = ",N[trace[pp.f],20]];
