<<Format.m
<<Optimize.m

LAngBegin[l_] :=(l*(l+1)*(l+2))/6+1;                                     
LMNIndex[l_,m_,n_] := LAngBegin[l+m+n]+(n*(2*(l+m+n)-n+3))/2+m;
LTD[l_] := l*(l+1)/2;

q[s_,t_]:= (r^s)*Exp[I*t*phi] LegendreP[s,t,Cos[theta]]/(s+t)! ;

prim[l_,m_,n_] := D[D[D[Exp[-Zeta*(x*x+y*y+z*z)],{z,n}],{y,m}],{x,l}];

ParmList   = {0};
ParmValue  = {0};

filename=StringJoin["HGToSP_L0to2_new.F90"];


MinEll = 0;
MaxEll = 2;

Do[Do[

   ncount=0;
   st = LTD[s]+t;

   rcc[st]=0;
   rss[st]=0;

   Print["s = ",s," t = ",t," st = ",st];

   Do[Do[Do[

      lmn = l+m+n;

      If[(l+m+n) == s,

	 tmp = ((-1)^(s))*prim[l,m,n];
	 tmp = Limit[tmp,z -> r*Cos[theta]];
	 tmp = Limit[tmp,y -> r*Sin[theta]*Sin[phi]];
	 tmp = Limit[tmp,x -> r*Sin[theta]*Cos[phi]];

	 Factor[tmp];

         tmp = ((Zeta/Pi)^(3/2))*Integrate[tmp*q[s,t]*r*r*Sin[theta],{theta,0,Pi},
               {phi,0,2*Pi},{r,0,Infinity}];


	 tmpreal=N[Re[tmp],16];
	 tmpimag=N[Im[tmp],16];
    	 
	 ldx=LMNIndex[l,m,n];
 	 ncount = ncount+1;

	 Print[" l=",l," m=",m," n=",n];
 	 Print[" ncount = ",ncount," ldx = ",ldx];

	 If[tmpreal != 0,
	   CC = ToString[Abs[tmpreal]];
	   CC = StringDrop[CC,{2,2}];
	   CC = StringInsert[CC,"N",1];
	   CC = ToExpression[CC];
	   If[Count[ParmList,CC] == 0,
	     ParmList  = Append[ParmList,CC];
	     ParmValue = Append[ParmValue,Abs[tmpreal]];
	   ];
	   If[tmpreal < 0,
	     If[tmpreal == -1.,
                rcc[st] = rcc[st]-Coefs[ldx];
    	     ];
	     If[tmpreal != -1.,
  	        rcc[st] = rcc[st]-CC*Coefs[ldx];
	     ];
	   ];
	   If[tmpreal > 0,
	     If[tmpreal == 1.,
  	        rcc[st] = rcc[st]+Coefs[ldx];
	     ];
	     If[tmpreal != 1.,
  	        rcc[st] = rcc[st]+CC*Coefs[ldx];
	     ];
	   ];
         ];
	 If[tmpimag != 0,
	   CC = ToString[Abs[tmpimag]];
	   CC = StringDrop[CC,{2,2}];
	   CC = StringInsert[CC,"N",1];
	   CC = ToExpression[CC];
	   If[Count[ParmList,CC] == 0,
	     ParmList  = Append[ParmList,CC];
	     ParmValue = Append[ParmValue,Abs[tmpimag]];
	   ];
	   If[tmpimag < 0,
	     If[tmpimag == -1.,
  	        rss[st] = rss[st]-Coefs[ldx];
	     ];
	     If[tmpimag != -1.,
  	        rss[st] = rss[st]-CC*Coefs[ldx];
	     ];
	   ];
	   If[tmpimag > 0,
	     If[tmpimag == 1.,
                rss[st] = rss[st]+Coefs[ldx];
	     ];
	     If[tmpimag != 1.,
  	        rss[st] = rss[st]+CC*Coefs[ldx];
	     ];
	   ];
         ];

    ];

   ,{n,0,MaxEll-m-l}],{m,0,MaxEll-l}],{l,0,MaxEll}];

   rcc[st] = Z*rcc[st];
   rss[st] = Z*rss[st];

   If[rcc[st] == 0, rcc[st] = Zero];
   If[rss[st] == 0, rss[st] = Zero];

,{t,0,s}],{s,MinEll,MaxEll}];



Print[ParmList];
Print[ParmValue];


Print[" Openning ",filename];
OpenWrite[filename];

WriteString[filename,"!\n"];

Do[

  WriteString[filename,"    REAL(DOUBLE),PARAMETER           :: ",ParmList[[n]]," = ",ParmValue[[n]],"\n"];

,{n,3,Length[ParmList]}];

WriteString[filename,"!\n"];


Do[

   Print["ELL = ",ell];

   If[ell == 0,
     WriteString[filename,"    IF(L == ",ell,") THEN \n"]; 
   ];
   If[ell != 0,
     WriteString[filename,"    ELSEIF(L == ",ell,") THEN \n"]; 
   ];

   Do[Do[

     st = LTD[s]+t;
     WriteString[filename,"       RC(",st,") = "];
     Write[filename,FortranForm[rcc[st]]];
     WriteString[filename,"       RS(",st,") = "];
     Write[filename,FortranForm[rss[st]]];

   ,{t,0,s}],{s,0,ell}];

  
   If[ell == MaxEll,
     WriteString[filename,"    ELSE \n"]; 
     WriteString[filename,"       WRITE(*,*) 'L = ',L,' Bad indexing in HGToSP' \n"];      
     WriteString[filename,"       STOP \n"];
     WriteString[filename,"    ENDIF \n"]; 
   ];

,{ell,MinEll,MaxEll}];

Close[filename];
