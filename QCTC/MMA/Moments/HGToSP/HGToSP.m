<<Format.m
<<Optimize.m


LAngBegin[l_] :=(l*(l+1)*(l+2))/6+1;                                     
LMNIndex[l_,m_,n_] := LAngBegin[l+m+n]+(n*(2*(l+m+n)-n+3))/2+m;
LTD[l_] := l*(l+1)/2;

p=Sqrt[px^2+py^2+pz^2];
costmp=pz/p;
sintmp=Sqrt[1-costmp*costmp];

q[s_,t_]:= (p^s)*Exp[I*t*ArcTan[px/py]] LegendreP[s,t,costmp]/(s+t)! ;

ParmList   = {0};
ParmValue  = {0};

filename=StringJoin["HGToSP_L5.F90"];

MinEll=5;
MaxEll=5;

Do[Do[

      ncount=0;
      st = LTD[s]+t;

      rcc[st]=0;
      rss[st]=0;

      Print["s = ",s," t = ",t," st = ",st];

      Do[ Do[ Do[ 

	 lmn = l+m+n;

         If[lmn == s,
 
	      tmp = D[D[D[q[s,t],{px,l}],{py,m}],{pz,n}];
	      tmp = Limit[tmp,px->0];
	      tmp = Limit[tmp,py->0];
	      tmp = Limit[tmp,pz->0];

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

      rcc[st] = PiZeta*rcc[st];
      rss[st] = PiZeta*rss[st];

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
     WriteString[filename,"    IF(LQ == ",ell,") THEN \n"]; 
   ];
   If[ell != 0,
     WriteString[filename,"    ELSEIF(LQ == ",ell,") THEN \n"]; 
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
     WriteString[filename,"       WRITE(*,*) 'LQ = ',LQ,' Bad indexing in HGToSP' \n"];      
     WriteString[filename,"       STOP \n"];
     WriteString[filename,"    ENDIF \n"]; 
   ];

,{ell,MinEll,MaxEll}];

Close[filename];


