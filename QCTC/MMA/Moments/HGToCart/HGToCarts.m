<<Format.m
<<Optimize.m

LAngBegin[l_] :=(l*(l+1)*(l+2))/6+1;                                     
LMNIndex[l_,m_,n_] := LAngBegin[l+m+n]+(n*(2*(l+m+n)-n+3))/2+m;
LTD[l_] := l*(l+1)/2;

Fact[ll_,mm_,nn_,l_,m_,n_] := D[D[D[(((x-RX)^ll)*((y-RY)^mm)*((z-RZ)^nn)/(ll!*mm!*nn!)),{z,n}],{y,m}],{x,l}];

prim = Exp[-Zeta*(x^2+y^2+z^2)]

CoefsRule={Coefs[ldx_]:>ToExpression[StringJoin["CoefsLPAR",ToString[ldx],"RPAR"]]};
CoefsRuleLMN={Coefs[l_,m_,n_]:>ToExpression[StringJoin["CoefsLPAR",ToString[l],ToString[m],ToString[n],"RPAR"]]};

filename=StringJoin["HGToCarts_L3.F90"];


MaxL   = 3;

Do[

Do[Do[Do[

   ncount=0;

   lldx=LMNIndex[ll,mm,nn];

   QCarts[lldx,ell]=0;

   If[(ll+mm+nn) < MaxL+1,

   Do[Do[Do[
      
      If[l < ll+1,If[m < mm+1,If[n < nn+1,

             tmp = ((Zeta/Pi)^(3/2))*Integrate[Fact[ll,mm,nn,l,m,n]*prim,
	       {x,-Infinity,Infinity},{y,-Infinity,Infinity},{z,-Infinity,Infinity}];

             tmp = Expand[tmp];
    	 
             ldx=LMNIndex[l,m,n];

	     Phase = (-1)^(l+m+n);

	     QCarts[lldx,ell] = QCarts[lldx,ell]+tmp*Phase*Coefs[ldx] /.CoefsRule;

      ];];];

   ,{n,0,ell}],{m,0,ell}],{l,0,ell}];


   QCarts[lldx,ell] = PiZeta*Expand[QCarts[lldx,ell]];

   If[QCarts[lldx,ell] == 0, QCarts[lldx,ell] = Zero];

   Print["ll = ",ll," mm = ",mm," nn = ",nn];

   Print["QCarts[",lldx,",",ell,"] = ",QCarts[lldx,ell]];

   ];


,{nn,0,MaxL}],{mm,0,MaxL}],{ll,0,MaxL}];


,{ell,0,MaxL}]

Print[" Openning ",filename];
OpenWrite[filename];

SetOptions[FortranAssign,AssignReplace->{"LPAR"->"(","RPAR"->")"},AssignBreak->{71,"\n"}];

Do[

   Print["ell = ",ell];

   If[ell == 0,
     WriteString[filename,"    IF(LP == ",ell,") THEN \n"]; 
   ];
   If[ell != 0,
     If[ell != MaxL,
       WriteString[filename,"    ELSEIF(LP == ",ell,") THEN \n"]; 
     ];
   ];
   If[ell == MaxL,
     WriteString[filename,"    ELSEIF(LP >= ",ell,") THEN \n"]; 
   ];

   Do[Do[Do[

     If[(l+m+n) < MaxL+1,

     ldx=LMNIndex[l,m,n];

     rhs =  QCarts[ldx,ell];

     lhs =  QCarts[l,m,n];
 
     Write[filename,FortranAssign[Evaluate[lhs],rhs]];

     ];

   ,{n,0,MaxL}],{m,0,MaxL}],{l,0,MaxL}];

  
   If[ell == MaxL,
     WriteString[filename,"    ENDIF \n"]; 
   ];

,{ell,0,MaxL}];

Close[filename];


