<<Format.m
<<Optimize.m

LAngBegin[l_] :=(l*(l+1)*(l+2))/6+1;                                     
LMNIndex[l_,m_,n_] := LAngBegin[l+m+n]+(n*(2*(l+m+n)-n+3))/2+m;
LTD[l_] := l*(l+1)/2;

Fact[ll_,mm_,nn_] := ((x-RX)^ll)*((y-RY)^mm)*((z-RZ)^nn)/(ll!*mm!*nn!)

prim[l_,m_,n_] := ((-1)^(l+m+n))*D[D[D[Exp[-Zeta*(x^2+y^2+z^2)],{z,n}],{y,m}],{x,l}];


filename=StringJoin["HGToCarts.form"];


MaxL   = 2;
MaxEll = 0;



Do[

Do[Do[Do[

   ncount=0;

   lldx=LMNIndex[ll,mm,nn];

   QCarts[lldx,ell]=0;

   Do[Do[Do[
      
      If[l < ll+1,If[m < mm+1,If[n < nn+1,

             tmp = ((Zeta/Pi)^(3/2))*Integrate[Fact[ll,mm,nn]*prim[l,m,n],
	       {x,-Infinity,Infinity},{y,-Infinity,Infinity},{z,-Infinity,Infinity}];

             tmp = Expand[tmp];
             tmp = Factor[tmp];

	     tmp = Limit[tmp,Zeta->Infinity];

	     tmp = N[tmp,16];

	     ldx=LMNIndex[l,m,n];
    	 
	     QCarts[lldx,ell] = QCarts[lldx,ell]+tmp*Coefs[ldx];

      ];];];

   ,{n,0,ell-m-l}],{m,0,ell-l}],{l,0,ell}];


   QCarts[lldx,ell] = PiZeta*QCarts[lldx,ell];

   If[QCarts[lldx,ell] == 0, QCarts[lldx,ell] = Zero];

   Print["ll = ",ll," mm = ",mm," nn = ",nn];

   Print["QCarts[",lldx,",",ell,"] = ",QCarts[lldx,ell]];


,{nn,0,MaxL}],{mm,0,MaxL}],{ll,0,MaxL}];


,{ell,0,MaxEll}]

Print[" Openning ",filename];
OpenWrite[filename];


Do[

   Print["ell = ",ell];

   If[ell == 0,
     WriteString[filename,"    IF(L == ",ell,") THEN \n"]; 
   ];
   If[ell != 0,
     WriteString[filename,"    ELSEIF(L == ",ell,") THEN \n"]; 
   ];

   Do[Do[Do[

     ldx=LMNIndex[l,m,n];

     rhs =  QCarts[ldx,ell];

     lhs =  QCarts[l,m,n];


     WriteString[filename,"!\n"];   
     WriteString[filename,"          QCarts[",l,",",m,",",n,"] = ",rhs,"\n"];

   ,{n,0,MaxL}],{m,0,MaxL}],{l,0,MaxL}];

  
   If[ell == MaxEll,
     WriteString[filename,"    ELSE \n"]; 
     WriteString[filename,"       WRITE(*,*) 'L = ',L,' Bad indexing in HGToCarts' \n"];      
     WriteString[filename,"       STOP \n"];
     WriteString[filename,"    ENDIF \n"]; 
   ];

,{ell,0,MaxEll}];

Close[filename];


