Needs[ "HeuristicSearch`Master`" ]

LBegin[L_]:=(L*(L+1)*(L+2))/6+1;
LEnd[L_]:=LBegin[L]+L*(L+1)/2+L;
lmndex[L_,M_,N_]:=LBegin[L+M+N]+N*(2*(L+M+N)-N+3)/2+M;
LHGTF[L_]:=(L+1)*(L+2)*(L+3)/6;                                    

LTot=2;

len=LHGTF[LTot];

Do[Do[                  
     
      lmnA=lmndex[0,0,n  ];
      lmnB=lmndex[0,0,n-1];
      lmnC=lmndex[0,0,n-2];

      R[lmnA,j]=R[lmnB,j+1]*PQz+R[lmnC,j+1]*(n-1);

,{j,0,LTot-n}],{n,1,LTot}];


Do[Do[Do[

      lmnA=lmndex[0,m,n  ];
      lmnB=lmndex[0,m-1,n];
      lmnC=lmndex[0,m-2,n];

      R[lmnA,j]=R[lmnB,j+1]*PQy+R[lmnC,j+1]*(m-1);

,{j,0,LTot-m-n}],{m,1,LTot-n}],{n,0,LTot}];

Do[Do[Do[Do[

      lmnA=lmndex[l,m,n  ];
      lmnB=lmndex[l-1,m,n];
      lmnC=lmndex[l-2,m,n];

      R[lmnA,j]=R[lmnB,j+1]*PQx+R[lmnC,j+1]*(l-1); 

,{j,0,LTot-l-m-n}],{l,1,LTot-m-n}],{m,0,LTot-n}],{n,0,LTot}];


poly=Expand[Sum[ hgrho[lmn]*R[lmn,0] ,{lmn,len}]];

               ToExpression[input]

Do[ 
string=StringJoin["hgX",ToString[lmn]];Print[string];
hgrho[lmn]=ToExpression[string];

Do[
string=StringJoin["RX",ToString[lmn],"X",ToString[j]];Print[string];
R[lmn,j]=ToExpression[string];
,{j,0,LTot}];

 ,{lmn,len}];
Print[poly];

Timing[HillClimbing[poly, EvaluationFunction -> PolynomialEvaluationCost,
SuccessorFunction -> PolynomialSuccessorFunction]]

