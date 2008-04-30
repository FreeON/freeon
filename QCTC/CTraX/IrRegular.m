
IrRegular[FileName_,Ell_,Tag_]:=Block[{},

(*  LOCAL DEFINITIONS *)
FactOlm0[m_]:=(-1)^m (2 m-1)!! ;
FactOlm2[m_,l_]:= 1/((l+m)*(l-m));
FactMlm2[m_,l_]:= (l+m-1)*(l-m-1);

WS[StringJoin["SUBROUTINE IrRegular",ToString[Tag],"(PQx,PQy,PQz,Cpq,Spq)"]];      
lensp=ToString[LSP[LHG+LLSP]];
WS[           "   IMPLICIT NONE"];
WS[StringJoin["   REAL*8 Cpq(0:",lensp,"),Spq(0:",lensp,"),LegendreP(0:",lensp,")"]];
WS[StringJoin["   REAL*8 Sine(0:",ToString[Ell],"),Cosine(0:",ToString[Ell],"),CoFact(2:",ToString[Ell],")"]];
WS[           "   REAL*8 PQx,PQy,PQz,PQ2,pqx2,pqy2,pqxy,PQ"];
WS[           "   REAL*8 Zero,One,Two,TT,TwoC,OneOvPQ"];
WS[           "   REAL*8 CoTan,OneOvPQxy,RS,SQ,PQToThMnsL"];
WS[           "   PARAMETER (Two=2D0,One=1D0,Zero=0D0)"];

If[Ell==0,WS["   Cpq(0)=One/DSQRT(PQ2)"];Return[]];

WS["   PQx2=PQx*PQx"];
WS["   PQy2=PQy*PQy"];
WS["   PQ=DSQRT(PQx2+PQy2+PQz*PQz)"];
WS["   OneOvPQ=One/PQ"];
WS["   CoTan=PQz*OneOvPQ"];      

WC["   Compute Cosine and Sine by recursion "];

WS["   Sine(0)  =Zero"];
WS["   Cosine(0)=One "];
WS["   PQxy=DSQRT(PQx2+PQy2) "];

WS["   IF(PQxy/=Zero)THEN  "];
WS["      OneOvPQxy=One/PQxy "];

(* Here is where to fix the ArcTan(Y/X)-> ArcTan(X/Y)! *)

WS["      Sine(1)  =PQx*OneOvPQxy  "];
WS["      Cosine(1)=PQy*OneOvPQxy "];
WS["   ELSE "];
WS["      Sine(1)  =0.70710678118654752D0 "];
WS["      Cosine(1)=0.70710678118654752D0 "];
WS["   ENDIF "];
If[Ell>1,WS["   TwoC=Two*Cosine(1) "]];

Do[
   ms=ToString[m]; 
   ms1=ToString[m-1]; 
   ms2=ToString[m-2]; 
   WS[StringJoin["   Sine(",ms,")=TwoC*Sine(",ms1,")-Sine(",ms2,")"]];
   WS[StringJoin["   Cosine(",ms,")=TwoC*Cosine(",ms1,")-Cosine(",ms2,") "]];
,{m,2,Ell}]; 

WC["   ----------------------------------------------------------------"];
WC["   Associated Legendre Polynomials by recursion"];
WC[""];
WC["   M=0,Ell; Compute LegendreP(M,M)"];
WS["   RS=One"];
WS["   Sq=SQRT(ABS(One-CoTan**2)) "];
Do[ 
   mdx=ToString[SPDex[m,m]];
   WS[StringJoin["   LegendreP(",mdx,")=",FF[FactOlm0[m]],"*RS"]];
   If[m!=Ell,WS["   RS=RS*Sq"]];
,{m,0,Ell}];

WC["   ----------------------------------------------------------------"];
WC["   M=0,Ell-1; Compute LegendreP(M+1,M)"];
Do[ 
   mdx=ToString[SPDex[m,m]];
   md1=ToString[SPDex[m+1,m]];
   WS[StringJoin["   LegendreP(",md1,")=",FF[2*m+1],"*CoTan*LegendreP(",mdx,")"]];
,{m,0,Ell-1}];
WC["   ----------------------------------------------------------------"];
Do[
   ldx=ToString[l];
   WS[StringJoin["   CoFact(",ldx,")=",FF[2*l-1],"*CoTan"]];
  ,{l,2,Ell}];
WC["   ----------------------------------------------------------------"];
WC["   M=0,Ell-2; L=M+2,Ell; Compute LegendreP(L,M)"];
Do[
   Do[ 
      ldx=ToString[l];
      lmdex=ToString[SPDex[l,m]];
      lmdm1=ToString[SPDex[l-1,m]];
      lmdm2=ToString[SPDex[l-2,m]];
      WS[StringJoin["   LegendreP(",lmdex,")=CoFact(",ldx,")*LegendreP(",lmdm1,")",FF[-FactMlm2[m,l]],"*LegendreP(",lmdm2,")"]];
     ,{l,m+2,Ell}]
  ,{m,0,Ell-2}];
WC["   ----------------------------------------------------------------"];
WC["   Compute Cpq(L,M) and Spq(L,M)"];
WS["   PQToThMnsL=OneOvPQ"];
Do[  
   Do[ 
      mdx=ToString[m];
      lmdx=ToString[SPDex[l,m]];
      WS[StringJoin["   TT=PQToThMnsL*LegendreP(",lmdx,")"]];
      WS[StringJoin["   Spq(",lmdx,")=TT*Sine(",mdx,")"]];
      WS[StringJoin["   Cpq(",lmdx,")=TT*Cosine(",mdx,")"]];
     ,{m,0,l}];
      If[l<Ell,WS["   PQToThMnsL=PQToThMnsL*OneOvPQ"]];
  ,{l,0,Ell}];

WS["RETURN"];
WS["END"];


]
