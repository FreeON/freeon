(* 
  This MMA code is part of the MondoSCF suite of 
  linear scaling electronic structure codes.  

  Matt Challacombe
  Los Alamos National Laboratory
  Copywrite 2000, The University of California
*)

(* GET THE MAX ANGULAR SYMMETRY TO BE USED *)

MondoHome=Environment["MONDO_HOME"];
If[MondoHome==$FAILED,
   Print["COULD NOT FIND $MONDO_HOME! CHECK YOUR .cshrc "];
   Abort[];
  ];
EllFile = StringJoin[MondoHome,"/Includes/Ell.m"];
Get[EllFile];

(* INDEXING OF ARRAY ELEMENTS *)

LBegin[L_]:=(L*(L+1)*(L+2))/6+1;
LMNDex[L_,M_,N_]:=LBegin[L+M+N]+N*(2*(L+M+N)-N+3)/2+M
SPDex[L_,M_]:=Binomial[L+1,2]+M;
LSP[L_]:=L*(L+3)/2;   

(* THE IRREGULAR SPHERICAL HARMONICS *)
(*
p=Sqrt[px^2+py^2+pz^2];
cs=pz/p;
sn=Sqrt[1-cs^2];
qc[l_,m_]:= p^l Cos[ m ArcTan[py,px]] LegendreP[l,m,cs]/(l+m)! ;
qs[l_,m_]:= p^l Sin[ m ArcTan[py,px]] LegendreP[l,m,cs]/(l+m)! ;
*)
ClearAll[p,px,py,pz,cos,sin,e,y,q];

p=Sqrt[px^2+py^2+pz^2]; cos=pz/p;
qc[n_,m_]:= p^n Cos[ m ArcTan[py,px]] LegendreP[n,m,cos]/(n+m)! ;
qs[n_,m_]:= p^n Sin[ m ArcTan[py,px]] LegendreP[n,m,cos]/(n+m)! ;

(* TRANSFORMATION FUNCTIONS FROM HG TO SP *)

HGToSPC[L_,M_,l_,m_,n_]:=Module[{LMN,DC},
                                 DC=D[D[D[qc[L,M],{px,l}],{py,m}],{pz,n}];
                                 DC=Limit[DC,px->0]; 
                                 DC=Limit[DC,py->0]; 
                                 DC=Limit[DC,pz->0]; 
                                 LMN=LMNDex[l,m,n];
                                 DC*PiZ*HGCo[LMN]
                                ];

HGToSPS[L_,M_,l_,m_,n_]:=Module[{LMN,DS},
                                 DS=D[D[D[qs[L,M],{px,l}],{py,m}],{pz,n}];
                                 DS=Limit[DS,px->0]; 
                                 DS=Limit[DS,py->0]; 
                                 DS=Limit[DS,pz->0];                                     
                                 LMN=LMNDex[l,m,n];
                                 DS*PiZ*HGCo[LMN]
                                ];

(* GENERATE THE EXPLICIT TRANSFORMATIONS *)

Do[Do[
      LM=SPDex[L,M];
      SPC[LM]=0;
      SPS[LM]=0;
      Do[Do[Do[  
               SPC[LM]=SPC[LM]+HGToSPC[L,M,l,m,n];
               SPS[LM]=SPS[LM]+HGToSPS[L,M,l,m,n];
      ,{n,0,L-m-l}],{m,0,L-l}],{l,0,L}];
      SPC[LM]=Factor[Simplify[SPC[LM]]];
      SPS[LM]=Factor[Simplify[SPS[LM]]];
,{M,0,L}],{L,0,HGEll}];

(* LOAD OPTIMIZING AND FORMATING ROUTINES FOR MMA AND SET ASSOCIATED OPTIONS *)

Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];
Get[StringJoin[MondoHome,"/MMA/Format.m"]];
Get[StringJoin[MondoHome,"/MMA/Optimize.m"]];

SetOptions[FortranAssign,AssignOptimize->False,AssignMaxSize->500,AssignBreak->{132," & \n          "},
                         AssignTemporary->{W,Array},AssignTemporary->{W,Array},AssignToArray->{Cp,Sp,HGCo}];

SetOptions[OpenWrite, PageWidth -> 200];

(* PUT THE TRANSFORMATIONS TO FILE *)

FileName=StringJoin[MondoHome,"/QCTC/HGToSP.Inc"];
Print[" Openned ",FileName];
OpenWrite[FileName];

Do[
   WriteString[FileName,"     CASE(",L,") \n"];
   Do[
      ls=ToString[l];
      Write[FileName,FortranAssign[c,SPC[l],AssignReplace->{"c"->StringJoin["C(",ls,")"]}]];
      Write[FileName,FortranAssign[s,SPS[l],AssignReplace->{"s"->StringJoin["S(",ls,")"]}]];
     ,{l,0,LSP[L]}];
  ,{L,0,HGEll}];
Close[FileName];          
Print[" Closed ",FileName];





