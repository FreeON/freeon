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

(* LOAD OPTIMIZING AND FORMATING ROUTINES FOR MMA AND SET ASSOCIATED OPTIONS *)

Get[StringJoin[MondoHome,"/MMA/Approximations.m"]];
Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];
Get[StringJoin[MondoHome,"/MMA/Format.m"]];
Get[StringJoin[MondoHome,"/MMA/Optimize.m"]];

FF[x_] := ToString[FixedNumberForm[SetPrecision[N[x,32],32], 16, 2]];

SetOptions[OpenWrite, PageWidth -> 200];

(* Define the Incomplete Gamma Function *)

Unprotect[D,Derivative,Function];
SetAttributes[F,Listable];
F[m_,T_] := Gamma[1/2 + m, 0, T]/(2*T^(1/2 + m))/;T!=0
F[m_,T_] := 1/(2 m + 1)/;T==0;
AssymptoticF[n_,T_]:=Abs[2 n-1]!!/(2 (2 T)^n) Sqrt[Pi/T];
PAC[m_,T_]:=N[Abs[F[m,T]-AssymptoticF[m,T]],30];
D[ F[m_,x_] , x_ ] := - F[m+1,x];
Derivative[0, 1][F][m_, x_ ] := - F[m+1,x];
Derivative[0, n_][F][m_, x_ ] := (-1)^n F[m+n,x];
Protect[D,Derivative,Function];

(* *)

WP = 50;

SetOptions[FindRoot, WorkingPrecision->WP, MaxIterations -> 1000, DampingFactor -> 100];
PACFind[m_, Err_] := x /. FindRoot[PAC[m, x] == Err, {x, {10^(-13), 100}}];

(*
HGEll=0; SPEll=0;
   *)

NPts = 2000;
ClearAll[list];
MinAcc=10^(3);
MaxAcc=10^(-16);

Do[

   TMin[L] = SetPrecision[PACFind[L, MinAcc], WP];
   TMax[L] = SetPrecision[PACFind[L, MaxAcc], WP];

   LTMin = Log[10, TMin[L]];
   LTMax = Log[10, TMax[L]];

   Delta = (LTMax - LTMin)/(NPts - 1);

   g[i_] := SetPrecision[10^(LTMin + i*Delta), WP];
   h[i_] := SetPrecision[-Log[PAC[L, g[i]]], WP];
 
   hb = h[0]; 
   he = h[NPts - 1];
 
   list = Table[{h[i], g[i]}, {i, 0, NPts - 1}];

   f = Interpolation[list];
   
   q[L] = RationalInterpolation[f[x], {x, 20, 20}, {x, hb, he}];

(*
   qtmp = SetPrecision[q[L], 14];
   Plot[qtmp - f[x], {x, hb, he}, PlotRange -> All];
   Plot[{qtmp, f[x]}, {x, hb, hb + 2}, PlotRange -> All];
*)
,{L,0,HGEll+SPEll}]

FileName=StringJoin[MondoHome,"/QCTC/PFunk.Inc"];
Print[" Openned ",FileName];
OpenWrite[FileName];

(*SetOptions[Optimize,OptimizePower->Binary,OptimizeVariable->{U,Array}];*)

SetOptions[FortranAssign,AssignOptimize->False,AssignMaxSize->300,
                         AssignBreak->{132," & \n          "},
                         AssignTemporary->{W,Array},
                         AssignIndent->"            "];

WriteString[FileName,StringJoin["     MinAcc=",FF[MinAcc],"\n"]];
WriteString[FileName,StringJoin["     MaxAcc=",FF[MaxAcc],"\n"]];
WriteString[FileName,"     SELECT CASE(Ell)\n"];

Do[
   WriteString[FileName,"     CASE(",L,") \n"];
   WriteString[FileName,"        IF(Ack>MinAcc)THEN\n"];
   WriteString[FileName,StringJoin["           PFunk=",FF[TMin[L]],"\n"]];
   WriteString[FileName,"        ELSEIF(Ack<MaxAcc)THEN\n"];
   WriteString[FileName,StringJoin["           PFunk=",FF[TMax[L]],"\n"]];
   WriteString[FileName,"        ELSE\n"];
   Write[FileName,FortranAssign[PFunk,Horner[q[L]]]];
   WriteString[FileName,"        ENDIF\n"]
  ,{L,0,HGEll+SPEll}];
WriteString[FileName,"     END SELECT\n"];
Close[FileName];          

Print[" Closed ",FileName];
