(* 
!------------------------------------------------------------------------------
!--  This code is part of the MondoSCF suite of programs for linear scaling 
!    electronic structure theory and ab initio molecular dynamics.
!
!--  Copyright (c) 2001, the Regents of the University of California.  
!    This SOFTWARE has been authored by an employee or employees of the 
!    University of California, operator of the Los Alamos National Laboratory 
!    under Contract No. W-7405-ENG-36 with the U.S. Department of Energy.  
!    The U.S. Government has rights to use, reproduce, and distribute this 
!    SOFTWARE.  The public may copy, distribute, prepare derivative works 
!    and publicly display this SOFTWARE without charge, provided that this 
!    Notice and any statement of authorship are reproduced on all copies.  
!    Neither the Government nor the University makes any warranty, express 
!    or implied, or assumes any liability or responsibility for the use of 
!    this SOFTWARE.  If SOFTWARE is modified to produce derivative works, 
!    such modified SOFTWARE should be clearly marked, so as not to confuse 
!    it with the version available from LANL.  The return of derivative works
!    to the primary author for integration and general release is encouraged. 
!    The first publication realized with the use of MondoSCF shall be
!    considered a joint work.  Publication of the results will appear
!    under the joint authorship of the researchers nominated by their
!    respective institutions. In future publications of work performed
!    with MondoSCF, the use of the software shall be properly acknowledged,
!    e.g. in the form "These calculations have been performed using MondoSCF, 
!    a suite of programs for linear scaling electronic structure theory and
!    ab initio molecular dynamics", and given appropriate citation.  
!------------------------------------------------------------------------------
     Author: Matt Challacombe
     COMPUTE FUNCTIONS THAT RETURN THE ARGUMENT T TO THE GAMMA FUNCTIONS F[m,T]
     THAT RESULT FROM USING THE THE MULTIPOLE APPROXIMATION TO WITHIN A SPECIFIED 
     ERROR:  THESE FUNCTIONS GIVE T (THE PAC) FOR A GIVEN PENETRATION ERROR.
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

FileName="PFunk.Inc";
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
