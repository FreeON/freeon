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
!    Author: Matt Challacombe
!    WRITE EXPLICIT CODE FOR COULOMB INTEGRALS OVER HERMITE GAUSSIANS
!------------------------------------------------------------------------------
*)
(* GET THE MAX ANGULAR SYMMETRY TO BE USED *)

MondoHome=Environment["MONDO_HOME"];
If[MondoHome==$FAILED,
   Print["COULD NOT FIND $MONDO_HOME! CHECK YOUR .cshrc "];
   Abort[];
  ];
EllFile = StringJoin[MondoHome,"/Includes/Ell.m"];
Get[EllFile];

(* INDEXING OF HG ARRAY ELEMENTS *)

LBegin[L_]:=(L*(L+1)*(L+2))/6+1;
LMNDex[L_,M_,N_]:=LBegin[L+M+N]+N*(2*(L+M+N)-N+3)/2+M;
HGLen[L_]:=(L+1)*(L+2)*(L+3)/6;                                    

(* LOAD OPTIMIZING AND FORMATING ROUTINES FOR MMA AND SET ASSOCIATED OPTIONS *)

Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];
Get[StringJoin[MondoHome,"/MMA/Format.m"]];
Get[StringJoin[MondoHome,"/MMA/Optimize.m"]];

SetOptions[Optimize,OptimizeCoefficients->True,
                    OptimizeTimes -> True, 
                    OptimizePlus -> True, 
                    OptimizePower -> Binary,
                    OptimizeVariable->{U,Array}]

SetOptions[FortranAssign,AssignOptimize->False,
                         AssignMaxSize->500,AssignBreak->{132," & \n          "},
                         AssignTemporary->{W,Array},
                         AssignToArray->{AuxR,Co,r,RR}];

SetOptions[OpenWrite, PageWidth -> 200];

(* GLOBAL FORMATING FUNCTIONS *)

FF[x_] := ToString[FixedNumberForm[SetPrecision[N[x,32],32], 16, 2]];
WS[String_]:=WriteString[FileName,"      ",String,"\n"];
WC[String_]:=WriteString[FileName,"!     ",String,"\n"];

<<Gammas.m;
<<RecThreeT.m;
<<HGContract.m;

(* WRITE CODE FOR AUXILIARY INTEGEGRALS AND 3 TERM RECURENCE RELATIONS *)

FileName="HGTraX.Inc"; 
Print[" Openned ",FileName];
OpenWrite[FileName];

GammasPrelim;

WS["SELECT CASE(Ell)"]
Do[
   ldx=ToString[l];
   WS[StringJoin["CASE(",ldx,")"]];
   Gammas[l];   
   ThreeTermRR[l];
  ,{l,0,2*HGEll-1}]
WS["END SELECT"]

(* WRITE CODE FOR HG KET CONTRACTION *)

WS["SELECT CASE(LCode)"]
Do[
   Do[
      lcode=lbra*100+lket;
      ldx=ToString[lcode];
      WS[StringJoin["CASE(",ldx,")"]];
      HGContract[lbra,lket];
     ,{lket,0,HGEll-1}];
  ,{lbra,0,HGEll}];
WS["END SELECT"]

Close[FileName];
Print[" Closed ",FileName];
