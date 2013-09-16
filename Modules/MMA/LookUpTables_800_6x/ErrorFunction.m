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
!    WRITE CHEBYSHEV INTERPOLATES FOR QUICK ERROR FUNCTION LOOKUP
!------------------------------------------------------------------------------
*)

(* GET THE MONDO INFO *)

MondoHome=Environment["FREEON_HOME"];
If[MondoHome==$FAILED,
   Print["COULD NOT FIND $FREEON_HOME! CHECK YOUR .cshrc "];
   Abort[];
  ];

(* LOAD FIXEDNUBERFORM *)

Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];

<<Rules.m;

(*                                          *)

WP = 30;
NTerms =  4;
n = NTerms ;
switch = 5.6;
NInterps = 1200;

FunctionList[x_] := {Erf[x]};

FuncList          = {"Erf"};

NFunctions=1;

ClearAll[Delta];

Delta=N[switch/(2 NInterps+2),30 ];

MaxError=DTest[13,switch,WP];

(* WRITE ARRAY DIMENSIONS TO FILE *)

FileName="ErrorFunction.Inc";
Print[" Openned ",FileName];
OpenWrite[FileName];

WriteString[FileName,"      REAL(DOUBLE), PARAMETER :: Erf_Grid  =",FF[1/(2*Delta)],"\n"];
WriteString[FileName,"      REAL(DOUBLE), PARAMETER :: Erf_Switch=",FF[switch],"\n"];
WriteString[FileName,StringJoin["      REAL(DOUBLE),DIMENSION(0:",ToString[NInterps],") :: & \n"]];
(*                                   *)
FunkStr="";
Do[
   Funk=FuncList[[i]];   
   TotStr=StringJoin["                 ",Funk,"_0,",Funk,"_1,",Funk,"_2,",Funk,"_3,",Funk,"_4, & \n"];
   FunkStr=StringJoin[FunkStr,TotStr]
  ,{i,NFunctions-1}];
Funk=FuncList[[NFunctions]];   
TotStr=StringJoin["                 ",Funk,"_0,",Funk,"_1,",Funk,"_2,",Funk,"_3,",Funk,"_4"];
FunkStr=StringJoin[FunkStr,TotStr];
WriteString[FileName,FunkStr,"\n"];
(*                                   *)
FunkStr="      INTEGER :: ";
Do[FunkStr=StringJoin[FunkStr,"I",FuncList[[i]],","],{i,NFunctions-1}];
FunkStr=StringJoin[FunkStr,"I",FuncList[[NFunctions]]];
WriteString[FileName,FunkStr,"\n"];

Do[
   WriteString[FileName,"      INCLUDE '",FuncList[[i]],".Inc' \n"];
  ,{i,NFunctions}];

Close[FileName];

(* WRITE ARRAY VALUES TO INDIVIDUAL FILES *)

<<BurnTables.m;

Close[FileName];
