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
!    WRITE CHEBYSHEV INTERPOLATES FOR QUICK GAMMA FUNCTION LOOKUP
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

NFunctions=2*HGEll+2;

(* LOAD FIXEDNUBERFORM *)

Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];

(*                                          *)

<<Rules.m;

(* WRITE ASSYMPTOTIC COEFFICIENTS TO FILE *)

FileName="GammaAssymptotics.Inc";
Print[" Openned ",FileName];
OpenWrite[FileName];
WriteString[FileName,"      REAL(DOUBLE), PARAMETER :: GammAss(0:",ToString[NFunctions-1],")= (/& \n"];
Do[ 
WriteString[FileName,"                                ",FF[Abs[2 n - 1]!!/(2 (2)^n) Sqrt[Pi]],", & \n"];
   ,{n,0,NFunctions-2}];
n=NFunctions-1
WriteString[FileName,"                                ",FF[Abs[2 n - 1]!!/(2 (2)^n) Sqrt[Pi]],"/) \n"];
WriteString[FileName,"\n"];
Close[FileName];

Abort[];

(*                                          *)

WP = 30;
NTerms =  6;
n = NTerms ;
switch = 26;
NInterps = 800;

FunctionList[x_] := {F[0,x],F[1,x],F[2,x],F[3,x],F[4,x],F[5,x],F[6,x],F[7,x],F[8,x],
			      F[9,x],F[10,x],F[11,x],F[12,x],F[13,x],F[14,x],F[15,x],F[16,x]}; 

(*
FunctionList[x_] := {Exp[x]*F[0,x],Exp[x]*F[1,x],Exp[x]*F[2,x],Exp[x]*F[3,x],Exp[x]*F[4,x],Exp[x]*F[5,x],
                     Exp[x]*F[6,x],Exp[x]*F[7,x],Exp[x]*F[8,x],Exp[x]*F[9,x],Exp[x]*F[10,x],Exp[x]*F[11,x],
                     Exp[x]*F[12,x],Exp[x]*F[13,x],Exp[x]*F[14,x],Exp[x]*F[15,x],Exp[x]*F[16,x]}; 
*)


FuncList          = {"F0","F1","F2","F3","F4","F5","F6","F7","F8",  
                      "F9","F10","F11","F12","F13","F14","F15","F16"}; 

ClearAll[Delta];

Delta=N[switch/(2 NInterps+2),30 ];

MaxError=DTest[2*HGEll+2,switch,WP];

Print["MaxError = ",MaxError];
Abort[];

(* WRITE ARRAY DIMENSIONS TO FILE *)

FileName="GammaDimensions.Inc";
Print[" Openned ",FileName];
OpenWrite[FileName];
WriteString[FileName,"      INTEGER,      PARAMETER :: Gamma_Mesh  =",ToString[NInterps],"\n"];
WriteString[FileName,"\n"];
Close[FileName];

FileName="GammaDimensions_77.Inc";
Print[" Openned ",FileName];
OpenWrite[FileName];
WriteString[FileName,"      INTEGER Gamma_Mesh \n"];
WriteString[FileName,"      PARAMETER (Gamma_Mesh  =",ToString[NInterps],")\n"];
WriteString[FileName,"\n"];
Close[FileName];

FileName="GammaGrid.Inc";
Print[" Openned ",FileName];
OpenWrite[FileName];
WriteString[FileName,"      REAL(DOUBLE), PARAMETER :: Gamma_Grid  =",FF[1/(2*Delta)],"\n"];
WriteString[FileName,"      REAL(DOUBLE), PARAMETER :: Gamma_Switch=",FF[switch],"\n"];
WriteString[FileName,"\n"];
Close[FileName];


FileName="GammaGrid_77.Inc";
Print[" Openned ",FileName];
OpenWrite[FileName];
WriteString[FileName,"      REAL*8 Gamma_Grid,Gamma_Switch \n"];
WriteString[FileName,"      PARAMETER (Gamma_Grid =",FF[1/(2*Delta)],") \n"];
WriteString[FileName,"      PARAMETER (Gamma_Switch=",FF[switch],")\n"];
WriteString[FileName,"\n"];
Close[FileName];

(* WRITE ASSYMPTOTIC COEFFICIENTS TO FILE *)

FileName="GammaAssymptotics.Inc";
Print[" Openned ",FileName];
OpenWrite[FileName];
WriteString[FileName,"      REAL(DOUBLE), PARAMETER :: GammAss(0:",ToString[NFunctions-1],")= (/& \n"];
Do[ 
WriteString[FileName,"                                ",FF[Abs[2 n - 1]!!/(2 (2)^n) Sqrt[Pi]],", & \n"];
   ,{n,0,NFunctions-2}];
n=NFunctions-1
WriteString[FileName,"                                ",FF[Abs[2 n - 1]!!/(2 (2)^n) Sqrt[Pi]],"/) \n"];
WriteString[FileName,"\n"];
Close[FileName];



(* WRITE ARRAY VALUES TO INDIVIDUAL FILES *)

<<BurnTables.m;


