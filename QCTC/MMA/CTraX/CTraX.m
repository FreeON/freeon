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
!    WRITE EXPLICIT CODE FOR MULTIPOLE GENERATION AND CONTRACTION
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

(* INDEXING OF SP ARRAY ELEMENTS: NOTE WE START FROM 0 NOT 1*)

SPDex[L_,M_]:=Binomial[L+1,2]+M;
LSP[L_]:=L*(L+3)/2;   

(* LOAD OPTIMIZING AND FORMATING ROUTINES FOR MMA AND SET ASSOCIATED OPTIONS *)

Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];
Get[StringJoin[MondoHome,"/MMA/Format.m"]];
Get[StringJoin[MondoHome,"/MMA/Optimize.m"]];

SetOptions[FortranAssign,AssignOptimize->False,AssignMaxSize->500,
                         AssignBreak->{132," & \n          "},AssignTemporary->{W,Array}]
SetOptions[Optimize,OptimizeVar)iable->{R,Sequence},OptimizeFunction->False]
SetOptions[OpenWrite, PageWidth -> 800];

(* GLOBAL FORMATING FUNCTIONS *)

FF[x_] := ToString[FixedNumberForm[SetPrecision[N[x,32],32], 16, 2]];
WS[String_]:=WriteString[FileName,"      ",String,"\n"];
WC[String_]:=WriteString[FileName,"!     ",String,"\n"];

(* LOAD MODULE FOR WRITING EXPLICIT SOURCE FOR COMPUTATION OF IRREGULAR HARMONICS *)

bkbreak=7;

<<IrRegular.m; 
<<Contract.m;

FileName="IrRegulars.Inc"; 
Print[" Openned ",FileName];
OpenWrite[FileName];
WS["SELECT CASE(Ell)"]
Do[
   ldx=ToString[l];
   WS[StringJoin["CASE(",ldx,")"]];
   IrRegular[FileName,l];
  ,{l,0,HGEll+SPEll}]
WS["END SELECT"]
Close[FileName];
Print[" Closed ",FileName];

Do[
Do[
   LCode=100*LHG+LSP;
   If[LHG+LSP>bkbreak, 
      FileName=StringJoin["CTraX",ToString[LCode],".F90"]; 
      Print[" Openned ",FileName];
      OpenWrite[FileName];
      WS[StringJoin["SUBROUTINE CTraX",ToString[LCode],"(Q)"]];
      WS[           "   USE PoleNodeType"];
      WS[           "   USE PoleGlobals"];
      WS[           "   USE ERIGlobals"];
      WS[           "   TYPE(PoleNode) :: Q"];
      Contract[FileName,LHG,LSP];
      WS[StringJoin["END SUBROUTINE CTraX",ToString[LCode]]];
      Close[FileName];
      Print[" Closed ",FileName];
    ];
,{LSP,0,SPEll}]
,{LHG,0,HGEll}];

FileName="CTraX.Inc"; 
Print[" Openned ",FileName];
OpenWrite[FileName];
WS["SELECT CASE(LCode)"]
Do[Do[
      ldx=ToString[LHG*100+LSP];
      WS[StringJoin["CASE(",ldx,")"]];
      If[LHG+LSP<=bkbreak,
         Contract[FileName,LHG,LSP];  
        ,
         WS[StringJoin["   CALL CTraX",ldx,"(Q)"]];
        ];
,{LSP,0,SPEll}]
,{LHG,0,HGEll}];
WS["END SELECT"]
Close[FileName];
Print[" Closed ",FileName];

