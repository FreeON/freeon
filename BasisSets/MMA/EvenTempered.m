
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
!    GENERATE UNCONTRACTED REGULARIZED-EVEN-TEMPERED BASIS SETS
!------------------------------------------------------------------------------
*)

MondoHome=Environment["MONDO_HOME"];
If[MondoHome==$FAILED,
   Print["COULD NOT FIND $MONDO_HOME! CHECK YOUR .cshrc "];
   Abort[];
  ];

FileName="RET.bas";
OpenWrite[FileName];
Print[" Oppened ",FileName];

Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];
Get[StringJoin[MondoHome,"/MMA/Format.m"]];
Get[StringJoin[MondoHome,"/MMA/Optimize.m"]];

(* GLOBAL FORMATING FUNCTIONS *)

FF[x_] := ToString[FixedNumberForm[SetPrecision[N[x,32],6], 6, 2]];

(* RET Parameters *)

<<RETParameters.m

(* Defs for exponents *)

beta[s_,H_,n_]:=Exp[Exp[b[s,H] Log[n]+bp[s,H]]];
alpha[s_,H_,n_]:=Exp[a[s,H]*Log[beta[s,H,n]-1]+ap[s,H]];
 
(* Defs for writing exponents *)
 
PunchSet[n_,A_,s_,FileName_]:=Block[{z},
                                    If[NumberQ[alpha[s,A,n]],
                                       Do[z=alpha[s,A,n] *beta[s,A,n]^k;
                                          WriteString[FileName," ",s,"  1  1.00 \n    ",FF[z],"     1.00  \n"];
                                         ,{k,1,n}];
                                      ,
                                       Print["Not even ",s," type functions defined for ",A," aborting... "];
                                       Abort[];
                                      ];
                                   ];
(* WRITE BASIS FOR H *)
 
WriteString[FileName," ",H,"   0\n"];
PunchSet[6,H,s,FileName];
(* polarization *)
(*
WriteString[FileName," P   1  1.00\n       0.75000000      1.00000000\n"];
*)
WriteString[FileName," ****\n"];
 
(* WRITE BASIS FOR 0 *)
 
WriteString[FileName," ",O,"   0\n"];
PunchSet[14,O,s,FileName];
PunchSet[6,O,p,FileName];
(*
WriteString[FileName," D   1  1.00\n       1.29200000      1.00000000\n"];
 *)
WriteString[FileName," ****\n"];

(* WRITE BASIS FOR C *)
 
WriteString[FileName," ",C,"   0\n"];
PunchSet[6,C,s,FileName];
PunchSet[3,C,p,FileName];
(*
WriteString[FileName," D   1  1.00\n       1.29200000      1.00000000\n"];
 *)
WriteString[FileName," ****\n"];
 
Close[FileName];
Print[" Closed ",FileName];










