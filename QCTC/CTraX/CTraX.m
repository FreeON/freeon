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

MondoHome=Environment["FREEON_HOME"];
If[MondoHome==$FAILED,
   Print["COULD NOT FIND $FREEON_HOME! CHECK YOUR .cshrc "];
   Abort[];
  ];
EllFile = StringJoin[MondoHome,"/Includes/Ell.m"];
Get[EllFile];
(* Override of expansion length *)

(* INDEXING OF SP ARRAY ELEMENTS: NOTE WE START FROM 0 NOT 1*)

SPDex[L_,M_]:=Binomial[L+1,2]+M;
LSP[L_]:=L*(L+3)/2;   

(* LOAD OPTIMIZING AND FORMATING ROUTINES FOR MMA AND SET ASSOCIATED OPTIONS *)

Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];
Get[StringJoin[MondoHome,"/MMA/Format.m"]];
Get[StringJoin[MondoHome,"/MMA/Optimize.m"]];

SetOptions[FortranAssign,AssignOptimize->False,AssignMaxSize->500,AssignPrecision -> Infinity,
                         AssignBreak->{132," & \n          "},AssignTemporary->{W,Array}]
SetOptions[Optimize,OptimizeVariable->{R,Sequence},OptimizeFunction->False]
SetOptions[OpenWrite, PageWidth -> 800];

(* GLOBAL FORMATING FUNCTIONS *)

FF[x_] := ToString[FixedNumberForm[SetPrecision[N[x,32],32], 16, 2]];
FF8[x_] := ToString[FixedNumberForm[SetPrecision[N[x,32],32], 6, 2]];
WS[String_]:=WriteString[FileName,"      ",String,"\n"];
WC[String_]:=WriteString[FileName,"!     ",String,"\n"];

(* LOAD MODULE FOR WRITING EXPLICIT SOURCE FOR COMPUTATION OF IRREGULAR HARMONICS *)

<<IrRegular.m; 
<<Contract.m;

Do[
   FileName=StringJoin["CTraX",ToString[LLSP],".F"]; 
   Print[" Openned ",FileName];
   OpenWrite[FileName];
   Do[
      LCode=100*LHG+LLSP;
      SubRoutineList=Contract[FileName,LCode,LHG,LLSP]; 
      IrRegular[FileName,LHG+LLSP,LCode];
      WS[StringJoin["SUBROUTINE CTraX",ToString[LCode],"(PQx,PQy,PQz,Cq,Sq,Cp,Sp)"]];
      WS[           "   IMPLICIT NONE"];
      WS[           "   REAL*8 PQx,PQy,PQz"];
      LDXP=ToString[LSP[LHG]];
      LDXQ=ToString[LSP[LLSP]];
      WS[StringJoin["   REAL*8 Cq(0:",LDXQ,"),Sq(0:",LDXQ,"),Cp(0:",LDXP,"),Sp(0:",LDXP,")"]];
      lensp=ToString[LSP[LHG+LLSP]];
      WS[StringJoin["   REAL*8 Cpq(0:",lensp,"),Spq(0:",lensp,")"]];
      WS[StringJoin["   CALL IrRegular",ToString[LCode],"(PQx,PQy,PQz,Cpq,Spq)"]];      
      Do[WS[SubRoutineList[[i]]];,{i,Length[SubRoutineList]}];
      WS[StringJoin["RETURN"]]
      WS[StringJoin["END"]]
   ,{LHG,0,Min[HGEll,LLSP-2]}];
   Close[FileName];
   Print[" Closed ",FileName];
,{LLSP,4,HGEll+14,2}]

(* NOW WRITE THE IRREGULARS FORTRAN77 FILE  *)

(* NOW WRITE THE F90 INCLUDES FILE FOR CALLS TO CTRAX *)

FileName="CTraX.Inc"; 
Print[" Openned ",FileName];
OpenWrite[FileName];
Print["hgell = ",HGEll];
Do[Do[
      ldx=ToString[LHG*100+LLSP];
      WS[StringJoin["CASE(",ldx,")"]];
      WS["   DO N=1,NFar"];
      WS["      Q=>Far(N)%P"];
      WS["      PQx=QC%Prim%P(1)-Q%Pole%Center(1)"];
      WS["      PQy=QC%Prim%P(2)-Q%Pole%Center(2)"];
      WS["      PQz=QC%Prim%P(3)-Q%Pole%Center(3)"];
      WS[StringJoin["      CALL CTraX",ldx,"(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))"]];
      WS["   ENDDO"];   
	 ,{LHG,0,Min[HGEll,LLSP-2]}];
      ,{LLSP,4,HGEll+14,2}];
   WS["CASE DEFAULT"];
WS["   STOP ' Unsuported expansion orders in TreeWalk.F90'"];
   WS["END SELECT"];
Close[FileName];
Print[" Closed ",FileName];

