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
!    WRITE EXPLICIT CODE FOR CONVERSION OF HERMITE GAUSSIAN MULTIPOLE TO
!    SPHERICALS
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

(* GENERATE THE EXPLICIT TRANSFORMATIONS *)

f[R_,Ell_]:= 1/(R-Delta)*(Delta/R)^(Ell+1);

R:=Sqrt[r r];

HGEll=4;


Do[Do[ 
      Err=D[D[ f[R,Ell-L],{r,M}],{r,L}];
      MACErr[M,L]= Simplify[PowerExpand[Simplify[Together[Err]]]];
,{M,0,HGEll+1}],{L,0,HGEll}];

Do[Do[ 
      error[M,L]=Sum[Sum[ABS[BraCo[m] KetCo[l] MACErr[m,l]],{m,0,M}],{l,0,L}];
,{M,0,HGEll+1}],{L,0,HGEll}];

(* LOAD OPTIMIZING AND FORMATING ROUTINES FOR MMA AND SET ASSOCIATED OPTIONS *)

Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];
Get[StringJoin[MondoHome,"/MMA/Format.m"]];
Get[StringJoin[MondoHome,"/MMA/Optimize.m"]];



SetOptions[Optimize,OptimizePower->Binary,OptimizeFunction->True,OptimizeVariable->{o,Array}];

SetOptions[FortranAssign,AssignOptimize->True,AssignMaxSize->400,AssignPrecision->Infinity,AssignBreak->{132," & \n          "},AssignTemporary->{w,Array}];

SetOptions[OpenWrite, PageWidth -> 200];

(* PUT THE TRANSFORMATIONS TO FILE *)

FileName="MACErrBnd4.Inc";
Print[" Openned ",FileName];
OpenWrite[FileName];
Do[Do[
      WriteString[FileName,"     CASE(",100*M+L,") \n"];
      LHSList={MACError};
      RHSList={error[M,L]};
      Write[FileName,FortranAssign[Evaluate[LHSList],RHSList,AssignReplace->{" "->"","Delta"->"Q%MAC%Delta","Ell"->"Q%Pole%Ell","KetCo"->"Q%MAC%O","BraCo"->"QC%MAC%O"}]];
     ,{M,0,HGEll+1}];
  ,{L,0,HGEll}];

Close[FileName];          
Print[" Closed ",FileName];
