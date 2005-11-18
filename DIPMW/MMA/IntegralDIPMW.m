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
(* INDEXING OF ARRAY ELEMENTS *)
MaxEll=2;
MaxOrder=3;
LBegin[L_]:=(L*(L+1)*(L+2))/6+1;
LMNDex[L_,M_,N_]:=LBegin[L+M+N]+N*(2*(L+M+N)-N+3)/2+M;
IJK[i_,j_,k_,O_] = 1+i+(O+1)*j+(O+1)*(O+1)*k;

(* LOAD OPTIMIZING AND FORMATING ROUTINES FOR MMA AND SET ASSOCIATED OPTIONS *)

Get["FixedNumberForm.m"];
Get["Format.m"];
Get["Optimize.m"];

SetOptions[FortranAssign,AssignOptimize->False,AssignPrecision->20,AssignMaxSize->500,AssignBreak->{180," & \n "},
           AssignTemporary->{W,Array},AssignTemporary->{W,Array},AssignToArray->{RhoSum,HGCo}];

SetOptions[OpenWrite, PageWidth -> 200];

(* PUT THE TRANSFORMATIONS TO FILE *)


FileName="DivRho_"
FileName=StringJoin[FileName,ToString[MaxEll],"_",ToString[MaxOrder],".Inc"];
Print[" Openned ",FileName];
OpenWrite[FileName];

Rho=0;
Expon = Exp[-Zeta*(RQx*RQx+RQy*RQy+RQz*RQz)];
LamX[0]=1;
LamY[0]=1;
LamZ[0]=1;
LamX[1]=TwoZeta*RQx;
LamY[1]=TwoZeta*RQy;
LamZ[1]=TwoZeta*RQz;
Do[
  LamX[l]=TwoZeta*(RQx*LamX[l-1]-(l-1)*LamX[l-2]);
  LamY[l]=TwoZeta*(RQy*LamY[l-1]-(l-1)*LamY[l-2]);
  LamZ[l]=TwoZeta*(RQz*LamZ[l-1]-(l-1)*LamZ[l-2]);
,{l,2,MaxEll}];
Do[
  Do[
    Do[
      LMN = LMNDex[l,m,n];
      Print["LMN = ",LMN,"  ",l,m,n];
      HGKet = LamX[l]*LamY[m]*LamZ[n];
      Rho = Rho+HGCo[LMN]*HGKet;
    ,{n,0,MaxEll-m-l}];
  ,{m,0,MaxEll-l}];
,{l,0,MaxEll}];
Rho = Expand[Factor[Simplify[Rho]]];
Rho = Rho*Expon;
Do[
  Do[
    Do[
      ijk = IJK[i,j,k,MaxOrder];
      Print["IJK = ",i,j,k];
      tmp = D[Rho,{RQx,i},{RQy,j},{RQz,k}];
      tmp = tmp/Expon;
      tmp = Factor[Simplify[Expand[tmp]]];
      tmp = tmp*Xpt+RhoSum[ijk];
      Write[FileName,FortranAssign[RXX,tmp,AssignReplace->{"RXX"->StringJoin["RhoSum(",ToString[ijk],") "],
                                                           "HGCo"->"Node%Co"}]];
    ,{i,0,MaxOrder}];
  ,{j,0,MaxOrder}];
,{k,0,MaxOrder}];

Close[FileName];          
Print[" Closed ",FileName];






