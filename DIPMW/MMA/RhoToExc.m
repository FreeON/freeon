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

MaxOrder=3;
IJK[i_,j_,k_,O_] = 1+i+(O+1)*j+(O+1)*(O+1)*k;

(* LOAD OPTIMIZING AND FORMATING ROUTINES FOR MMA AND SET ASSOCIATED OPTIONS *)

Get["FixedNumberForm.m"];
Get["Format.m"];
Get["Optimize.m"];

SetOptions[FortranAssign,AssignOptimize->False,AssignPrecision->20,AssignMaxSize->500,AssignBreak->{160," & \n "},
           AssignTemporary->{W,Array},AssignTemporary->{W,Array},AssignToArray->{RhoXX}];
SetOptions[OpenWrite, PageWidth -> 200];

(* PUT THE TRANSFORMATIONS TO FILE *)


FileName="RhoToExc_"
FileName=StringJoin[FileName,ToString[MaxOrder],".Inc"];
Print[" Openned ",FileName];
OpenWrite[FileName];

RhoXX[1]=0;
Subs={};
Do[
  Do[
    Do[
      ijk = IJK[i,j,k,MaxOrder];
      Elem  = D[rho[x,y,z],{x,i},{y,j},{z,k}]->RhoXX[ijk];
      Tmp   = Append[Subs,Elem];
      Subs  = Tmp;
    ,{i,0,MaxOrder}];
  ,{j,0,MaxOrder}];
,{k,0,MaxOrder}];
Do[
  Do[
    Do[
      ijk = IJK[i,j,k,MaxOrder];
      Print["IJK = ",i,j,k];
      Do[  
         Tmp = D[rho[x,y,z]^Q/Q!,{x,i},{y,j},{z,k}];
         Tmp = Tmp /. Subs;
         Tmp = Factor[Simplify[Tmp]];
	 Write[FileName,FortranAssign[RT,Tmp,AssignReplace->{"RT"->StringJoin["RhoTmp",ToString[Q]]}]];
      ,{Q,1,i+j+k}];

      If[i+j+k==0,ExcTmp = d0Exc];
      If[i+j+k==1,ExcTmp = d1Exc*RhoTmp1];
      If[i+j+k==2,ExcTmp = d1Exc*RhoTmp1+d2Exc*RhoTmp2];
      If[i+j+k==3,ExcTmp = d1Exc*RhoTmp1+d2Exc*RhoTmp2+d3Exc*RhoTmp3];
      If[i+j+k==4,ExcTmp = d1Exc*RhoTmp1+d2Exc*RhoTmp2+d3Exc*RhoTmp3+d4Exc*RhoTmp4];
      If[i+j+k==5,ExcTmp = d1Exc*RhoTmp1+d2Exc*RhoTmp2+d3Exc*RhoTmp3+d4Exc*RhoTmp4+d5Exc*RhoTmp5];
      If[i+j+k==6,ExcTmp = d1Exc*RhoTmp1+d2Exc*RhoTmp2+d3Exc*RhoTmp3+d4Exc*RhoTmp4+d5Exc*RhoTmp5+d6Exc*RhoTmp6];
      If[i+j+k==7,ExcTmp = d1Exc*RhoTmp1+d2Exc*RhoTmp2+d3Exc*RhoTmp3+d4Exc*RhoTmp4+d5Exc*RhoTmp5+d6Exc*RhoTmp6+d7Exc*RhoTmp7];
      If[i+j+k==8,ExcTmp = d1Exc*RhoTmp1+d2Exc*RhoTmp2+d3Exc*RhoTmp3+d4Exc*RhoTmp4+d5Exc*RhoTmp5+d6Exc*RhoTmp6+d7Exc*RhoTmp7+d8Exc*RhoTmp8];
      If[i+j+k==9,ExcTmp = d1Exc*RhoTmp1+d2Exc*RhoTmp2+d3Exc*RhoTmp3+d4Exc*RhoTmp4+d5Exc*RhoTmp5+d6Exc*RhoTmp6+d7Exc*RhoTmp7+d8Exc*RhoTmp8+d9Exc*RhoTmp9];

      If[i+j+k==0,VxcTmp = d0Vxc];
      If[i+j+k==1,VxcTmp = d1Vxc*RhoTmp1];
      If[i+j+k==2,VxcTmp = d1Vxc*RhoTmp1+d2Vxc*RhoTmp2];
      If[i+j+k==3,VxcTmp = d1Vxc*RhoTmp1+d2Vxc*RhoTmp2+d3Vxc*RhoTmp3];
      If[i+j+k==4,VxcTmp = d1Vxc*RhoTmp1+d2Vxc*RhoTmp2+d3Vxc*RhoTmp3+d4Vxc*RhoTmp4];
      If[i+j+k==5,VxcTmp = d1Vxc*RhoTmp1+d2Vxc*RhoTmp2+d3Vxc*RhoTmp3+d4Vxc*RhoTmp4+d5Vxc*RhoTmp5];
      If[i+j+k==6,VxcTmp = d1Vxc*RhoTmp1+d2Vxc*RhoTmp2+d3Vxc*RhoTmp3+d4Vxc*RhoTmp4+d5Vxc*RhoTmp5+d6Vxc*RhoTmp6];
      If[i+j+k==7,VxcTmp = d1Vxc*RhoTmp1+d2Vxc*RhoTmp2+d3Vxc*RhoTmp3+d4Vxc*RhoTmp4+d5Vxc*RhoTmp5+d6Vxc*RhoTmp6+d7Vxc*RhoTmp7];
      If[i+j+k==8,VxcTmp = d1Vxc*RhoTmp1+d2Vxc*RhoTmp2+d3Vxc*RhoTmp3+d4Vxc*RhoTmp4+d5Vxc*RhoTmp5+d6Vxc*RhoTmp6+d7Vxc*RhoTmp7+d8Vxc*RhoTmp8];
      If[i+j+k==9,VxcTmp = d1Vxc*RhoTmp1+d2Vxc*RhoTmp2+d3Vxc*RhoTmp3+d4Vxc*RhoTmp4+d5Vxc*RhoTmp5+d6Vxc*RhoTmp6+d7Vxc*RhoTmp7+d8Vxc*RhoTmp8+d9Vxc*RhoTmp9];


      Write[FileName,FortranAssign[EEE,ExcTmp,AssignReplace->{"EEE"->StringJoin["Exc(",ToString[ijk],") "]}]];
      Write[FileName,FortranAssign[VVV,VxcTmp,AssignReplace->{"VVV"->StringJoin["Vxc(",ToString[ijk],") "]}]];

    ,{i,0,MaxOrder}];
  ,{j,0,MaxOrder}];
,{k,0,MaxOrder}];

Close[FileName];          
Print[" Closed ",FileName];






