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
!    WRITE EXPLICIT CODE FOR COMPUTATION OF THE ECP ANGULAR INTEGRALS OF THE 
!    FIRST KIND, OMEGA(LAMBDA,I,J,K)
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

Ell=(HGEll-1)/2;

Z[l_,m_,x_,y_,z_]:=Module[{prefactor,em},
    em=Abs[m];
    prefactor=Sqrt[(2 l+1)(l-em)!/(2  Pi(l+em)!)];
    dp=D[LegendreP[l,u],{u,em}]/.u\[Rule]z;
    If[m>0, midfact=prefactor*(1-I+(1+I))/2  ((x+ I y)^em+(x-I  y)^em)/2];
    If[m\[Equal]0,midfact=Sqrt[(2 l+1)/(4  Pi)]] ;
    If[m<0, midfact=prefactor*(1-I-(1+I))/2  ((x+ I y)^em-(x-I  y)^em)/2];
    midfact*dp]

Do[Do[
    Zxyz[lambda,mu]=0;
    Do[Do[Do[
              v=(1/(l! m! n!)) D[ Z[lambda,mu,x,y,z],{x,l},{y,m},{z,n}];
              v=v/.{x->0,y->0,z->0}; 
              If[v!=0,Zxyz[lambda,mu]=Zxyz[lambda,mu]+v*X^l Y^m Z^n;];
          ,{l,0,lambda-n-m}],{m,0,lambda-n}],{n,0,lambda}];
,{mu,-lambda,lambda}],{lambda,0,Ell}]

EvenZero[i_]:=If[Mod[i,2]==0,1,0];
FactorialList={X^a_*Y^b_*Z^c_:>EvenZero[a]EvenZero[b]EvenZero[c](a-1)!!(b-1)!!(c-1)!!/(a+b+c+1)!!};

(* HERE IS THE ANGULAR INTEGRALS *)

omega[lambda_,i_,j_,k_]:=Module[{sum,kpart,apart},
      sum=0;
       Do[
          kpart=Zxyz[lambda,mu]/.{X->Kx,Y->Ky,Z->Kz};
          apart=4 Pi Expand[Apart[(Zxyz[lambda,mu]*X^Eye Y^Jay Z^Kay)]];
          apart=apart/.FactorialList;
          apart=apart/.{Eye->i,Jay->j,Kay->k};
          sum=sum+kpart;
         ,{mu,-lambda,lambda}];
       N[Factor[Expand[Simplify[sum]]]]]

(* LOAD OPTIMIZING AND FORMATING ROUTINES FOR MMA AND SET ASSOCIATED OPTIONS *)

Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];
Get[StringJoin[MondoHome,"/MMA/Format.m"]];
Get[StringJoin[MondoHome,"/MMA/Optimize.m"]];

SetOptions[FortranAssign,AssignOptimize->False,AssignMaxSize->300,AssignBreak->{132," & \n          "},
                         AssignTemporary->{W,Array},AssignTemporary->{W,Array}];

SetOptions[OpenWrite, PageWidth -> 200];

(* PUT THE TRANSFORMATIONS TO FILE *)

FileName="Omega1.Inc";
Print[" Openned ",FileName];
OpenWrite[FileName];

Do[WriteString[FileName,"     CASE(",lambda,") \n"];
   Write[FileName,FortranAssign[Omega1,omega[lambda,0,0,0]]];
,{lambda,0,Ell}]

Close[FileName];          
Print[" Closed ",FileName];


