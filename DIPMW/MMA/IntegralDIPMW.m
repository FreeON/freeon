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
!    Author: C. J. Tymczak 
!    Inregrate Exp[-bet*x]*x^n from A to AB to B
!------------------------------------------------------------------------------
*)

(*---*)
NAsy      = 12;
MaxMaxEll = 0;
Dorder=0;

(* Initialize *)

(* Dorder=0 *)
WveL[x_,0,0] = 1+x;
WveR[x_,0,0] = 1-x;

(* Dorder=1 *)
WveL[x_,1,0] = 1-3*x*x-2*x^3;
WveR[x_,1,0] = 1-3*x*x+2*x^3;
WveL[x_,1,1] = x*(1+2*x+x*x);
WveR[x_,1,1] = x*(1-2*x+x*x);

(* Dorder=2 *)
WveL[x_,2,0] = (1+10*x^3+15*x^4+6*x^5);
WveR[x_,2,0] = (1-10*x^3+15*x^4-6*x^5);
WveL[x_,2,1] = x*(1-6*x^2-8*x^3-3*x^4);
WveR[x_,2,1] = x*(1-6*x^2+8*x^3-3*x^4);
WveL[x_,2,2] = x*x*(1/2+(3/2)*x+(3/2)*x^2+(1/2)*x^3);
WveR[x_,2,2] = x*x*(1/2-(3/2)*x+(3/2)*x^2-(1/2)*x^3);


(* Dorder=3 *)
WveL[x_,3,0] = (1 - 35*x^4 - 84*x^5 - 70*x^6 - 20*x^7);
WveR[x_,3,0] = (1 - 35*x^4 + 84*x^5 - 70*x^6 + 20*x^7);
WveL[x_,3,1] = x*(1 + 20*x^3 + 45*x^4 + 36*x^5 + 10*x^6);
WveR[x_,3,1] = x*(1 - 20*x^3 + 45*x^4 - 36*x^5 + 10*x^6);
WveL[x_,3,2] = x*x*((1/2) - 5*x^2 - 10*x^3 - (15/2)*x^4 - 2*x^5);
WveR[x_,3,2] = x*x*((1/2) - 5*x^2 + 10*x^3 - (15/2)*x^4 + 2*x^5);
WveL[x_,3,3] = x*x*x*((1/6) + (2/3)*x + x^2 + (2/3)*x^3 + (1/6)*x^4);
WveR[x_,3,3] = x*x*x*((1/6) - (2/3)*x + x^2 - (2/3)*x^3 + (1/6)*x^4);


(*--Ell stuff--*)
FFun[x_,0]=1;
FFun[x_,1]=Expand[2*(RZeta1-x)];
Do[
   FFun[x_,L] = Expand[2*( (RZeta1-x)*FFun[x,L-1] - (L-1)*FFun[x,L-2] )] ;
,{L,2,MaxMaxEll}];

(* LOAD OPTIMIZING AND FORMATING ROUTINES FOR MMA AND SET ASSOCIATED OPTIONS *)

Get["FixedNumberForm.m"];
Get["Format.m"];
Get["Optimize.m"];

SetOptions[FortranAssign,AssignOptimize->False,AssignPrecision->16,AssignMaxSize->500,AssignBreak->{180," & \n "}];
SetOptions[OpenWrite, PageWidth -> 200];

(* PUT THE TRANSFORMATIONS TO FILE *)

Do[

FileName="IntWavelet_";
FileName=StringJoin[FileName,ToString[Dorder],"_",ToString[MaxEll],".Inc"];
Print[" Openned ",FileName];
OpenWrite[FileName];

S1 = StringJoin["IntD",ToString[Dorder],"L"];
If[Dorder==0,
  MinPX = 2;
  SS = StringJoin["    ",S1,"0 = One/(2.D0*SqZeta*XZeta1)"];
];
If[Dorder==1,
  MinPX = 4;
  Do[
    SS = StringJoin["    ","XZeta",ToString[P]," = XZeta",ToString[P-1],"*XZeta1"];
    WriteString[FileName,SS,"\n"];
  ,{P,2,MinPX-1}];
  SS = StringJoin["    ",S1,"0 = One/(8.D0*SqZeta*XZeta3)"];
];
If[Dorder==2,
  MinPX = 6;
  Do[
    SS = StringJoin["    ","XZeta",ToString[P]," = XZeta",ToString[P-1],"*XZeta1"];
    WriteString[FileName,SS,"\n"];
  ,{P,2,MinPX-1}];
  SS = StringJoin["    ",S1,"0 = One/(32.D0*SqZeta*XZeta5)"];
];
If[Dorder==3,
  MinPX = 8;
  Do[
    SS = StringJoin["    ","XZeta",ToString[P]," = XZeta",ToString[P-1],"*XZeta1"];
    WriteString[FileName,SS,"\n"];
  ,{P,2,MinPX-1}];
  SS = StringJoin["    ",S1,"0 = One/(128.D0*SqZeta*XZeta7)"];
];
WriteString[FileName,SS,"\n"];
Do[
  SS = StringJoin["    ",S1,ToString[Ell]," = ",S1,ToString[Ell-1],"*SqZeta\n"];
  WriteString[FileName,SS];
,{Ell,1,MaxEll}];
WriteString[FileName,"    SELECT CASE(SType)\n"];

Do[
If[SType==-1,
   WriteString[FileName,"    CASE(-1)\n"]
];
If[SType==0,
   WriteString[FileName,"    CASE(0)\n"]
];
If[SType==1,
   WriteString[FileName,"    CASE(1)\n"]
];

Do[
   Do[
      Print["Ell=",Ell," DD=",DD];
      If[SType==-1,
         Tmp1 = 0;
         Tmp2 = ((2*XZeta1)^(2*Dorder+1))*Integrate[FFun[x,Ell]*WveR[x/XZeta1,Dorder,DD]*Exp[-(x-RZeta1)^2],{x,0,XZeta1}];
      ];
      If[SType==1,
         Tmp1= ((2*XZeta1)^(2*Dorder+1))*Integrate[FFun[x,Ell]*WveL[x/XZeta1,Dorder,DD]*Exp[-(x-RZeta1)^2],{x,-XZeta1,0}];
         Tmp2= 0;
      ];	
      If[SType==0,
         Tmp1= ((2*XZeta1)^(2*Dorder+1))*Integrate[FFun[x,Ell]*WveL[x/XZeta1,Dorder,DD]*Exp[-(x-RZeta1)^2],{x,-XZeta1,0}];
         Tmp2= ((2*XZeta1)^(2*Dorder+1))*Integrate[FFun[x,Ell]*WveR[x/XZeta1,Dorder,DD]*Exp[-(x-RZeta1)^2],{x,0, XZeta1}];
      ];
      Tmp = Expand[Tmp1+Tmp2];
      Tmp = Tmp /. Erf[RZeta1-XZeta1]->ErfRmX0;
      Tmp = Tmp /. Erf[RZeta1+XZeta1]->ErfRpX0;
      Tmp = Tmp /. Erf[RZeta1]->Erf0;
      Tmp = Tmp /. Exp[-RZeta1^2]->Exp0;
      Tmp = Tmp /. Exp[-(RZeta1-XZeta1)^2]->ExpRmX0;
      Tmp = Tmp /. Exp[-(RZeta1+XZeta1)^2]->ExpRpX0;
      Tmp = Tmp /. Exp[4*RZeta1*XZeta1-(RZeta1+XZeta1)^2] -> ExpRmX0;
      Tmp = Tmp /. Exp[-(RZeta1+XZeta1)^2+XZeta1*(2*RZeta1+XZeta1)] -> Exp0;
      Tmp = Tmp /. Exp[-RZeta1^2+2*RZeta1*XZeta1-XZeta1^2] -> ExpRmX0;
      Tmp = Tmp /. Exp[-2*RZeta1^2 + (RZeta1 - XZeta1)^2 + 2*RZeta1*XZeta1 - XZeta1^2] -> Exp0;
      Tmp = Tmp /. Exp[-RZeta1^2 - 2*RZeta1*XZeta1 - 2*XZeta1^2 + 2*XZeta1*(RZeta1 +XZeta1)] -> Exp0;
      Tmp = Tmp /. Exp[-RZeta1^2 - 2*RZeta1*XZeta1 - XZeta1^2]->ExpRpX0;
      Tmp = Expand[Tmp];
      Print["Int = ",Tmp];
      IntInt[DD,Ell] = Tmp;
   ,{Ell,0,MaxEll}];
,{DD,0,Dorder}];

(*    Determine the Maxium RZeta and X Zeta Power *)

MaxPX = 1;
MaxPR = 1;
Do[
  Do[
      MaxPX = Max[MaxPX,Exponent[IntInt[DD,Ell],XZeta1]];
      MaxPR = Max[MaxPR,Exponent[IntInt[DD,Ell],RZeta1]];
   ,{Ell,0,MaxEll}];
,{DD,0,Dorder}];
Do[
  SS = StringJoin["        ","XZeta",ToString[P]," = XZeta",ToString[P-1],"*XZeta1"];
  WriteString[FileName,SS,"\n"];
,{P,MinPX,MaxPX}];
Do[ 
  SS = StringJoin["        ","RZeta",ToString[P]," = RZeta",ToString[P-1],"*RZeta1"];
  WriteString[FileName,SS,"\n"];
,{P,2,MaxPR}];

(*  *) 

Do[
  Do[
      Print["Ell=",Ell," DD=",DD];
      Fac = ToExpression[StringJoin["IntD",ToString[Dorder],"L",ToString[Ell]]];

      Tmp = IntInt[DD,Ell];
      Tmp = Tmp /. RZeta1^2 -> RZeta2;
      Tmp = Tmp /. XZeta1^2 -> XZeta2;
      Tmp = Tmp /. RZeta1^3 -> RZeta3;
      Tmp = Tmp /. XZeta1^3 -> XZeta3;
      Tmp = Tmp /. RZeta1^4 -> RZeta4;
      Tmp = Tmp /. XZeta1^4 -> XZeta4;
      Tmp = Tmp /. RZeta1^5 -> RZeta5;
      Tmp = Tmp /. XZeta1^5 -> XZeta5;
      Tmp = Tmp /. RZeta1^6 -> RZeta6;
      Tmp = Tmp /. XZeta1^6 -> XZeta6;
      Tmp = Tmp /. RZeta1^7 -> RZeta7;
      Tmp = Tmp /. XZeta1^7 -> XZeta7;
      Tmp = Tmp /. RZeta1^8 -> RZeta8;
      Tmp = Tmp /. XZeta1^8 -> XZeta8;
      Tmp = Tmp /. RZeta1^9 -> RZeta9;
      Tmp = Tmp /. XZeta1^9 -> XZeta9;
      Tmp = Tmp /. RZeta1^10 -> RZeta10;
      Tmp = Tmp /. XZeta1^10 -> XZeta10;
      Tmp = Tmp /. RZeta1^11 -> RZeta11;
      Tmp = Tmp /. XZeta1^11 -> XZeta11;
      Tmp = Tmp /. RZeta1^12 -> RZeta12;
      Tmp = Tmp /. XZeta1^12 -> XZeta12;
      Tmp = Tmp /. RZeta1^13 -> RZeta13;
      Tmp = Tmp /. XZeta1^13 -> XZeta13;
      Tmp = Tmp /. RZeta1^14 -> RZeta14;
      Tmp = Tmp /. XZeta1^14 -> XZeta14;
      Tmp = Tmp /. RZeta1^15 -> RZeta15;
      Tmp = Tmp /. XZeta1^15 -> XZeta15;
      Tmp = Tmp /. RZeta1^16 -> RZeta16;
      Tmp = Tmp /. XZeta1^16 -> XZeta16;
      Tmp = Tmp /. RZeta1^17 -> RZeta17;
      Tmp = Tmp /. XZeta1^17 -> XZeta17;
      Tmp = Tmp /. RZeta1^18 -> RZeta18;
      Tmp = Tmp /. XZeta1^18 -> XZeta18;
      Tmp = Tmp /. RZeta1^19 -> RZeta19;
      Tmp = Tmp /. XZeta1^19 -> XZeta19;
      Tmp = Tmp /. RZeta1^20 -> RZeta20;
      Tmp = Tmp /. XZeta1^20 -> XZeta20;

      TmpFun=Fac*Tmp;
      Write[FileName,FortranAssign[WV,TmpFun,AssignReplace->{"WV"->StringJoin["WV(",ToString[DD],",",ToString[Ell],")"]} ]];

   ,{Ell,0,MaxEll}];
,{DD,0,Dorder}];



,{SType,-1,1}];
WriteString[FileName,"    END SELECT\n"]

Close[FileName];          
Print[" Closed ",FileName];

,{MaxEll,0,MaxMaxEll}];


(* PUT THE TRANSFORMATIONS TO FILE *)

Do[

FileName="IntWavelet_Asy_";
FileName=StringJoin[FileName,ToString[Dorder],"_",ToString[MaxEll],".Inc"];
Print[" Openned ",FileName];
OpenWrite[FileName];

S1 = StringJoin["IntD",ToString[Dorder],"L"];
SS = StringJoin["    ",S1,"0 = Exp0/SqZeta"];
WriteString[FileName,SS,"\n"];
Do[
  SS = StringJoin["    ",S1,ToString[Ell]," = ",S1,ToString[Ell-1],"*SqZeta\n"];
  WriteString[FileName,SS];
,{Ell,1,MaxEll}];

WriteString[FileName,"    SELECT CASE(SType)\n"]
Do[
If[SType==-1,
   WriteString[FileName,"    CASE(-1)\n"]
];
If[SType==0,
   WriteString[FileName,"    CASE(0)\n"]
];
If[SType==1,
   WriteString[FileName,"    CASE(1)\n"]
];
Do[
   Do[
      Print["Ell=",Ell," DD=",DD];
      FFunTmp = Expand[Normal[Series[FFun[x,Ell]*Exp[-(x-RZeta1)^2],{x,0,NAsy}]]/Exp[-RZeta1*RZeta1]];
      If[SType==-1,
         Tmp1= 0;
         Tmp2= (2*XZeta1)^(2*Dorder+1)*Integrate[FFunTmp*WveR[x/XZeta1,Dorder,DD],{x,0, XZeta1}];
      ];
      If[SType==1,
         Tmp1= (2*XZeta1)^(2*Dorder+1)*Integrate[FFunTmp*WveL[x/XZeta1,Dorder,DD],{x,-XZeta1,0}];
         Tmp2= 0;
      ];	
      If[SType==0,
         Tmp1= (2*XZeta1)^(2*Dorder+1)*Integrate[FFunTmp*WveL[x/XZeta1,Dorder,DD],{x,-XZeta1,0}];
         Tmp2= (2*XZeta1)^(2*Dorder+1)*Integrate[FFunTmp*WveR[x/XZeta1,Dorder,DD],{x,0, XZeta1}];
      ];
      Tmp = Expand[Tmp1+Tmp2];
      Print["Int = ",Tmp];
      IntInt[DD,Ell] = Tmp;

      NN = Exponent[Tmp,XZeta1];
      Error[RZeta1_] = Coefficient[Tmp,XZeta1^NN];
      TmpTmp = Maximize[Error[x]*Exp[-x*x],x];
      MaxCoef = N[TmpTmp[[1]],32];
      MinDX = (10^-14/MaxCoef)^(1/(NN+1));
      Print["MinDX = ",MinDX];

   ,{Ell,0,MaxEll}];
,{DD,0,Dorder}];

(*    Determine the Maxium RZeta and X Zeta Power *)

MaxPX = 1;
MaxPR = 1;
Do[
  Do[
      MaxPX = Max[MaxPX,Exponent[IntInt[DD,Ell],XZeta1]];
      MaxPR = Max[MaxPR,Exponent[IntInt[DD,Ell],RZeta1]];
   ,{Ell,0,MaxEll}];
,{DD,0,Dorder}];
Do[
  SS = StringJoin["        ","XZeta",ToString[P]," = XZeta",ToString[P-1],"*XZeta1"];
  WriteString[FileName,SS,"\n"];
,{P,2,MaxPX}];
Do[ 
  SS = StringJoin["        ","RZeta",ToString[P]," = RZeta",ToString[P-1],"*RZeta1"];
  WriteString[FileName,SS,"\n"];
,{P,2,MaxPR}];

Do[
  Do[

      Fac = ToExpression[StringJoin["IntD",ToString[Dorder],"L",ToString[Ell]]];

      Tmp = Tmp /. RZeta1^2 -> RZeta2;
      Tmp = Tmp /. XZeta1^2 -> XZeta2;
      Tmp = Tmp /. RZeta1^3 -> RZeta3;
      Tmp = Tmp /. XZeta1^3 -> XZeta3;
      Tmp = Tmp /. RZeta1^4 -> RZeta4;
      Tmp = Tmp /. XZeta1^4 -> XZeta4;
      Tmp = Tmp /. RZeta1^5 -> RZeta5;
      Tmp = Tmp /. XZeta1^5 -> XZeta5;
      Tmp = Tmp /. RZeta1^6 -> RZeta6;
      Tmp = Tmp /. XZeta1^6 -> XZeta6;
      Tmp = Tmp /. RZeta1^7 -> RZeta7;
      Tmp = Tmp /. XZeta1^7 -> XZeta7;
      Tmp = Tmp /. RZeta1^8 -> RZeta8;
      Tmp = Tmp /. XZeta1^8 -> XZeta8;
      Tmp = Tmp /. RZeta1^9 -> RZeta9;
      Tmp = Tmp /. XZeta1^9 -> XZeta9;
      Tmp = Tmp /. RZeta1^10 -> RZeta10;
      Tmp = Tmp /. XZeta1^10 -> XZeta10;
      Tmp = Tmp /. RZeta1^11 -> RZeta11;
      Tmp = Tmp /. XZeta1^11 -> XZeta11;
      Tmp = Tmp /. RZeta1^12 -> RZeta12;
      Tmp = Tmp /. XZeta1^12 -> XZeta12;
      Tmp = Tmp /. RZeta1^13 -> RZeta13;
      Tmp = Tmp /. XZeta1^13 -> XZeta13;
      Tmp = Tmp /. RZeta1^14 -> RZeta14;
      Tmp = Tmp /. XZeta1^14 -> XZeta14;
      Tmp = Tmp /. RZeta1^15 -> RZeta15;
      Tmp = Tmp /. XZeta1^15 -> XZeta15;
      Tmp = Tmp /. RZeta1^16 -> RZeta16;
      Tmp = Tmp /. XZeta1^16 -> XZeta16;
      Tmp = Tmp /. RZeta1^17 -> RZeta17;
      Tmp = Tmp /. XZeta1^17 -> XZeta17;
      Tmp = Tmp /. RZeta1^18 -> RZeta18;
      Tmp = Tmp /. XZeta1^18 -> XZeta18;
      Tmp = Tmp /. RZeta1^19 -> RZeta19;
      Tmp = Tmp /. XZeta1^19 -> XZeta19;
      Tmp = Tmp /. RZeta1^20 -> RZeta20;
      Tmp = Tmp /. XZeta1^20 -> XZeta20;


      TmpFun=Fac*Tmp;
      Write[FileName,FortranAssign[WV,TmpFun,AssignReplace->{"WV"->StringJoin["WV(",ToString[DD],",",ToString[Ell],")"]} ]];


   ,{Ell,0,MaxEll}];
,{DD,0,Dorder}];
,{SType,-1,1}];
WriteString[FileName,"    END SELECT\n"]

Close[FileName];          
Print[" Closed ",FileName];

,{MaxEll,0,MaxMaxEll}];







