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
!    WRITE EXPLICIT CODE FOR CONTRACTION OF SP C AND S TENSORS 
!------------------------------------------------------------------------------
*)
Contract[FileName_,EllP_,EllQ_]:=Block[{},

Print[" EllP = ",EllP," EllQ = ",EllQ];

SRules={Cp[ m_,l_]:>(-1)^(m) Cp[Abs[m],l]/;m<0,Sp[m_,l_]:>(-1)^(m+1)   Sp[Abs[m],l]/;m<0,
        Cpq[m_,l_]:>(-1)^(m)Cpq[Abs[m],l]/;m<0,Spq[m_,l_]:> (-1)^(m+1) Spq[Abs[m],l]/;m<0,
        Cq[ m_,l_]:>(-1)^(m) Cq[Abs[m],l]/;m<0,Sq[m_,l_]:> (-1)^(m+1)  Sq[Abs[m],l]/;m<0,
        Sp[ m_,l_]:>0/;m==0,Spq[m_,l_]:>0/;m==0,Sq[m_,l_]:>0/;m==0};

TRules={Cq[m_,l_]:>Cq[SPDex[l,m]],Sq[m_,l_]:>Sq[SPDex[l,m]],
        Cpq[m_,l_]:>Cpq[SPDex[l,m]],Spq[m_,l_]:>Spq[SPDex[l,m]]};

TmpA=0;
Do[Do[TmpA=TmpA+(-1)^(lp)*Sum[Sum[Cp[mp,lp]*Cpq[mp+mq,lp+lq]*Cq[mq,lq] 
                                 +Sp[mp,lp]*Spq[mp+mq,lp+lq]*Cq[mq,lq]
                                 -Sp[mp,lp]*Cpq[mp+mq,lp+lq]*Sq[mq,lq]
                                 +Cp[mp,lp]*Spq[mp+mq,lp+lq]*Sq[mq,lq]
                          ,{mq,-lq,lq}],{mp,-lp,lp}];
,{lq,0,EllQ}],{lp,0,EllP}];

(* EXPLOIT REDUNDANCIES BY APPLYING SRULES *)

Tmp=Expand[TmpA]/.SRules;
Tmp=Expand[Tmp]//.{2 x_:>Two x,-2 x_:>-Two x};
Tmp=Expand[Tmp];

(* EXPRESION FOR WHICH NO REDUNDANCIES EXISTED *)

Tmp1=Tmp/.Two->0;

(* EXPRESION THAT CONTAINS FACTORS OF TWO *)

Tmp2=Expand[Tmp-Tmp1];
Tmp2=Tmp2/.Two->1;

(* EXTRACT FACTORS OF ONE AND TWO AND CONVERT TO SP INDEXING *)

Do[Do[ 
       CF1[m,l]=Coefficient[Tmp1,Cp[m,l]]/.TRules;
       SF1[m,l]=Coefficient[Tmp1,Sp[m,l]]/.TRules;
       CF2[m,l]=Coefficient[Tmp2,Cp[m,l]]/.TRules;
       SF2[m,l]=Coefficient[Tmp2,Sp[m,l]]/.TRules;
,{m,0,l}],{l,0,EllP}];

(* SIMPLIFY AND COLLECT FACTORS *)

CList={};
SList={};
CollectQ=Flatten[Table[Table[ {Cq[m,l],Sq[m,l]},{m,0,l}],{l,0,EllQ}]]/.TRules;
Do[Do[ CF2[m,l]=Simplify[Collect[CF2[m,l],CollectQ]];
       SF2[m,l]=Simplify[Collect[SF2[m,l],CollectQ]];
,{m,0,l}],{l,0,EllP}];

Do[Do[
      ldx=ToString[SPDex[l,m]];
      WC[StringJoin["  -------------------------------- Lp=",ToString[l],", Mp = ",ToString[m],"-------------------------------------"]];
      RList={" "->"","Cq("->"Q%C(","Sq("->"Q%S(" };
      AList={Cq,Sq,Cpq,Spq,Cp,Sp};

      If[TrueQ[CF1[m,l]==0],C1=0,
         Write[FileName,FortranAssign[COne,CF1[m,l],AssignToArray->AList,AssignReplace->RList]];C1=1];
      If[TrueQ[SF1[m,l]==0],S1=0,
         Write[FileName,FortranAssign[SOne,SF1[m,l],AssignToArray->AList,AssignReplace->RList]];S1=1];
      If[TrueQ[CF2[m,l]==0],C2=0,
         Write[FileName,FortranAssign[CTwo,CF2[m,l],AssignToArray->AList,AssignReplace->RList]];C2=Two];
      If[TrueQ[SF2[m,l]==0],S2=0,
         Write[FileName,FortranAssign[STwo,SF2[m,l],AssignToArray->AList,AssignReplace->RList]];S2=Two];

      RList={" "->"","c"->StringJoin["SPKetC(",ldx,")"],"s"->StringJoin["SPKetS(",ldx,")"]};

      If[TrueQ[C1+C2==0],,Write[FileName,FortranAssign[c,c+C1*COne+C2*CTwo,AssignToArray->AList,AssignReplace->RList]]];
      If[TrueQ[S1+S2==0],,Write[FileName,FortranAssign[s,s+C1*SOne+S2*STwo,AssignToArray->AList,AssignReplace->RList]]];

,{m,0,l}],{l,0,EllP}];

]
