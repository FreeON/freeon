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
!    MODULE TO WRITE SORTED ALMLOF STYLE MCMURCHIE-DAVIDSON 3 TERM RECURENCE
!    RELATIONS FOR THE HERMITE GAUSSIAN R[L,M,N] TENSORS W/ L+M+N<=LTot 
!------------------------------------------------------------------------------
*)
ThreeTermRR[LTot_]:=Module[{R,Tmp,Dex,RStr,K,L,M,n},
                           ClearAll[R,Tmp,Dex,RStr];
                           Do[
                              K=0;
                              Do[
                                 LDX1=LMNDex[L,0,0];
                                 LDX2=LMNDex[L-1,0,0];
                                 LDX3=LMNDex[L-2,0,0];
                                 K=K+1;
                                 Dex[K]=LDX1;
                                 Tmp[K]=PQx*r[LDX2]+(L-1)*r[LDX3];
                                 RStr[K]=StringJoin["r(",ToString[LDX1],")"];
                                 Do[
                                    MDX1=LMNDex[M,L,0];
                                    MDX2=LMNDex[M-1,L,0];
                                    MDX3=LMNDex[M-2,L,0];
                                    K=K+1;
                                    Dex[K]=MDX1;
                                    Tmp[K]=PQx*r[MDX2]+(M-1)*r[MDX3]; 
                                    RStr[K]=StringJoin["r(",ToString[MDX1],")"];
                                    Do[
                                       NDX1=LMNDex[n,M,L];
                                       NDX2=LMNDex[n-1,M,L];
                                       NDX3=LMNDex[n-2,M,L];
                                       K=K+1;
                                       Dex[K]=NDX1;
                                       Tmp[K]=PQx*r[NDX2]+(n-1)*r[NDX3];
                                       RStr[K]=StringJoin["r(",ToString[NDX1],")"];
                                      ,{n,LTot-J-L-M,1,-1}];
                                    MDX1=LMNDex[0,M,L];
                                    MDX2=LMNDex[0,M-1,L];
                                    MDX3=LMNDex[0,M-2,L];
                                    K=K+1;
                                    Dex[K]=MDX1;
                                    Tmp[K]=PQy*r[MDX2]+(M-1)*r[MDX3]; 
                                    RStr[K]=StringJoin["r(",ToString[MDX1],")"];
                                    MDX1=LMNDex[M,0,L];
                                    MDX2=LMNDex[M-1,0,L];
                                    MDX3=LMNDex[M-2,0,L];
                                    K=K+1;
                                    Dex[K]=MDX1;
                                    Tmp[K]=PQx*r[MDX2]+(M-1)*r[MDX3]; 
                                    RStr[K]=StringJoin["r(",ToString[MDX1],")"];
                                   ,{M,LTot-J-L,1,-1}];
                                 LDX1=LMNDex[0,L,0];
                                 LDX2=LMNDex[0,L-1,0];
                                 LDX3=LMNDex[0,L-2,0];
                                 K=K+1;
                                 Dex[K]=LDX1;
                                 Tmp[K]=PQy*r[LDX2]+(L-1)*r[LDX3];
                                 RStr[K]=StringJoin["r(",ToString[LDX1],")"];
                                 LDX1=LMNDex[0,0,L];
                                 LDX2=LMNDex[0,0,L-1];
                                 LDX3=LMNDex[0,0,L-2];
                                 K=K+1;
                                 Dex[K]=LDX1;
                                 Tmp[K]=PQz*r[LDX2]+(L-1)*r[LDX3];
                                 RStr[K]=StringJoin["r(",ToString[LDX1],")"];
                                ,{L,LTot-J,1,-1}];
                                KList=Table[{Dex[I],I},{I,K}];
                             (* Sorting on this index introduces a bug...
                                KList=Sort[KList]; *) 
                                Do[
                                   Q=KList[[I,2]];
                                    Write[FileName,FortranAssign[rr,Tmp[Q],AssignReplace->{"rr"->RStr[Q]}]];
                                  ,{I,K}];
                                 Tmp2=AuxR[J];
                                 Write[FileName,FortranAssign[rr,Tmp2,AssignReplace->{"rr"->"r(1)"}]];
                             ,{J,LTot,0,-1}];
];
