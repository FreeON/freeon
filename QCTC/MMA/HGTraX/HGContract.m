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
!    WRITE EXPLICIT CODE FOR CONTRACTION OF MD R TENSORS WITH 
!    HERMITE GAUSSIAN DENSITY COEFFICIENTS
!------------------------------------------------------------------------------
*)
HGContract[LBra_,LKet_]:=Module[{KList,PKet,BraDex,KetDex},
                                 Do[Do[Do[   
                                     BraDex=LMNDex[LB,MB,NB]; 
                                     PKet[BraDex]=0;
                                     Do[Do[Do[
                                              BraDex=LMNDex[LB,MB,NB];
                                              KetDex=LMNDex[LK,MK,NK];
                                              MDRDex=LMNDex[LK+LB,MK+MB,NK+NB];
                                              PKet[BraDex]=PKet[BraDex]+r[MDRDex]*Co[KetDex];
                                       ,{NK,0,LKet-LK-MK}],{MK,0,LKet-LK}],{LK,0,LKet}];
                              ,{NB,0,LBra-LB-MB}],{MB,0,LBra-LB}],{LB,0,LBra}];
                              KList=ExpandAll[Table[PKet[k],{k,HGLen[LBra]}]];
                              Do[
                                 KStr=ToString[K];
                                 KetStr=StringJoin["HGKet(",KStr,")=HGKet(",KStr,") +"];
                                 CoStr="Q%Co";
                                 Write[FileName,FortranAssign[pk,KList[[K]],
                                                              AssignReplace->{"pk ="->KetStr,"Co"->CoStr}]
                                      ];
                                ,{K,HGLen[LBra]}];
];