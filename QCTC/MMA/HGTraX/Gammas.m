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
!    WRITE EXPLICIT CODE FOR ASSOCIATED INCOMPLETE GAMMA FUNCTIONS USING 
!    DOWNWARD RECURRENCE FROM LTOT+2 INTERPOLATED FUNCTION AND INTRINSIC EXP
!------------------------------------------------------------------------------
*)
GammasPrelim:=Module[{},
                      WS["J=AINT(T*Gamma_Grid)"];                                              
];

FStr[L_]:=Module[{LSt},
                 LSt=ToString[L];
                 StringJoin["(F",LSt,"_0(J)+T*(F",LSt,"_1(J)+T*(F",LSt,"_2(J)+T*(F",LSt,"_3(J)+T*F",LSt,"_4(J)))))"]];

Gammas[LTot_]:=Module[{LSt,LSt2,GSt,GSt1,GFlt},
                      If[LTot==0,WS[StringJoin["  AuxR(0)=Upq*",FStr[0]]]];
                      If[LTot==1,
                         WS[StringJoin["  AuxR(0)=Upq*",FStr[0]]];
                         WS[StringJoin["  AuxR(1)=-Two*Omega*Upq*",FStr[1]]];
                        ];
                      If[LTot>1,
                         LSt=ToString[LTot];
                         WS["  o1=Upq"];
                         WS["  o2=-Two*Omega"];
                         WS["  ET=EXP(-T)"];
                         WS["  TwoT=Two*T"];
                         WS[StringJoin["  G(",LSt,")=",FStr[LTot]]];
                         Do[GSt=ToString[L-1];
                            GSt1=ToString[L];
                            GFlt=FF[1/(2 *(L-1)+1)];
                            WS[StringJoin["  G(",GSt,")=",GFlt,"*(TwoT*G(",GSt1,")+ET)"]];
                           ,{L,LTot,2,-1}];
                         WS["  G(0)=TwoT*G(1)+ET"];
                         Do[GSt=ToString[L];
                            WS[StringJoin["  AuxR(",GSt,")=o1*G(",GSt,")"]];
                            If[L<LTot,WS[StringJoin["  o1=o2*o1"]]];
                           ,{L,0,LTot}];
                        ]
]







(*

FStr[L_]:=Module[{LSt},
                 LSt=ToString[L];
                 StringJoin["(F",LSt,"_0(J)+T*(F",LSt,"_1(J)+T*(F",LSt,"_2(J)+T*(F",LSt,"_3(J)+T*F",LSt,"_4(J)))))"]];

Gammas[LTot_]:=Module[{LSt,LSt2,GSt,GSt1,GFlt},
                      If[LTot==0,WS[StringJoin["  AuxR(0)=Upq*",FStr[0]]]];
                      If[LTot==1,
                         WS[StringJoin["  AuxR(0)=Upq*",FStr[0]]];
                         WS[StringJoin["  AuxR(1)=-Two*Omega*Upq*",FStr[1]]];
                        ];
                      If[LTot>1,
                         LSt=ToString[LTot];
                         WS["  o1=Upq"];
                         WS["  o2=-Two*Omega"];
                         WS["  ET=EXP(-T)"];
                         WS["  TwoT=Two*T"];
                         WS[StringJoin["  AuxR(",LSt,")=",FStr[LTot]]];
                         Do[GSt=ToString[L-1];
                            GSt1=ToString[L];
                            GFlt=FF[1/(2 *(L-1)+1)];
                            WS[StringJoin["  AuxR(",GSt,")=",GFlt,"*(TwoT*AuxR(",GSt1,")+ET)"]];
                           ,{L,LTot,2,-1}];
                         WS["  AuxR(0)=TwoT*AuxR(1)+ET"];
                         Do[GSt=ToString[L];
                            WS[StringJoin["  AuxR(",GSt,")=o1*AuxR(",GSt,")"]];
                            If[L<LTot,WS[StringJoin["  o1=o2*o1"]]];
                           ,{L,0,LTot}];
                        ]
]


*)
