ThreeTermRR[LBra_,LKet_]:=Module[{R,LTot},
                           ClearAll[R];
                           LTot=LBra+LKet;
                           Print["LTot = ",LTot];
                           Do[
                              R[0,0,0,J]=AuxR[J];
                             ,{J,0,LTot}];
                           Do[ 
                              Do[ 
                                  R[0,0,N,J]=R[0,0,N-1,J+1]*PQz+(N-1)R[0,0,N-2,J+1];
                                ,{J,0,LTot-N}];
                             ,{N,1,LTot}];
                           Do[ 
                              Do[
                                 Do[ 
                                     R[0,M,N,J]=R[0,M-1,N,J+1]*PQy+(M-1)R[0,M-2,N,J+1];
                                   ,{J,0,LTot-N-M}];
                                 ,{M,1,LTot}];
                              ,{N,0,LTot}];

                           Do[
                              Do[ 
                                 Do[
                                    Do[ 
                                        R[L,M,N,J]=R[L-1,M,N,J+1]*PQx+(L-1)R[L-2,M,N,J+1];
                                      ,{J,0,LTot-N-M-L}];
                                   ,{L,1,LTot-M-N}];
                                 ,{M,0,LTot-N}];
                              ,{N,0,LTot}];                              

                              Do[Do[Do[   
                                       BraDex=LMNDex[LB,MB,NB]; 
                                       PKet[BraDex]=0;
                              Do[Do[Do[
                                       BraDex=LMNDex[LB,MB,NB];
                                       KetDex=LMNDex[LK,MK,NK];
                                       Print["NK = ",NK,"MK = ",MK,"LK = ",LK," KetDex = ",KetDex];
                                       PKet[BraDex]=PKet[BraDex]+R[LB+LK,MB+MK,NB+NK,0]*Co[KetDex];
                              ,{NK,0,LKet-LK-MK}],{MK,0,LKet-LK}],{LK,0,LKet}];
                              ,{NB,0,LBra-LB-MB}],{MB,0,LBra-LB}],{LB,0,LBra}];

                             KList=Simplify[ExpandAll[Table[PKet[k],{k,HGLen[LBra]}]]];

Print[KList];
                             FortranAssign[Ket,KList]]