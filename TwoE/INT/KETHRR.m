(*======================================--  KET HRR --============================================*)

PunchKetHRRFront[Subroutine_,Subname_,kc_,lc_]:=Block[{WS,BKType,LenBra,LenKet},

WriteString[Subroutine,StringJoin["   SUBROUTINE ",Subname,"(LB,HRR) \n"]];
kmin = Classes[[kc, 1]]; kmax = Classes[[kc, 2]];
lmin = Classes[[lc, 1]]; lmax = Classes[[lc, 2]];
LenK=LEnd[Classes[[kc,2]]]-LBegin[Classes[[kc,1]]]+1;
LenL=LEnd[Classes[[lc,2]]];
KEnd=LEnd[kmax+lmax];

KetMax=KEnd;

If[kc==2&&lc==2,KetMax=KetMax+LenK+LenL];
If[kc==2&&lc!=2,KetMax=KetMax+1];
If[lc==2&&kc!=2,KetMax=KetMax+LenK];

WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];
WS["USE DerivedTypes"];
WS["USE VScratchB"];
WS["USE GlobalScalars"];
WS["IMPLICIT REAL(DOUBLE) (W)"]; 
WS["INTEGER :: LB"];							      
WS[StringJoin["REAL(DOUBLE) :: HRR(1:LB,",ToString[KetMax],",",ToString[LenL],")"]];

];

(* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  *)

HRRKet[adex_,c_List,d_List]:=Module[{AB,CD,pa,pb,pc,pd,MaxEll,a1,b1,c1,d1,one,two,ddex,cdex},
                              pc=Position[c,Max[c]][[1, 1]];
                              pd=Position[d,Max[d]][[1, 1]];
			      (* Exit condition 1 *)
                              If[ c[[pc]] < 0 || d[[pd]] < 0 ,Return[0];]; 
			      (* Exit condition 2 *)
                              If[ d[[pd]]==0,
                                 cdex=LMNDex[c[[1]],c[[2]],c[[3]]];
                                 ddex=LMNDex[d[[1]],d[[2]],d[[3]]];
                                 BarExp=ToExpression[StringJoin["HRR[",ToString[adex],",",ToString[cdex],",",ToString[ddex],"]"]];
                                 Return[BarExp]
                                ];
			      (* Recursion *)
                              If[pd==1, one = {1, 0, 0}; CD = CDx; ];
                              If[pd==2, one = {0, 1, 0}; CD = CDy; ];
                              If[pd==3, one = {0, 0, 1}; CD = CDz; ];
                              c1=c+one;
                              d1=d-one;                
                              Return[HRRKet[adex,c1,d1]+CD HRRKet[adex,c,d1]]
                             ];

(* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  *)

PunchHRRKet:=Module[{oList,IList,Kount,a,b,c,d},
Do[Do[
 If[IntegralClass[Classes[[kc]]]>=IntegralClass[Classes[[lc]]] ,
 kmin = Classes[[kc, 1]]; kmax = Classes[[kc, 2]];
 lmin = Classes[[lc, 1]]; lmax = Classes[[lc, 2]];
 LenK=LEnd[Classes[[kc,2]]]-LBegin[Classes[[kc,1]]]+1;
 LenK1=LBegin[Classes[[kc,1]]]-1;
 LenL=LEnd[Classes[[lc,2]]]-LBegin[Classes[[lc,1]]]+1;
 KEnd=LEnd[kmax+lmax];
 KetMax=KEnd;
 If[kc==2&&lc==2,KetMax=KetMax+LenK+LenL];
 If[kc==2&&lc!=2,KetMax=KetMax+1];
 If[lc==2&&kc!=2,KetMax=KetMax+LenK];
 IList={};
 oList={"+"->"+  & \n        "};
 Kount = 0;
 Do[Do[Do[Do[
             Kount = Kount + 1;
             IList=Append[IList,ToExpression[StringJoin["KK",ToString[k],"xx",ToString[l],"BB"]]];
             oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->"!"];
                c = {lx[k], my[k], nz[k]};
                d = {lx[l], my[l], nz[l]};
                hrr=Normy[c]*Normy[d]*HRRKet[1,c,d];
                (* Ket addressing scheme for auxiliary SP integrals *)
                (*  1   LenK  1    LenL *)
	        (*  SS  SFk   SSl  SFl  *)
                If[ kc != 2 && lc == 2 && ll == 0 ,hrr=hrr/.{HRR[a_,b_,c_]:>HRR[a,KEnd+b-LenK1,c]}];
                If[ kc == 2 && lc != 2 && kl == 0 ,hrr=hrr/.{HRR[a_,b_,c_]:>HRR[a,KEnd+b,c]}];
                If[ kc == 2 && lc == 2 && kl >=  0 && ll == 0 ,hrr=hrr/.{HRR[a_,b_,c_]:>HRR[a,KEnd+b,c]}];
                If[ kc == 2 && lc == 2 && kl == 0 && ll > 0 ,hrr=hrr/.{HRR[a_,b_,c_]:>HRR[a,KEnd+LenK+b,c]}];
                If[hrr=!=HRR[1,k,l],
		   hrr=hrr/.{HRR[i_,k_,l_]:>HRR[abc,k,l]};
                   Kount = Kount + 1;
                   IList=Append[IList,Horner[hrr]];
                   HRRAddress=StringJoin["1:LB,",ToString[k],",",ToString[l]];
                   oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["HRR(",HRRAddress,")"]];
                  ]; 
 ,{k,LBegin[kl],LEnd[kl]}]
 ,{l,LBegin[ll],LEnd[ll]}]
 ,{kl,kmax,kmin,-1}]
 ,{ll,lmax,lmin,-1}];
  oList=Append[oList,{" "->"","xx"->",","BB"->")","KK"->"|","In"->"INTGRL","abc"->"1:LB"}];
  oList=Flatten[oList];
  If[Kount>0,
     KClass=ToString[IntegralClass[Classes[[kc]]]];
     LClass=ToString[IntegralClass[Classes[[lc]]]];
     HRRSubName=StringJoin["KetHRR",KClass,LClass];
     MakeList=Append[MakeList,StringJoin[HRRSubName,".o"]];
     RelsList=Append[RelsList,StringJoin[HRRSubName,".x"]];
     HRRSubroutine=StringJoin[HRRSubName,".F90"];
     OpenWrite[HRRSubroutine];
     PunchKetHRRFront[HRRSubroutine,HRRSubName,kc,lc];
     Write[HRRSubroutine,FortranAssign[o,IList,AssignReplace->oList]];
     WriteString[HRRSubroutine,StringJoin["END SUBROUTINE ",HRRSubName]];			       
     Close[HRRSubroutine];
    ];
 ]
 ,{lc,LC}],{kc,LC}]
 ];

(* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  *)

PunchHRRKetCalls[FileName_,ic_,jc_,kc_,lc_]:=Module[{oList,IList,Kount,a,b,c,d},
 imin = Classes[[ic, 1]]; imax = Classes[[ic, 2]];
 jmin = Classes[[jc, 1]]; jmax = Classes[[jc, 2]];
 kmin = Classes[[kc, 1]]; kmax = Classes[[kc, 2]];
 lmin = Classes[[lc, 1]]; lmax = Classes[[lc, 2]];
 LenI=LEnd[Classes[[ic,2]]]-LBegin[Classes[[ic,1]]]+1;
 LenJ=LEnd[Classes[[jc,2]]]-LBegin[Classes[[jc,1]]]+1;
 LenK=LEnd[Classes[[kc,2]]]-LBegin[Classes[[kc,1]]]+1;
 LenK1=LBegin[Classes[[kc,1]]]-1;
 LenL=LEnd[Classes[[lc,2]]]-LBegin[Classes[[lc,1]]]+1;
 BEnd=LEnd[imax+jmax];
 KEnd=LEnd[kmax+lmax];
 BraMax=BEnd;
 If[ic==2&&jc==2,BraMax=BraMax+LenI+LenJ];
 If[ic==2&&jc!=2,BraMax=BraMax+1];
 If[jc==2&&ic!=2,BraMax=BraMax+LenI];
 KetMax=KEnd;
 If[kc==2&&lc==2,KetMax=KetMax+LenK+LenL];
 If[kc==2&&lc!=2,KetMax=KetMax+1];
 If[lc==2&&kc!=2,KetMax=KetMax+LenK];
 
 KClass=ToString[IntegralClass[Classes[[kc]]]];
 LClass=ToString[IntegralClass[Classes[[lc]]]];
 HRRSubName=StringJoin["KetHRR",KClass,LClass];
 WriteString[FileName,StringJoin["      ! Generating (",CType[IntegralClass[{imin,imax}]],",0|",
                                                        CType[IntegralClass[{kmin,kmax}]],",",
                                                        CType[IntegralClass[{lmin,lmax}]],")^(0) \n"]];
 WriteString[FileName,StringJoin["      CALL ",HRRSubName,"(",ToString[BraMax],",HRR) \n"]]; 
];

