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
If[lmax!=0 || kc==2 || (kmax>1 && Normy[{1,1,0}]!=1) ,
 WriteString[FileName,StringJoin["      ! Generating (",CType[IntegralClass[{imin,imax}]],",0|",
                                                        CType[IntegralClass[{kmin,kmax}]],",",
                                                        CType[IntegralClass[{lmin,lmax}]],")^(0) \n"]];
 WriteString[FileName,StringJoin["      CALL ",HRRSubName,"(",ToString[BraMax],",HRR) \n"]]; 

   ,
 WriteString[FileName,StringJoin["      ! No need to generate (",CType[IntegralClass[{imin,imax}]],",0|",
                                                        CType[IntegralClass[{kmin,kmax}]],",",
                                                        CType[IntegralClass[{lmin,lmax}]],")^(0) \n"]];
]
];

