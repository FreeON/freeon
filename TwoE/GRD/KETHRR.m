
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

 If[Classes[[ic]]=={0,1},ISP=TRUE,ISP=FALSE];
 If[Classes[[jc]]=={0,1},JSP=TRUE,JSP=FALSE];
 If[Classes[[kc]]=={0,1},KSP=TRUE,KSP=FALSE];
 If[Classes[[lc]]=={0,1},LSP=TRUE,LSP=FALSE];

 If[ic==2&&jc==2,BraMax=BraMax+LenI+LenJ];
 If[ic==2&&jc!=2,BraMax=BraMax+1];
 If[jc==2&&ic!=2,BraMax=BraMax+LenI];

 KetMax=KEnd;

 If[kc==2&&lc==2,KetMax=KetMax+LenK+LenL];
 If[kc==2&&lc!=2,KetMax=KetMax+1];
 If[lc==2&&kc!=2,KetMax=KetMax+LenK];
 
 KC1=Classes[[kc]];

 If[kc==2,
    KC1[[2]]=KC1[[2]]+1;
    kmin1=kmin;
   ,
    KC1=KC1+1;
    kmin1=kmax+1;
    ];

 KClass=ToString[IntegralClass[Classes[[kc]]]];
 KClass1=ToString[IntegralClass[KC1]];
 LClass=ToString[IntegralClass[Classes[[lc]]]];


 ijmax=imax+jmax;
 ijmax1=ijmax+1;
 If[ic==2||jc==2,
    ijmin=0;
    ijmin1=0;
    ,
    ijmin=ijmax;
    ijmin1=ijmax1;
   ];

 HRRSubName=StringJoin["KetHRR",KClass,LClass];



						    WriteString[FileName,StringJoin["      ! Generating (",CType[IntegralClass[{ijmin,ijmax}]],",0|",
                                                        CType[IntegralClass[{kmin,kmax}]],",",
                                                        CType[IntegralClass[{lmin,lmax}]],")\n"]]; 

 WriteString[FileName,StringJoin["      CALL ",HRRSubName,"(",ToString[BraMax],",HRR) \n"]]; 

			  			    WriteString[FileName,StringJoin["      ! Generating (",CType[IntegralClass[{ijmin1,ijmax1}]],",0|",
                                                        CType[IntegralClass[{kmin,kmax}]],",",
                                                        CType[IntegralClass[{lmin,lmax}]],")^a\n"]]; 
 
 WriteString[FileName,StringJoin["      CALL ",HRRSubName,"(",ToString[BraMax],",HRRA) \n"]]; 

						    WriteString[FileName,StringJoin["      ! Generating (",CType[IntegralClass[{ijmin1,ijmax1}]],",0|",
                                                        CType[IntegralClass[{kmin,kmax}]],",",
                                                        CType[IntegralClass[{lmin,lmax}]],")^b\n"]]; 

 WriteString[FileName,StringJoin["      CALL ",HRRSubName,"(",ToString[BraMax],",HRRB) \n"]]; 


						    WriteString[FileName,StringJoin["      ! Generating (",CType[IntegralClass[{ijmin,ijmax}]],",0|",
                                                        CType[IntegralClass[{kmin1,kmax+1}]],",",
                                                        CType[IntegralClass[{lmin,lmax}]],")^c\n"]]; 
 HRRSubName=StringJoin["KetHRR",KClass1,LClass];

 WriteString[FileName,StringJoin["      CALL ",HRRSubName,"(",ToString[BraMax],",HRRC) \n"]]; 

];

