(*======================================--  BRA HRR --============================================*)

PunchBraHRRFront[Subroutine_,Subname_,ic_,jc_]:=Block[{WS,BKType,LenBra,LenKet},
 WriteString[Subroutine,StringJoin["   SUBROUTINE ",Subname,"(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) \n"]];

 imin = Classes[[ic, 1]]; imax = Classes[[ic, 2]];
 jmin = Classes[[jc, 1]]; jmax = Classes[[jc, 2]];

 LenI=LEnd[Classes[[ic,2]]]-LBegin[Classes[[ic,1]]];
 LenJ=LEnd[Classes[[jc,2]]]-LBegin[Classes[[jc,1]]];

 BEnd=LEnd[imax+jmax];

 BraMax=BEnd;

 WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];
 WS["USE DerivedTypes"];
 WS["USE VScratchB"];
 WS["USE GlobalScalars"];
 WS["IMPLICIT REAL(DOUBLE) (W)"]; 
 WS["INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet"];
 WS[StringJoin["REAL(DOUBLE)  :: HRR(*)"]];
 WS[StringJoin["REAL(DOUBLE)  :: INTGRL(*)"]];
];

(* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  *)

HRRBra[a_List,b_List]:=Module[{AB,CD,pa,pb,pc,pd,MaxEll,a1,b1,c1,d1,one,two,adex,cdex},
                              pa=Position[a,Max[b]][[1, 1]];
                              pb=Position[b,Max[b]][[1, 1]];
			      (* Exit condition 1 *)
                              If[ a[[pa]] < 0 || b[[pb]] < 0 ,Return[0];]; 
			      (* Exit condition 2 *)
                              If[ b[[pb]]==0 ,
                                 adex=LMNDex[a[[1]],a[[2]],a[[3]]];
                                 BarExp=ToExpression[StringJoin["HRR[",ToString[adex],"]"]];
                                 Return[BarExp]
                                ];
			      (* Recursion *)
                              If[pb==1, one = {1, 0, 0}; AB = ABx; ];
                              If[pb==2, one = {0, 1, 0}; AB = ABy; ];
                              If[pb==3, one = {0, 0, 1}; AB = ABz; ];
                              a1 = a + one;
                              b1 = b - one;                
                              Return[ HRRBra[a1, b1] + AB  HRRBra[a, b1]]
                            ];

(* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  *)

PunchHRRBraCalls[FileName_,ic_,jc_,kc_,lc_]:=Module[{oList,IList,Kount,a,b,c,d},

 imin = Classes[[ic, 1]]; imax = Classes[[ic, 2]];
 jmin = Classes[[jc, 1]]; jmax = Classes[[jc, 2]];
 kmin = Classes[[kc, 1]]; kmax = Classes[[kc, 2]];
 lmin = Classes[[lc, 1]]; lmax = Classes[[lc, 2]];
 LenI=LEnd[Classes[[ic,2]]]-LBegin[Classes[[ic,1]]];
 LenJ=LEnd[Classes[[jc,2]]]-LBegin[Classes[[jc,1]]];
 BEnd=LEnd[imax+jmax];

 IClass=ToString[IntegralClass[Classes[[ic]]]];
 JClass=ToString[IntegralClass[Classes[[jc]]]];
 HRRSubName=StringJoin["BraHRR",IClass,JClass];

 WriteString[FileName,StringJoin["      ! Generating (",CType[IntegralClass[{imin,imax}]],",",
                                                        CType[IntegralClass[{jmin,jmax}]],"|",
                                                        CType[IntegralClass[{kmin,kmax}]],",",
                                                        CType[IntegralClass[{lmin,lmax}]],")^(0) \n"]];

 WriteString[FileName,StringJoin["      DO L=",ToString[LBegin[lmin]],",",ToString[LEnd[lmax]],"\n"]];
 WriteString[FileName,StringJoin["         DO K=",ToString[LBegin[kmin]],",",ToString[LEnd[kmax]],"\n"]];
 WriteString[FileName,StringJoin["            CDOffSet=(OC+K-",ToString[LBegin[kmin]],")*LDC+(OD+L-",ToString[LBegin[lmin]],")*LDD \n"]];
 WriteString[FileName,StringJoin["            CALL ",HRRSubName,"(OA,OB,LDA,LDB,CDOffSet,HRR(1,",ToString[K],",",ToString[L],"),INTGRL) \n"]];
 WriteString[FileName,"          ENDDO \n"];
 WriteString[FileName,"      ENDDO \n "];

];

(* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  *)

PunchHRRBra:=Module[{oList,IList,Kount,a,b,c,d},

Do[Do[

 imin = Classes[[ic, 1]]; imax = Classes[[ic, 2]];
 jmin = Classes[[jc, 1]]; jmax = Classes[[jc, 2]];
 LenI=LEnd[Classes[[ic,2]]]-LBegin[Classes[[ic,1]]];
 LenI1=LBegin[Classes[[ic,2]]]-1;
 LenJ=LEnd[Classes[[jc,2]]]-LBegin[Classes[[jc,1]]];
 BEnd=LEnd[imax+jmax];

 IClass=ToString[IntegralClass[Classes[[ic]]]];
 JClass=ToString[IntegralClass[Classes[[jc]]]];
 HRRSubName=StringJoin["BraHRR",IClass,JClass];

 IList={};
 oList={};
 Kount = 0;
 Do[Do[Do[Do[
             a = {lx[i], my[i], nz[i]};
             b = {lx[j], my[j], nz[j]};

	     hrr=Normy[a]*Normy[b]*HRRBra[a,b];

             (* Bra addressing scheme for auxiliary SP integrals *)
             (*  1   LenI  1    LenJ *)
	     (*  SS  SFi   SSj  SFj  *)

             If[ ic == 2 && jc != 2 && il == 0 ,hrr=hrr/.{HRR[a_]:>HRR[BEnd+a]}];
             If[ ic != 2 && jc == 2 && jl == 0 ,hrr=hrr/.{HRR[a_]:>HRR[BEnd+a-LenI1]}];
             If[ ic == 2 && jc == 2 && il >= 0 && jl == 0 ,hrr=hrr/.{HRR[a_]:>HRR[BEnd+a]}];
             If[ ic == 2 && jc == 2 && il == 0 && jl >  0 ,hrr=hrr/.{HRR[a_]:>HRR[BEnd+LenI+1+a]}];
             If[hrr=!=HRR[i,k,l],

                OffSetString=StringJoin["(OA+",ToString[i-LBegin[imin]],")*LDA+(OB+",ToString[j-LBegin[jmin]],")*LDB+CDOffSet"];
                Kount = Kount + 1;

                IList=Append[IList,ToExpression[StringJoin["BB",ToString[i],"xx",ToString[j],"KK"]]];
                oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["OffSet=",OffSetString," !"]];
                Kount = Kount + 1;
                IList=Append[IList,Horner[hrr]];
                HRRAddress=StringJoin[ToString[i],",",ToString[j]];
                oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["INTGRL(OffSet)"]];
		];

 ,{i,LBegin[il],LEnd[il]}]
 ,{j,LBegin[jl],LEnd[jl]}]
 ,{il,imin,imax}]
 ,{jl,jmin,jmax}];

  oList=Append[oList,{" "->"","xx"->",","BB"->"(","KK"->"|","In"->"INTGRL","+"->"+  & \n        "}];
  oList=Flatten[oList];

 If[Kount>0,
    MakeList=Append[MakeList,StringJoin[HRRSubName,".o"]];
    RelsList=Append[RelsList,StringJoin[HRRSubName,".x"]];
    HRRSubroutine=StringJoin[HRRSubName,".F90"];
    OpenWrite[HRRSubroutine];
    PunchBraHRRFront[HRRSubroutine,HRRSubName,ic,jc];
    Write[HRRSubroutine,FortranAssign[o,IList,AssignReplace->oList]];
    WriteString[HRRSubroutine,StringJoin["END SUBROUTINE ",HRRSubName]];			       
    Close[HRRSubroutine];
   ];

,{ic,LC}],{jc,LC}];

];
