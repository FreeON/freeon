
(* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  *)

PunchVRRCalls[FileName_,BraEll_,KetEll_]:=Module[{spaces,CLenBra,CLenKet,CLenEhm},

WS[String_]:=WriteString[FileName,"            ",String,"\n"];

CLenBra=ToString[LEnd[BraEll+1]];
CLenKet=ToString[LEnd[KetEll+1]];
CLenEhm=ToString[BraEll+KetEll+1];
Do[Do[
If[iell+kell>0 ,

   VRRSubName=StringJoin["VRR",CType[IntegralClass[{iell,iell}]],"0",CType[IntegralClass[{kell,kell}]],"0"];
   Do[ 

      ms=ToString[m];
      m1=ToString[m-1];

      WS[StringJoin["! Generating (",CType[IntegralClass[{iell,iell}]],"0|", \
                                     CType[IntegralClass[{kell,kell}]],"0)^(",m1,")"]]; 

      (*
      If[iell==0&&kell==1,
         WS[StringJoin["VRR(1,2,",m1,")=QCx*VRR(1,1,",m1,")+WQx*VRR(1,1,",ms,")"]];
         WS[StringJoin["VRR(1,3,",m1,")=QCy*VRR(1,1,",m1,")+WQy*VRR(1,1,",ms,")"]];
         WS[StringJoin["VRR(1,4,",m1,")=QCz*VRR(1,1,",m1,")+WQz*VRR(1,1,",ms,")"]];
        ];

      If[iell==1&&kell==0,
	 WS[StringJoin["VRR(2,1,",m1,")=PAx*VRR(1,1,",m1,")+WPx*VRR(1,1,",ms,") " ]];
	 WS[StringJoin["VRR(3,1,",m1,")=PAy*VRR(1,1,",m1,")+WPy*VRR(1,1,",ms,") " ]];
	 WS[StringJoin["VRR(4,1,",m1,")=PAz*VRR(1,1,",m1,")+WPz*VRR(1,1,",ms,") " ]];
        ];

      If[iell==1&&kell==1,
         WS[StringJoin["VRR(2,2,",m1,")=QCx*VRR(2,1,",m1,")+HfxZpE*VRR(1,1,",ms,")+WQx*VRR(2,1,",ms,") "]];
         WS[StringJoin["VRR(2,3,",m1,")=QCy*VRR(2,1,",m1,")+WQy*VRR(2,1,",ms,") "]];
         WS[StringJoin["VRR(2,4,",m1,")=QCz*VRR(2,1,",m1,")+WQz*VRR(2,1,",ms,") "]];
         WS[StringJoin["VRR(3,2,",m1,")=QCx*VRR(3,1,",m1,")+WQx*VRR(3,1,",ms,") "]];
         WS[StringJoin["VRR(3,3,",m1,")=QCy*VRR(3,1,",m1,")+HfxZpE*VRR(1,1,",ms,")+WQy*VRR(3,1,",ms,") "]];
         WS[StringJoin["VRR(3,4,",m1,")=QCz*VRR(3,1,",m1,")+WQz*VRR(3,1,",ms,") "]];
         WS[StringJoin["VRR(4,2,",m1,")=QCx*VRR(4,1,",m1,")+WQx*VRR(4,1,",ms,") "]];
         WS[StringJoin["VRR(4,3,",m1,")=QCy*VRR(4,1,",m1,")+WQy*VRR(4,1,",ms,") "]];
         WS[StringJoin["VRR(4,4,",m1,")=QCz*VRR(4,1,",m1,")+HfxZpE*VRR(1,1,",ms,")+WQz*VRR(4,1,",ms,") "]];
        ];

      If[iell==2&&kell==0,
         WS[StringJoin["VRR(5,1,",m1,")=PAx*VRR(2,1,",m1,")+r1x2Z*(VRR(1,1,",m1,")-ExZpE*VRR(1,1,",ms,"))+WPx*VRR(2,1,",ms,")"]];
         WS[StringJoin["VRR(6,1,",m1,")=PAx*VRR(3,1,",m1,")+WPx*VRR(3,1,",ms,")"]];
         WS[StringJoin["VRR(7,1,",m1,")=PAy*VRR(3,1,",m1,")+r1x2Z*(VRR(1,1,",m1,")-ExZpE*VRR(1,1,",ms,"))+WPy*VRR(3,1,",ms,")"]];
         WS[StringJoin["VRR(8,1,",m1,")=PAx*VRR(4,1,",m1,")+WPx*VRR(4,1,",ms,")"]];
         WS[StringJoin["VRR(9,1,",m1,")=PAy*VRR(4,1,",m1,")+WPy*VRR(4,1,",ms,")"]];
         WS[StringJoin["VRR(10,1,",m1,")=PAz*VRR(4,1,",m1,")+r1x2Z*(VRR(1,1,",m1,")-ExZpE*VRR(1,1,",ms,"))+WPz*VRR(4,1,",ms,")"]];
        ];

      If[(iell>2&&kell==0) || (iell==0 && kell>1) || (iell>1 && kell>=1) || (iell >=1 && kell>1 ),
         WS[StringJoin["CALL ",VRRSubName,"(",CLenBra,",",CLenKet,",VRR(1,1,",ToString[m-1],"),VRR(1,1,",ToString[m],"))"]]; 
        ]
       *)
      WS[StringJoin["CALL ",VRRSubName,"(",CLenBra,",",CLenKet,",VRR(1,1,",ToString[m-1],"),VRR(1,1,",ToString[m],"))"]]; 

     ,{m,BraEll+KetEll-(iell+kell)+2,1,-1}]; 

];
,{iell,0,BraEll+1}];
,{kell,0,KetEll+1}];
];

(*  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  *)

(*======================--  VRR 2 HRR CONTRACTION STEP (A FUCKING BITCH IF SP) --===============================*)

PunchContractCalls[Subroutine_,ic_,jc_,kc_,lc_]:=Block[{WS,BKString,bra,ket},
      WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];
      WS["      ! Contracting ... "] ;

      imax = Classes[[ic, 2]];
      jmax = Classes[[jc, 2]];
      kmax = Classes[[kc, 2]];
      lmax = Classes[[lc, 2]];

      BLen=LEnd[imax+jmax+1];
      KLen=LEnd[kmax+lmax+1];

      LMax=imax+jmax+kmax+lmax;

      IClass=ToString[IntegralClass[Classes[[ic]]]];
      JClass=ToString[IntegralClass[Classes[[jc]]]];
      KClass=ToString[IntegralClass[Classes[[kc]]]];
      LClass=ToString[IntegralClass[Classes[[lc]]]];

      ContractName=StringJoin["CNTRCT",IClass,JClass,KClass,LClass];

      WS[StringJoin["      CALL DBLAXPY( ",ToString[BLen*KLen],",HRR(1,1,1),       VRR(1,1,0)) "]];   
      WS[StringJoin["      CALL DBLAXPZY(",ToString[BLen*KLen],",HRRA(1,1,1),Alpha,VRR(1,1,0)) "]];   
      WS[StringJoin["      CALL DBLAXPZY(",ToString[BLen*KLen],",HRRB(1,1,1),Beta, VRR(1,1,0)) "]];   
      WS[StringJoin["      CALL DBLAXPZY(",ToString[BLen*KLen],",HRRC(1,1,1),Gamma,VRR(1,1,0)) "]];   
						       (*
      If[ic==2 || jc==2 || kc==2 || lc==2 || LMax < 6  ,
          WS[StringJoin["      CALL ",ContractName,"(VRR,HRR)"]];
	 ,
          WS[StringJoin["      CALL DBLAXPY(",ToString[BEnd*KEnd],",HRR(1,1,1),VRR(1,1,0)) "]];   
        ];
							*)
 ];
PunchVRRContract[Subroutine_,ic_,jc_,kc_,lc_,LBra_,LKet_]:=Block[{WS,BKString,bra,ket},

Return[];

      imax = Classes[[ic, 2]];
      jmax = Classes[[jc, 2]];
      kmax = Classes[[kc, 2]];
      lmax = Classes[[lc, 2]];
      LMax=imax+jmax+kmax+lmax;
      If[ic==2||jc==2||kc==2||lc==2|| LMax < 6  ,

      WS[String_]:=WriteString[Subroutine,"    ",String,"\n"];
      IClass=ToString[IntegralClass[Classes[[ic]]]];
      JClass=ToString[IntegralClass[Classes[[jc]]]];
      KClass=ToString[IntegralClass[Classes[[kc]]]];
      LClass=ToString[IntegralClass[Classes[[lc]]]];
      LenI=LEnd[Classes[[ic,2]]]-LBegin[Classes[[ic,1]]]+1;
      LenJ=LEnd[Classes[[jc,2]]]-LBegin[Classes[[jc,1]]]+1;
      LenK=LEnd[Classes[[kc,2]]]-LBegin[Classes[[kc,1]]]+1;
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
      ContractName=StringJoin["CNTRCT",IClass,JClass,KClass,LClass];
      WS[StringJoin["SUBROUTINE ",ContractName,"(VRR,HRR)"]];
      WS["  USE DerivedTypes"];
      WS[StringJoin["  REAL(DOUBLE)  :: VRR(",ToString[LEnd[LBra]],",",ToString[LEnd[LKet]],",0:",ToString[LBra+LKet],")"]];
      WS[StringJoin["  REAL(DOUBLE)  :: HRR(",ToString[BraMax],",",ToString[KetMax],",",ToString[LEnd[lmax]],")"]];

      BEnd=LEnd[LBra];
      KEnd=LEnd[LKet];
      LenI=LEnd[Classes[[ic,2]]]-LBegin[Classes[[ic,1]]];
      LenI1=LBegin[Classes[[ic,2]]]-1;
      LenJ=LEnd[Classes[[jc,2]]]-LBegin[Classes[[jc,1]]];
      LenK=LEnd[Classes[[kc,2]]]-LBegin[Classes[[kc,1]]];
      LenK1=LBegin[Classes[[kc,2]]]-1;
      LenL=LEnd[Classes[[lc,2]]]-LBegin[Classes[[lc,1]]];

      Do[
         lb=lx[bra];
         mb=my[bra];
         nb=nz[bra];
         braell=lb+mb+nb;
         Do[
	    lk=lx[ket];
            mk=my[ket];
            nk=nz[ket];
            ketell=lk+mk+nk;
	    BraAd1=bra;BraAd2=bra;
            KetAd1=ket;KetAd2=ket;
            BraCo1="";BraCo2="";
            KetCo1="";KetCo2="";

            (* Bra addressing scheme for auxiliary SP integrals *)
            (* OFF SET:        1            LenI |  1    LenJ   *)
            (* TARGET CLASSS:  SiS or SiSj  FiS  |  SSj   SFj   *)
            (*                |   Address 1      |  Address 2   *)

            BK00=StringJoin[ToString[bra],",",ToString[ket],","];
	    WS[StringJoin["  HRR(",BK00,"1)=HRR(",BK00,"1)+VRR(",BK00,"0)"]];
            (* ----------------------------------------------------------------------------------------*)
            If[ic==2 && jc!=2 && braell <= Classes[[jc,2]] , BraCo1="SpFnB*" ; BraAd1=BEnd+bra];
            If[kc==2 && lc!=2 && ketell <= Classes[[lc,2]] , KetCo1="SpFnK*" ; KetAd1=KEnd+ket];

            If[ic!=2 && jc==2 && braell == Classes[[ic,1]] , BraCo2="FnSpB*" ; BraAd2=BEnd+bra-LenI1];
            If[kc!=2 && lc==2 && ketell == Classes[[kc,1]] , KetCo2="FnSpK*" ; KetAd2=KEnd+ket-LenK1];

	    If[ BraCo1!=""&&KetCo1=="",
                BK11=StringJoin[ToString[BraAd1],",",ToString[KetAd1],","];
	    WS[StringJoin["  HRR(",BK11,"1)=HRR(",BK11,"1)+",BraCo1,KetCo1,"VRR(",BK00,"0)"]];
              ];

            If[BraCo1==""&&KetCo1!="",
                BK11=StringJoin[ToString[BraAd1],",",ToString[KetAd1],","];
	    WS[StringJoin["  HRR(",BK11,"1)=HRR(",BK11,"1)+",BraCo1,KetCo1,"VRR(",BK00,"0)"]];
              ];

	    If[BraCo2!=""&&KetCo2=="",
                BK22=StringJoin[ToString[BraAd2],",",ToString[KetAd2],","];
	    WS[StringJoin["  HRR(",BK22,"1)=HRR(",BK22,"1)+",BraCo2,KetCo2,"VRR(",BK00,"0)"]];
              ];

	    If[BraCo2==""&&KetCo2!="",
                BK22=StringJoin[ToString[BraAd2],",",ToString[KetAd2],","];
	    WS[StringJoin["  HRR(",BK22,"1)=HRR(",BK22,"1)+",BraCo2,KetCo2,"VRR(",BK00,"0)"]];
              ];

	    (* ----------------------------------------------------------------------------------------*)
	    If[ BraCo2 != "" && KetCo2 != "" || BraCo1 != "" && KetCo1 != ""  ,
                BK12=StringJoin[ToString[BraAd1],",",ToString[KetAd2],","];
                BK21=StringJoin[ToString[BraAd2],",",ToString[KetAd1],","];
                WS[StringJoin["  HRR(",BK12,"1)=HRR(",BK12,"1)+",BraCo1,KetCo2,"VRR(",BK00,"0)"]]
                WS[StringJoin["  HRR(",BK21,"1)=HRR(",BK21,"1)+",BraCo2,KetCo1,"VRR(",BK00,"0)"]]
              ];

	    (* ----------------------------------------------------------------------------------------*)
            BraCo1="";
            BraCo2="";
            KetCo1="";
            KetCo2="";
	    (* ----------------------------------------------------------------------------------------*)
            If[ic==2 && jc==2 ,
               If[braell == 0               , BraCo1="SpSpB*" ; BraAd1=BEnd+bra;
                                              BraCo2="SpFnB*" ; BraAd2=BEnd+LenI+bra+1];
               If[braell == Classes[[jc,2]] , BraCo1="FnSpB*" ; BraAd1=BEnd+bra;
                                              BraCo2="SpFnB*" ; BraAd2=BEnd+LenI+bra+1];
               ];
            If[kc==2 && lc==2,
               If[ketell == 0               , KetCo1="SpSpK*" ; KetAd1=KEnd+ket;
                                              KetCo2="SpFnK*" ; KetAd2=KEnd+LenK+ket+1];
               If[ketell == Classes[[lc,2]] , KetCo1="FnSpK*" ; KetAd1=KEnd+ket;
                                              KetCo2="SpFnK*" ; KetAd2=KEnd+LenK+ket+1];
	      ];
	    If[ BraCo1 != "" && KetCo2 != "" , 
                BK12=StringJoin[ToString[BraAd1],",",ToString[KetAd2],","];
                WS[StringJoin["  HRR(",BK12,"1)=HRR(",BK12,"1)+",BraCo1,KetCo2,"VRR(",BK00,"0)"]]
              ];
	    If[ BraCo1 != "" && KetCo1 != "" , 
                BK11=StringJoin[ToString[BraAd1],",",ToString[KetAd1],","];
                WS[StringJoin["  HRR(",BK11,"1)=HRR(",BK11,"1)+",BraCo1,KetCo1,"VRR(",BK00,"0)"]]
              ];
	    If[ BraCo2 != "" && KetCo1 != "" , 
                BK21=StringJoin[ToString[BraAd2],",",ToString[KetAd1],","];
                WS[StringJoin["  HRR(",BK21,"1)=HRR(",BK21,"1)+",BraCo2,KetCo1,"VRR(",BK00,"0)"]]
              ];
	    If[ BraCo2 != "" && KetCo2 != "" , 
                BK22=StringJoin[ToString[BraAd2],",",ToString[KetAd2],","];
                WS[StringJoin["  HRR(",BK22,"1)=HRR(",BK22,"1)+",BraCo2,KetCo2,"VRR(",BK00,"0)"]]
              ];
	    (* ----------------------------------------------------------------------------------------*)
        ,{ket,1,LEnd[LKet]}]
        ,{bra,1,LEnd[LBra]}]; 

   WS[StringJoin["END SUBROUTINE ",ContractName]];
]

];

