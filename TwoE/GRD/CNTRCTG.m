
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

      If[iell==0&&kell==2,
         WS[StringJoin["VRR(1,5,",m1,")=r1x2E*VRR(1,1,",m1,")+QCx*VRR(1,2,",m1,")-r1x2E*ZxZpE*VRR(1,1,",ms,")+WQx*VRR(1,2,",ms,")"]];
         WS[StringJoin["VRR(1,6,",m1,")=QCx*VRR(1,3,",m1,")+WQx*VRR(1,3,",ms,")"]];
         WS[StringJoin["VRR(1,7,",m1,")=r1x2E*VRR(1,1,",m1,")+QCy*VRR(1,3,",m1,")-r1x2E*ZxZpE*VRR(1,1,",ms,")+WQy*VRR(1,3,",ms,")"]];
         WS[StringJoin["VRR(1,8,",m1,")=QCx*VRR(1,4,",m1,")+WQx*VRR(1,4,",ms,")"]];
         WS[StringJoin["VRR(1,9,",m1,")=QCy*VRR(1,4,",m1,")+WQy*VRR(1,4,",ms,")"]];
         WS[StringJoin["VRR(1,10,",m1,")=r1x2E*VRR(1,1,",m1,")+QCz*VRR(1,4,",m1,")-r1x2E*ZxZpE*VRR(1,1,",ms,")+WQz*VRR(1,4,",ms,")"]];
	 ];

      If[(iell>2&&kell==0) || (iell==0 && kell>2) || (iell>1 && kell>=1) || (iell >=1 && kell>1 ),
         WS[StringJoin["CALL ",VRRSubName,"(",CLenBra,",",CLenKet,",VRR(1,1,",ToString[m-1],"),VRR(1,1,",ToString[m],"))"]]; 
        ]

      (*
      WS[StringJoin["CALL ",VRRSubName,"(",CLenBra,",",CLenKet,",VRR(1,1,",ToString[m-1],"),VRR(1,1,",ToString[m],"))"]]; 
       *)

     ,{m,BraEll+KetEll-(iell+kell)+2,1,-1}]; 

];
,{iell,0,BraEll+1}];
,{kell,0,KetEll+1}];
];

(*======================--  VRR 2 HRR CONTRACTION STEP (A FUCKING BITCH IF SP) --===============================*)

PunchContractCalls[Subroutine_,ic_,jc_,kc_,lc_]:=Block[{WS,BKString,bra,ket},
      WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];
      WS["      ! Contracting ... "] ;

      imax = Classes[[ic, 2]];
      jmax = Classes[[jc, 2]];
      kmax = Classes[[kc, 2]];
      lmax = Classes[[lc, 2]];
      LMax=imax+jmax+kmax+lmax;
      IClass=ToString[IntegralClass[Classes[[ic]]]];
      JClass=ToString[IntegralClass[Classes[[jc]]]];
      KClass=ToString[IntegralClass[Classes[[kc]]]];
      LClass=ToString[IntegralClass[Classes[[lc]]]];
      ContractName=StringJoin["CNTRCTG",IClass,JClass,KClass,LClass];
      (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
      If[DoStress==0,
        WS[StringJoin["      CALL ",ContractName,"(VRR,HRR,Alpha,HRRA,Beta,HRRB,Gamma,HRRC,HRRS(1,1,1,1),PQIJ(1))"]];
      ,
        WS[StringJoin["      CALL ",ContractName,"(VRR,HRR,Alpha,HRRA,Beta,HRRB,Gamma,HRRC)"]];
      ];
      (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
 ];

PunchVRRContract[Subroutine_,ic_,jc_,kc_,lc_,LBra_,LKet_]:=Block[{WS,BKString,bra,ket},

      imax = Classes[[ic, 2]];
      jmax = Classes[[jc, 2]];
      kmax = Classes[[kc, 2]];
      lmax = Classes[[lc, 2]];
      LMax=imax+jmax+kmax+lmax;

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
   
      BEnd1=LEnd[imax+jmax+1];
      KEnd1=LEnd[kmax+lmax+1];
      BraMax1=BEnd1;
      KetMax1=KEnd1;
      LenBra1=LEnd[LBra+1];
      LenKet1=LEnd[LKet+1];

      LenI=LEnd[Classes[[ic,2]]]-LBegin[Classes[[ic,1]]];
      LenI1=LBegin[Classes[[ic,2]]]-1;
      LenJ=LEnd[Classes[[jc,2]]]-LBegin[Classes[[jc,1]]];
      LenK=LEnd[Classes[[kc,2]]]-LBegin[Classes[[kc,1]]];
      LenK1=LBegin[Classes[[kc,2]]]-1;
      LenL=LEnd[Classes[[lc,2]]]-LBegin[Classes[[lc,1]]];

      ContractName=StringJoin["CNTRCTG",IClass,JClass,KClass,LClass];
      (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
      If[DoStress==0,
        WS[StringJoin["SUBROUTINE ",ContractName,"(VRR,HRR,Alpha,HRRA,Beta,HRRB,Gamma,HRRC,HRRS,PQIJ)"]];
      ,
        WS[StringJoin["SUBROUTINE ",ContractName,"(VRR,HRR,Alpha,HRRA,Beta,HRRB,Gamma,HRRC)"]];
      ];
      (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
      WS["  USE DerivedTypes"];
      WS["  USE VScratchB"];
      (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
      If[DoStress==0,
        WS["  IMPLICIT NONE"];
      ];
      (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
      WS["  INTEGER :: K"];
      WS["  REAL(DOUBLE)  :: Alpha,Beta,Gamma"];
      WS[StringJoin["  REAL(DOUBLE), DIMENSION(",ToString[BraMax],",",ToString[KetMax],",",ToString[LEnd[lmax]],") :: HRR "]];
      WS[StringJoin["  REAL(DOUBLE), DIMENSION(",ToString[BraMax1],",",ToString[KetMax],",",ToString[LEnd[lmax]],") :: HRRA,HRRB "]];
      WS[StringJoin["  REAL(DOUBLE), DIMENSION(",ToString[BraMax],",",ToString[KetMax1],",",ToString[LEnd[lmax]],") :: HRRC "]];
      WS[StringJoin["  REAL(DOUBLE)  :: VRR(",ToString[LenBra1],",",ToString[LenKet1],",0:",ToString[LBra+LKet+1],")"]];
      (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
      If[DoStress==0,
        WS[StringJoin["  REAL(DOUBLE)  :: HRRS(",ToString[BraMax],",",ToString[KetMax],",",ToString[LEnd[lmax]],",9),PQIJ(9)"]];
        WS[StringJoin["  INTEGER :: IJ,K"]];
      ];
      (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)

      If[LKet>=2,
         WS[StringJoin["  DO K=1,",ToString[LEnd[LKet]]]];
         Do[Do[
                BK00=StringJoin[ToString[bra],",K,"];
                If[LB<=LBra,
                   WS[StringJoin["     HRR(",BK00,"1)=HRR(",BK00,"1)+VRR(",BK00,"0)"]];
      (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
                   If[DoStress==0,
                     WS[StringJoin["     DO IJ=1,9"]];
                     WS[StringJoin["       HRRS(",BK00,"1,IJ)=HRRS(",BK00,"1,IJ)+PQIJ(IJ)*VRR(",BK00,"1)"]];
                     WS[StringJoin["     ENDDO"]];
                   ];
      (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
                   WS[StringJoin["     HRRC(",BK00,"1)=HRRC(",BK00,"1)+Gamma*VRR(",BK00,"0)"]];
                  ];
                WS[StringJoin["     HRRA(",BK00,"1)=HRRA(",BK00,"1)+Alpha*VRR(",BK00,"0)"]];
                WS[StringJoin["     HRRB(",BK00,"1)=HRRB(",BK00,"1)+Beta*VRR(",BK00,"0)"]];
              ,{bra,LBegin[LB],LEnd[LB]}],{LB,0,LBra+1}]; 
         WS["  ENDDO"];
         WS[StringJoin["  DO K=",ToString[LBegin[LKet+1]],",",ToString[LEnd[LKet+1]]]];
         Do[Do[
               BK00=StringJoin[ToString[bra],",K,"];
               WS[StringJoin["     HRRC(",BK00,"1)=HRRC(",BK00,"1)+Gamma*VRR(",BK00,"0)"]];
	      ,{bra,LBegin[LB],LEnd[LB]}],{LB,0,LBra}]; 
         WS["  ENDDO"];
	 ];
    (* ==========================================================================================*)

      Do[Do[
         lb=lx[bra];
         mb=my[bra];
         nb=nz[bra];
         braell=lb+mb+nb;
         Do[Do[
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

            If[LKet<2,

               If[bra<=BraMax && ket<=KetMax,
                  WS[StringJoin["  HRR(",BK00,"1)=HRR(",BK00,"1)+VRR(",BK00,"0)"]];
      (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
                  If[DoStress==0,
                    WS[StringJoin["  DO IJ=1,9"]];
                    WS[StringJoin["    HRRS(",BK00,"1,IJ)=HRRS(",BK00,"1,IJ)+PQIJ(IJ)*VRR(",BK00,"1)"]];
                    WS[StringJoin["  ENDDO"]];
                  ];
      (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
                 ];

               If[bra<=BraMax1 && ket<=KetMax,
                  WS[StringJoin["  HRRA(",BK00,"1)=HRRA(",BK00,"1)+Alpha*VRR(",BK00,"0)"]];
                  WS[StringJoin["  HRRB(",BK00,"1)=HRRB(",BK00,"1)+Beta*VRR(",BK00,"0)"]];
                 ];

                If[bra<=BraMax && ket<=KetMax1,
                  WS[StringJoin["  HRRC(",BK00,"1)=HRRC(",BK00,"1)+Gamma*VRR(",BK00,"0)"]];
                 ];

              ]; 


	    ,{ket,LBegin[LK],LEnd[LK]}],{LK,0,LKet+1}]
            ,{bra,LBegin[LB],LEnd[LB]}],{LB,0,LBra+1}]; 

   WS[StringJoin["END SUBROUTINE ",ContractName]];

];

