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
!    Author: Matt Challacombe and Valery "wheels" Weber
!    WRITE EXPLICIT CODE FOR COMPUTATION OF THE ECP ANGULAR INTEGRALS OF THE 
!    FIRST KIND, OMEGA(LAMBDA,I,J,K)
!------------------------------------------------------------------------------
*)

(* GET THE MAX ANGULAR SYMMETRY TO BE USED *)

MondoHome=Environment["MONDO_HOME"];
If[MondoHome==$FAILED,
   Print["COULD NOT FIND $MONDO_HOME! CHECK YOUR .cshrc "]; 
   Abort[];
  ,
   Print[" Using ",MondoHome," as your home directory"];
  ];
Print[" Using ",MondoHome," as your home directory"];
EllFile = StringJoin[MondoHome,"/Includes/Ell.m"];
Get[EllFile];

(***************** SOME PRELIMINARY INDEXING DEFINITIONS ***************)

LBegin[L_]:=(L*(L+1)*(L+2))/6+1;
LEnd[L_]:=LBegin[L]+L*(L+1)/2+L;
LMNDex[L_, M_, N_] := LBegin[L + M + N] + N*(2*(L + M + N) - N + 3)/2 + M;

IntegralClass[Ell_List] := Ell[[2]]*(Ell[[2]] + 1)/2 + Ell[[1]] + 1;

Normy[LMN_List]:=Module[{Fct,X,Y,L,M,N},
			Return[1]; 
                        L=LMN[[1]];
                        M=LMN[[2]];
                        N=LMN[[3]];
                        Fct[n_]:=(2 n-1)!!;
                        X=Fct[L+M+N];
                        Y=Fct[L]*Fct[M]*Fct[N];
			Return[Sqrt[X/Y]]];

(* Minimal 
   Classes = { {0,0},{1,1}} 
 *)

Classes = { {0,0},{0,1},{1,1},{2,2},{3,3}};

(* Maximal    Classes = { {0,0},{0,1},{1,1},{2,2},{3,3}}  *)

CType[1]  = "s";
CType[2]  = "sp";
CType[3]  = "p";
CType[6]  = "d";
CType[10] = "f";
CType[15] = "g";
CType[21] = "h";
CType[28] = "i";
CType[36] = "j";
CType[45] = "k";

LC=Length[Classes];
MxEll = Max[Classes];

Print[" LC = ",LC];

Do[Do[Do[ 
         lmn=LMNDex[l, m, n];
         lx[lmn] = l;
         my[lmn] = m;
         nz[lmn] = n;
,{l,0,2 MxEll-m-n}]
,{m,0,2 MxEll-n}]
,{n,0,2 MxEll}];


(**************** THE VRR RELATIONS: THX OS ! ************************)



(*=============================== 
                                  LOAD OPTIMIZING AND FORMATING ROUTINES FOR MMA AND SET ASSOCIATED OPTIONS
                                                                                                            =============================== *)

Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];
Get[StringJoin[MondoHome,"/MMA/Format.m"]];
Get[StringJoin[MondoHome,"/MMA/Optimize.m"]];

FF[x_] := ToString[FixedNumberForm[SetPrecision[N[x,32],32], 16, 2]];

SetOptions[FortranAssign,AssignOptimize->False,AssignMaxSize->5000,AssignBreak->{5000,"         "},AssignIndent->"      ",AssignTemporary->{W,Sequence}];
SetOptions[Optimize,OptimizeVariable->{V,Array},OptimizeTimes->True,OptimizePlus->True,OptimizeCoefficients->True,OptimizeFunction->False]; 
SetOptions[OpenWrite, PageWidth -> 200];

SetAttributes[o,NHoldAll];
SetAttributes[AuxR,NHoldAll];
SetAttributes[VRR0,NHoldAll];
SetAttributes[VRR1,NHoldAll];
SetAttributes[VRR2,NHoldAll];
SetAttributes[HRR,NHoldAll];
SetAttributes[SSSS,NHoldAll];
SetAttributes[MBarN,NHoldAll];

(*======================================--  ROUTINE FOR ASSOCIATED INCOMPLETE GAMMA FUNCTIONS F_M(T)  --============================================*)

PunchGammas[Subroutine_,LTot_]:=Block[{WS,FStr,Gammas},

WS[String_]:=WriteString[Subroutine,"            ",String,"\n"];

FStr[L_]:=Module[{LSt},
                 LSt=ToString[L];
                 StringJoin["(F",LSt,"_0(L)+T*(F",LSt,"_1(L)+T*(F",LSt,"_2(L)+T*(F",LSt,"_3(L)+T*F",LSt,"_4(L)))))"]];


Gammas[EllTot_]:=Module[{LSt,LSt2,GSt,GSt1,GFlt,GammAss},

              GammAss[n_] := FF[ Abs[2 n - 1]!!/(2 (2 )^n) Sqrt[Pi] ];

                                WS["IF(T<Gamma_Switch)THEN"];
                                WS["  L=AINT(T*Gamma_Grid)"];			
         If[EllTot==0,
                     WS[StringJoin["  VRR(1,1,0)=Upq*",FStr[0]]];
            ];

         If[EllTot==1,
                     WS[StringJoin["  VRR(1,1,0)=Upq*",FStr[0]]];
                     WS[StringJoin["  VRR(1,1,1)=Upq*",FStr[1]]];
            ];
          If[EllTot>1,
             LSt=ToString[EllTot];
                     WS["  ET=EXP(-T)"];
                     WS["  TwoT=Two*T"];
          WS[StringJoin["  W",LSt,"=",FStr[LTot]]];
                         Do[GSt=ToString[L-1];
                            GSt1=ToString[L];
                            GFlt=FF[1/(2 *(L-1)+1)];
          WS[StringJoin["  W",GSt,"=",GFlt,"*(TwoT*W",GSt1,"+ET)"]];
                           ,{L,EllTot,2,-1}];
                     WS["  W0=TwoT*W1+ET"];
                         Do[GSt=ToString[L];
          WS[StringJoin["  VRR(1,1,",GSt,")=Upq*W",GSt,""]];
                           ,{L,0,EllTot}];
            ]; 

                     WS["ELSE"];
                     WS["  InvT=One/T"];
                     WS["  SqInvT=DSQRT(InvT)"];
                         Do[GSt=ToString[L];

           WS[StringJoin["  VRR(1,1,",GSt,")=",GammAss[L],"*Upq*SqInvT"]];

If[L<LTot,WS[StringJoin["  SqInvT=SqInvT*InvT"]]];
                           ,{L,0,EllTot}];

                      WS["ENDIF"];

      ];

Gammas[LTot];
];

(*======================================--  THE OBARA-SAIKA VERTICAL RECURENCE RELATIONS  --============================================*)

VRR2[a_List,c_List,m_,as_,cs_]:=Module[{p,q,PA,QC,WP,WQ,one,two,a1,a2,c1,c2,Ai1,Ci1,CiO2zn,AiO2zn,lmna,lmnc}, 

pmin=Position[a,Min[a]][[1, 1]];
qmin=Position[c,Min[c]][[1, 1]];
If[ a[[pmin]] < 0 || c[[qmin]] < 0, Return[0] ];                                        

TotA=Sum[a[[i]],{i,3}];
TotC=Sum[c[[i]],{i,3}];
TotL=TotA+TotC;

If[TotA<=as && TotC<=cs,
   lmna=LMNDex[a[[1]],a[[2]],a[[3]]];                 
   lmnc=LMNDex[c[[1]],c[[2]],c[[3]]];
   Return[ToExpression[StringJoin[ "VRR",ToString[m],"[",ToString[lmna],",",ToString[lmnc],"]"]]]
  ];

p=Position[a,Max[a]][[1, 1]];
q=Position[c,Max[c]][[1, 1]];
If[a[[p]]>c[[q]],
   (* Taking down ell on bra side *)
   If[p == 1, one = {1, 0, 0}; two={2,0,0}; PA=PAx; WP=WPx; ];
   If[p == 2, one = {0, 1, 0}; two={0,2,0}; PA=PAy; WP=WPy; ];
   If[p == 3, one = {0, 0, 1}; two={0,0,2}; PA=PAz; WP=WPz; ];
   Ai1=(a[[p]]-1)*r1x2Z;   (* (a[i]-1)/(2 zeta)   *)
   CiO2zn=c[[p]]*HfxZpE;   (* c[i]/(2(zeta+eta)) *)
   a1 = a - one;
   a2 = a - two;
   c1 = c - one;                                  
   Return[PA*VRR[a1,c,m,as,cs]+WP*VRR[a1,c,m+1,as,cs]+Ai1*(VRR[a2,c,m,as,cs]-ExZpE*VRR[a2,c,m+1,as,cs]) +CiO2zn*VRR[a1,c1,m+1,as,cs]]
  ,(* Taking down ell on ket side *)
   If[q == 1, one = {1, 0, 0}; two={2,0,0}; QC=QCx; WQ=WQx;];
   If[q == 2, one = {0, 1, 0}; two={0,2,0}; QC=QCy; WQ=WQy;];
   If[q == 3, one = {0, 0, 1}; two={0,0,2}; QC=QCz; WQ=WQz;];
   Ci1=(c[[q]]-1)*r1x2E ; (* (c[i]-1)/(2 eta)   *)
   AiO2zn=a[[q]]*HfxZpE; (* a[i]/(2(zeta+eta)) *)
   c1 = c - one;
   c2 = c - two;
   a1 = a - one;  
   Return[QC*VRR[a,c1,m,as,cs]+WQ*VRR[a,c1,m+1,as,cs]+Ci1*(VRR[a,c2,m,as,cs]-ZxZpE*VRR[a,c2,m+1,as,cs])+AiO2zn*VRR[a1,c1,m+1,as,cs]]
  ];
];

VRR[a_List,c_List,m_,Ell_]:=Module[{p,q,PA,QC,WP,WQ,one,two,a1,a2,c1,c2,Ai1,Ci1,CiO2zn,AiO2zn,lmna,lmnc}, 

pmin=Position[a,Min[a]][[1, 1]];
qmin=Position[c,Min[c]][[1, 1]];
If[ a[[pmin]] < 0 || c[[qmin]] < 0, Return[0] ];                                        

TotA=Sum[a[[i]],{i,3}];
TotC=Sum[c[[i]],{i,3}];
TotL=TotA+TotC;

If[TotL<=Ell,
   lmna=LMNDex[a[[1]],a[[2]],a[[3]]];                 
   lmnc=LMNDex[c[[1]],c[[2]],c[[3]]];
   Return[ToExpression[StringJoin[ "VRR",ToString[m],"[",ToString[lmna],",",ToString[lmnc],"]"]]]
  ];

p=Position[a,Max[a]][[1, 1]];
q=Position[c,Max[c]][[1, 1]];
If[a[[p]]>c[[q]],
   (* Taking down ell on bra side *)
   If[p == 1, one = {1, 0, 0}; two={2,0,0}; PA=PAx; WP=WPx; ];
   If[p == 2, one = {0, 1, 0}; two={0,2,0}; PA=PAy; WP=WPy; ];
   If[p == 3, one = {0, 0, 1}; two={0,0,2}; PA=PAz; WP=WPz; ];
   Ai1=(a[[p]]-1)*r1x2Z;   (* (a[i]-1)/(2 zeta)   *)
   CiO2zn=c[[p]]*HfxZpE;   (* c[i]/(2(zeta+eta)) *)
   a1 = a - one;
   a2 = a - two;
   c1 = c - one;                                  
   Return[PA*VRR[a1,c,m,Ell]+WP*VRR[a1,c,m+1,Ell]+Ai1*(VRR[a2,c,m,Ell]-ExZpE*VRR[a2,c,m+1,Ell]) +CiO2zn*VRR[a1,c1,m+1,Ell]]
  ,(* Taking down ell on ket side *)
   If[q == 1, one = {1, 0, 0}; two={2,0,0}; QC=QCx; WQ=WQx;];
   If[q == 2, one = {0, 1, 0}; two={0,2,0}; QC=QCy; WQ=WQy;];
   If[q == 3, one = {0, 0, 1}; two={0,0,2}; QC=QCz; WQ=WQz;];
   Ci1=(c[[q]]-1)*r1x2E ; (* (c[i]-1)/(2 eta)   *)
   AiO2zn=a[[q]]*HfxZpE; (* a[i]/(2(zeta+eta)) *)
   c1 = c - one;
   c2 = c - two;
   a1 = a - one;  
   Return[QC*VRR[a,c1,m,Ell]+WQ*VRR[a,c1,m+1,Ell]+Ci1*(VRR[a,c2,m,Ell]-ZxZpE*VRR[a,c2,m+1,Ell])+AiO2zn*VRR[a1,c1,m+1,Ell]]
  ];
];

(* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  *)

PunchVRRCalls[FileName_,BraEll_,KetEll_]:=Module[{spaces,CLenBra,CLenKet,CLenEhm},
spaces="            ";	 
CLenBra=ToString[LEnd[BraEll]];
CLenKet=ToString[LEnd[KetEll]];
CLenEhm=ToString[BraEll+KetEll];
Do[Do[
If[iell+kell>0,
   VRRSubName=StringJoin["VRR",CType[IntegralClass[{iell,iell}]],"0",CType[IntegralClass[{kell,kell}]],"0"];
   Do[
      WriteString[FileName,StringJoin[spaces,"! Generating (",CType[IntegralClass[{iell,iell}]],"0|", \
                                                              CType[IntegralClass[{kell,kell}]],"0)^(",ToString[m-1],") \n"]]; 
      WriteString[FileName,StringJoin[spaces,"CALL ",VRRSubName,"(",CLenBra,",",CLenKet,",VRR(1,1,",ToString[m-1],"),VRR(1,1,",ToString[m],")) \n"]]; 
     ,{m,BraEll+KetEll-(iell+kell)+1,1,-1}];
  ]
,{iell,0,BraEll}];
,{kell,0,KetEll}];
];

(* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  *)

PunchVRRFront[Subroutine_,Subname_]:=Block[{WS,BKType,LenBra,LenKet},
WriteString[Subroutine,StringJoin["   SUBROUTINE ",Subname,"(LB,LK,VRR0,VRR1) \n"]];
WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];
WS["USE DerivedTypes"];
WS["USE VScratchB"];
WS["USE GlobalScalars"];
WS["IMPLICIT REAL(DOUBLE) (W)"]; 
WS["INTEGER :: LB,LK"];
WS[StringJoin["REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1"]];
];

(* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  *)

PunchVRRClass[BraEll_,KetEll_]:=Module[{oList,IList,Kount,a,c,CollectList},

spaces="            ";	 
TStart=TimeUsed[];
TStartTotal=TStart;
(* CollectList=Flatten[Table[VRR[i,j,m],{i,BraEll},{j,KetEll},{m,BraEll+KetEll}]]; *)
MaxSameEll=Min[BraEll,KetEll];
Do[Do[
If[iell+kell>0,
   TStart=TimeUsed[];
   as=Max[0,iell-1];
   cs=Max[0,kell-1];   
   TStart=TimeUsed[];oList={" "->""};IList={};Kount = 0;
   Do[Do[
         Kount=Kount+1;
         a={lx[i], my[i], nz[i]};
         c={lx[k], my[k], nz[k]};
         vrr=VRR[a,c,0, kell+iell-1];

(*       vrrList=Coefficient[vrr,CollectList];
         vrr=Sum[co[vrrList[[i]]] CollectList[[i]],{i,Length[CollectList]}];*)
         IList=Append[IList,vrr];
         MBarString=StringJoin["VRR0(",ToString[i],",",ToString[k],")"];
         oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->MBarString];
    ,{k,LBegin[kell],LEnd[kell]}]
    ,{i,LBegin[iell],LEnd[iell]}];

   VRRSubName=StringJoin["VRR",CType[IntegralClass[{iell,iell}]],"0",CType[IntegralClass[{kell,kell}]],"0"];
   MakeList=Append[MakeList,StringJoin[VRRSubName,".o"]];
   RelsList=Append[RelsList,StringJoin[VRRSubName,".x"]];
   VRRSubroutine=StringJoin[VRRSubName,".F90"];
   OpenWrite[VRRSubroutine];
   PunchVRRFront[VRRSubroutine,VRRSubName];
   Write[VRRSubroutine,FortranAssign[o,IList,AssignReplace->oList]];
   WriteString[VRRSubroutine,StringJoin["END SUBROUTINE ",VRRSubName]];			       
   Close[VRRSubroutine];
   Print["Done with ",VRRSubName," Time = ",TimeUsed[]-TStart];
  ]
,{iell,0,BraEll}];
,{kell,0,KetEll}];
];

(*======================================--  VRR 2 HRR CONTRACTION STEP (A FUCKING BITCH IF SP) --==========================================*)

PunchVRRContract[Subroutine_,ic_,jc_,kc_,lc_,LBra_,LKet_]:=Block[{WS,BKString,bra,ket},

      WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];

      WS["      ! Contracting ... "] ;

      BEnd=LEnd[LBra];
      KEnd=LEnd[LKet];
      LenI=LEnd[Classes[[ic,2]]]-LBegin[Classes[[ic,1]]];
      LenI1=LBegin[Classes[[ic,2]]]-1;
      LenJ=LEnd[Classes[[jc,2]]]-LBegin[Classes[[jc,1]]];
      LenK=LEnd[Classes[[kc,2]]]-LBegin[Classes[[kc,1]]];
      LenK1=LBegin[Classes[[kc,2]]]-1;
      LenL=LEnd[Classes[[lc,2]]]-LBegin[Classes[[lc,1]]];

      WS[StringJoin["      CALL DBLAXPY(",ToString[BEnd*KEnd],",HRR(1,1,1),VRR(1,1,0))"]];


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
	    (*
            WS[StringJoin["      HRR(",BK00,"1)=HRR(",BK00,"1)+VRR(",BK00,"0)"]];
	     *)

	    (* ----------------------------------------------------------------------------------------*)
            If[ic==2 && jc!=2 && braell <= Classes[[jc,2]] , BraCo1="SpFnB*" ; BraAd1=BEnd+bra];
            If[kc==2 && lc!=2 && ketell <= Classes[[lc,2]] , KetCo1="SpFnK*" ; KetAd1=KEnd+ket];

            If[ic!=2 && jc==2 && braell == Classes[[ic,1]] , BraCo2="FnSpB*" ; BraAd2=BEnd+bra-LenI1];
            If[kc!=2 && lc==2 && ketell == Classes[[kc,1]] , KetCo2="FnSpK*" ; KetAd2=KEnd+ket-LenK1];

	    If[ BraCo1!=""&&KetCo1=="",
                BK11=StringJoin[ToString[BraAd1],",",ToString[KetAd1],","];
	    WS[StringJoin["      HRR(",BK11,"1)=HRR(",BK11,"1)+",BraCo1,KetCo1,"VRR(",BK00,"0)"]];
              ];

            If[BraCo1==""&&KetCo1!="",
                BK11=StringJoin[ToString[BraAd1],",",ToString[KetAd1],","];
	    WS[StringJoin["      HRR(",BK11,"1)=HRR(",BK11,"1)+",BraCo1,KetCo1,"VRR(",BK00,"0)"]];
              ];

	    If[BraCo2!=""&&KetCo2=="",
                BK22=StringJoin[ToString[BraAd2],",",ToString[KetAd2],","];
	    WS[StringJoin["      HRR(",BK22,"1)=HRR(",BK22,"1)+",BraCo2,KetCo2,"VRR(",BK00,"0)"]];
              ];

	    If[BraCo2==""&&KetCo2!="",
                BK22=StringJoin[ToString[BraAd2],",",ToString[KetAd2],","];
	    WS[StringJoin["      HRR(",BK22,"1)=HRR(",BK22,"1)+",BraCo2,KetCo2,"VRR(",BK00,"0)"]];
              ];

	    (* ----------------------------------------------------------------------------------------*)
	    If[ BraCo2 != "" && KetCo2 != "" || BraCo1 != "" && KetCo1 != ""  ,
                BK12=StringJoin[ToString[BraAd1],",",ToString[KetAd2],","];
                BK21=StringJoin[ToString[BraAd2],",",ToString[KetAd1],","];
                WS[StringJoin["      HRR(",BK12,"1)=HRR(",BK12,"1)+",BraCo1,KetCo2,"VRR(",BK00,"0)"]]
                WS[StringJoin["      HRR(",BK21,"1)=HRR(",BK21,"1)+",BraCo2,KetCo1,"VRR(",BK00,"0)"]]
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
                WS[StringJoin["      HRR(",BK12,"1)=HRR(",BK12,"1)+",BraCo1,KetCo2,"VRR(",BK00,"0)"]]
              ];
	    If[ BraCo1 != "" && KetCo1 != "" , 
                BK11=StringJoin[ToString[BraAd1],",",ToString[KetAd1],","];
                WS[StringJoin["      HRR(",BK11,"1)=HRR(",BK11,"1)+",BraCo1,KetCo1,"VRR(",BK00,"0)"]]
              ];
	    If[ BraCo2 != "" && KetCo1 != "" , 
                BK21=StringJoin[ToString[BraAd2],",",ToString[KetAd1],","];
                WS[StringJoin["      HRR(",BK21,"1)=HRR(",BK21,"1)+",BraCo2,KetCo1,"VRR(",BK00,"0)"]]
              ];
	    If[ BraCo2 != "" && KetCo2 != "" , 
                BK22=StringJoin[ToString[BraAd2],",",ToString[KetAd2],","];
                WS[StringJoin["      HRR(",BK22,"1)=HRR(",BK22,"1)+",BraCo2,KetCo2,"VRR(",BK00,"0)"]]
              ];
	    (* ----------------------------------------------------------------------------------------*)
        ,{ket,1,LEnd[LKet]}]
        ,{bra,1,LEnd[LBra]}]; 
      WS["   ENDDO ! (M0| loop"];            					 
      WS["ENDDO ! |N0) loop"];            					 
];


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

PunchHRRBra[FileName_,ic_,jc_,kc_,lc_]:=Module[{oList,IList,Kount,a,b,c,d},

 imin = Classes[[ic, 1]]; imax = Classes[[ic, 2]];
 jmin = Classes[[jc, 1]]; jmax = Classes[[jc, 2]];
 kmin = Classes[[kc, 1]]; kmax = Classes[[kc, 2]];
 lmin = Classes[[lc, 1]]; lmax = Classes[[lc, 2]];
 LenI=LEnd[Classes[[ic,2]]]-LBegin[Classes[[ic,1]]];
 LenJ=LEnd[Classes[[jc,2]]]-LBegin[Classes[[jc,1]]];
 BEnd=LEnd[imax+jmax];

 IClass=ToString[IntegralClass[Classes[[ic]]]];
 JClass=ToString[IntegralClass[Classes[[jc]]]];
 HRRSubName=StringJoin["HRR_Bra_",IClass,JClass];

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
             If[ ic != 2 && jc == 2 && jl == 0 ,hrr=hrr/.{HRR[a_]:>HRR[BEnd+a-LenI+1]}];
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

  oList=Append[oList,{" "->"","xx"->",","BB"->"(","KK"->"|","In"->"INTGRL"}];
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
];

(*======================================--  KET HRR --============================================*)

PunchKetHRRFront[Subroutine_,Subname_,ic_,jc_,kc_,lc_]:=Block[{WS,BKType,LenBra,LenKet},

           WriteString[Subroutine,StringJoin["   SUBROUTINE ",Subname,"(HRR) \n"]];

           imin = Classes[[ic, 1]]; imax = Classes[[ic, 2]];
           jmin = Classes[[jc, 1]]; jmax = Classes[[jc, 2]];
           kmin = Classes[[kc, 1]]; kmax = Classes[[kc, 2]];
           lmin = Classes[[lc, 1]]; lmax = Classes[[lc, 2]];
           LenI=LEnd[Classes[[ic,2]]]-LBegin[Classes[[ic,1]]]+1;
           LenJ=LEnd[Classes[[jc,2]]]-LBegin[Classes[[jc,1]]]+1;
           LenK=LEnd[Classes[[kc,2]]]-LBegin[Classes[[kc,1]]]+1;
           LenL=LEnd[Classes[[lc,2]]];
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

           WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];
	   WS["USE DerivedTypes"];
	   WS["USE VScratchB"];
	   WS["USE GlobalScalars"];
           WS["IMPLICIT REAL(DOUBLE) (W)"]; 
           WS[StringJoin["REAL(DOUBLE)  :: HRR(",ToString[BraMax],",",ToString[KetMax],",",ToString[LenL],")"]];

];

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

PunchHRRKet[FileName_,ic_,jc_,kc_,lc_]:=Module[{oList,IList,Kount,a,b,c,d},

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
                  
 IList={};
 oList={"+"->"+  & \n        "};

 Kount = 0;
 Do[Do[Do[Do[
             Kount = Kount + 1;
	     
             IList=Append[IList,ToExpression[StringJoin["KK",ToString[k],"xx",ToString[l],"BB"]]];
             oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->"!"];
	     (*             Do[ *)
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
                   HRRAddress=StringJoin["1:",ToString[BraMax],",",ToString[k],",",ToString[l]];
                   oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["HRR(",HRRAddress,")"]];
                  ]; 
		(*		,{i,1,BraMax}] *)
 ,{k,LBegin[kl],LEnd[kl]}]
 ,{l,LBegin[ll],LEnd[ll]}]
 ,{kl,kmax,kmin,-1}]
 ,{ll,lmax,lmin,-1}];

  oList=Append[oList,{" "->"","xx"->",","BB"->")","KK"->"|","In"->"INTGRL","abc"->StringJoin["1:",ToString[BraMax]]}];
  oList=Flatten[oList];

			If[Kount>0,
			   IClass=ToString[IntegralClass[Classes[[ic]]]];
			   JClass=ToString[IntegralClass[Classes[[jc]]]];
			   KClass=ToString[IntegralClass[Classes[[kc]]]];
			   LClass=ToString[IntegralClass[Classes[[lc]]]];
			   HRRSubName=StringJoin["HRR_Ket_",IClass,JClass,KClass,LClass];
			   MakeList=Append[MakeList,StringJoin[HRRSubName,".o"]];
			   RelsList=Append[RelsList,StringJoin[HRRSubName,".x"]];
			   HRRSubroutine=StringJoin[HRRSubName,".F90"];
			   OpenWrite[HRRSubroutine];
			   WriteString[FileName,StringJoin["      ! Generating (",CType[IntegralClass[{imin,imax}]],",0|",
                                                                                  CType[IntegralClass[{kmin,kmax}]],",",
                                                                                  CType[IntegralClass[{lmin,lmax}]],")^(0) \n"]];

                           WriteString[FileName,StringJoin["      CALL ",HRRSubName,"(HRR) \n"]]; 
                           PunchKetHRRFront[HRRSubroutine,HRRSubName,ic,jc,kc,lc];
                           Write[HRRSubroutine,FortranAssign[o,IList,AssignReplace->oList]];
                           WriteString[HRRSubroutine,StringJoin["END SUBROUTINE ",HRRSubName]];			       
                           Close[HRRSubroutine];
                          ];
 ];


PunchFront[Subroutine_,ic_,jc_,kc_,lc_,IJKL_]:=Block[{WS,LBra,LKet,BKType,LenBra,LenKet},

	   imin = Classes[[ic, 1]]; imax = Classes[[ic, 2]];
	   jmin = Classes[[jc, 1]]; jmax = Classes[[jc, 2]];
	   kmin = Classes[[kc, 1]]; kmax = Classes[[kc, 2]];
           lmin = Classes[[lc, 1]]; lmax = Classes[[lc, 2]];

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

           LBra=imax+jmax;
           LKet=kmax+lmax;
	   BKType=100*LBra+LKet;					   

           WriteString[Subroutine,StringJoin["   SUBROUTINE IntB",ToString[IJKL],"(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, & \n", \
                                             "                              OA,LDA,OB,LDB,OC,LDC,OD,LDD,PBC,INTGRL) \n"]];

	   WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];

	   WS["USE DerivedTypes"];
	   WS["USE VScratchB"];
	   WS["USE GlobalScalars"];
           (*WS["USE ONX2DataType"];*)
           WS["USE ShellPairStruct"];
	   If[LBra+LKet==1,
              WS["USE GammaF0"];
              WS["USE GammaF1"];,
              WS[StringJoin["USE GammaF",ToString[LBra+LKet]]]];
           WS["IMPLICIT REAL(DOUBLE) (W)"]; 
           WS["INTEGER        :: LBra,LKet,CDOffSet"];
           WS["REAL(DOUBLE)   :: PrmBufB(8,LBra),PrmBufK(8,LKet)"];
	   WS["TYPE(SmallAtomInfo) :: ACInfo,BDInfo"];
           WS["TYPE(PBCInfo) :: PBC"];
	   LenBra=LEnd[LBra];
           LenKet=LEnd[LKet];
           WS[StringJoin["REAL(DOUBLE)  :: INTGRL(*)"]];
           WS["REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz"];
           WS["REAL(DOUBLE)  :: Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz"];
           WS["REAL(DOUBLE)  :: PQx,PQy,PQz,FPQx,FPQy,FPQz"];
           WS["REAL(DOUBLE)  :: Zeta,Eta,Omega,Up,Uq,Upq"];
           WS["REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT"];
           If[ic==2||jc==2,						     
              WS["REAL(DOUBLE)  :: SpSpB,FnSpB,SpFnB"]];
           If[kc==2||lc==2,						     
              WS["REAL(DOUBLE)  :: SpSpK,FnSpK,SpFnK"]];


           WS[StringJoin["REAL(DOUBLE)  :: VRR(",ToString[LEnd[LBra]],",",ToString[LEnd[LKet]],",0:",ToString[LBra+LKet],")"]];
           WS[StringJoin["REAL(DOUBLE)  :: HRR(",ToString[BraMax],",",ToString[KetMax],",",ToString[LEnd[lmax]],")"]];

           WS["INTEGER       :: OffSet,OA,LDA,OB,LDB,OC,LDC,OD,LDD,I,J,K,L"];
           WS["EXTERNAL InitDbl"];						     
           WS[StringJoin["CALL InitDbl(",ToString[BraMax],"*",ToString[KetMax],",HRR(1,1,1))"]];

           WS["Ax=ACInfo%Atm1X"];
           WS["Ay=ACInfo%Atm1Y"];
           WS["Az=ACInfo%Atm1Z"];
           WS["Bx=ACInfo%Atm2X"];
           WS["By=ACInfo%Atm2Y"];
           WS["Bz=ACInfo%Atm2Z"];
           WS["Cx=BDInfo%Atm1X"];
           WS["Cy=BDInfo%Atm1Y"];
           WS["Cz=BDInfo%Atm1Z"];
           WS["Dx=BDInfo%Atm2X"];
           WS["Dy=BDInfo%Atm2Y"];
           WS["Dz=BDInfo%Atm2Z"];
           WS["ABx=Ax-Bx"];
           WS["ABy=Ay-By"];
           WS["ABz=Az-Bz"];
           WS["CDx=Cx-Dx"];
           WS["CDy=Cy-Dy"];
           WS["CDz=Cz-Dz"];
           WS["DO J=1,LKet ! K^2 VRR |N0) loop "];
           WS["   Eta=PrmBufK(1,J)"];
           WS["   Qx=PrmBufK(2,J)"];
           WS["   Qy=PrmBufK(3,J)"];
           WS["   Qz=PrmBufK(4,J)"];
           WS["   Uq=PrmBufK(5,J)"];
           If[kc==2&&lc!=2,WS["   SpFnK=PrmBufK(6,J)"]];
           If[kc!=2&&lc==2,WS["   FnSpK=PrmBufK(6,J)"]];
           If[kc==2&&lc==2,WS["   SpSpK=PrmBufK(6,J)"];
                           WS["   FnSpK=PrmBufK(7,J)"];
                           WS["   SpFnK=PrmBufK(8,J)"]];
           WS["   QCx=Qx-Cx"];
           WS["   QCy=Qy-Cy"];
           WS["   QCz=Qz-Cz"];
           WS["   DO K=1,LBra ! K^2 VRR (M0| loop "];
           WS["      Zeta=PrmBufB(1,K)"];
           WS["      Px=PrmBufB(2,K)"];
           WS["      Py=PrmBufB(3,K)"];
           WS["      Pz=PrmBufB(4,K)"];
           WS["      Up=PrmBufB(5,K)"];
           If[ic==2&&jc!=2,WS["      SpFnB=PrmBufB(6,K)"]];
           If[ic!=2&&jc==2,WS["      FnSpB=PrmBufB(6,K)"]];
           If[ic==2&&jc==2,WS["      SpSpB=PrmBufB(6,K)"];
                           WS["      FnSpB=PrmBufB(7,K)"];
                           WS["      SpFnB=PrmBufB(8,K)"]];
           WS["      r1xZpE=One/(Zeta+Eta)"];
	   WS["      Upq=SQRT(r1xZpE)*Up*Uq"];					   
           WS["      HfxZpE=Half/(Zeta+Eta)"];
           WS["      r1x2E=Half/Eta"];
           WS["      r1x2Z=Half/Zeta"];
           WS["      ExZpE=Eta*r1xZpE"];
           WS["      ZxZpE=Zeta*r1xZpE"];
           WS["      Omega=Eta*Zeta*r1xZpE"];
           WS["      PAx=Px-Ax"];
           WS["      PAy=Py-Ay"];
           WS["      PAz=Pz-Az"];
           WS["      PQx=Px-Qx"];
           WS["      PQy=Py-Qy"];
           WS["      PQz=Pz-Qz"];
           WS["      ! Begin Minimum Image Convention "];
           WS["      FPQx = PQx*PBC%InvBoxSh%D(1,1)+PQy*PBC%InvBoxSh%D(1,2)+PQz*PBC%InvBoxSh%D(1,3)"];
           WS["      FPQy = PQy*PBC%InvBoxSh%D(2,2)+PQz*PBC%InvBoxSh%D(2,3)"];
           WS["      FPQz = PQz*PBC%InvBoxSh%D(3,3)"];
           WS["      IF(PBC%AutoW%I(1)==1)FPQx=FPQx-ANINT(FPQx-SIGN(1D-15,FPQx))"];
           WS["      IF(PBC%AutoW%I(2)==1)FPQy=FPQy-ANINT(FPQy-SIGN(1D-15,FPQy))"];
           WS["      IF(PBC%AutoW%I(3)==1)FPQz=FPQz-ANINT(FPQz-SIGN(1D-15,FPQz))"];
           WS["      PQx=FPQx*PBC%BoxShape%D(1,1)+FPQy*PBC%BoxShape%D(1,2)+FPQz*PBC%BoxShape%D(1,3)"];
           WS["      PQy=FPQy*PBC%BoxShape%D(2,2)+FPQz*PBC%BoxShape%D(2,3)"];
           WS["      PQz=FPQz*PBC%BoxShape%D(3,3)"];
           WS["      ! End MIC"];
           WS["      WPx = -Eta*PQx*r1xZpE"];
           WS["      WPy = -Eta*PQy*r1xZpE"];
           WS["      WPz = -Eta*PQz*r1xZpE"];
           WS["      WQx = Zeta*PQx*r1xZpE"];
           WS["      WQy = Zeta*PQy*r1xZpE"];
           WS["      WQz = Zeta*PQz*r1xZpE"];
           WS["      T=Omega*(PQx*PQx+PQy*PQy+PQz*PQz)"];
           PunchGammas[Subroutine,LBra+LKet];
];

MakeList={};
RelsList={};

PunchVRRClass[MxEll*2,MxEll*2];

IncludeFile="ERIInclude.Inc";
OpenWrite[IncludeFile];
Print[" Openned ",IncludeFile];
WSI[String_]:=WriteString[IncludeFile,"   ",String,"\n"];

WSI["SELECT CASE(IntType)"];

Do[Do[Do[Do[


   If[IntegralClass[Classes[[ic]]]>=IntegralClass[Classes[[jc]]]&& \
      IntegralClass[Classes[[kc]]]>=IntegralClass[Classes[[lc]]]&& \
      IntegralClass[Classes[[ic]]]*100+IntegralClass[Classes[[jc]]]>= \
      IntegralClass[Classes[[kc]]]*100+IntegralClass[Classes[[lc]]],


            CommentLine=StringJoin["(",CType[IntegralClass[Classes[[ic]]]]," ", \
                                       CType[IntegralClass[Classes[[jc]]]],"|", \
                                       CType[IntegralClass[Classes[[kc]]]]," ", \
                                       CType[IntegralClass[Classes[[lc]]]],")"];
	    CommentLine=StringJoin["! ---------------------------------------------------------- \n", \
                                   "! COMPUTES THE INTEGRAL CLASS ",CommentLine," \n", \
                                   "! ---------------------------------------------------------- \n"];                                             

	    imin = Classes[[ic, 1]]; imax = Classes[[ic, 2]];
	    jmin = Classes[[jc, 1]]; jmax = Classes[[jc, 2]];
	    kmin = Classes[[kc, 1]]; kmax = Classes[[kc, 2]];
	    lmin = Classes[[lc, 1]]; lmax = Classes[[lc, 2]];


          ijklFlag=1000000*IntegralClass[Classes[[ic]]] \
                  +10000*IntegralClass[Classes[[jc]]] \
	          +100*IntegralClass[Classes[[kc]]] \
	          +IntegralClass[Classes[[lc]]];


            ijklType=1000*IntegralClass[Classes[[ic]]] \
                      +100*IntegralClass[Classes[[jc]]] \
                        +10*IntegralClass[Classes[[kc]]] \
	                    +IntegralClass[Classes[[lc]]];


             WSI[StringJoin["CASE(",ToString[ijklFlag],")"]];
	     WSI[StringJoin["CALL IntB",ToString[ijklType],"(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & \n", 
                            "                        BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & \n",
                            "                        ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & \n",
                            "                        OffSet%A  ,1              , & \n",
                            "                        OffSet%C-1,NBFA           , & \n",
                            "                        OffSet%B-1,NBFA*NBFC      , & \n",
                            "                        OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) \n"]];

Print["ijklType=",ijklType," i=",IntegralClass[Classes[[ic]]]," j=",IntegralClass[Classes[[jc]]]," k=",IntegralClass[Classes[[kc]]]," l=",IntegralClass[Classes[[lc]]]];

	   Subroutine=StringJoin["IntB",ToString[ijklType],".F90"];
	   OpenWrite[Subroutine];
	   Print[" Openned ",Subroutine];
	   WS[String_]:=WriteString[Subroutine,"   ",String,"\n"];
	   WriteString[Subroutine,CommentLine]; 
	   MakeList=Append[MakeList,StringJoin["IntB",ToString[ijklType],".o"]];
           RelsList=Append[RelsList,StringJoin["IntB",ToString[ijklType],".x"]]; 
           BraEll=imax+jmax;
           KetEll=kmax+lmax;
           PunchFront[Subroutine,ic,jc,kc,lc,ijklType]; 
           PunchVRRCalls[Subroutine,BraEll,KetEll];
           PunchVRRContract[Subroutine,ic,jc,kc,lc,BraEll,KetEll];
           PunchHRRKet[Subroutine,ic,jc,kc,lc]; 
	   PunchHRRBra[Subroutine,ic,jc,kc,lc]; 
           WS[StringJoin["END SUBROUTINE IntB",ToString[ijklType]]];
           Close[Subroutine];
           Print[" Closed ",Subroutine];
      ]; 


,{ic,1,LC}]
,{jc,1,LC}]
,{kc,1,LC}]
,{lc,1,LC}];

WSI["CASE DEFAULT"];
WSI[" ! STOP 'MISSED AN INTEGRAL' "];
WSI["END SELECT "];

Close[IncludeFile];

TAB=FromCharacterCode[9];

MakeList=Union[MakeList];
RelsList=Union[RelsList];
MakeString=""
RelsString=""
LM=Length[MakeList];
LR=Length[RelsList];
Do[MakeString=StringJoin[MakeString,MakeList[[i]],"\\\n"],{i,LM-1}]
MakeString=StringJoin[MakeString,MakeList[[LM]],"\n"];
Do[RelsString=StringJoin[RelsString,RelsList[[i]],"\\\n"],{i,LR-1}]
RelsString=StringJoin[RelsString,RelsList[[LR]],"\n"];
MakeList=MakeString;
RelsList=RelsString;

MakeList=StringJoin["IntObjs=",MakeList];
RelsList=StringJoin["IntRels=",RelsList]; 

(**************** Print out the Makefile ************************)

Makefile="Makefile"

OpenWrite[Makefile];

WriteString[Makefile,"include $(MONDO_COMPILER)\n"];
WriteString[Makefile,"include $(MONDO_HOME)/Includes/Suffixes\n"];
WriteString[Makefile,"include $(MONDO_HOME)/Includes/RemoveAll\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"CPPMISC =\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"EXTRA_INCLUDES=\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"TwoELibs=$(F90_LIBRARIES) -lTwoE $(LAPACK)\n"] 
WriteString[Makefile,"#\n"];
WriteString[Makefile,"MiscObjs=ShellPairStruct.o VScratchB.o GetAtomPair.o InitDbl.o \n"];
WriteString[Makefile,"#\n"];
Print[MakeList];
Print[RelsList];
WriteString[Makefile,MakeList];
WriteString[Makefile,"#\n"];
WriteString[Makefile,RelsList];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"all:    TwoE\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"clean:  CTwoE\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"purge:clean\n",TAB,"rm -f $(MONDO_LIB)/libTwoE.a\n",TAB,"rm -f $(REMOVESRCE)\n"];
WriteString[Makefile,"ppurge:purge\n",TAB,"rm -f IntB*.F90 HRR*.F90 VRR*.F90 \n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"TwoE:$(MiscObjs) $(IntObjs)\n",TAB,"$(AR) $(ARFLAGS) $(MONDO_LIB)/libTwoE.a $(?:.F90=.o)\n",TAB,"$(RANLIB) $(MONDO_LIB)/libTwoE.a\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"CTwoE:\n",TAB,"rm -f $(REMOVEMISC) $(REMOVEMODS)\n",TAB,"rm -f \#*\n",TAB,"rm -f *~\n",TAB,"ln -s /dev/null core\n",TAB,"ls -l\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"N4K: TwoE N4K.o ERIInclude.Inc\n",TAB,"$(LD) $(LDFLAGS) N4K N4K.o $(TwoELibs)\n"];
Close[Makefile];


