(*======================================--  BRA HRR --============================================*)

PunchBraHRRFront[Subroutine_,Subname_,ic_,jc_]:=Block[{WS,BKType,LenBra,LenKet},
 WriteString[Subroutine,StringJoin["   SUBROUTINE ",Subname,"(OA,OB,LDA,LDB,NINT,CDOffSet,HRR,HRRA,HRRB,GRADIENT) \n"]];

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
 WS["IMPLICIT NONE"];
 WS["INTEGER       :: OA,OB,LDA,LDB,NINT,CDOffSet,OffSet"];
 WS[StringJoin["REAL(DOUBLE)  :: HRR(*),HRRA(*),HRRB(*)"]];
 WS[StringJoin["REAL(DOUBLE)  :: GRADIENT(NINT,12)"]];
];

(* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  *)

HRRBra[Label_String,a_List,b_List]:=Module[{AB,CD,pa,pb,pc,pd,MaxEll,a1,b1,c1,d1,one,two,adex,cdex},
                              pa=Position[a,Max[b]][[1, 1]];
                              pb=Position[b,Max[b]][[1, 1]];
			      (* Exit condition 1 *)
                              If[ a[[pa]] < 0 || b[[pb]] < 0 ,Return[0];]; 
			      (* Exit condition 2 *)
                              If[ b[[pb]]==0 ,
                                 adex=LMNDex[a[[1]],a[[2]],a[[3]]];
                                 BarExp=ToExpression[StringJoin["HRR",Label,"[",ToString[adex],"]"]];
                                 Return[BarExp]
                                ];
			      (* Recursion *)
                              If[pb==1, one = {1, 0, 0}; AB = ABx; ];
                              If[pb==2, one = {0, 1, 0}; AB = ABy; ];
                              If[pb==3, one = {0, 0, 1}; AB = ABz; ];
                              a1 = a + one;
                              b1 = b - one;                
                              Return[ HRRBra[Label,a1, b1] + AB  HRRBra[Label,a, b1]]
                            ];

(* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  *)

PunchHRRBraCalls[FileName_,ic_,jc_,kc_,lc_]:=Module[{oList,IList,Kount,a,b,c,d,WS},

 imin = Classes[[ic, 1]]; imax = Classes[[ic, 2]];
 jmin = Classes[[jc, 1]]; jmax = Classes[[jc, 2]];
 kmin = Classes[[kc, 1]]; kmax = Classes[[kc, 2]];
 lmin = Classes[[lc, 1]]; lmax = Classes[[lc, 2]];
 LenI=LEnd[Classes[[ic,2]]]-LBegin[Classes[[ic,1]]];
 LenJ=LEnd[Classes[[jc,2]]]-LBegin[Classes[[jc,1]]];
 BEnd=LEnd[imax+jmax];
 HRRSubName=StringJoin["BraHRR",IClass,JClass];

 IClass=ToString[IntegralClass[Classes[[ic]]]];
 JClass=ToString[IntegralClass[Classes[[jc]]]];

 WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];

 WS[StringJoin["DO L=",ToString[LBegin[lmin]],",",ToString[LEnd[lmax]]]];
 Do[ 
    KS=ToString[k];
    WS[StringJoin["\n         !K = ",KS]];
    WS[StringJoin["   CDOffSet=(OC+",KS,"-",ToString[LBegin[kmin]],")*LDC+(OD+L-",ToString[LBegin[lmin]],")*LDD"]];

    HRRAddress=StringJoin["1,",KS,",",ToString[L]];

    WS[StringJoin["   ! Generating (",CType[IntegralClass[{imin,imax}]],"',",
                                      CType[IntegralClass[{jmin,jmax}]],"|",ToString[k],",L)  and ",
                                  "(",CType[IntegralClass[{imin,imax}]],",",
                                      CType[IntegralClass[{jmin,jmax}]],"'|",ToString[k],",L)"]];
    (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
    If[DoStress==0, 
      WS[StringJoin["   CALL ",HRRSubName,"ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(",HRRAddress,"),&\n",
      "                      HRRA(",HRRAddress,"),HRRB(",HRRAddress,"),GRADIENTS(1,1),FP(1),&\n",
      "                      STRESS(1,1))"]];
    ,
      WS[StringJoin["   CALL ",HRRSubName,"ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR(",HRRAddress,"),&\n                          ",
                                        "HRRA(",HRRAddress,"),HRRB(",HRRAddress,"),GRADIENTS(1,1))"]];
    ];
    (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
    LK=lx[k];
    MK=my[k];
    NK=nz[k];
    KPX=ToString[LMNDex[LK+1,MK,NK]];
    KPY=ToString[LMNDex[LK,MK+1,NK]];
    KPZ=ToString[LMNDex[LK,MK,NK+1]];
    KMX=ToString[LMNDex[LK-1,MK,NK]];
    KMY=ToString[LMNDex[LK,MK-1,NK]];
    KMZ=ToString[LMNDex[LK,MK,NK-1]];
    BS=ToString[BEnd];

    WS[StringJoin["   ! Generating (",CType[IntegralClass[{imin,imax}]],",",
                                      CType[IntegralClass[{jmin,jmax}]],"|",ToString[k],"_x,L)  and ",
                                  "(",CType[IntegralClass[{imin,imax}]],",",
                                      CType[IntegralClass[{jmin,jmax}]],"|",ToString[k],",L_x)"]];
    If[ LK==0 , 
      (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
      If[DoStress==0, 
        WS[StringJoin["   CALL ",HRRSubName,"cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,&\n",
        "                      HRRC(1,",KPX,",L),GRADIENTS(1,1),FP(1),STRESS(1,1))"]];
      ,
        WS[StringJoin["   CALL ",HRRSubName,"cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,HRRC(1,",KPX,",L),GRADIENTS(1,1))"]];
      ];
      (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
    ,
       WS[StringJoin["   HRRTmp(1:",BS,")=HRRC(1:",BS,",",KPX,",L)-",ToString[LK],"D0*HRR(1:",BS,",",KMX,",L)"]];
      (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
      If[DoStress==0, 
        WS[StringJoin["   CALL ",HRRSubName,"cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,&\n",
        "                      HRRTmp,GRADIENTS(1,1),FP(1),STRESS(1,1))"]]; 
      ,
        WS[StringJoin["   CALL ",HRRSubName,"cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,0,HRRTmp,GRADIENTS(1,1))"]]; 
      ];
      (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
    ];


    WS[StringJoin["   ! Generating (",CType[IntegralClass[{imin,imax}]],",",
                                      CType[IntegralClass[{jmin,jmax}]],"|",ToString[k],"_y,L)  and ",
                                  "(",CType[IntegralClass[{imin,imax}]],",",
                                      CType[IntegralClass[{jmin,jmax}]],"|",ToString[k],",L_y)"]];
    If[ MK==0 , 
      (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
      If[DoStress==0, 
        WS[StringJoin["   CALL ",HRRSubName,"cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,&\n",
        "                      HRRC(1,",KPY,",L),GRADIENTS(1,1),FP(1),STRESS(1,1))"]];
      ,
        WS[StringJoin["   CALL ",HRRSubName,"cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRC(1,",KPY,",L),GRADIENTS(1,1))"]];
      ];
      (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
    ,
        WS[StringJoin["   HRRTmp(1:",BS,")=HRRC(1:",BS,",",KPY,",L)-",ToString[MK],"D0*HRR(1:",BS,",",KMY,",L)"]];
      (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
      If[DoStress==0, 
        WS[StringJoin["   CALL ",HRRSubName,"cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,&\n",
        "                      HRRTmp,GRADIENTS(1,1),FP(1),STRESS(1,1))"]]; 
      ,
        WS[StringJoin["   CALL ",HRRSubName,"cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,1,HRRTmp,GRADIENTS(1,1))"]]; 
      ];
      (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
    ];


    WS[StringJoin["   ! Generating (",CType[IntegralClass[{imin,imax}]],",",
                                      CType[IntegralClass[{jmin,jmax}]],"|",ToString[k],"_z,L)  and ",
                                  "(",CType[IntegralClass[{imin,imax}]],",",
                                      CType[IntegralClass[{jmin,jmax}]],"|",ToString[k],",L_z)"]];
    If[ NK==0 , 
      (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
      If[DoStress==0, 
        WS[StringJoin["   CALL ",HRRSubName,"cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,&\n", 
        "                      HRRC(1,",KPZ,",L),GRADIENTS(1,1),FP(1),STRESS(1,1))"]];
      ,
        WS[StringJoin["   CALL ",HRRSubName,"cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRC(1,",KPZ,",L),GRADIENTS(1,1))"]];
      ];
      (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
    ,
        WS[StringJoin["   HRRTmp(1:",BS,")=HRRC(1:",BS,",",KPZ,",L)-",ToString[NK],"D0*HRR(1:",BS,",",KMZ,",L)"]];
      (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
      If[DoStress==0, 
        WS[StringJoin["   CALL ",HRRSubName,"cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,&\n", 
        "                      HRRTmp,GRADIENTS(1,1),FP(1),STRESS(1,1))"]]; 
      ,
        WS[StringJoin["   CALL ",HRRSubName,"cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,2,HRRTmp,GRADIENTS(1,1))"]]; 
      ];
      (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
    ];


 ,{k,LBegin[kmin],LEnd[kmax]}];
 WriteString[FileName,"      ENDDO \n "];

 (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
 If[DoStress==0, 
   WriteString[FileName,StringJoin["     ! Stress: Generating (",CType[IntegralClass[{imin,imax}]],",",
                                                       CType[IntegralClass[{jmin,jmax}]],"|",
                                                       CType[IntegralClass[{kmin,kmax}]],",",
                                                       CType[IntegralClass[{lmin,lmax}]],")^(1) \n"]];
   WriteString[FileName,           "      DO J=1,3\n"];
   WriteString[FileName,           "      DO I=1,3\n"];
(*
   WriteString[FileName,           "      IF(PBC%AutoW%I(I).EQ.1.AND.PBC%AutoW%I(J).EQ.1) THEN \n"];
*)
   WriteString[FileName,           "      IJ=3*(J-1)+I\n"];
   WriteString[FileName,StringJoin["        DO L=",ToString[LBegin[lmin]],",",ToString[LEnd[lmax]],"\n"]];
   WriteString[FileName,StringJoin["          DO K=",ToString[LBegin[kmin]],",",ToString[LEnd[kmax]],"\n"]];
   WriteString[FileName,StringJoin["            CDOffSet=(OC+K-",ToString[LBegin[kmin]],")*LDC+(OD+L-",ToString[LBegin[lmin]],")*LDD \n"]];
   WriteString[FileName,StringJoin["            CALL ",HRRSubName,"(OA,OB,LDA,LDB,CDOffSet,HRRS(1,",ToString[K],",",ToString[L],",IJ),STRESS(1,IJ))\n"]];
   WriteString[FileName,           "          ENDDO \n"];
   WriteString[FileName,           "        ENDDO \n"];
(*
   WriteString[FileName,           "      ENDIF \n"];
*)
   WriteString[FileName,           "      ENDDO \n"];
   WriteString[FileName,           "      ENDDO \n"];
 ];
 (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)

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

 CartS={"x","y","z"};

 (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)

 IList={};
 oList={};
 Kount = 0;
 Do[Do[Do[Do[
             a = {lx[i], my[i], nz[i]};
             b = {lx[j], my[j], nz[j]};

             OffSetString=StringJoin["(OA+",ToString[i-LBegin[imin]],")*LDA+(OB+",ToString[j-LBegin[jmin]],")*LDB+CDOffSet"];
             Kount = Kount + 1;
             IList=Append[IList,ZP];
             oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["OffSet=",OffSetString," !"]];

             Do[
                plus={0,0,0}; plus[[cart]]=+1;
                mnus={0,0,0}; mnus[[cart]]=-1;


                Kount = Kount + 1;
                IList=Append[IList,ToExpression[StringJoin["BB",ToString[i],"DV",CartS[[cart]],"XX",ToString[j],"KK"]]];
                oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->"!"];
 
                (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
                If[DoStress==0, 

                  hrr=Normy[a]*Normy[b]*(HRRBra["A",a+plus,b]- a[[cart]]  HRRBra["",a+mnus,b]);
                  Kount = Kount + 1;
                  IList=Append[IList,Horner[hrr]];
                  HRRAddress=StringJoin[ToString[i],",",ToString[j]];
                  oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["DUM"]];


                  (* Set Gradient *)

                  Kount = Kount + 1;
                  tmpp=DUM+GRADIENT[OffSet,cart+GOA-1];
                  IList=Append[IList,tmpp];
                  oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["GRADIENT(OffSet,",ToString[cart+GOA-1],")"]];

                  (* Set Stress *)

		  Kount = Kount + 1;
		  (*tmpp=DUM*FP[1]+STRESS[OffSet,1+3*(cart-1)];*)
                  tmpp=DUM*FP[1]+STRESS[OffSet,0+cart];
                  IList=Append[IList,tmpp];
                  (*oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,",ToString[1+3*(cart-1)],")"]];*)
		  oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,",ToString[0+cart],")"]];

                  Kount = Kount + 1;
                  (*tmpp=DUM*FP[2]+STRESS[OffSet,2+3*(cart-1)];*)
                  tmpp=DUM*FP[2]+STRESS[OffSet,3+cart];
                  IList=Append[IList,tmpp];
                  (*oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,",ToString[2+3*(cart-1)],")"]];*) 
                  oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,",ToString[3+cart],")"]];

                  Kount = Kount + 1;

                  (*tmpp=DUM*FP[3]+STRESS[OffSet,3+3*(cart-1)];*)
                  tmpp=DUM*FP[3]+STRESS[OffSet,6+cart];
                  IList=Append[IList,tmpp];
                  (*oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,",ToString[3+3*(cart-1)],")"]];*)  
                  oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,",ToString[6+cart],")"]];  

                ,
                  hrr=Normy[a]*Normy[b]*(HRRBra["A",a+plus,b]- a[[cart]]  HRRBra["",a+mnus,b])+GRADIENT[OffSet,cart+GOA-1];
                  Kount = Kount + 1;
                  IList=Append[IList,Horner[hrr]];
                  HRRAddress=StringJoin[ToString[i],",",ToString[j]];
                  oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["GRADIENT(OffSet,",ToString[cart+GOA-1],")"]];
                ];
                (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)

                Kount = Kount + 1;
                IList=Append[IList,ToExpression[StringJoin["BB",ToString[i],"XX",ToString[j],"DV",CartS[[cart]],"KK"]]];
                oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->"!"];
                (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
                If[DoStress==0, 

                  hrr=Normy[a]*Normy[b]*(HRRBra["B",a,b+plus]- b[[cart]]  HRRBra["",a,b+mnus]);
                  Kount = Kount + 1;
                  IList=Append[IList,Horner[hrr]];
                  HRRAddress=StringJoin[ToString[i],",",ToString[j]];
                  oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["DUM"]];

                  (* Set Gradient *)

                  Kount = Kount + 1;
                  tmpp=DUM+GRADIENT[OffSet,cart+GOB-1];
                  IList=Append[IList,tmpp];
                  oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["GRADIENT(OffSet,",ToString[cart+GOB-1],")"]];

                  (* Set Stress *)

                  Kount = Kount + 1;
                  (*tmpp=DUM*FP[7]+STRESS[OffSet,1+3*(cart-1)];*)
                  tmpp=DUM*FP[7]+STRESS[OffSet,0+cart];
                  IList=Append[IList,tmpp];
                  (*oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,",ToString[1+3*(cart-1)],")"]];*)  
                  oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,",ToString[0+cart],")"]];  

                  Kount = Kount + 1;
                  (*tmpp=DUM*FP[8]+STRESS[OffSet,2+3*(cart-1)];*)
                  tmpp=DUM*FP[8]+STRESS[OffSet,3+cart];
                  IList=Append[IList,tmpp];
                  (*oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,",ToString[2+3*(cart-1)],")"]];*)    
                  oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,",ToString[3+cart],")"]];    

                  Kount = Kount + 1;
                  (*tmpp=DUM*FP[9]+STRESS[OffSet,3+3*(cart-1)];*)
                  tmpp=DUM*FP[9]+STRESS[OffSet,6+cart];
                  IList=Append[IList,tmpp];
                  (*oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,",ToString[3+3*(cart-1)],")"]];*)   
                  oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,",ToString[6+cart],")"]];   

                ,
                  hrr=Normy[a]*Normy[b]*(HRRBra["B",a,b+plus]- b[[cart]]  HRRBra["",a,b+mnus])+GRADIENT[OffSet,cart+GOB-1];
                  Kount = Kount + 1;
                  IList=Append[IList,Horner[hrr]];
                  HRRAddress=StringJoin[ToString[i],",",ToString[j]];
                  oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["GRADIENT(OffSet,",ToString[cart+GOB-1],")"]];
                ];
                (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)

             ,{cart,1,3}]
 ,{i,LBegin[il],LEnd[il]}]
 ,{j,LBegin[jl],LEnd[jl]}]
 ,{il,imin,imax}]
 ,{jl,jmin,jmax}];

  oList=Append[oList,{" "->"","DV"->"_","ZP"->"","XX"->",","BB"->"(","KK"->"|","+"->"+&\n                         ","CC"->")"}];
  oList=Flatten[oList];

  MakeList=Append[MakeList,StringJoin[HRRSubName,".o"]];
  RelsList=Append[RelsList,StringJoin[HRRSubName,".x"]];
  HRRSubroutine=StringJoin[HRRSubName,".F90"];

  OpenWrite[HRRSubroutine];
  WSS[String_]:=WriteString[HRRSubroutine,"    ",String,"\n"];

  (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
  If[DoStress==0,
    WSS[StringJoin["SUBROUTINE ",HRRSubName, \
        "ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT,FP,STRESS)"]];
  ,
    WSS[StringJoin["SUBROUTINE ",HRRSubName,"ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT)"]];
  ];
  (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)

  WSS["  USE DerivedTypes"];
  WSS["  USE VScratchB"];
  WSS["  USE GlobalScalars"];
  WSS["  IMPLICIT NONE"];
  WSS["  INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,OffSet"];
  WSS[StringJoin["  REAL(DOUBLE)  :: HRR(*),HRRA(*),HRRB(*)"]];
  WSS[StringJoin["  REAL(DOUBLE)  :: GRADIENT(NINT,12)"]];

  (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
  If[DoStress==0,
    WSS[StringJoin["  REAL(DOUBLE)  :: STRESS(NINT,9),FP(9),DUM"]];
  ];
  (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)

  Write[HRRSubroutine,FortranAssign[o,IList,AssignReplace->oList]];
  WSS[StringJoin["END SUBROUTINE ",HRRSubName,"ab"]];			       

 (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)

 IList={};
 oList={};
 Kount = 0;
 Do[Do[Do[Do[
             a = {lx[i], my[i], nz[i]};
             b = {lx[j], my[j], nz[j]};

             OffSetString=StringJoin["(OA+",ToString[i-LBegin[imin]],")*LDA+(OB+",ToString[j-LBegin[jmin]],")*LDB+CDOffSet"];
             Kount = Kount + 1;
             IList=Append[IList,ZP];
             oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["OffSet=",OffSetString," !"]];

             Kount = Kount + 1;
             IList=Append[IList,ToExpression[StringJoin["BB",ToString[i],"XX",ToString[j],"KK"]]];
             oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->"!"];
 

  (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
             If[DoStress==0,

               hrr=Normy[a]*Normy[b]*HRRBra["",a,b];
               Kount = Kount + 1;
               IList=Append[IList,Horner[hrr]];
               HRRAddress=StringJoin[ToString[i],",",ToString[j]];
               oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["DUM"]];

               (* Set Gradient *)
               Kount = Kount + 1;

               tmpp=DUM+GRADIENT[OffSet,Cart+GOC];
               IList=Append[IList,tmpp];
               oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["GRADIENT(OffSet,",ToString[Cart+GOC],")"]];

               (* Set Stress *)

               Kount = Kount + 1;
               (*tmpp=DUM*FP[4]+STRESS[OffSet,1+3*Cart];*)
               tmpp=DUM*FP[4]+STRESS[OffSet,1+Cart];
               IList=Append[IList,tmpp];
               (*oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,1+3*Cart)"]];*)
               oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,1+Cart)"]];

               Kount = Kount + 1;
               (*tmpp=DUM*FP[5]+STRESS[OffSet,2+3*Cart];*)
               tmpp=DUM*FP[5]+STRESS[OffSet,4+Cart];
               IList=Append[IList,tmpp];
               (*oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,2+3*Cart)"]];*)
               oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,4+Cart)"]];

               Kount = Kount + 1;
               (*tmpp=DUM*FP[6]+STRESS[OffSet,3+3*Cart];*)
               tmpp=DUM*FP[6]+STRESS[OffSet,7+Cart];
               IList=Append[IList,tmpp];
               (*oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,3+3*Cart)"]];*)
               oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["STRESS(OffSet,7+Cart)"]];

             ,
               hrr=Normy[a]*Normy[b]*HRRBra["",a,b]+GRADIENT[OffSet,cart+GOC];
               Kount = Kount + 1;
               IList=Append[IList,Horner[hrr]];
               HRRAddress=StringJoin[ToString[i],",",ToString[j]];
               oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["GRADIENT(OffSet,",ToString[Cart+GOC],")"]]; 
             ];
  (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)

             (* Use translational symmetry *)
             hrr=QRS;
             Kount = Kount + 1;
             IList=Append[IList,Horner[hrr]];
             HRRAddress=StringJoin[ToString[i],",",ToString[j]];
             oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["GRADIENT(OffSet,",ToString[Cart+GOD],")"]];


 ,{i,LBegin[il],LEnd[il]}]
 ,{j,LBegin[jl],LEnd[jl]}]
 ,{il,imin,imax}]
 ,{jl,jmin,jmax}];

 spaces="                                ";
 oList=Append[oList,{" "->"","DV"->"_","ZP"->"","XX"->",","BB"->"(","KK"->"|","+"->"+&\n                                ",
     "QRS"->StringJoin["-GRADIENT(OffSet,Cart+GOA)&\n",spaces,
                       "-GRADIENT(OffSet,Cart+GOB)&\n",spaces,"-GRADIENT(OffSet,Cart+GOC)"]}];

 oList=Flatten[oList];

  (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
  If[DoStress==0,
    WSS[StringJoin["SUBROUTINE ",HRRSubName, \
        "cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT,FP,STRESS)"]];
  ,
    WSS[StringJoin["SUBROUTINE ",HRRSubName,"cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT)"]];
  ];
  (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)

  WSS["  USE DerivedTypes"];
  WSS["  USE VScratchB"];
  WSS["  USE GlobalScalars"];
  WSS["  IMPLICIT NONE"];
  WSS["  INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,Cart,CDOffSet,OffSet"];
  WSS[StringJoin["  REAL(DOUBLE)  :: HRR(*)"]];
  WSS[StringJoin["  REAL(DOUBLE)  :: GRADIENT(NINT,12)"]];

  (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
  If[DoStress==0,
    WSS[StringJoin["  REAL(DOUBLE)  :: STRESS(NINT,9),FP(9),DUM"]];
  ];
  (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)

  Write[HRRSubroutine,FortranAssign[o,IList,AssignReplace->oList]];
  WSS[StringJoin["END SUBROUTINE ",HRRSubName,"cd"]];			       
 Close[HRRSubroutine];


,{ic,LC}],{jc,LC}];

];
