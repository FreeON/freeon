
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


   (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
   (*  MIC-VRR MIC-VRR MIC-VRR MIC-VRR MIC-VRR MIC-VRR MIC-VRR MIC-VRR MIC-VRR MIC-VRR MIC-VRR MIC-VRR   *)
   (* Do we want Stress? *)
   If[DoStress==0,

   CLenBra=ToString[LEnd[BraEll]];
   CLenKet=ToString[LEnd[KetEll]];
   CLenEhm=ToString[BraEll+KetEll];

   CLenBra1=ToString[LEnd[BraEll+1]];
   CLenKet1=ToString[LEnd[KetEll+1]];
   CLenEhm1=ToString[BraEll+KetEll+1];

   If[BraEll+KetEll>0,WS["!MAY BE BETTER TO PUT WHAT FOLLOWS IN A DO LOOP!"];];

   (* Start loop over ixyz=x,y,z *)
   Do[

      (* Punch IF-loop *)
      If[ixyz==1,WS["IF(PBC%AutoW%I(1).EQ.1) THEN"];];
      If[ixyz==2,WS["IF(PBC%AutoW%I(2).EQ.1) THEN"];];
      If[ixyz==3,WS["IF(PBC%AutoW%I(3).EQ.1) THEN"];];

      (* Punch Initialization *)
      Do[
         GSt=ToString[L];
         If[L>0,
            If[ixyz==1,WS[StringJoin["VRRS(1,1,",ToString[L-1],",1)=PQx*VRR(1,1,",GSt,")"]];];
            If[ixyz==2,WS[StringJoin["VRRS(1,1,",ToString[L-1],",2)=PQy*VRR(1,1,",GSt,")"]];];
            If[ixyz==3,WS[StringJoin["VRRS(1,1,",ToString[L-1],",3)=PQz*VRR(1,1,",GSt,")"]];];
         ];
     ,{L,0,BraEll+KetEll+1}];

      (* Punch MIC-VRR *)
    Do[Do[
      If[iell+kell>0,
        MVRRSubName=StringJoin["MVRR",CType[IntegralClass[{iell,iell}]],"0",CType[IntegralClass[{kell,kell}]],"0"];
        Do[
      
           ms=ToString[m];
           m1=ToString[m-1];
           p1=ToString[m+1];

           WS[StringJoin["! MIC-VRR: Generating [",CType[IntegralClass[{iell,iell}]],"0|", \
                                                   CType[IntegralClass[{kell,kell}]],"0]^(",m1,")"]];
      
           If[iell==0&&kell==1,
             If[ixyz==1,
             (* i=x *)
             WS[StringJoin["VRRS(1,2,",m1,",1)=QCx*VRRS(1,1,",m1,",1)+WQx*VRRS(1,1,",ms,",1)-r1x2E*VRR(1,1,",ms,")"]];
             WS[StringJoin["VRRS(1,3,",m1,",1)=QCy*VRRS(1,1,",m1,",1)+WQy*VRRS(1,1,",ms,",1) "]];
             WS[StringJoin["VRRS(1,4,",m1,",1)=QCz*VRRS(1,1,",m1,",1)+WQz*VRRS(1,1,",ms,",1) "]];
             ];
             If[ixyz==2,
             (* i=y *)
             WS[StringJoin["VRRS(1,2,",m1,",2)=QCx*VRRS(1,1,",m1,",2)+WQx*VRRS(1,1,",ms,",2)"]];
             WS[StringJoin["VRRS(1,3,",m1,",2)=QCy*VRRS(1,1,",m1,",2)+WQy*VRRS(1,1,",ms,",2)-r1x2E*VRR(1,1,",ms,")"]];
             WS[StringJoin["VRRS(1,4,",m1,",2)=QCz*VRRS(1,1,",m1,",2)+WQz*VRRS(1,1,",ms,",2)"]];
             ];
             If[ixyz==3,
             (* i=z *)
             WS[StringJoin["VRRS(1,2,",m1,",3)=QCx*VRRS(1,1,",m1,",3)+WQx*VRRS(1,1,",ms,",3)"]];
             WS[StringJoin["VRRS(1,3,",m1,",3)=QCy*VRRS(1,1,",m1,",3)+WQy*VRRS(1,1,",ms,",3)"]];
             WS[StringJoin["VRRS(1,4,",m1,",3)=QCz*VRRS(1,1,",m1,",3)+WQz*VRRS(1,1,",ms,",3)-r1x2E*VRR(1,1,",ms,")"]];
             ];
           ];
        
           If[iell==1&&kell==0,
             If[ixyz==1,
             (* i=x *)
             WS[StringJoin["VRRS(2,1,",m1,",1)=PAx*VRRS(1,1,",m1,",1)+WPx*VRRS(1,1,",ms,",1)+r1x2Z*VRR(1,1,",ms,")"]];
             WS[StringJoin["VRRS(3,1,",m1,",1)=PAy*VRRS(1,1,",m1,",1)+WPy*VRRS(1,1,",ms,",1) "]];
             WS[StringJoin["VRRS(4,1,",m1,",1)=PAz*VRRS(1,1,",m1,",1)+WPz*VRRS(1,1,",ms,",1) "]];
             ];
             If[ixyz==2,
             (* i=y *)
             WS[StringJoin["VRRS(2,1,",m1,",2)=PAx*VRRS(1,1,",m1,",2)+WPx*VRRS(1,1,",ms,",2)"]];
             WS[StringJoin["VRRS(3,1,",m1,",2)=PAy*VRRS(1,1,",m1,",2)+WPy*VRRS(1,1,",ms,",2)+r1x2Z*VRR(1,1,",ms,")"]];
             WS[StringJoin["VRRS(4,1,",m1,",2)=PAz*VRRS(1,1,",m1,",2)+WPz*VRRS(1,1,",ms,",2)"]];
             ];
             If[ixyz==3,
             (* i=z *)
             WS[StringJoin["VRRS(2,1,",m1,",3)=PAx*VRRS(1,1,",m1,",3)+WPx*VRRS(1,1,",ms,",3)"]];
             WS[StringJoin["VRRS(3,1,",m1,",3)=PAy*VRRS(1,1,",m1,",3)+WPy*VRRS(1,1,",ms,",3)"]];
             WS[StringJoin["VRRS(4,1,",m1,",3)=PAz*VRRS(1,1,",m1,",3)+WPz*VRRS(1,1,",ms,",3)+r1x2Z*VRR(1,1,",ms,")"]];
             ];
           ];
        
           If[iell==1&&kell==1,
             If[ixyz==1,
             (* i=x *)
             WS[StringJoin["VRRS(2,2,",m1,",1)=QCx*VRRS(2,1,",m1,",1)+WQx*VRRS(2,1,",ms,",1)+HfxZpE*VRRS(1,1,",ms,",1)",\
                           "-r1x2E *VRR(2,1,",ms,")"]];
	     WS[StringJoin["VRRS(2,3,",m1,",1)=QCy*VRRS(2,1,",m1,",1)+WQy*VRRS(2,1,",ms,",1)"]];
             WS[StringJoin["VRRS(2,4,",m1,",1)=QCz*VRRS(2,1,",m1,",1)+WQz*VRRS(2,1,",ms,",1)"]];
             WS[StringJoin["VRRS(3,2,",m1,",1)=QCx*VRRS(3,1,",m1,",1)+WQx*VRRS(3,1,",ms,",1)-r1x2E *VRR(3,1,",ms,")"]];
             WS[StringJoin["VRRS(3,3,",m1,",1)=QCy*VRRS(3,1,",m1,",1)+WQy*VRRS(3,1,",ms,",1)+HfxZpE*VRRS(1,1,",ms,",1)"]];
             WS[StringJoin["VRRS(3,4,",m1,",1)=QCz*VRRS(3,1,",m1,",1)+WQz*VRRS(3,1,",ms,",1)"]];
             WS[StringJoin["VRRS(4,2,",m1,",1)=QCx*VRRS(4,1,",m1,",1)+WQx*VRRS(4,1,",ms,",1)-r1x2E *VRR(4,1,",ms,")"]];
             WS[StringJoin["VRRS(4,3,",m1,",1)=QCy*VRRS(4,1,",m1,",1)+WQy*VRRS(4,1,",ms,",1)"]];
             WS[StringJoin["VRRS(4,4,",m1,",1)=QCz*VRRS(4,1,",m1,",1)+WQz*VRRS(4,1,",ms,",1)+HfxZpE*VRRS(1,1,",ms,",1)"]];
             ];
             If[ixyz==2,
             (* i=y *)
             WS[StringJoin["VRRS(2,2,",m1,",2)=QCx*VRRS(2,1,",m1,",2)+WQx*VRRS(2,1,",ms,",2)+HfxZpE*VRRS(1,1,",ms,",2)"]];
	     WS[StringJoin["VRRS(2,3,",m1,",2)=QCy*VRRS(2,1,",m1,",2)+WQy*VRRS(2,1,",ms,",2)-r1x2E *VRR(2,1,",ms,")"]];
             WS[StringJoin["VRRS(2,4,",m1,",2)=QCz*VRRS(2,1,",m1,",2)+WQz*VRRS(2,1,",ms,",2)"]];
             WS[StringJoin["VRRS(3,2,",m1,",2)=QCx*VRRS(3,1,",m1,",2)+WQx*VRRS(3,1,",ms,",2)"]];
             WS[StringJoin["VRRS(3,3,",m1,",2)=QCy*VRRS(3,1,",m1,",2)+WQy*VRRS(3,1,",ms,",2)+HfxZpE*VRRS(1,1,",ms,",2)",\
                           "-r1x2E *VRR(3,1,",ms,")"]];
             WS[StringJoin["VRRS(3,4,",m1,",2)=QCz*VRRS(3,1,",m1,",2)+WQz*VRRS(3,1,",ms,",2)"]];
             WS[StringJoin["VRRS(4,2,",m1,",2)=QCx*VRRS(4,1,",m1,",2)+WQx*VRRS(4,1,",ms,",2)"]];
             WS[StringJoin["VRRS(4,3,",m1,",2)=QCy*VRRS(4,1,",m1,",2)+WQy*VRRS(4,1,",ms,",2)-r1x2E *VRR(4,1,",ms,")"]];
             WS[StringJoin["VRRS(4,4,",m1,",2)=QCz*VRRS(4,1,",m1,",2)+WQz*VRRS(4,1,",ms,",2)+HfxZpE*VRRS(1,1,",ms,",2)"]];
             ];
             If[ixyz==3,
             (* i=z *)
             WS[StringJoin["VRRS(2,2,",m1,",3)=QCx*VRRS(2,1,",m1,",3)+WQx*VRRS(2,1,",ms,",3)+HfxZpE*VRRS(1,1,",ms,",3)"]];
	     WS[StringJoin["VRRS(2,3,",m1,",3)=QCy*VRRS(2,1,",m1,",3)+WQy*VRRS(2,1,",ms,",3)"]];
             WS[StringJoin["VRRS(2,4,",m1,",3)=QCz*VRRS(2,1,",m1,",3)+WQz*VRRS(2,1,",ms,",3)-r1x2E *VRR(2,1,",ms,")"]];
             WS[StringJoin["VRRS(3,2,",m1,",3)=QCx*VRRS(3,1,",m1,",3)+WQx*VRRS(3,1,",ms,",3)"]];
             WS[StringJoin["VRRS(3,3,",m1,",3)=QCy*VRRS(3,1,",m1,",3)+WQy*VRRS(3,1,",ms,",3)+HfxZpE*VRRS(1,1,",ms,",3)"]];
             WS[StringJoin["VRRS(3,4,",m1,",3)=QCz*VRRS(3,1,",m1,",3)+WQz*VRRS(3,1,",ms,",3)-r1x2E *VRR(3,1,",ms,")"]];
             WS[StringJoin["VRRS(4,2,",m1,",3)=QCx*VRRS(4,1,",m1,",3)+WQx*VRRS(4,1,",ms,",3)"]];
             WS[StringJoin["VRRS(4,3,",m1,",3)=QCy*VRRS(4,1,",m1,",3)+WQy*VRRS(4,1,",ms,",3)"]];
             WS[StringJoin["VRRS(4,4,",m1,",3)=QCz*VRRS(4,1,",m1,",3)+WQz*VRRS(4,1,",ms,",3)+HfxZpE*VRRS(1,1,",ms,",3)",\
                           "-r1x2E *VRR(4,1,",ms,")"]];
             ];
           ];
        
           If[iell==2&&kell==0,
             If[ixyz==1,
             (* i=x *)
             WS[StringJoin["VRRS( 5,1,",m1,",1)=PAx*VRRS(2,1,",m1,",1)",\
                                              "+r1x2Z*(VRRS(1,1,",m1,",1)-ExZpE*VRRS(1,1,",ms,",1))",\
                                              "+WPx*VRRS(2,1,",ms,",1)",\
                                              "+r1x2Z*VRR(2,1,",ms,")"]];
             WS[StringJoin["VRRS( 6,1,",m1,",1)=PAx*VRRS(3,1,",m1,",1)",\
                                              "+WPx*VRRS(3,1,",ms,",1)",\
                                              "+r1x2Z*VRR(3,1,",ms,")"]];
             WS[StringJoin["VRRS( 7,1,",m1,",1)=PAy*VRRS(3,1,",m1,",1)",\
                                              "+r1x2Z*(VRRS(1,1,",m1,",1)-ExZpE*VRRS(1,1,",ms,",1))",\
                                              "+WPy*VRRS(3,1,",ms,",1)"]];
             WS[StringJoin["VRRS( 8,1,",m1,",1)=PAx*VRRS(4,1,",m1,",1)",\
                                              "+WPx*VRRS(4,1,",ms,",1)",\
                                              "+r1x2Z*VRR(4,1,",ms,")"]];
             WS[StringJoin["VRRS( 9,1,",m1,",1)=PAy*VRRS(4,1,",m1,",1)",\
                                              "+WPy*VRRS(4,1,",ms,",1)"]];
             WS[StringJoin["VRRS(10,1,",m1,",1)=PAz*VRRS(4,1,",m1,",1)",\
                                              "+r1x2Z*(VRRS(1,1,",m1,",1)-ExZpE*VRRS(1,1,",ms,",1))",\
                                              "+WPz*VRRS(4,1,",ms,",1)"]];
             ];
             If[ixyz==2,
             (* i=y *)
             WS[StringJoin["VRRS( 5,1,",m1,",2)=PAx*VRRS(2,1,",m1,",2)",\
                                              "+r1x2Z*(VRRS(1,1,",m1,",2)-ExZpE*VRRS(1,1,",ms,",2))",\
                                              "+WPx*VRRS(2,1,",ms,",2)"]];
             WS[StringJoin["VRRS( 6,1,",m1,",2)=PAx*VRRS(3,1,",m1,",2)",\
                                              "+WPx*VRRS(3,1,",ms,",2)"]];
             WS[StringJoin["VRRS( 7,1,",m1,",2)=PAy*VRRS(3,1,",m1,",2)",\
                                              "+r1x2Z*(VRRS(1,1,",m1,",2)-ExZpE*VRRS(1,1,",ms,",2))",\
                                              "+WPy*VRRS(3,1,",ms,",2)",\
                                              "+r1x2Z*VRR(3,1,",ms,")"]];
             WS[StringJoin["VRRS( 8,1,",m1,",2)=PAx*VRRS(4,1,",m1,",2)",\
                                              "+WPx*VRRS(4,1,",ms,",2)"]];
             WS[StringJoin["VRRS( 9,1,",m1,",2)=PAy*VRRS(4,1,",m1,",2)",\
                                              "+WPy*VRRS(4,1,",ms,",2)",\
                                              "+r1x2Z*VRR(4,1,",ms,")"]];
             WS[StringJoin["VRRS(10,1,",m1,",2)=PAz*VRRS(4,1,",m1,",2)",\
                                              "+r1x2Z*(VRRS(1,1,",m1,",2)-ExZpE*VRRS(1,1,",ms,",2))",\
                                              "+WPz*VRRS(4,1,",ms,",2)"]];
             ];
             If[ixyz==3,
             (* i=z *)
             WS[StringJoin["VRRS( 5,1,",m1,",3)=PAx*VRRS(2,1,",m1,",3)",\
                                              "+r1x2Z*(VRRS(1,1,",m1,",3)-ExZpE*VRRS(1,1,",ms,",3))",\
                                              "+WPx*VRRS(2,1,",ms,",3)"]];
             WS[StringJoin["VRRS( 6,1,",m1,",3)=PAx*VRRS(3,1,",m1,",3)",\
                                              "+WPx*VRRS(3,1,",ms,",3)"]];
             WS[StringJoin["VRRS( 7,1,",m1,",3)=PAy*VRRS(3,1,",m1,",3)",\
                                              "+r1x2Z*(VRRS(1,1,",m1,",3)-ExZpE*VRRS(1,1,",ms,",3))",\
                                              "+WPy*VRRS(3,1,",ms,",3)"]];
             WS[StringJoin["VRRS( 8,1,",m1,",3)=PAx*VRRS(4,1,",m1,",3)",\
                                              "+WPx*VRRS(4,1,",ms,",3)"]];
             WS[StringJoin["VRRS( 9,1,",m1,",3)=PAy*VRRS(4,1,",m1,",3)",\
                                              "+WPy*VRRS(4,1,",ms,",3)"]];
             WS[StringJoin["VRRS(10,1,",m1,",3)=PAz*VRRS(4,1,",m1,",3)",\
                                              "+r1x2Z*(VRRS(1,1,",m1,",3)-ExZpE*VRRS(1,1,",ms,",3))",\
                                              "+WPz*VRRS(4,1,",ms,",3)",\
                                              "+r1x2Z*VRR(4,1,",ms,")"]];
             ];
           ];
                                                                                                                            
           If[iell==0&&kell==2,
             If[ixyz==1,
             (* i=x *)
             WS[StringJoin["VRRS(1, 5,",m1,",1)=QCx*VRRS(1,2,",m1,",1)",\
                                              "+r1x2E*(VRRS(1,1,",m1,",1)-ZxZpE*VRRS(1,1,",ms,",1))",\
                                              "+WQx*VRRS(1,2,",ms,",1)",\
                                              "-r1x2E*VRR(1,2,",ms,")"]];
             WS[StringJoin["VRRS(1, 6,",m1,",1)=QCx*VRRS(1,3,",m1,",1)",\
                                              "+WQx*VRRS(1,3,",ms,",1)",\
                                              "-r1x2E*VRR(1,3,",ms,")"]];
             WS[StringJoin["VRRS(1, 7,",m1,",1)=QCy*VRRS(1,3,",m1,",1)",\
                                              "+r1x2E*(VRRS(1,1,",m1,",1)-ZxZpE*VRRS(1,1,",ms,",1))",\
                                              "+WQy*VRRS(1,3,",ms,",1)"]];
             WS[StringJoin["VRRS(1, 8,",m1,",1)=QCx*VRRS(1,4,",m1,",1)",\
                                              "+WQx*VRRS(1,4,",ms,",1)",\
                                              "-r1x2E*VRR(1,4,",ms,")"]];
             WS[StringJoin["VRRS(1, 9,",m1,",1)=QCy*VRRS(1,4,",m1,",1)",\
                                              "+WQy*VRRS(1,4,",ms,",1)"]];
             WS[StringJoin["VRRS(1,10,",m1,",1)=QCz*VRRS(1,4,",m1,",1)",\
                                              "+r1x2E*(VRRS(1,1,",m1,",1)-ZxZpE*VRRS(1,1,",ms,",1))",\
                                              "+WQz*VRRS(1,4,",ms,",1)"]];
             ];
             If[ixyz==2,
             (* i=y *)
             WS[StringJoin["VRRS(1, 5,",m1,",2)=QCx*VRRS(1,2,",m1,",2)",\
                                              "+r1x2E*(VRRS(1,1,",m1,",2)-ZxZpE*VRRS(1,1,",ms,",2))",\
                                              "+WQx*VRRS(1,2,",ms,",2)"]];
             WS[StringJoin["VRRS(1, 6,",m1,",2)=QCx*VRRS(1,3,",m1,",2)",\
                                              "+WQx*VRRS(1,3,",ms,",2)"]];
             WS[StringJoin["VRRS(1, 7,",m1,",2)=QCy*VRRS(1,3,",m1,",2)",\
                                              "+r1x2E*(VRRS(1,1,",m1,",2)-ZxZpE*VRRS(1,1,",ms,",2))",\
                                              "+WQy*VRRS(1,3,",ms,",2)",\
                                              "-r1x2E*VRR(1,3,",ms,")"]];
             WS[StringJoin["VRRS(1, 8,",m1,",2)=QCx*VRRS(1,4,",m1,",2)",\
                                              "+WQx*VRRS(1,4,",ms,",2)"]];
             WS[StringJoin["VRRS(1, 9,",m1,",2)=QCy*VRRS(1,4,",m1,",2)",\
                                              "+WQy*VRRS(1,4,",ms,",2)",\
                                              "-r1x2E*VRR(1,4,",ms,")"]];
             WS[StringJoin["VRRS(1,10,",m1,",2)=QCz*VRRS(1,4,",m1,",2)",\
                                              "+r1x2E*(VRRS(1,1,",m1,",2)-ZxZpE*VRRS(1,1,",ms,",2))",\
                                              "+WQz*VRRS(1,4,",ms,",2)"]];
             ];
             If[ixyz==3,
             (* i=z *)
             WS[StringJoin["VRRS(1, 5,",m1,",3)=QCx*VRRS(1,2,",m1,",3)",\
                                              "+r1x2E*(VRRS(1,1,",m1,",3)-ZxZpE*VRRS(1,1,",ms,",3))",\
                                              "+WQx*VRRS(1,2,",ms,",3)"]];
             WS[StringJoin["VRRS(1, 6,",m1,",3)=QCx*VRRS(1,3,",m1,",3)",\
                                              "+WQx*VRRS(1,3,",ms,",3)"]];
             WS[StringJoin["VRRS(1, 7,",m1,",3)=QCy*VRRS(1,3,",m1,",3)",\
                                              "+r1x2E*(VRRS(1,1,",m1,",3)-ZxZpE*VRRS(1,1,",ms,",3))",\
                                              "+WQy*VRRS(1,3,",ms,",3)"]];
             WS[StringJoin["VRRS(1, 8,",m1,",3)=QCx*VRRS(1,4,",m1,",3)",\
                                              "+WQx*VRRS(1,4,",ms,",3)"]];
             WS[StringJoin["VRRS(1, 9,",m1,",3)=QCy*VRRS(1,4,",m1,",3)",\
                                              "+WQy*VRRS(1,4,",ms,",3)"]];
             WS[StringJoin["VRRS(1,10,",m1,",3)=QCz*VRRS(1,4,",m1,",3)",\
                                              "+r1x2E*(VRRS(1,1,",m1,",3)-ZxZpE*VRRS(1,1,",ms,",3))",\
                                              "+WQz*VRRS(1,4,",ms,",3)",\
                                              "-r1x2E*VRR(1,4,",ms,")"]];
             ];
           ];
(*
           If[iell==2&&kell==1,
             If[ixyz==1,
             (* i=x *)
             (* xx,x,y,z *)
             WS[StringJoin["VRRS(5,2,",m1,",1)=QCx*VRRS(5,1,",m1,",1)+WQx*VRRS(5,1,",ms,",1)",\
                                              "+2d0*HfxZpE*VRRS(2,1,",ms,",1)",\
                                              "-r1x2E*VRR(5,1,",ms,")"]];
             WS[StringJoin["VRRS(5,3,",m1,",1)=QCy*VRRS(5,1,",m1,",1)+WQy*VRRS(5,1,",ms,",1)"]];
             WS[StringJoin["VRRS(5,4,",m1,",1)=QCz*VRRS(5,1,",m1,",1)+WQz*VRRS(5,1,",ms,",1)"]];
             (* yx,x,y,z *)
             WS[StringJoin["VRRS(6,2,",m1,",1)=QCx*VRRS(6,1,",m1,",1)+WQx*VRRS(6,1,",ms,",1)",\
                                              "+HfxZpE*VRRS(3,1,",ms,",1)",\
                                              "-r1x2E*VRR(6,1,",ms,")"]];
             WS[StringJoin["VRRS(6,3,",m1,",1)=QCy*VRRS(6,1,",m1,",1)+WQy*VRRS(6,1,",ms,",1)",\
                                              "+HfxZpE*VRRS(2,1,",ms,",1)"]];
             WS[StringJoin["VRRS(6,4,",m1,",1)=QCz*VRRS(6,1,",m1,",1)+WQz*VRRS(6,1,",ms,",1)"]];
             (* yy,x,y,z *)
             WS[StringJoin["VRRS(7,2,",m1,",1)=QCx*VRRS(7,1,",m1,",1)+WQx*VRRS(7,1,",ms,",1)",\
                                              "-r1x2E*VRR(7,1,",ms,")"]];
             WS[StringJoin["VRRS(7,3,",m1,",1)=QCy*VRRS(7,1,",m1,",1)+WQy*VRRS(7,1,",ms,",1)",\
                                              "+2d0*HfxZpE*VRRS(3,1,",ms,",1)"]];
             WS[StringJoin["VRRS(7,4,",m1,",1)=QCz*VRRS(7,1,",m1,",1)+WQz*VRRS(7,1,",ms,",1)"]];
             (* zx,x,y,z *)
             WS[StringJoin["VRRS(8,2,",m1,",1)=QCx*VRRS(8,1,",m1,",1)+WQx*VRRS(8,1,",ms,",1)",\
                                              "+HfxZpE*VRRS(4,1,",ms,",1)",\
                                              "-r1x2E*VRR(8,1,",ms,")"]];
             WS[StringJoin["VRRS(8,3,",m1,",1)=QCy*VRRS(8,1,",m1,",1)+WQy*VRRS(8,1,",ms,",1)"]];
             WS[StringJoin["VRRS(8,4,",m1,",1)=QCz*VRRS(8,1,",m1,",1)+WQz*VRRS(8,1,",ms,",1)",\
                                              "+HfxZpE*VRRS(2,1,",ms,",1)"]];
             (* zy,x,y,z *)
             WS[StringJoin["VRRS(9,2,",m1,",1)=QCx*VRRS(9,1,",m1,",1)+WQx*VRRS(9,1,",ms,",1)",\
                                              "-r1x2E*VRR(9,1,",ms,")"]];
             WS[StringJoin["VRRS(9,3,",m1,",1)=QCy*VRRS(9,1,",m1,",1)+WQy*VRRS(9,1,",ms,",1)",\
                                              "+HfxZpE*VRRS(4,1,",ms,",1)"]];
             WS[StringJoin["VRRS(9,4,",m1,",1)=QCz*VRRS(9,1,",m1,",1)+WQz*VRRS(9,1,",ms,",1)",\
                                              "+HfxZpE*VRRS(3,1,",ms,",1)"]];
             (* zz,x,y,z *)
             WS[StringJoin["VRRS(10,2,",m1,",1)=QCx*VRRS(10,1,",m1,",1)+WQx*VRRS(10,1,",ms,",1)",\
                                              "-r1x2E*VRR(10,1,",ms,")"]];
             WS[StringJoin["VRRS(10,3,",m1,",1)=QCy*VRRS(10,1,",m1,",1)+WQy*VRRS(10,1,",ms,",1)"]];
             WS[StringJoin["VRRS(10,4,",m1,",1)=QCz*VRRS(10,1,",m1,",1)+WQz*VRRS(10,1,",ms,",1)",\
                                              "+2d0*HfxZpE*VRRS(4,1,",ms,",1)"]];
             ];
             If[ixyz==2,
             (* i=y *)
             (* xx,x,y,z *)
             WS[StringJoin["VRRS(5,2,",m1,",2)=QCx*VRRS(5,1,",m1,",2)+WQx*VRRS(5,1,",ms,",2)",\
                                              "+2d0*HfxZpE*VRRS(2,1,",ms,",2)"]];
             WS[StringJoin["VRRS(5,3,",m1,",2)=QCy*VRRS(5,1,",m1,",2)+WQy*VRRS(5,1,",ms,",2)",\
                                              "-r1x2E*VRR(5,1,",ms,")"]];
             WS[StringJoin["VRRS(5,4,",m1,",2)=QCz*VRRS(5,1,",m1,",2)+WQz*VRRS(5,1,",ms,",2)"]];
             (* yx,x,y,z *)
             WS[StringJoin["VRRS(6,2,",m1,",2)=QCx*VRRS(6,1,",m1,",2)+WQx*VRRS(6,1,",ms,",2)",\
                                              "+HfxZpE*VRRS(3,1,",ms,",2)"]];
             WS[StringJoin["VRRS(6,3,",m1,",2)=QCy*VRRS(6,1,",m1,",2)+WQy*VRRS(6,1,",ms,",2)",\
                                              "+HfxZpE*VRRS(2,1,",ms,",2)",\
                                              "-r1x2E*VRR(6,1,",ms,")"]];
             WS[StringJoin["VRRS(6,4,",m1,",2)=QCz*VRRS(6,1,",m1,",2)+WQz*VRRS(6,1,",ms,",2)"]];
             (* yy,x,y,z *)
             WS[StringJoin["VRRS(7,2,",m1,",2)=QCx*VRRS(7,1,",m1,",2)+WQx*VRRS(7,1,",ms,",2)"]];
             WS[StringJoin["VRRS(7,3,",m1,",2)=QCy*VRRS(7,1,",m1,",2)+WQy*VRRS(7,1,",ms,",2)",\
                                              "+2d0*HfxZpE*VRRS(3,1,",ms,",2)",\
                                              "-r1x2E*VRR(7,1,",ms,")"]];
             WS[StringJoin["VRRS(7,4,",m1,",2)=QCz*VRRS(7,1,",m1,",2)+WQz*VRRS(7,1,",ms,",2)"]];
             (* zx,x,y,z *)
             WS[StringJoin["VRRS(8,2,",m1,",2)=QCx*VRRS(8,1,",m1,",2)+WQx*VRRS(8,1,",ms,",2)",\
                                              "+HfxZpE*VRRS(4,1,",ms,",2)"]];
             WS[StringJoin["VRRS(8,3,",m1,",2)=QCy*VRRS(8,1,",m1,",2)+WQy*VRRS(8,1,",ms,",2)",\
                                              "-r1x2E*VRR(8,1,",ms,")"]];
             WS[StringJoin["VRRS(8,4,",m1,",2)=QCz*VRRS(8,1,",m1,",2)+WQz*VRRS(8,1,",ms,",2)",\
                                              "+HfxZpE*VRRS(2,1,",ms,",2)"]];
             (* zy,x,y,z *)
             WS[StringJoin["VRRS(9,2,",m1,",2)=QCx*VRRS(9,1,",m1,",2)+WQx*VRRS(9,1,",ms,",2)"]];
             WS[StringJoin["VRRS(9,3,",m1,",2)=QCy*VRRS(9,1,",m1,",2)+WQy*VRRS(9,1,",ms,",2)",\
                                              "+HfxZpE*VRRS(4,1,",ms,",2)",\
                                              "-r1x2E*VRR(9,1,",ms,")"]];
             WS[StringJoin["VRRS(9,4,",m1,",2)=QCz*VRRS(9,1,",m1,",2)+WQz*VRRS(9,1,",ms,",2)",\
                                              "+HfxZpE*VRRS(3,1,",ms,",2)"]];
             (* zz,x,y,z *)
             WS[StringJoin["VRRS(10,2,",m1,",2)=QCx*VRRS(10,1,",m1,",2)+WQx*VRRS(10,1,",ms,",2)"]];
             WS[StringJoin["VRRS(10,3,",m1,",2)=QCy*VRRS(10,1,",m1,",2)+WQy*VRRS(10,1,",ms,",2)",\
                                              "-r1x2E*VRR(10,1,",ms,")"]];
             WS[StringJoin["VRRS(10,4,",m1,",2)=QCz*VRRS(10,1,",m1,",2)+WQz*VRRS(10,1,",ms,",2)",\
                                              "+2d0*HfxZpE*VRRS(4,1,",ms,",2)"]];
             ];
             If[ixyz==3,
             (* i=z *)
             (* xx,x,y,z *)
             WS[StringJoin["VRRS(5,2,",m1,",3)=QCx*VRRS(5,1,",m1,",3)+WQx*VRRS(5,1,",ms,",3)",\
                                              "+2d0*HfxZpE*VRRS(2,1,",ms,",3)"]];
             WS[StringJoin["VRRS(5,3,",m1,",3)=QCy*VRRS(5,1,",m1,",3)+WQy*VRRS(5,1,",ms,",3)"]];
             WS[StringJoin["VRRS(5,4,",m1,",3)=QCz*VRRS(5,1,",m1,",3)+WQz*VRRS(5,1,",ms,",3)",\
                                              "-r1x2E*VRR(5,1,",ms,")"]];
             (* yx,x,y,z *)
             WS[StringJoin["VRRS(6,2,",m1,",3)=QCx*VRRS(6,1,",m1,",3)+WQx*VRRS(6,1,",ms,",3)",\
                                              "+HfxZpE*VRRS(3,1,",ms,",3)"]];
             WS[StringJoin["VRRS(6,3,",m1,",3)=QCy*VRRS(6,1,",m1,",3)+WQy*VRRS(6,1,",ms,",3)",\
                                              "+HfxZpE*VRRS(2,1,",ms,",3)"]];
             WS[StringJoin["VRRS(6,4,",m1,",3)=QCz*VRRS(6,1,",m1,",3)+WQz*VRRS(6,1,",ms,",3)",\
                                              "-r1x2E*VRR(6,1,",ms,")"]];
             (* yy,x,y,z *)
             WS[StringJoin["VRRS(7,2,",m1,",3)=QCx*VRRS(7,1,",m1,",3)+WQx*VRRS(7,1,",ms,",3)"]];
             WS[StringJoin["VRRS(7,3,",m1,",3)=QCy*VRRS(7,1,",m1,",3)+WQy*VRRS(7,1,",ms,",3)",\
                                              "+2d0*HfxZpE*VRRS(3,1,",ms,",3)"]];
             WS[StringJoin["VRRS(7,4,",m1,",3)=QCz*VRRS(7,1,",m1,",3)+WQz*VRRS(7,1,",ms,",3)",\
                                              "-r1x2E*VRR(7,1,",ms,")"]];
             (* zx,x,y,z *)
             WS[StringJoin["VRRS(8,2,",m1,",3)=QCx*VRRS(8,1,",m1,",3)+WQx*VRRS(8,1,",ms,",3)",\
                                              "+HfxZpE*VRRS(4,1,",ms,",3)"]];
             WS[StringJoin["VRRS(8,3,",m1,",3)=QCy*VRRS(8,1,",m1,",3)+WQy*VRRS(8,1,",ms,",3)"]];
             WS[StringJoin["VRRS(8,4,",m1,",3)=QCz*VRRS(8,1,",m1,",3)+WQz*VRRS(8,1,",ms,",3)",\
                                              "+HfxZpE*VRRS(2,1,",ms,",3)",\
                                              "-r1x2E*VRR(8,1,",ms,")"]];
             (* zy,x,y,z *)
             WS[StringJoin["VRRS(9,2,",m1,",3)=QCx*VRRS(9,1,",m1,",3)+WQx*VRRS(9,1,",ms,",3)"]];
             WS[StringJoin["VRRS(9,3,",m1,",3)=QCy*VRRS(9,1,",m1,",3)+WQy*VRRS(9,1,",ms,",3)",\
                                              "+HfxZpE*VRRS(4,1,",ms,",3)"]];
             WS[StringJoin["VRRS(9,4,",m1,",3)=QCz*VRRS(9,1,",m1,",3)+WQz*VRRS(9,1,",ms,",3)",\
                                              "+HfxZpE*VRRS(3,1,",ms,",3)",\
                                              "-r1x2E*VRR(9,1,",ms,")"]];
             (* zz,x,y,z *)
             WS[StringJoin["VRRS(10,2,",m1,",3)=QCx*VRRS(10,1,",m1,",3)+WQx*VRRS(10,1,",ms,",3)"]];
             WS[StringJoin["VRRS(10,3,",m1,",3)=QCy*VRRS(10,1,",m1,",3)+WQy*VRRS(10,1,",ms,",3)"]];
             WS[StringJoin["VRRS(10,4,",m1,",3)=QCz*VRRS(10,1,",m1,",3)+WQz*VRRS(10,1,",ms,",3)",\
                                              "+2d0*HfxZpE*VRRS(4,1,",ms,",3)",\
                                              "-r1x2E*VRR(10,1,",ms,")"]];
             ];
           ];

           If[iell==1&&kell==2,
             If[ixyz==1,
             (* i=x *)
             (* x,y,z,xx *)
             WS[StringJoin["VRRS(2,5,",m1,",1)=PAx*VRRS(1,5,",m1,",1)+WPx*VRRS(1,5,",ms,",1)",\
                                              "+2d0*HfxZpE*VRRS(1,2,",ms,",1)",\
                                              "+r1x2Z*VRR(1,5,",ms,")"]];
             WS[StringJoin["VRRS(3,5,",m1,",1)=PAy*VRRS(1,5,",m1,",1)+WPy*VRRS(1,5,",ms,",1)"]];
             WS[StringJoin["VRRS(4,5,",m1,",1)=PAz*VRRS(1,5,",m1,",1)+WPz*VRRS(1,5,",ms,",1)"]];
             (* x,y,z,yx *)
             WS[StringJoin["VRRS(2,6,",m1,",1)=PAx*VRRS(1,6,",m1,",1)+WPx*VRRS(1,6,",ms,",1)",\
                                              "+HfxZpE*VRRS(1,3,",ms,",1)",\
                                              "+r1x2Z*VRR(1,6,",ms,")"]];
             WS[StringJoin["VRRS(3,6,",m1,",1)=PAy*VRRS(1,6,",m1,",1)+WPy*VRRS(1,6,",ms,",1)",\
                                              "+HfxZpE*VRRS(1,2,",ms,",1)"]];
             WS[StringJoin["VRRS(4,6,",m1,",1)=PAz*VRRS(1,6,",m1,",1)+WPz*VRRS(1,6,",ms,",1)"]];
             (* x,y,z,yy *)
             WS[StringJoin["VRRS(2,7,",m1,",1)=PAx*VRRS(1,7,",m1,",1)+WPx*VRRS(1,7,",ms,",1)",\
                                              "+r1x2Z*VRR(1,7,",ms,")"]];
             WS[StringJoin["VRRS(3,7,",m1,",1)=PAy*VRRS(1,7,",m1,",1)+WPy*VRRS(1,7,",ms,",1)",\
                                              "+2d0*HfxZpE*VRRS(1,3,",ms,",1)"]];
             WS[StringJoin["VRRS(4,7,",m1,",1)=PAz*VRRS(1,7,",m1,",1)+WPz*VRRS(1,7,",ms,",1)"]];
             (* x,y,z,zx *)
             WS[StringJoin["VRRS(2,8,",m1,",1)=PAx*VRRS(1,8,",m1,",1)+WPx*VRRS(1,8,",ms,",1)",\
                                              "+HfxZpE*VRRS(1,4,",ms,",1)",\
                                              "+r1x2Z*VRR(1,8,",ms,")"]];
             WS[StringJoin["VRRS(3,8,",m1,",1)=PAy*VRRS(1,8,",m1,",1)+WPy*VRRS(1,8,",ms,",1)"]];
             WS[StringJoin["VRRS(4,8,",m1,",1)=PAz*VRRS(1,8,",m1,",1)+WPz*VRRS(1,8,",ms,",1)",\
                                              "+HfxZpE*VRRS(1,2,",ms,",1)"]];
             (* x,y,z,zy *)
             WS[StringJoin["VRRS(2,9,",m1,",1)=PAx*VRRS(1,9,",m1,",1)+WPx*VRRS(1,9,",ms,",1)",\
                                              "+r1x2Z*VRR(1,9,",ms,")"]];
             WS[StringJoin["VRRS(3,9,",m1,",1)=PAy*VRRS(1,9,",m1,",1)+WPy*VRRS(1,9,",ms,",1)",\
                                              "+HfxZpE*VRRS(1,4,",ms,",1)"]];
             WS[StringJoin["VRRS(4,9,",m1,",1)=PAz*VRRS(1,9,",m1,",1)+WPz*VRRS(1,9,",ms,",1)",\
                                              "+HfxZpE*VRRS(1,3,",ms,",1)"]];
             (* x,y,z,zz *)
             WS[StringJoin["VRRS(2,10,",m1,",1)=PAx*VRRS(1,10,",m1,",1)+WPx*VRRS(1,10,",ms,",1)",\
                                              "+r1x2Z*VRR(1,10,",ms,")"]];
             WS[StringJoin["VRRS(3,10,",m1,",1)=PAy*VRRS(1,10,",m1,",1)+WPy*VRRS(1,10,",ms,",1)"]];
             WS[StringJoin["VRRS(4,10,",m1,",1)=PAz*VRRS(1,10,",m1,",1)+WPz*VRRS(1,10,",ms,",1)",\
                                              "+2d0*HfxZpE*VRRS(1,4,",ms,",1)"]];
             ];
             If[ixyz==2,
             (* i=y *)
             (* x,y,z,xx *)
             WS[StringJoin["VRRS(2,5,",m1,",2)=PAx*VRRS(1,5,",m1,",2)+WPx*VRRS(1,5,",ms,",2)",\
                                              "+2d0*HfxZpE*VRRS(1,2,",ms,",2)"]];
             WS[StringJoin["VRRS(3,5,",m1,",2)=PAy*VRRS(1,5,",m1,",2)+WPy*VRRS(1,5,",ms,",2)",\
                                              "+r1x2Z*VRR(1,5,",ms,")"]];
             WS[StringJoin["VRRS(4,5,",m1,",2)=PAz*VRRS(1,5,",m1,",2)+WPz*VRRS(1,5,",ms,",2)"]];
             (* x,y,z,yx *)
             WS[StringJoin["VRRS(2,6,",m1,",2)=PAx*VRRS(1,6,",m1,",2)+WPx*VRRS(1,6,",ms,",2)",\
                                              "+HfxZpE*VRRS(1,3,",ms,",2)"]];
             WS[StringJoin["VRRS(3,6,",m1,",2)=PAy*VRRS(1,6,",m1,",2)+WPy*VRRS(1,6,",ms,",2)",\
                                              "+HfxZpE*VRRS(1,2,",ms,",2)",\
                                              "+r1x2Z*VRR(1,6,",ms,")"]];
             WS[StringJoin["VRRS(4,6,",m1,",2)=PAz*VRRS(1,6,",m1,",2)+WPz*VRRS(1,6,",ms,",2)"]];
             (* x,y,z,yy *)
             WS[StringJoin["VRRS(2,7,",m1,",2)=PAx*VRRS(1,7,",m1,",2)+WPx*VRRS(1,7,",ms,",2)"]];
             WS[StringJoin["VRRS(3,7,",m1,",2)=PAy*VRRS(1,7,",m1,",2)+WPy*VRRS(1,7,",ms,",2)",\
                                              "+2d0*HfxZpE*VRRS(1,3,",ms,",2)",\
                                              "+r1x2Z*VRR(1,7,",ms,")"]];
             WS[StringJoin["VRRS(4,7,",m1,",2)=PAz*VRRS(1,7,",m1,",2)+WPz*VRRS(1,7,",ms,",2)"]];
             (* x,y,z,zx *)
             WS[StringJoin["VRRS(2,8,",m1,",2)=PAx*VRRS(1,8,",m1,",2)+WPx*VRRS(1,8,",ms,",2)",\
                                              "+HfxZpE*VRRS(1,4,",ms,",2)"]];
             WS[StringJoin["VRRS(3,8,",m1,",2)=PAy*VRRS(1,8,",m1,",2)+WPy*VRRS(1,8,",ms,",2)",\
                                              "+r1x2Z*VRR(1,8,",ms,")"]];
             WS[StringJoin["VRRS(4,8,",m1,",2)=PAz*VRRS(1,8,",m1,",2)+WPz*VRRS(1,8,",ms,",2)",\
                                              "+HfxZpE*VRRS(1,2,",ms,",2)"]];
             (* x,y,z,zy *)
             WS[StringJoin["VRRS(2,9,",m1,",2)=PAx*VRRS(1,9,",m1,",2)+WPx*VRRS(1,9,",ms,",2)"]];
             WS[StringJoin["VRRS(3,9,",m1,",2)=PAy*VRRS(1,9,",m1,",2)+WPy*VRRS(1,9,",ms,",2)",\
                                              "+HfxZpE*VRRS(1,4,",ms,",2)",\
                                              "+r1x2Z*VRR(1,9,",ms,")"]];
             WS[StringJoin["VRRS(4,9,",m1,",2)=PAz*VRRS(1,9,",m1,",2)+WPz*VRRS(1,9,",ms,",2)",\
                                              "+HfxZpE*VRRS(1,3,",ms,",2)"]];
             (* x,y,z,zz *)
             WS[StringJoin["VRRS(2,10,",m1,",2)=PAx*VRRS(1,10,",m1,",2)+WPx*VRRS(1,10,",ms,",2)"]];
             WS[StringJoin["VRRS(3,10,",m1,",2)=PAy*VRRS(1,10,",m1,",2)+WPy*VRRS(1,10,",ms,",2)",\
                                              "+r1x2Z*VRR(1,10,",ms,")"]];
             WS[StringJoin["VRRS(4,10,",m1,",2)=PAz*VRRS(1,10,",m1,",2)+WPz*VRRS(1,10,",ms,",2)",\
                                              "+2d0*HfxZpE*VRRS(1,4,",ms,",2)"]];
             ];
             If[ixyz==3,
             (* i=z *)
             (* x,y,z,xx *)
             WS[StringJoin["VRRS(2,5,",m1,",3)=PAx*VRRS(1,5,",m1,",3)+WPx*VRRS(1,5,",ms,",3)",\
                                              "+2d0*HfxZpE*VRRS(1,2,",ms,",3)"]];
             WS[StringJoin["VRRS(3,5,",m1,",3)=PAy*VRRS(1,5,",m1,",3)+WPy*VRRS(1,5,",ms,",3)"]];
             WS[StringJoin["VRRS(4,5,",m1,",3)=PAz*VRRS(1,5,",m1,",3)+WPz*VRRS(1,5,",ms,",3)",\
                                              "+r1x2Z*VRR(1,5,",ms,")"]];
             (* x,y,z,yx *)
             WS[StringJoin["VRRS(2,6,",m1,",3)=PAx*VRRS(1,6,",m1,",3)+WPx*VRRS(1,6,",ms,",3)",\
                                              "+HfxZpE*VRRS(1,3,",ms,",3)"]];
             WS[StringJoin["VRRS(3,6,",m1,",3)=PAy*VRRS(1,6,",m1,",3)+WPy*VRRS(1,6,",ms,",3)",\
                                              "+HfxZpE*VRRS(1,2,",ms,",3)"]];
             WS[StringJoin["VRRS(4,6,",m1,",3)=PAz*VRRS(1,6,",m1,",3)+WPz*VRRS(1,6,",ms,",3)",\
                                              "+r1x2Z*VRR(1,6,",ms,")"]];
             (* x,y,z,yy *)
             WS[StringJoin["VRRS(2,7,",m1,",3)=PAx*VRRS(1,7,",m1,",3)+WPx*VRRS(1,7,",ms,",3)"]];
             WS[StringJoin["VRRS(3,7,",m1,",3)=PAy*VRRS(1,7,",m1,",3)+WPy*VRRS(1,7,",ms,",3)",\
                                              "+2d0*HfxZpE*VRRS(1,3,",ms,",3)"]];
             WS[StringJoin["VRRS(4,7,",m1,",3)=PAz*VRRS(1,7,",m1,",3)+WPz*VRRS(1,7,",ms,",3)",\
                                              "+r1x2Z*VRR(1,7,",ms,")"]];
             (* x,y,z,zx *)
             WS[StringJoin["VRRS(2,8,",m1,",3)=PAx*VRRS(1,8,",m1,",3)+WPx*VRRS(1,8,",ms,",3)",\
                                              "+HfxZpE*VRRS(1,4,",ms,",3)"]];
             WS[StringJoin["VRRS(3,8,",m1,",3)=PAy*VRRS(1,8,",m1,",3)+WPy*VRRS(1,8,",ms,",3)"]];
             WS[StringJoin["VRRS(4,8,",m1,",3)=PAz*VRRS(1,8,",m1,",3)+WPz*VRRS(1,8,",ms,",3)",\
                                              "+HfxZpE*VRRS(1,2,",ms,",3)",\
                                              "+r1x2Z*VRR(1,8,",ms,")"]];
             (* x,y,z,zy *)
             WS[StringJoin["VRRS(2,9,",m1,",3)=PAx*VRRS(1,9,",m1,",3)+WPx*VRRS(1,9,",ms,",3)"]];
             WS[StringJoin["VRRS(3,9,",m1,",3)=PAy*VRRS(1,9,",m1,",3)+WPy*VRRS(1,9,",ms,",3)",\
                                              "+HfxZpE*VRRS(1,4,",ms,",3)"]];
             WS[StringJoin["VRRS(4,9,",m1,",3)=PAz*VRRS(1,9,",m1,",3)+WPz*VRRS(1,9,",ms,",3)",\
                                              "+HfxZpE*VRRS(1,3,",ms,",3)",\
                                              "+r1x2Z*VRR(1,9,",ms,")"]];
             (* x,y,z,zz *)
             WS[StringJoin["VRRS(2,10,",m1,",3)=PAx*VRRS(1,10,",m1,",3)+WPx*VRRS(1,10,",ms,",3)"]];
             WS[StringJoin["VRRS(3,10,",m1,",3)=PAy*VRRS(1,10,",m1,",3)+WPy*VRRS(1,10,",ms,",3)"]];
             WS[StringJoin["VRRS(4,10,",m1,",3)=PAz*VRRS(1,10,",m1,",3)+WPz*VRRS(1,10,",ms,",3)",\
                                              "+2d0*HfxZpE*VRRS(1,4,",ms,",3)",\
                                              "+r1x2Z*VRR(1,10,",ms,")"]];
             ];
           ];


           If[iell==2&&kell==2,
             If[ixyz==1,
             (* i=x *)
             (* xx,xx *)
             WS[StringJoin["VRRS(5,5,",m1,",1)=QCx*VRRS(5,2,",m1,",1)+WQx*VRRS(5,2,",ms,",1)",\
                                    "+r1x2E*(VRRS(5,1,",m1,",1)-ZxZpE*VRRS(5,1,",ms,",1)) &"]];
             WS[StringJoin["             +2d0*HfxZpE*VRRS(2,2,",ms,",1)",\
                                         "-r1x2E*VRR(5,2,",ms,")"]];
             (* xx,yx *)
             WS[StringJoin["VRRS(5,6,",m1,",1)=QCy*VRRS(5,2,",m1,",1)+WQy*VRRS(5,2,",ms,",1)"]];
             (* xx,yy *)
             WS[StringJoin["VRRS(5,7,",m1,",1)=QCy*VRRS(5,3,",m1,",1)+WQy*VRRS(5,3,",ms,",1)",\
                                    "+r1x2E*(VRRS(5,1,",m1,",1)-ZxZpE*VRRS(5,1,",ms,",1))"]];
             (* xx,zx *)
             WS[StringJoin["VRRS(5,8,",m1,",1)=QCz*VRRS(5,2,",m1,",1)+WQz*VRRS(5,2,",ms,",1)"]];
             (* xx,zy *)
             WS[StringJoin["VRRS(5,9,",m1,",1)=QCz*VRRS(5,3,",m1,",1)+WQz*VRRS(5,3,",ms,",1)"]];
             (* xx,zz *)
             WS[StringJoin["VRRS(5,10,",m1,",1)=QCz*VRRS(5,4,",m1,",1)+WQz*VRRS(5,4,",ms,",1)",\
                                    "+r1x2E*(VRRS(5,1,",m1,",1)-ZxZpE*VRRS(5,1,",ms,",1))"]];
             (* yx,xx *)
             WS[StringJoin["VRRS(6,5,",m1,",1)=QCx*VRRS(6,2,",m1,",1)+WQx*VRRS(6,2,",ms,",1)",\
                                    "+r1x2E*(VRRS(6,1,",m1,",1)-ZxZpE*VRRS(6,1,",ms,",1)) &"]];
             WS[StringJoin["             +HfxZpE*VRRS(3,2,",ms,",1)",\
                                         "-r1x2E*VRR(6,2,",ms,")"]];
             (* yx,yx *)
             WS[StringJoin["VRRS(6,6,",m1,",1)=QCy*VRRS(6,2,",m1,",1)+WQy*VRRS(6,2,",ms,",1)",\
                                        "+HfxZpE*VRRS(2,2,",ms,",1)"]];
             (* yx,yy *)
             WS[StringJoin["VRRS(6,7,",m1,",1)=QCy*VRRS(6,3,",m1,",1)+WQy*VRRS(6,3,",ms,",1)",\
                                    "+r1x2E*(VRRS(6,1,",m1,",1)-ZxZpE*VRRS(6,1,",ms,",1)) &"]];
             WS[StringJoin["             +HfxZpE*VRRS(2,3,",ms,",1)"]];
             (* yx,zx *)
             WS[StringJoin["VRRS(6,8,",m1,",1)=QCz*VRRS(6,2,",m1,",1)+WQz*VRRS(6,2,",ms,",1)"]];
             (* yx,zy *)
             WS[StringJoin["VRRS(6,9,",m1,",1)=QCz*VRRS(6,3,",m1,",1)+WQz*VRRS(6,3,",ms,",1)"]];
             (* yx,zz *)
             WS[StringJoin["VRRS(6,10,",m1,",1)=QCz*VRRS(6,4,",m1,",1)+WQz*VRRS(6,4,",ms,",1)",\
                                    "+r1x2E*(VRRS(6,1,",m1,",1)-ZxZpE*VRRS(6,1,",ms,",1))"]];
             (* yy,xx *)
             WS[StringJoin["VRRS(7,5,",m1,",1)=QCx*VRRS(7,2,",m1,",1)+WQx*VRRS(7,2,",ms,",1)",\
                                    "+r1x2E*(VRRS(7,1,",m1,",1)-ZxZpE*VRRS(7,1,",ms,",1)) &"]];
             WS[StringJoin["             -r1x2E*VRR(7,2,",ms,")"]];
             (* yy,yx *)
             WS[StringJoin["VRRS(7,6,",m1,",1)=QCy*VRRS(7,2,",m1,",1)+WQy*VRRS(7,2,",ms,",1)",\
                                   "+2d0*HfxZpE*VRRS(3,2,",ms,",1)"]];
             (* yy,yy *)
             WS[StringJoin["VRRS(7,7,",m1,",1)=QCy*VRRS(7,3,",m1,",1)+WQy*VRRS(7,3,",ms,",1)",\
                                    "+r1x2E*(VRRS(7,1,",m1,",1)-ZxZpE*VRRS(7,1,",ms,",1)) &"]];
             WS[StringJoin["             +2d0*HfxZpE*VRRS(3,3,",ms,",1)"]];
             (* yy,zx *)
             WS[StringJoin["VRRS(7,8,",m1,",1)=QCz*VRRS(7,2,",m1,",1)+WQz*VRRS(7,2,",ms,",1)"]];
             (* yy,zy *)
             WS[StringJoin["VRRS(7,9,",m1,",1)=QCz*VRRS(7,3,",m1,",1)+WQz*VRRS(7,3,",ms,",1)"]];
             (* yy,zz *)
             WS[StringJoin["VRRS(7,10,",m1,",1)=QCz*VRRS(7,4,",m1,",1)+WQz*VRRS(7,4,",ms,",1)",\
                                    "+r1x2E*(VRRS(7,1,",m1,",1)-ZxZpE*VRRS(7,1,",ms,",1))"]];
             (* zx,xx *)
             WS[StringJoin["VRRS(8,5,",m1,",1)=QCx*VRRS(8,2,",m1,",1)+WQx*VRRS(8,2,",ms,",1)",\
                                    "+r1x2E*(VRRS(8,1,",m1,",1)-ZxZpE*VRRS(8,1,",ms,",1)) &"]];
             WS[StringJoin["             +HfxZpE*VRRS(4,2,",ms,",1)",\
                                        "-r1x2E*VRR(8,2,",ms,")"]];
             (* zx,yx *)
             WS[StringJoin["VRRS(8,6,",m1,",1)=QCy*VRRS(8,2,",m1,",1)+WQy*VRRS(8,2,",ms,",1)"]];
             (* zx,yy *)
             WS[StringJoin["VRRS(8,7,",m1,",1)=QCy*VRRS(8,3,",m1,",1)+WQy*VRRS(8,3,",ms,",1)",\
                                    "+r1x2E*(VRRS(8,1,",m1,",1)-ZxZpE*VRRS(8,1,",ms,",1)) "]];
             (* zx,zx *)
             WS[StringJoin["VRRS(8,8,",m1,",1)=QCz*VRRS(8,2,",m1,",1)+WQz*VRRS(8,2,",ms,",1)",\
                                        "+HfxZpE*VRRS(2,2,",ms,",1)"]];
             (* zx,zy *)
             WS[StringJoin["VRRS(8,9,",m1,",1)=QCz*VRRS(8,3,",m1,",1)+WQz*VRRS(8,3,",ms,",1)",\
                                        "+HfxZpE*VRRS(2,3,",ms,",1)"]];
             (* zx,zz *)
             WS[StringJoin["VRRS(8,10,",m1,",1)=QCz*VRRS(8,4,",m1,",1)+WQz*VRRS(8,4,",ms,",1)",\
                                    "+r1x2E*(VRRS(8,1,",m1,",1)-ZxZpE*VRRS(8,1,",ms,",1)) &"]];
             WS[StringJoin["             +HfxZpE*VRRS(2,4,",ms,",1)"]];

             (* zy,xx *)
             WS[StringJoin["VRRS(9,5,",m1,",1)=QCx*VRRS(9,2,",m1,",1)+WQx*VRRS(9,2,",ms,",1)",\
                                    "+r1x2E*(VRRS(9,1,",m1,",1)-ZxZpE*VRRS(9,1,",ms,",1)) &"]];
             WS[StringJoin["            -r1x2E*VRR(9,2,",ms,")"]];
             (* zy,yx *)
             WS[StringJoin["VRRS(9,6,",m1,",1)=QCy*VRRS(9,2,",m1,",1)+WQy*VRRS(9,2,",ms,",1)",\
                                        "+HfxZpE*VRRS(4,2,",ms,",1)"]];
             (* zy,yy *)
             WS[StringJoin["VRRS(9,7,",m1,",1)=QCy*VRRS(9,3,",m1,",1)+WQy*VRRS(9,3,",ms,",1)",\
                                    "+r1x2E*(VRRS(9,1,",m1,",1)-ZxZpE*VRRS(9,1,",ms,",1)) &"]];
             WS[StringJoin["            +HfxZpE*VRRS(4,3,",ms,",1)"]];
             (* zy,zx *)
             WS[StringJoin["VRRS(9,8,",m1,",1)=QCz*VRRS(9,2,",m1,",1)+WQz*VRRS(9,2,",ms,",1)",\
                                        "+HfxZpE*VRRS(3,2,",ms,",1)"]];
             (* zy,zy *)
             WS[StringJoin["VRRS(9,9,",m1,",1)=QCz*VRRS(9,3,",m1,",1)+WQz*VRRS(9,3,",ms,",1)",\
                                        "+HfxZpE*VRRS(3,3,",ms,",1)"]];
             (* zy,zz *)
             WS[StringJoin["VRRS(9,10,",m1,",1)=QCz*VRRS(9,4,",m1,",1)+WQz*VRRS(9,4,",ms,",1)",\
                                    "+r1x2E*(VRRS(9,1,",m1,",1)-ZxZpE*VRRS(9,1,",ms,",1)) &"]];
             WS[StringJoin["             +HfxZpE*VRRS(3,4,",ms,",1)"]];
             (* zz,xx *)
             WS[StringJoin["VRRS(10,5,",m1,",1)=QCx*VRRS(10,2,",m1,",1)+WQx*VRRS(10,2,",ms,",1)",\
                                    "+r1x2E*(VRRS(10,1,",m1,",1)-ZxZpE*VRRS(10,1,",ms,",1)) &"]];
             WS[StringJoin["            -r1x2E*VRR(10,2,",ms,")"]];
             (* zz,yx *)
             WS[StringJoin["VRRS(10,6,",m1,",1)=QCy*VRRS(10,2,",m1,",1)+WQy*VRRS(10,2,",ms,",1)"]];
             (* zz,yy *)
             WS[StringJoin["VRRS(10,7,",m1,",1)=QCy*VRRS(10,3,",m1,",1)+WQy*VRRS(10,3,",ms,",1)",\
                                    "+r1x2E*(VRRS(10,1,",m1,",1)-ZxZpE*VRRS(10,1,",ms,",1)) "]];
             (* zz,zx *)
             WS[StringJoin["VRRS(10,8,",m1,",1)=QCz*VRRS(10,2,",m1,",1)+WQz*VRRS(10,2,",ms,",1)",\
                                        "+2d0*HfxZpE*VRRS(4,2,",ms,",1)"]];
             (* zz,zy *)
             WS[StringJoin["VRRS(10,9,",m1,",1)=QCz*VRRS(10,3,",m1,",1)+WQz*VRRS(10,3,",ms,",1)",\
                                        "+2d0*HfxZpE*VRRS(4,3,",ms,",1)"]];
             (* zz,zz *)
             WS[StringJoin["VRRS(10,10,",m1,",1)=QCz*VRRS(10,4,",m1,",1)+WQz*VRRS(10,4,",ms,",1)",\
                                    "+r1x2E*(VRRS(10,1,",m1,",1)-ZxZpE*VRRS(10,1,",ms,",1)) &"]];
             WS[StringJoin["             +2d0*HfxZpE*VRRS(4,4,",ms,",1)"]];
             ];
             If[ixyz==2,
             (* i=y *)
             (* xx,xx *)
             WS[StringJoin["VRRS(5,5,",m1,",2)=QCx*VRRS(5,2,",m1,",2)+WQx*VRRS(5,2,",ms,",2)",\
                                    "+r1x2E*(VRRS(5,1,",m1,",2)-ZxZpE*VRRS(5,1,",ms,",2)) &"]];
             WS[StringJoin["             +2d0*HfxZpE*VRRS(2,2,",ms,",2)"]];
             (* xx,yx *)
             WS[StringJoin["VRRS(5,6,",m1,",2)=QCy*VRRS(5,2,",m1,",2)+WQy*VRRS(5,2,",ms,",2)",\
                                         "-r1x2E*VRR(5,2,",ms,")"]];
             (* xx,yy *)
             WS[StringJoin["VRRS(5,7,",m1,",2)=QCy*VRRS(5,3,",m1,",2)+WQy*VRRS(5,3,",ms,",2)",\
                                    "+r1x2E*(VRRS(5,1,",m1,",2)-ZxZpE*VRRS(5,1,",ms,",2)) &"]];
             WS[StringJoin["             -r1x2E*VRR(5,3,",ms,")"]];
             (* xx,zx *)
             WS[StringJoin["VRRS(5,8,",m1,",2)=QCz*VRRS(5,2,",m1,",2)+WQz*VRRS(5,2,",ms,",2)"]];
             (* xx,zy *)
             WS[StringJoin["VRRS(5,9,",m1,",2)=QCz*VRRS(5,3,",m1,",2)+WQz*VRRS(5,3,",ms,",2)"]];
             (* xx,zz *)
             WS[StringJoin["VRRS(5,10,",m1,",2)=QCz*VRRS(5,4,",m1,",2)+WQz*VRRS(5,4,",ms,",2)",\
                                    "+r1x2E*(VRRS(5,1,",m1,",2)-ZxZpE*VRRS(5,1,",ms,",2))"]];
             (* yx,xx *)
             WS[StringJoin["VRRS(6,5,",m1,",2)=QCx*VRRS(6,2,",m1,",2)+WQx*VRRS(6,2,",ms,",2)",\
                                    "+r1x2E*(VRRS(6,1,",m1,",2)-ZxZpE*VRRS(6,1,",ms,",2)) &"]];
             WS[StringJoin["             +HfxZpE*VRRS(3,2,",ms,",2)"]]
             (* yx,yx *)
             WS[StringJoin["VRRS(6,6,",m1,",2)=QCy*VRRS(6,2,",m1,",2)+WQy*VRRS(6,2,",ms,",2)",\
                                        "+HfxZpE*VRRS(2,2,",ms,",2) &"]];
             WS[StringJoin["             -r1x2E*VRR(6,2,",ms,")"]];
             (* yx,yy *)
             WS[StringJoin["VRRS(6,7,",m1,",2)=QCy*VRRS(6,3,",m1,",2)+WQy*VRRS(6,3,",ms,",2)",\
                                    "+r1x2E*(VRRS(6,1,",m1,",2)-ZxZpE*VRRS(6,1,",ms,",2)) &"]];
             WS[StringJoin["             +HfxZpE*VRRS(2,3,",ms,",2)",\
                                         "-r1x2E*VRR(6,3,",ms,")"]];
             (* yx,zx *)
             WS[StringJoin["VRRS(6,8,",m1,",2)=QCz*VRRS(6,2,",m1,",2)+WQz*VRRS(6,2,",ms,",2)"]];
             (* yx,zy *)
             WS[StringJoin["VRRS(6,9,",m1,",2)=QCz*VRRS(6,3,",m1,",2)+WQz*VRRS(6,3,",ms,",2)"]];
             (* yx,zz *)
             WS[StringJoin["VRRS(6,10,",m1,",2)=QCz*VRRS(6,4,",m1,",2)+WQz*VRRS(6,4,",ms,",2)",\
                                    "+r1x2E*(VRRS(6,1,",m1,",2)-ZxZpE*VRRS(6,1,",ms,",2))"]];
             (* yy,xx *)
             WS[StringJoin["VRRS(7,5,",m1,",2)=QCx*VRRS(7,2,",m1,",2)+WQx*VRRS(7,2,",ms,",2)",\
                                    "+r1x2E*(VRRS(7,1,",m1,",2)-ZxZpE*VRRS(7,1,",ms,",2))"]];
             (* yy,yx *)
             WS[StringJoin["VRRS(7,6,",m1,",2)=QCy*VRRS(7,2,",m1,",2)+WQy*VRRS(7,2,",ms,",2)",\
                                   "+2d0*HfxZpE*VRRS(3,2,",ms,",2) &"]];
             WS[StringJoin["             -r1x2E*VRR(7,2,",ms,")"]];
             (* yy,yy *)
             WS[StringJoin["VRRS(7,7,",m1,",2)=QCy*VRRS(7,3,",m1,",2)+WQy*VRRS(7,3,",ms,",2)",\
                                    "+r1x2E*(VRRS(7,1,",m1,",2)-ZxZpE*VRRS(7,1,",ms,",2)) &"]];
             WS[StringJoin["             +2d0*HfxZpE*VRRS(3,3,",ms,",2)",\
                                        "-r1x2E*VRR(7,3,",ms,")"]];
             (* yy,zx *)
             WS[StringJoin["VRRS(7,8,",m1,",2)=QCz*VRRS(7,2,",m1,",2)+WQz*VRRS(7,2,",ms,",2)"]];
             (* yy,zy *)
             WS[StringJoin["VRRS(7,9,",m1,",2)=QCz*VRRS(7,3,",m1,",2)+WQz*VRRS(7,3,",ms,",2)"]];
             (* yy,zz *)
             WS[StringJoin["VRRS(7,10,",m1,",2)=QCz*VRRS(7,4,",m1,",2)+WQz*VRRS(7,4,",ms,",2)",\
                                    "+r1x2E*(VRRS(7,1,",m1,",2)-ZxZpE*VRRS(7,1,",ms,",2))"]];
             (* zx,xx *)
             WS[StringJoin["VRRS(8,5,",m1,",2)=QCx*VRRS(8,2,",m1,",2)+WQx*VRRS(8,2,",ms,",2)",\
                                    "+r1x2E*(VRRS(8,1,",m1,",2)-ZxZpE*VRRS(8,1,",ms,",2)) &"]];
             WS[StringJoin["             +HfxZpE*VRRS(4,2,",ms,",2)"]];
             (* zx,yx *)
             WS[StringJoin["VRRS(8,6,",m1,",2)=QCy*VRRS(8,2,",m1,",2)+WQy*VRRS(8,2,",ms,",2)",\
                                        "-r1x2E*VRR(8,2,",ms,")"]];
             (* zx,yy *)
             WS[StringJoin["VRRS(8,7,",m1,",2)=QCy*VRRS(8,3,",m1,",2)+WQy*VRRS(8,3,",ms,",2)",\
                                    "+r1x2E*(VRRS(8,1,",m1,",2)-ZxZpE*VRRS(8,1,",ms,",2)) &"]];
             WS[StringJoin["             -r1x2E*VRR(8,3,",ms,")"]];
             (* zx,zx *)
             WS[StringJoin["VRRS(8,8,",m1,",2)=QCz*VRRS(8,2,",m1,",2)+WQz*VRRS(8,2,",ms,",2)",\
                                        "+HfxZpE*VRRS(2,2,",ms,",2)"]];
             (* zx,zy *)
             WS[StringJoin["VRRS(8,9,",m1,",2)=QCz*VRRS(8,3,",m1,",2)+WQz*VRRS(8,3,",ms,",2)",\
                                        "+HfxZpE*VRRS(2,3,",ms,",2)"]];
             (* zx,zz *)
             WS[StringJoin["VRRS(8,10,",m1,",2)=QCz*VRRS(8,4,",m1,",2)+WQz*VRRS(8,4,",ms,",2)",\
                                    "+r1x2E*(VRRS(8,1,",m1,",2)-ZxZpE*VRRS(8,1,",ms,",2)) &"]];
             WS[StringJoin["             +HfxZpE*VRRS(2,4,",ms,",2)"]];
             (* zy,xx *)
             WS[StringJoin["VRRS(9,5,",m1,",2)=QCx*VRRS(9,2,",m1,",2)+WQx*VRRS(9,2,",ms,",2)",\
                                    "+r1x2E*(VRRS(9,1,",m1,",2)-ZxZpE*VRRS(9,1,",ms,",2))"]];
             (* zy,yx *)
             WS[StringJoin["VRRS(9,6,",m1,",2)=QCy*VRRS(9,2,",m1,",2)+WQy*VRRS(9,2,",ms,",2)",\
                                        "+HfxZpE*VRRS(4,2,",ms,",2)",\
                                        "-r1x2E*VRR(9,2,",ms,")"]];
             (* zy,yy *)
             WS[StringJoin["VRRS(9,7,",m1,",2)=QCy*VRRS(9,3,",m1,",2)+WQy*VRRS(9,3,",ms,",2)",\
                                    "+r1x2E*(VRRS(9,1,",m1,",2)-ZxZpE*VRRS(9,1,",ms,",2)) &"]];
             WS[StringJoin["            +HfxZpE*VRRS(4,3,",ms,",2)",\
                                        "-r1x2E*VRR(9,3,",ms,")"]];
             (* zy,zx *)
             WS[StringJoin["VRRS(9,8,",m1,",2)=QCz*VRRS(9,2,",m1,",2)+WQz*VRRS(9,2,",ms,",2)",\
                                        "+HfxZpE*VRRS(3,2,",ms,",2)"]];
             (* zy,zy *)
             WS[StringJoin["VRRS(9,9,",m1,",2)=QCz*VRRS(9,3,",m1,",2)+WQz*VRRS(9,3,",ms,",2)",\
                                        "+HfxZpE*VRRS(3,3,",ms,",2)"]];
             (* zy,zz *)
             WS[StringJoin["VRRS(9,10,",m1,",2)=QCz*VRRS(9,4,",m1,",2)+WQz*VRRS(9,4,",ms,",2)",\
                                    "+r1x2E*(VRRS(9,1,",m1,",2)-ZxZpE*VRRS(9,1,",ms,",2)) &"]];
             WS[StringJoin["             +HfxZpE*VRRS(3,4,",ms,",2)"]];
             (* zz,xx *)
             WS[StringJoin["VRRS(10,5,",m1,",2)=QCx*VRRS(10,2,",m1,",2)+WQx*VRRS(10,2,",ms,",2)",\
                                    "+r1x2E*(VRRS(10,1,",m1,",2)-ZxZpE*VRRS(10,1,",ms,",2))"]];
             (* zz,yx *)
             WS[StringJoin["VRRS(10,6,",m1,",2)=QCy*VRRS(10,2,",m1,",2)+WQy*VRRS(10,2,",ms,",2)",\
                                         "-r1x2E*VRR(10,2,",ms,")"]];
             (* zz,yy *)
             WS[StringJoin["VRRS(10,7,",m1,",2)=QCy*VRRS(10,3,",m1,",2)+WQy*VRRS(10,3,",ms,",2)",\
                                    "+r1x2E*(VRRS(10,1,",m1,",2)-ZxZpE*VRRS(10,1,",ms,",2)) &"]];
             WS[StringJoin["            -r1x2E*VRR(10,3,",ms,")"]];
             (* zz,zx *)
             WS[StringJoin["VRRS(10,8,",m1,",2)=QCz*VRRS(10,2,",m1,",2)+WQz*VRRS(10,2,",ms,",2)",\
                                        "+2d0*HfxZpE*VRRS(4,2,",ms,",2)"]];
             (* zz,zy *)
             WS[StringJoin["VRRS(10,9,",m1,",2)=QCz*VRRS(10,3,",m1,",2)+WQz*VRRS(10,3,",ms,",2)",\
                                        "+2d0*HfxZpE*VRRS(4,3,",ms,",2)"]];
             (* zz,zz *)
             WS[StringJoin["VRRS(10,10,",m1,",2)=QCz*VRRS(10,4,",m1,",2)+WQz*VRRS(10,4,",ms,",2)",\
                                    "+r1x2E*(VRRS(10,1,",m1,",2)-ZxZpE*VRRS(10,1,",ms,",2)) &"]];
             WS[StringJoin["             +2d0*HfxZpE*VRRS(4,4,",ms,",2)"]];
             ];
             If[ixyz==3,
             (* i=z *)
             (* xx,xx *)
             WS[StringJoin["VRRS(5,5,",m1,",3)=QCx*VRRS(5,2,",m1,",3)+WQx*VRRS(5,2,",ms,",3)",\
                                    "+r1x2E*(VRRS(5,1,",m1,",3)-ZxZpE*VRRS(5,1,",ms,",3)) &"]];
             WS[StringJoin["             +2d0*HfxZpE*VRRS(2,2,",ms,",3)"]];
             (* xx,yx *)
             WS[StringJoin["VRRS(5,6,",m1,",3)=QCy*VRRS(5,2,",m1,",3)+WQy*VRRS(5,2,",ms,",3)"]];
             (* xx,yy *)
             WS[StringJoin["VRRS(5,7,",m1,",3)=QCy*VRRS(5,3,",m1,",3)+WQy*VRRS(5,3,",ms,",3)",\
                                    "+r1x2E*(VRRS(5,1,",m1,",3)-ZxZpE*VRRS(5,1,",ms,",3))"]];
             (* xx,zx *)
             WS[StringJoin["VRRS(5,8,",m1,",3)=QCz*VRRS(5,2,",m1,",3)+WQz*VRRS(5,2,",ms,",3)",\
                                         "-r1x2E*VRR(5,2,",ms,")"]];
             (* xx,zy *)
             WS[StringJoin["VRRS(5,9,",m1,",3)=QCz*VRRS(5,3,",m1,",3)+WQz*VRRS(5,3,",ms,",3)",\
                                         "-r1x2E*VRR(5,3,",ms,")"]];
             (* xx,zz *)
             WS[StringJoin["VRRS(5,10,",m1,",3)=QCz*VRRS(5,4,",m1,",3)+WQz*VRRS(5,4,",ms,",3)",\
                                    "+r1x2E*(VRRS(5,1,",m1,",3)-ZxZpE*VRRS(5,1,",ms,",3)) &"]];
             WS[StringJoin["             -r1x2E*VRR(5,4,",ms,")"]];
             (* yx,xx *)
             WS[StringJoin["VRRS(6,5,",m1,",3)=QCx*VRRS(6,2,",m1,",3)+WQx*VRRS(6,2,",ms,",3)",\
                                    "+r1x2E*(VRRS(6,1,",m1,",3)-ZxZpE*VRRS(6,1,",ms,",3)) &"]];
             WS[StringJoin["             +HfxZpE*VRRS(3,2,",ms,",3)"]]
             (* yx,yx *)
             WS[StringJoin["VRRS(6,6,",m1,",3)=QCy*VRRS(6,2,",m1,",3)+WQy*VRRS(6,2,",ms,",3)",\
                                        "+HfxZpE*VRRS(2,2,",ms,",3)"]];
             (* yx,yy *)
             WS[StringJoin["VRRS(6,7,",m1,",3)=QCy*VRRS(6,3,",m1,",3)+WQy*VRRS(6,3,",ms,",3)",\
                                    "+r1x2E*(VRRS(6,1,",m1,",3)-ZxZpE*VRRS(6,1,",ms,",3)) &"]];
             WS[StringJoin["             +HfxZpE*VRRS(2,3,",ms,",3)"]];
             (* yx,zx *)
             WS[StringJoin["VRRS(6,8,",m1,",3)=QCz*VRRS(6,2,",m1,",3)+WQz*VRRS(6,2,",ms,",3)",\
                                         "-r1x2E*VRR(6,2,",ms,")"]];
             (* yx,zy *)
             WS[StringJoin["VRRS(6,9,",m1,",3)=QCz*VRRS(6,3,",m1,",3)+WQz*VRRS(6,3,",ms,",3)",\
                                         "-r1x2E*VRR(6,3,",ms,")"]];
             (* yx,zz *)
             WS[StringJoin["VRRS(6,10,",m1,",3)=QCz*VRRS(6,4,",m1,",3)+WQz*VRRS(6,4,",ms,",3)",\
                                    "+r1x2E*(VRRS(6,1,",m1,",3)-ZxZpE*VRRS(6,1,",ms,",3)) &"]];
             WS[StringJoin["             -r1x2E*VRR(6,4,",ms,")"]];
             (* yy,xx *)
             WS[StringJoin["VRRS(7,5,",m1,",3)=QCx*VRRS(7,2,",m1,",3)+WQx*VRRS(7,2,",ms,",3)",\
                                    "+r1x2E*(VRRS(7,1,",m1,",3)-ZxZpE*VRRS(7,1,",ms,",3))"]];
             (* yy,yx *)
             WS[StringJoin["VRRS(7,6,",m1,",3)=QCy*VRRS(7,2,",m1,",3)+WQy*VRRS(7,2,",ms,",3)",\
                                   "+2d0*HfxZpE*VRRS(3,2,",ms,",3)"]];
             (* yy,yy *)
             WS[StringJoin["VRRS(7,7,",m1,",3)=QCy*VRRS(7,3,",m1,",3)+WQy*VRRS(7,3,",ms,",3)",\
                                    "+r1x2E*(VRRS(7,1,",m1,",3)-ZxZpE*VRRS(7,1,",ms,",3)) &"]];
             WS[StringJoin["             +2d0*HfxZpE*VRRS(3,3,",ms,",3)"]];
             (* yy,zx *)
             WS[StringJoin["VRRS(7,8,",m1,",3)=QCz*VRRS(7,2,",m1,",3)+WQz*VRRS(7,2,",ms,",3)",\
                                         "-r1x2E*VRR(7,2,",ms,")"]];
             (* yy,zy *)
             WS[StringJoin["VRRS(7,9,",m1,",3)=QCz*VRRS(7,3,",m1,",3)+WQz*VRRS(7,3,",ms,",3)",\
                                         "-r1x2E*VRR(7,3,",ms,")"]];
             (* yy,zz *)
             WS[StringJoin["VRRS(7,10,",m1,",3)=QCz*VRRS(7,4,",m1,",3)+WQz*VRRS(7,4,",ms,",3)",\
                                    "+r1x2E*(VRRS(7,1,",m1,",3)-ZxZpE*VRRS(7,1,",ms,",3)) &"]];
             WS[StringJoin["             -r1x2E*VRR(7,4,",ms,")"]];
             (* zx,xx *)
             WS[StringJoin["VRRS(8,5,",m1,",3)=QCx*VRRS(8,2,",m1,",3)+WQx*VRRS(8,2,",ms,",3)",\
                                    "+r1x2E*(VRRS(8,1,",m1,",3)-ZxZpE*VRRS(8,1,",ms,",3)) &"]];
             WS[StringJoin["             +HfxZpE*VRRS(4,2,",ms,",3)"]];
             (* zx,yx *)
             WS[StringJoin["VRRS(8,6,",m1,",3)=QCy*VRRS(8,2,",m1,",3)+WQy*VRRS(8,2,",ms,",3)"]];
             (* zx,yy *)
             WS[StringJoin["VRRS(8,7,",m1,",3)=QCy*VRRS(8,3,",m1,",3)+WQy*VRRS(8,3,",ms,",3)",\
                                    "+r1x2E*(VRRS(8,1,",m1,",3)-ZxZpE*VRRS(8,1,",ms,",3)) "]];
             (* zx,zx *)
             WS[StringJoin["VRRS(8,8,",m1,",3)=QCz*VRRS(8,2,",m1,",3)+WQz*VRRS(8,2,",ms,",3)",\
                                        "+HfxZpE*VRRS(2,2,",ms,",3)",\
                                         "-r1x2E*VRR(8,2,",ms,")"]];
             (* zx,zy *)
             WS[StringJoin["VRRS(8,9,",m1,",3)=QCz*VRRS(8,3,",m1,",3)+WQz*VRRS(8,3,",ms,",3)",\
                                        "+HfxZpE*VRRS(2,3,",ms,",3)",\
                                         "-r1x2E*VRR(8,3,",ms,")"]];
             (* zx,zz *)
             WS[StringJoin["VRRS(8,10,",m1,",3)=QCz*VRRS(8,4,",m1,",3)+WQz*VRRS(8,4,",ms,",3)",\
                                    "+r1x2E*(VRRS(8,1,",m1,",3)-ZxZpE*VRRS(8,1,",ms,",3)) &"]];
             WS[StringJoin["             +HfxZpE*VRRS(2,4,",ms,",3)",\
                                         "-r1x2E*VRR(8,4,",ms,")"]];
             (* zy,xx *)
             WS[StringJoin["VRRS(9,5,",m1,",3)=QCx*VRRS(9,2,",m1,",3)+WQx*VRRS(9,2,",ms,",3)",\
                                    "+r1x2E*(VRRS(9,1,",m1,",3)-ZxZpE*VRRS(9,1,",ms,",3))"]];
             (* zy,yx *)
             WS[StringJoin["VRRS(9,6,",m1,",3)=QCy*VRRS(9,2,",m1,",3)+WQy*VRRS(9,2,",ms,",3)",\
                                        "+HfxZpE*VRRS(4,2,",ms,",3)"]];
             (* zy,yy *)
             WS[StringJoin["VRRS(9,7,",m1,",3)=QCy*VRRS(9,3,",m1,",3)+WQy*VRRS(9,3,",ms,",3)",\
                                    "+r1x2E*(VRRS(9,1,",m1,",3)-ZxZpE*VRRS(9,1,",ms,",3)) &"]];
             WS[StringJoin["            +HfxZpE*VRRS(4,3,",ms,",3)"]];
             (* zy,zx *)
             WS[StringJoin["VRRS(9,8,",m1,",3)=QCz*VRRS(9,2,",m1,",3)+WQz*VRRS(9,2,",ms,",3)",\
                                        "+HfxZpE*VRRS(3,2,",ms,",3) &"]];
             WS[StringJoin["             -r1x2E*VRR(9,2,",ms,")"]];
             (* zy,zy *)
             WS[StringJoin["VRRS(9,9,",m1,",3)=QCz*VRRS(9,3,",m1,",3)+WQz*VRRS(9,3,",ms,",3)",\
                                        "+HfxZpE*VRRS(3,3,",ms,",3)",\
                                         "-r1x2E*VRR(9,3,",ms,")"]];
             (* zy,zz *)
             WS[StringJoin["VRRS(9,10,",m1,",3)=QCz*VRRS(9,4,",m1,",3)+WQz*VRRS(9,4,",ms,",3)",\
                                    "+r1x2E*(VRRS(9,1,",m1,",3)-ZxZpE*VRRS(9,1,",ms,",3)) &"]];
             WS[StringJoin["             +HfxZpE*VRRS(3,4,",ms,",3)",\
                                         "-r1x2E*VRR(9,4,",ms,")"]];
             (* zz,xx *)
             WS[StringJoin["VRRS(10,5,",m1,",3)=QCx*VRRS(10,2,",m1,",3)+WQx*VRRS(10,2,",ms,",3)",\
                                    "+r1x2E*(VRRS(10,1,",m1,",3)-ZxZpE*VRRS(10,1,",ms,",3))"]];
             (* zz,yx *)
             WS[StringJoin["VRRS(10,6,",m1,",3)=QCy*VRRS(10,2,",m1,",3)+WQy*VRRS(10,2,",ms,",3)"]];
             (* zz,yy *)
             WS[StringJoin["VRRS(10,7,",m1,",3)=QCy*VRRS(10,3,",m1,",3)+WQy*VRRS(10,3,",ms,",3)",\
                                    "+r1x2E*(VRRS(10,1,",m1,",3)-ZxZpE*VRRS(10,1,",ms,",3))"]];
             (* zz,zx *)
             WS[StringJoin["VRRS(10,8,",m1,",3)=QCz*VRRS(10,2,",m1,",3)+WQz*VRRS(10,2,",ms,",3)",\
                                        "+2d0*HfxZpE*VRRS(4,2,",ms,",3) &"]];
             WS[StringJoin["             -r1x2E*VRR(10,2,",ms,")"]];
             (* zz,zy *)
             WS[StringJoin["VRRS(10,9,",m1,",3)=QCz*VRRS(10,3,",m1,",3)+WQz*VRRS(10,3,",ms,",3)",\
                                        "+2d0*HfxZpE*VRRS(4,3,",ms,",3) &"]];
             WS[StringJoin["             -r1x2E*VRR(10,3,",ms,")"]];
             (* zz,zz *)
             WS[StringJoin["VRRS(10,10,",m1,",3)=QCz*VRRS(10,4,",m1,",3)+WQz*VRRS(10,4,",ms,",3)",\
                                    "+r1x2E*(VRRS(10,1,",m1,",3)-ZxZpE*VRRS(10,1,",ms,",3)) &"]];
             WS[StringJoin["             +2d0*HfxZpE*VRRS(4,4,",ms,",3)",\
                                         "-r1x2E*VRR(10,4,",ms,")"]];
             ];
           ];

*)

           (* Punch general MVRR call *)
           If[(iell>2&&kell==0) || (iell==0 && kell>2) || (iell>1 && kell>=1) || (iell >=1 && kell>1 ),
             WS[StringJoin["CALL ",MVRRSubName,"(",ToString[ixyz],\
                           ",",CLenBra,",",CLenKet,\
                           ",VRRS(1,1,",ToString[m-1],",",ToString[ixyz],")",\
                           ",VRRS(1,1,",ToString[m],",",ToString[ixyz],")",\
                           ",",CLenBra1,",",CLenKet1,",VRR(1,1,",ToString[m],"))"]];
           ];

        ,{m,BraEll+KetEll-(iell+kell)+1,1,-1}];
      ];

   ,{iell,0,BraEll}];
   ,{kell,0,KetEll}];

   (* Punch ENDIF-loop *)
   WS["ENDIF"];

   (*End loop over ixyz=x,y,z*)
   ,{ixyz,1,3}];

   (* Do we want Stress? *)
   ];
   (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)

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
        WS[StringJoin["      CALL ",ContractName,"(VRR,HRR,Alpha,HRRA,Beta,HRRB,Gamma,HRRC, &"]];
        WS[StringJoin["               & VRRS,HRRS(1,1,1,1),PQJ(1),PBC%AutoW%I(1))"]];
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
        WS[StringJoin["SUBROUTINE ",ContractName,"(VRR,HRR,Alpha,HRRA,Beta,HRRB,Gamma,HRRC,VRRS,HRRS,PQJ,IW)"]];
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
        WS[StringJoin["  REAL(DOUBLE)  :: VRRS(",ToString[LEnd[LBra]],",",ToString[LEnd[LKet]],",0:",ToString[LBra+LKet],",3)"]];
        WS[StringJoin["  REAL(DOUBLE)  :: HRRS(",ToString[BraMax],",",ToString[KetMax],",",ToString[LEnd[lmax]],",9),PQJ(3)"]];
        WS[StringJoin["  INTEGER :: IJ,J,I,IW(3)"]];
      ];
      (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)

      If[LKet>=2,
         WS[StringJoin["  DO K=1,",ToString[LEnd[LKet]]]];
         Do[Do[
                BK00=StringJoin[ToString[bra],",K,"];
                If[LB<=LBra,
                   WS[StringJoin["     HRR(",BK00,"1)=HRR(",BK00,"1)+VRR(",BK00,"0)"]];
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
(* ------------------------------------------------------------------ *)
(* ------------------------------------------------------------------ *)
(* ------------------------------------------------------------------ *)
   (*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
   (* Do we want Stress? *)   
   If[DoStress==0,
      If[LKet>=2,
         WS["  DO J=1,3"];
         WS["  DO I=1,3"];
         WS["  IF(IW(J).EQ.1.AND.IW(I).EQ.1) THEN"];
         WS["  IJ=3*(J-1)+I"];
         WS[StringJoin["  DO K=1,",ToString[LEnd[LKet]]]];
         Do[Do[
                BK00=StringJoin[ToString[bra],",K,"];
                If[LB<=LBra,
                  WS[StringJoin["     HRRS(",BK00,"1,IJ)=HRRS(",BK00,"1,IJ)+PQJ(J)*VRRS(",BK00,"0,I)"]];
                ];
              ,{bra,LBegin[LB],LEnd[LB]}],{LB,0,LBra+1}]; 
           WS["  ENDDO !K"];
           WS["  ENDIF"];
           WS["  ENDDO !I"];
           WS["  ENDDO !J"];
         ,
         WS["  IJ=1"];
         WS["  DO J=1,3"];
         WS["  IF(IW(J).EQ.1) THEN"];
         Do[
           WS[StringJoin["    IF(IW(",ToString[ixyz],").EQ.1) THEN"]];
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
                     If[ixyz==1,WS[StringJoin["    HRRS(",BK00,"1,IJ)=HRRS(",BK00,"1,IJ)+PQJ(J)*VRRS(",BK00,"0,1)"]];];
                     If[ixyz==2,WS[StringJoin["    HRRS(",BK00,"1,IJ)=HRRS(",BK00,"1,IJ)+PQJ(J)*VRRS(",BK00,"0,2)"]];];
                     If[ixyz==3,WS[StringJoin["    HRRS(",BK00,"1,IJ)=HRRS(",BK00,"1,IJ)+PQJ(J)*VRRS(",BK00,"0,3)"]];];
                   ];
                 ]; 
               ,{ket,LBegin[LK],LEnd[LK]}],{LK,0,LKet+1}]
             ,{bra,LBegin[LB],LEnd[LB]}],{LB,0,LBra+1}]; 
             WS["    ENDIF"];
             WS["    IJ=IJ+1"];
           ,{ixyz,1,3}];
           WS["  ELSE"];
           WS["    IJ=IJ+3"];
           WS["  ENDIF"];
           WS["  ENDDO !J"];
         ];
   (* Do we want Stress? *)   
   ];
   (*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)

   WS[StringJoin["END SUBROUTINE ",ContractName]];
];

