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
			(*			Return[1]; *)
                        L=LMN[[1]];
                        M=LMN[[2]];
                        N=LMN[[3]];
                        Fct[n_]:=(2 n-1)!!;
                        X=Fct[L+M+N];
                        Y=Fct[L]*Fct[M]*Fct[N];
			Return[Sqrt[X/Y]]];

Classes = { {0,0},{0,1},{1,1},{2,2},{3,3},{4,4}};

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

Do[Do[Do[ 
         lmn=LMNDex[l, m, n];
         lx[lmn] = l;
         my[lmn] = m;
         nz[lmn] = n;
,{l,0,2 MxEll-m-n}]
,{m,0,2 MxEll-n}]
,{n,0,2 MxEll}];


(*==================== LOAD OPTIMIZING AND FORMATING ROUTINES FOR MMA AND SET ASSOCIATED OPTIONS====================== *)

Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];
Get[StringJoin[MondoHome,"/MMA/Format.m"]];
Get[StringJoin[MondoHome,"/MMA/Optimize.m"]];

FF[x_] := ToString[FixedNumberForm[SetPrecision[N[x,32],32], 16, 2]];

SetOptions[FortranAssign,AssignOptimize->True,AssignMaxSize->5000,
           AssignBreak->{5000,"         "},AssignIndent->"      ",AssignTemporary->{W,Sequence}];
SetOptions[Optimize,OptimizeVariable->{V,Array},OptimizeTimes->True,OptimizePlus->True,
           OptimizeCoefficients->True,OptimizeFunction->False]; 
SetOptions[OpenWrite, PageWidth -> 200];

SetAttributes[o,NHoldAll];
SetAttributes[AuxR,NHoldAll];
SetAttributes[VRR0,NHoldAll];
SetAttributes[VRR1,NHoldAll];


(*===========================--  THE OBARA-SAIKA VERTICAL RECURENCE RELATIONS  --=================================*)

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

PunchVRRFront[Subroutine_,Subname_]:=Block[{BKType,LenBra,LenKet,WS},
WriteString[Subroutine,StringJoin["   SUBROUTINE ",Subname,"(LB,LK,VRR0,VRR1) \n"]];
WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];
WS["USE DerivedTypes"];
WS["USE VScratchB"];
WS["USE GlobalScalars"];
WS["IMPLICIT REAL(DOUBLE) (W)"]; 
WS["INTEGER :: LB,LK"];
WS[StringJoin["REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1"]];
];

(*  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  *)

PunchVRRClass[BraEll_,KetEll_]:=Module[{oList,IList,Kount,a,c,CollectList},

spaces="            ";	 
TStart=TimeUsed[];
TStartTotal=TStart;
MaxSameEll=Min[BraEll,KetEll];
Do[Do[
If[iell+kell>0,
   TStart=TimeUsed[];
   as=Max[0,iell-1];
   cs=Max[0,kell-1];   
   TStart=TimeUsed[];oList={" "->""};IList={};Kount = 0;
   CollectList=Flatten[Table[{VRR0[i,j],VRR1[i,j]},{i,LEnd[BraEll]},{j,LEnd[KetEll]}]]; 
   Do[Do[
         Kount=Kount+1;
         a={lx[i], my[i], nz[i]};
         c={lx[k], my[k], nz[k]};
	 (*
         vrr=VRR[a,c,0, kell+iell-1];
*)
         vrr=Expand[VRR[a,c,0, kell+iell-1]];
         vrrList=Coefficient[vrr,CollectList];
         vrr=Factor[Sum[vrrList[[i]] CollectList[[i]],{i,Length[CollectList]}]];
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

MakeList={};
RelsList={};

PunchVRRClass[MxEll*2,MxEll*2];

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

MakeList=StringJoin["VRRObjs=",MakeList];
RelsList=StringJoin["ReleaseFiles=",RelsList]; 

(**************** Print out the Makefile ************************)

Makefile="Makefile"

OpenWrite[Makefile];

WriteString[Makefile,"include $(MONDO_COMPILER)\n"];
WriteString[Makefile,"include $(MONDO_HOME)/Includes/Suffixes\n"];
WriteString[Makefile,"include $(MONDO_HOME)/Includes/RemoveAll\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"CPPMISC =\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"EXTRA_INCLUDES=-I$(MONDO_HOME)/TwoE/MSC/\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"VRRLibs=$(F90_LIBRARIES) -lVRR $(LAPACK)\n"] 
WriteString[Makefile,"#\n"];
WriteString[Makefile,MakeList];
WriteString[Makefile,"#\n"];
WriteString[Makefile,RelsList];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"all:    VRR\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"clean:  CVRR\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"release: clean $(ReleaseFiles)\n",TAB,"rm -f *.F90 \n",];
WriteString[Makefile,"purge: clean\n",
                            TAB,"rm -f $(MONDO_LIB)/libVRR.a\n",
                            TAB,"rm -f $(REMOVESRCE)\n"];
WriteString[Makefile,"ppurge: purge\n",
                              TAB,"rm -f VRR*.F90 \n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"source:\n",
                             TAB,"math<VRR.m>VRR.out\n",
                             TAB,"rm VRR.out\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"VRR: $(VRRObjs)\n",
                           TAB,"$(AR) $(ARFLAGS) $(MONDO_LIB)/libVRR.a $(?:.F90=.o)\n",
                           TAB,"$(RANLIB) $(MONDO_LIB)/libVRR.a\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"CVRR:\n",
                           TAB,"rm -f $(REMOVEMISC) $(REMOVEMODS)\n",
                           TAB,"rm -f \#*\n",
                           TAB,"rm -f *~\n",
                           TAB,"ln -s /dev/null core\n",
                           TAB,"ls -l\n"];
Close[Makefile];


