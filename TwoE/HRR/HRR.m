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

Classes = {{0,0},{0,1},{0,2},{1,1},{2,2},{3,3},{4,4}}; 

CType[1]  = "s";
CType[2]  = "sp";
CType[3]  = "p";
CType[4]  = "spd";
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

SetOptions[FortranAssign,AssignOptimize->False,AssignMaxSize->5000,
           AssignBreak->{5000,"         "},AssignIndent->"      ",AssignTemporary->{W,Sequence}];
SetOptions[Optimize,OptimizeVariable->{V,Array},OptimizeTimes->True,OptimizePlus->True,
           OptimizeCoefficients->True,OptimizeFunction->False]; 
SetOptions[OpenWrite, PageWidth -> 200];

SetAttributes[o,NHoldAll];
SetAttributes[AuxR,NHoldAll];
SetAttributes[HRR,NHoldAll];

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
 oList={"+"->"+  & \n                        "};
 Kount = 0;
 Do[Do[
          Do[Do[
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
                   Kount = Kount + 1;
                   IList=Append[IList,ToExpression[StringJoin["KK",ToString[k],"xx",ToString[l],"BB"]]];
                   oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->"!"];
		   hrr=hrr/.{HRR[i_,k_,l_]:>HRR[abc,k,l]};
                   Kount = Kount + 1;
                   IList=Append[IList,Horner[hrr]];
                   HRRAddress=StringJoin["1:LB,",ToString[k],",",ToString[l]];
                   oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["HRR(",HRRAddress,")"]];
                  ]; 
 ,{k,LBegin[kl],LEnd[kl]}]
 ,{l,LBegin[ll],LEnd[ll]}]
 ,{kl,kmax,0,-1}]
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


(**************** CREATE THE SUBROUTINES  ************************)
MakeList={};
RelsList={};

PunchHRRKet;
(**************** DONE CREATING THE SUBROUTINES  ************************)

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

MakeList=StringJoin["HRRObjs=",MakeList];
RelsList=StringJoin["HRRRels=",RelsList]; 

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
WriteString[Makefile,"HRRLibs=$(F90_LIBRARIES) -lHRR $(LAPACK)\n"] 
WriteString[Makefile,"#\n"];
WriteString[Makefile,MakeList];
WriteString[Makefile,"#\n"];
WriteString[Makefile,RelsList];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"all:    HRR\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"clean:  CHRR\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"purge: clean\n",
                            TAB,"rm -f $(MONDO_LIB)/libHRR.a\n",
                            TAB,"rm -f $(REMOVESRCE)\n"];
WriteString[Makefile,"ppurge: purge\n",
                              TAB,"rm -f KetHRR*.F90 \n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"source:\n",
                             TAB,"math<HRR.m>HRR.out\n",
                             TAB,"rm HRR.out\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"HRR: $(HRRObjs)\n",
                           TAB,"$(AR) $(ARFLAGS) $(MONDO_LIB)/libHRR.a $(?:.F90=.o)\n",
                           TAB,"$(RANLIB) $(MONDO_LIB)/libHRR.a\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"CHRR:\n",
                           TAB,"rm -f $(REMOVEMISC) $(REMOVEMODS)\n",
                           TAB,"rm -f \#*\n",
                           TAB,"rm -f *~\n",
                           TAB,"ln -s /dev/null core\n",
                           TAB,"ls -l\n"];
Close[Makefile];


