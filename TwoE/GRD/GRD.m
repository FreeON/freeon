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

(*Classes = { {0,0},{0,1},{1,1},{2,2} };*)
Classes = { {0,0},{0,1},{1,1},{2,2},{3,3} };

(*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
(* DoStres, Yes=0, No!=0 *)
DoStress=1;
(*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)

CType[1]  = "s";
CType[2]  = "sp";
CType[3]  = "p";
CType[4]  = "spd";
CType[6]  = "d";
CType[7]  = "spdf";
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
,{l,0,2 MxEll-m-n}];
,{m,0,2 MxEll-n}];
,{n,0,2 MxEll}];


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
SetAttributes[HRR,NHoldAll];
SetAttributes[HRRA,NHoldAll];
SetAttributes[HRRB,NHoldAll];
SetAttributes[HRRC,NHoldAll];
SetAttributes[GRADIENT,NHoldAll];

(*> STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)
If[DoStress==0,
  SetAttributes[STRESS,NHoldAll];
  SetAttributes[PQIJ,NHoldAll];
  SetAttributes[FP,NHoldAll];
];
(*< STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS *)

Get["GAMMAS.m"];
Get["CNTRCTG.m"];
Get["BRAHRR.m"];
Get["KETHRR.m"];
Get["FRONTM.m"];

MakeList={};
RelsList={};

PunchHRRBra;

Do[Do[Do[Do[

   If[IntegralClass[Classes[[ic]]]>=IntegralClass[Classes[[jc]]]&& \
      IntegralClass[Classes[[kc]]]>=IntegralClass[Classes[[lc]]]&& \
      IntegralClass[Classes[[ic]]]*100+IntegralClass[Classes[[jc]]]>= \
      IntegralClass[Classes[[kc]]]*100+IntegralClass[Classes[[lc]]],


            CommentLine=StringJoin["(",CType[IntegralClass[Classes[[ic]]]]," ", \
                                       CType[IntegralClass[Classes[[jc]]]],"|", \
                                       CType[IntegralClass[Classes[[kc]]]]," ", \
                                       CType[IntegralClass[Classes[[lc]]]],")"];
	    CommentLine=StringJoin["! ---------------------------------------------------------- \n",\
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


             ijklType=ijklFlag;

	   Subroutine=StringJoin["dIntB",ToString[ijklType],".F90"];
	   OpenWrite[Subroutine];
	   Print[" Openned ",Subroutine];
	   WS[String_]:=WriteString[Subroutine,String,"\n"];
	   WriteString[Subroutine,CommentLine]; 
	   MakeList=Append[MakeList,StringJoin["dIntB",ToString[ijklType],".o"]];
           RelsList=Append[RelsList,StringJoin["dIntB",ToString[ijklType],".x"]]; 
           BraEll=imax+jmax;
           KetEll=kmax+lmax;
           PunchFront[Subroutine,ic,jc,kc,lc,ijklType]; 
           PunchVRRCalls[Subroutine,BraEll,KetEll];
           PunchContractCalls[Subroutine,ic,jc,kc,lc];
	   WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];
           WS["   ENDDO ! (M0| loop"];            					 
           WS["ENDDO ! |N0) loop"];            					 
           PunchHRRKetCalls[Subroutine,ic,jc,kc,lc]; 
	   PunchHRRBraCalls[Subroutine,ic,jc,kc,lc]; 
	   WS[String_]:=WriteString[Subroutine,"   ",String,"\n"];
           WS[StringJoin["END SUBROUTINE dIntB",ToString[ijklType]]];
           PunchVRRContract[Subroutine,ic,jc,kc,lc,BraEll,KetEll];

           Close[Subroutine];
           Print[" Closed ",Subroutine];
      ]; 

,{ic,1,LC}]
,{jc,1,LC}]
,{kc,1,LC}]
,{lc,1,LC}];

TAB=FromCharacterCode[9];

MakeList=Union[MakeList];
RelsList=Union[RelsList];
MakeString="";
RelsString="";
LM=Length[MakeList];
LR=Length[RelsList];
Do[MakeString=StringJoin[MakeString,MakeList[[i]],"\\\n"],{i,LM-1}];
MakeString=StringJoin[MakeString,MakeList[[LM]],"\n"];
Do[RelsString=StringJoin[RelsString,RelsList[[i]],"\\\n"],{i,LR-1}];
RelsString=StringJoin[RelsString,RelsList[[LR]],"\n"];
MakeList=MakeString;
RelsList=RelsString;

MakeList=StringJoin["IntObjs=",MakeList];
RelsList=StringJoin["ReleaseFiles=",RelsList]; 

(**************** Print out the Makefile ************************)

Makefile="Makefile";

OpenWrite[Makefile];

WriteString[Makefile,"include $(MONDO_COMPILER)\n"];
WriteString[Makefile,"include $(MONDO_HOME)/Includes/Suffixes\n"];
WriteString[Makefile,"include $(MONDO_HOME)/Includes/RemoveAll\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"CPPMISC =\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"EXTRA_INCLUDES=-I$(MONDO_HOME)/TwoE/MSC/\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"GRDLibs=$(F90_LIBRARIES) -lGRD $(LAPACK)\n"] 
WriteString[Makefile,"#\n"];
WriteString[Makefile,MakeList];
WriteString[Makefile,"#\n"];
WriteString[Makefile,RelsList];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"all:    GRD\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"clean:  CGRD\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"release: clean $(ReleaseFiles)\n",TAB,"rm -f *.F90\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"purge:clean\n",
                            TAB,"rm -f $(MONDO_LIB)/libGRD.a\n",
                            TAB,"rm -f $(REMOVESRCE)\n"];
WriteString[Makefile,"ppurge:purge\n",
                             TAB,"rm -f *.Inc dIntB*.F90 Ket*.F90 Bra*.F90 VRR*.F90 \n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"source:\n",
                             TAB,"math<GRD.m>GRD.out\n",
                             TAB,"math<DERIInterface.m>GRD.out\n",
                             TAB,"rm GRD.out\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"GRD:$(MiscObjs) $(IntObjs)\n",
                          TAB,"$(AR) $(ARFLAGS) $(MONDO_LIB)/libGRD.a $(?:.F90=.o)\n",
                          TAB,"$(RANLIB) $(MONDO_LIB)/libGRD.a\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"CGRD:\n",
                           TAB,"rm -f $(REMOVEMISC) $(REMOVEMODS)\n",
                           TAB,"rm -f \#*\n",
                           TAB,"rm -f *~\n",
                           TAB,"ln -s /dev/null core\n",
                           TAB,"ls -l\n"];
Close[Makefile];


