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
!    Author: Matt Challacombe
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

Classes = { {1,1}};
CType[1] = "S";
CType[2] = "SP";
CType[3] = "P";
CType[6] = "D";

LC=Length[Classes];
MxEll = Max[Classes];

Do[Do[Do[ 
         lmn=LMNDex[l, m, n];
         lx[lmn] = l;
         my[lmn] = m;
         nz[lmn] = n;
,{l,0,MxEll-m-n}]
,{m,0,MxEll-n}]
,{n,0,MxEll}];

Do[ Print[CType[IntegralClass[Classes[[i]]]]], {i, 1, Length[Classes]}];

SetAttributes[X,NHoldAll];

(**************** THE HRR RELATIONS: THX HGP! ************************)

KetHRR[a_List,b_List,c_List,d_List]:=Module[{p, CD, c1, d1,adex,cdex,c1dex}, 
                                            p=Position[d, Max[d]][[1, 1]];
                                            If[p == 1, one = {1, 0, 0}; CD = x1; ];
                                            If[p == 2, one = {0, 1, 0}; CD = y1; ];
                                            If[p == 3, one = {0, 0, 1}; CD = z1; ];
                                            c1=c+one;
                                            d1=d-one;                
                                            If[Max[d1] == 0,                                       
					       adex=LMNDex[a[[1]],a[[2]],a[[3]]];
					       cdex=LMNDex[c[[1]],c[[2]],c[[3]]];
					       c1dex=LMNDex[c1[[1]],c1[[2]],c1[[3]]];
					       v1=ToExpression[StringJoin["MBarNu",ToString[adex],"v",ToString[c1dex],"w"]];
                                               v2=ToExpression[StringJoin["MBarNu",ToString[adex],"v",ToString[cdex],"w"]];
                                               v1+CD*v2 
                                              ,
                                               Return[KetHRR[a,b,c1,d1]+CD KetHRR[a,b,c,d1]]
                                               ]
                                           ]

HRR[a_List,b_List,c_List,d_List]:=Module[{p,AB,a1,b1}, 
                                         p=Position[b,Max[b]][[1, 1]];
                                         If[p == 1, one = {1, 0, 0}; AB = x2; ];
                                         If[p == 2, one = {0, 1, 0}; AB = y2; ];
                                         If[p == 3, one = {0, 0, 1}; AB = z2; ];
                                         a1 = a + one;
                                         b1 = b - one;                
                                         If[Max[b1] == 0, 
                                            Return[ KetHRR[a1, b1, c, d] + AB  KetHRR[a, b1, c, d]],
                                            Return[HRR[a1, b1, c, d] + AB  HRR[a, b1, c, d]]
                                           ] 
                                        ];


(* LOAD OPTIMIZING AND FORMATING ROUTINES FOR MMA AND SET ASSOCIATED OPTIONS *)

Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];
Get[StringJoin[MondoHome,"/MMA/Format.m"]];
Get[StringJoin[MondoHome,"/MMA/Optimize.m"]];

SetOptions[Optimize,OptimizeVariable->{V,Array},OptimizeTimes->True];
SetOptions[FortranAssign,AssignOptimize->True,AssignMaxSize->400,AssignBreak->{132," & \n          "},AssignTemporary->{W,Array}];
SetOptions[OpenWrite, PageWidth -> 200];

(* PUT THE TRANSFORMATIONS TO FILE *)

PunchHRRClass[FileName_,ic_,jc_,kc_,lc_]:=Module[{oList,IList,Kount,a,b,c,d},
						 imin = Classes[[ic, 1]]; imax = Classes[[ic, 2]];
						 jmin = Classes[[jc, 1]]; jmax = Classes[[jc, 2]];
						 kmin = Classes[[kc, 1]]; kmax = Classes[[kc, 2]];
						 lmin = Classes[[lc, 1]]; lmax = Classes[[lc, 2]];
						 oList={" "->"","u"->"(","v"->",","w"->")","x1"->"CDx","y1"->"CDy","z1"->"CDz","x2"->"ABx","y2"->"ABy","z2"->"ABz"};
						 IList={};
						 Kount = 0;
                                                 Do[Do[Do[Do[
                                                 Do[Do[Do[Do[
                                                             Kount = Kount + 1;
                                                             a = {lx[i], my[i], nz[i]};
                                                             b = {lx[j], my[j], nz[j]};
                                                             c = {lx[k], my[k], nz[k]};
                                                             d = {lx[l], my[l], nz[l]};
                                                             IList=Append[IList,HRR[a,b,c,d]];
                                                             oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->StringJoin["I(",ToString[i],",",ToString[j],",",ToString[k],",",ToString[l],")"]];
                                                ,{i,LBegin[il],LEnd[il]}]
                                                ,{j,LBegin[jl],LEnd[jl]}]
                                                ,{k,LBegin[kl],LEnd[kl]}]
                                                ,{l,LBegin[ll],LEnd[ll]}]
                                                ,{il,imin,imax}]
                                                ,{jl,jmin,jmax}]
                                                ,{kl,kmin,kmax}]
                                                ,{ll,lmin,lmax}];
                                                Write[FileName,FortranAssign[o,IList,AssignReplace->oList]];
                                               ];

Do[Do[Do[Do[
	    imin = Classes[[ic, 1]]; imax = Classes[[ic, 2]];
	    jmin = Classes[[jc, 1]]; jmax = Classes[[jc, 2]];
	    kmin = Classes[[kc, 1]]; kmax = Classes[[kc, 2]];
	    lmin = Classes[[lc, 1]]; lmax = Classes[[lc, 2]];
            ijklType=1000*IntegralClass[Classes[[ic]]]
                     +100*IntegralClass[Classes[[jc]]]
                      +10*IntegralClass[Classes[[kc]]]
	                 +IntegralClass[Classes[[lc]]];
	   Subroutine=StringJoin["HRR",ToString[ijklType],".F90"];
	   WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];
	   Print[" Openned ",Subroutine];
	   OpenWrite[Subroutine];
	   WriteString[Subroutine,"      SUBROUTINE HRR",ToString[ijklType],"(AOff,BOff,COff,DOff,ABx,ABy,ABz,CDx,CDy,CDz) \n"];
	   WS["  USE TwoEGlobals  ! Avoiding copy in copy out since 1990"];
	   WS["  INTEGER :: AOff,BOff,COff,DOff"];
	   WS["  REAL(DOUBLE) :: ABx,ABy,ABz,CDx,CDy,CDz "];
           PunchHRRClass[Subroutine,ic,jc,kc,lc];
	   (*
	   Kab=Pcd (ac|bd} 
	   Kac=Pbd (ab|cd} 
	    *)
	   WS["! We will just do the re-indexing in a cheezy way, untill we get it right, then hardwire into explicit code "];
	   WS["! K(1:ABlk,1:CBlk)=P(1:BBlk,1:DBlk)*ERI(1:DBlk,1:BBlk,1:ABlk,1:CBlk)--Kac=Pbd (ab|cd)"];
	   WS[                    "         ACount=AOff"];  
	   WriteString[Subroutine,"         DO IA=",ToString[LBegin[imin]],",",ToString[LEnd[imax]]," \n"];
	   WS[                    "            ACount=ACount+1"];  
	   WS[                    "            BCount=BOff"];  
	   WriteString[Subroutine,"            DO IB=",ToString[LBegin[jmin]],",",ToString[LEnd[jmax]]," \n"];
	   WS[                    "               BCount=BCount+1"];  
	   WS[                    "               BCount=COff"];  
	   WriteString[Subroutine,"               DO IC=",ToString[LBegin[kmin]],",",ToString[LEnd[kmax]]," \n"];
	   WS[                    "                  CCount=CCount+1"];  
	   WS[                    "                  DCount=DOff"];  
	   WriteString[Subroutine,"                  DO ID=",ToString[LBegin[lmin]],",",ToString[LEnd[lmax]]," \n"];
	   WS[                    "                     DCount=DCount+1"];  
           WS[                    "                     ERI(DCount,BCount,ACount,CCount)=I(IA,IB,IC,ID)  "];
           WS["          ENDDO "];
           WS["       ENDDDO "];
           WS["    ENDDO "];
           WS["  ENDDO 	 "];

           WS[StringJoin["  END SUBROUTINE HRR",ToString[ijklType]]];
           Close[Subroutine];
           Print[" Closed ",Subroutine];
,{ic,1,LC}]
,{jc,1,LC}]
,{kc,1,LC}]
,{lc,1,LC}];

