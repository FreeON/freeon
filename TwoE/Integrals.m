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

Classes = {{0,0},{1,1} };
CType[1] = "S";
CType[2] = "SP";
CType[3] = "P";
CType[6] = "D";

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

VRR[a_List,c_List,m_]:=Module[{p,q,PA,QC,WP,WQ,one,two,a1,a2,c1,c2,Ai1,Ci1,CiO2zn,AiO2zn}, 
                              p=Position[a,Max[a]][[1, 1]];
                              q=Position[c,Max[c]][[1, 1]];
                              If[a[[p]]<0 || c[[p]]< 0,Return[0]];                                        
                              If[a[[p]]==0 && c[[p]]==0,Return[SSSS[m]]];                                        
			      If[a[[p]]>c[[q]],
			        (* Taking down ell on bra side *)
			         If[p == 1, one = {1, 0, 0}; two={2,0,0}; PA=PAx; WP=WPx; ];
			         If[p == 2, one = {0, 1, 0}; two={0,2,0}; PA=PAy; WP=WPy; ];
			         If[p == 3, one = {0, 0, 1}; two={0,0,2}; PA=PAz; WP=WPz; ];
				 Ai1=(a[[p]]-1)*r1x2Z;  (* (a[i]-1)/(2 zeta)   *)
				 CiO2zn=c[[p]]*HfxZpE; (* c[i]/(2(zeta+eta)) *)
                                 a1 = a - one;
                                 a2 = a - two;
			         c1 = c - one;                                  
                                 Return[PA*VRR[a1,c,m]+WP*VRR[a1,c,m+1]+Ai1*(VRR[a2,c,m]  \
                                        -ExZpE*VRR[a2,c,m+1]) +CiO2zn*VRR[a1,c1,m+1]]
			       ,(* Taking down ell on ket side *)
			         If[q == 1, one = {1, 0, 0}; two={2,0,0}; QC=QCx; WQ=WQx;];
			         If[q == 2, one = {0, 1, 0}; two={0,2,0}; QC=QCy; WQ=WQy;];
			         If[q == 3, one = {0, 0, 1}; two={0,0,2}; QC=QCz; WQ=WQz;];
                                 Ci1=(c[[q]]-1)*r1x2E ; (* (c[i]-1)/(2 eta)   *)
			         AiO2zn=a[[q]]*HfxZpE; (* a[i]/(2(zeta+eta)) *)
                                 c1 = c - one;
                                 c2 = c - two;
                                 a1 = a - one;  
                                 Return[QC*VRR[a,c1,m]+WQ*VRR[a,c1,m+1]+Ci1*(VRR[a,c2,m] \
                                        -ZxZpE*VRR[a,c2,m+1])+AiO2zn*VRR[a1,c1,m+1]]
                               ];
                             ];

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

FF[x_] := ToString[FixedNumberForm[SetPrecision[N[x,32],32], 16, 2]];

SetOptions[Optimize,OptimizeVariable->{V,Sequence},OptimizeTimes->True];
SetOptions[OpenWrite, PageWidth -> 200];

SetAttributes[o,NHoldAll];
SetAttributes[SSSS,NHoldAll];
SetAttributes[X,NHoldAll];

(* PUT THE TRANSFORMATIONS TO FILE *)

PunchVRRClass[FileName_,BraEll_,KetEll_]:=Module[{oList,IList,Kount,a,c},
						 oList={" "->""};
						 IList={};
                                                 sList={};
						 Kount = 0;
                                                 Do[Do[
                                                       Kount = Kount + 1;
                                                       a = {lx[i], my[i], nz[i]};
                                                       c = {lx[k], my[k], nz[k]};
                                                       IList=Append[IList,VRR[a,c,0]+o[Kount]];
						       MBarString=StringJoin["MBarN(",ToString[i],",",ToString[k],")"];
                                                       oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->MBarString];
                                                       sList=Append[sList,StringJoin["SSSS(",ToString[i-1+k-1],")"]-> \
                                                                          StringJoin["AuxR(",ToString[i-1+k-1],")"]];
                                                 ,{i,1,LEnd[BraEll]}],{k,1,LEnd[KetEll]}];
						 oList=Join[oList,sList];
                                                 Write[FileName,FortranAssign[o,IList,AssignReplace->oList]];
];                                                ;



PunchGammas[Subroutine_,LTot_]:=Block[{WS,FStr,Gammas},

WS[String_]:=WriteString[Subroutine,"            ",String,"\n"];

FStr[L_]:=Module[{LSt},
                 LSt=ToString[L];
                 StringJoin["(F",LSt,"_0(L)+T*(F",LSt,"_1(L)+T*(F",LSt,"_2(L)+T*(F",LSt,"_3(L)+T*F",LSt,"_4(L)))))"]];


Gammas[EllTot_]:=Module[{LSt,LSt2,GSt,GSt1,GFlt,GammAss},


              GammAss[n_] := FF[ Abs[2 n - 1]!!/(2 (2 )^n) Sqrt[Pi] ];


                                WS["IF(T<Gamma_Switch)THEN"];
                                WS["  L=AINT(T*Gamma_Grid)"];			
         If[EllTot==1,
	    Print[" EllTot= ",EllTot];
                     WS[StringJoin["   AuxR(0)=",FStr[0]]];
                     WS[StringJoin["   AuxR(1)=",FStr[1]]];

            ];
          If[EllTot>1,
             LSt=ToString[EllTot];

                     WS["  V1=Upq"];
                     WS["  V2=-Two*Omega"];
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
          WS[StringJoin["  AuxR(",GSt,")=V1*W",GSt,""]];
                            If[L<LTot,WS[StringJoin["  V1=V2*V1"]]];
                           ,{L,0,EllTot}];
            ]; 

                     WS["ELSE"];
                     WS["  InvT=One/T"];
                     WS["  SqInvT=Upq*DSQRT(InvT)"];
                         Do[GSt=ToString[L];

           WS[StringJoin["  AuxR(",GSt,")=",GammAss[L],"*SqInvT"]];

If[L<LTot,WS[StringJoin["  SqInvT=SqInvT*InvT"]]];
                           ,{L,0,EllTot}];

                      WS["ENDIF"];

      ];

Gammas[LTot];
];

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


PunchFront[Subroutine_,IMax_,JMax_,KMax_,LMax_,IJKL_]:=Block[{WS,LBra,LKet,BKType,LenBra,LenKet},
           LBra=IMax+JMax;
           LKet=KMax+LMax;
	   BKType=100*LBra+LKet;					   

(*

 *)
           WriteString[Subroutine,StringJoin["   SUBROUTINE Int",ToString[IJKL],"(Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz, &  \n",  
                                                           "                      Dx,Dy,Dz,ShlPrAC2,ShlPrBD2,I) \n"]];

	   WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];

	   WS["USE DerivedTypes"];
	   WS["USE GlobalScalars"];
	   If[LBra+LKet==1,
              WS["USE GammaF0"];
              WS["USE GammaF1"];,
              WS[StringJoin["USE GammaF",ToString[LBra+LKet]]]];

	   WS["USE ShellPairStruct"];

           WS["IMPLICIT REAL(DOUBLE) (V,W)"];

	   WS["TYPE(ShellPair), POINTER :: ShlPrAC2,ShlPrBD2 "];
	   LenBra=LEnd[LBra];
           LenKet=LEnd[LKet];
           WS[StringJoin["REAL(DOUBLE),DIMENSION(",ToString[LBra+LKet],") :: AuxR"]];
           WS[StringJoin["REAL(DOUBLE),DIMENSION(",ToString[LenBra],",",ToString[LenKet],") :: MBarN=0D0"]];


           WS[StringJoin["REAL(DOUBLE),DIMENSION(",ToString[LEnd[IMax]],",", \
                                                   ToString[LEnd[JMax]],",",
                                                   ToString[LEnd[KMax]],",",
                                                   ToString[LEnd[LMax]],") :: I"]];

           WS["REAL(DOUBLE)  :: Zeta,Eta,r1xZpE,HfxZpE,r1x2E,r1x2Z,ExZpE,ZxZpE,Omega,Up,Uq,Upq"];
           WS["REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz,Wx,Wy,Wz"];
           WS["REAL(DOUBLE)  :: QCx,QCy,QCz,PAx,PAy,PAz,PQx,PQy,PQz,WPx,WPy,WPz,WQx,WQy,WQz   "];
           WS["REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT"];

           WS["INTEGER       :: J,K,L"];

           WS["DO J=1,ShlPrBD2%L ! K^2 VRR |N0) loop "];

           WS["   Eta=ShlPrBD2%SP(1,J)"];
           WS["   Qx =ShlPrBD2%SP(2,J)"];
           WS["   Qy =ShlPrBD2%SP(3,J)"];
           WS["   Qz =ShlPrBD2%SP(4,J)"];
           WS["   Uq =ShlPrBD2%SP(5,J)"];

           WS["   QCx=Qx-Cx"];
           WS["   QCy=Qy-Cy"];
           WS["   QCz=Qz-Cz"];


           WS["   DO K=1,ShlPrAC2%L ! K^2 VRR (M0| loop "];

           WS["      Zeta=ShlPrAC2%SP(1,K)"];
           WS["      Px  =ShlPrAC2%SP(2,K)"];
           WS["      Py  =ShlPrAC2%SP(3,K)"];
           WS["      Pz  =ShlPrAC2%SP(4,K)"];
           WS["      Up  =ShlPrAC2%SP(5,K)"];

           WS["      r1xZpE=One/(Zeta+Eta)"];
	   WS["      Upq=SQRT(r1xZpE)*Up*Uq"];					   
           WS["      HfxZpE=Half/(Zeta+Eta)"];
           WS["      r1x2E=Half/Eta"];
           WS["      r1x2Z=Half/Zeta"];
           WS["      ExZpE=Eta*r1xZpE"];
           WS["      ZxZpE=Zeta*r1xZpE"];
           WS["      Omega=ExZpe+ZxZpE"];

           WS["      Wx=(Zeta*Px+Eta*Qx)*r1xZpE"];
           WS["      Wy=(Zeta*Py+Eta*Qy)*r1xZpE"];
           WS["      Wz=(Zeta*Pz+Eta*Qz)*r1xZpE"];

           WS["      PAx=Px-Ax"];
           WS["      PAy=Py-Ay"];
           WS["      PAz=Pz-Az"];



           WS["      PQx=Px-Qx"];
           WS["      PQy=Py-Qy"];
           WS["      PQz=Pz-Qz"];
           WS["      WPx=Wx-Px"];
           WS["      WPy=Wy-Py"];
           WS["      WPz=Wz-Pz"];
           WS["      WQx=Wx-Qx"];
           WS["      WQy=Wy-Qy"];
           WS["      WQz=Wz-Qz"];

           WS["      T=Omega*(PQx*PQx+PQy*PQy+PQz*PQz)"];

           PunchGammas[Subroutine,LBra+LKet];
];
PunchVRRBack[Subroutine_,BKType_]:=Block[{WS},
      WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];
      WS["   ENDDO ! (M0| loop"];            					 
      WS["ENDDO ! |N0) loop"];            					 
];


MakeList="TwoEObjs= \\ \n";
RelsList="TwoERels= \\ \n";

Do[Do[Do[Do[

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

            ijklType=1000*IntegralClass[Classes[[ic]]] \
                     +100*IntegralClass[Classes[[jc]]] \
                      +10*IntegralClass[Classes[[kc]]] \
	                 +IntegralClass[Classes[[lc]]];

	   Subroutine=StringJoin["Int",ToString[ijklType],".F90"];

	   OpenWrite[Subroutine];
	   Print[" Openned ",Subroutine];

	   WS[String_]:=WriteString[Subroutine,"   ",String,"\n"];

	   WriteString[Subroutine,CommentLine]; 

	   MakeList=StringJoin[MakeList,StringJoin["Int",ToString[ijklType],".o \\ \n"]];
           RelsList=StringJoin[RelsList,StringJoin["Int",ToString[ijklType],".x \\ \n"]]; 

           BraEll=imax+jmax;
           KetEll=kmax+lmax;

           PunchFront[Subroutine,imax,jmax,kmax,lmax,ijklType]; 

           SetOptions[FortranAssign,AssignOptimize->True,AssignMaxSize->200, \
           AssignBreak->{132," & \n          "},AssignIndent->"            ",\
           AssignTemporary->{W,Sequence}];


           PunchVRRClass[Subroutine,BraEll,KetEll];

           PunchVRRBack[Subroutine,BKType];


           SetOptions[FortranAssign,AssignOptimize->True,AssignMaxSize->200, \
           AssignBreak->{132," & \n          "},AssignIndent->"      ",\
           AssignTemporary->{W,Sequence}];

           WS["   ! HRR "];
           PunchHRRClass[Subroutine,ic,jc,kc,lc]; 


           WS[StringJoin["END SUBROUTINE Int",ToString[ijklType]]];

           Close[Subroutine];
           Print[" Closed ",Subroutine];

,{ic,1,LC}]
,{jc,1,LC}]
,{kc,1,LC}]
,{lc,1,LC}];


Print[MakeList];

Print[RelsList];






