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

(* Minimal *)
   Classes = { {0,0},{1,1},{2,2}} 



(*   Classes = { {0,0},{1,1},{2,2}}*)

(* Maximal 
   Classes = { {0,0},{0,1},{1,1},{2,2},{3,3}}
 *)

CType[1]  = "S";
CType[2]  = "SP";
CType[3]  = "P";
CType[6]  = "D";
CType[10] = "F";

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

                              pmin=Position[a,Min[a]][[1, 1]];
                              qmin=Position[c,Min[c]][[1, 1]];
                              If[ a[[pmin]] < 0 || c[[qmin]] < 0,Return[0]];                                        


                              p=Position[a,Max[a]][[1, 1]];
                              q=Position[c,Max[c]][[1, 1]];
        
                              If[a[[p]]==0 && c[[q]]==0,
                                 Return[ToExpression[StringJoin["AuxR",ToString[m]]]]];                                        

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
                                 Return[PA*VRR[a1,c,m]+WP*VRR[a1,c,m+1]+Ai1*(VRR[a2,c,m]-ExZpE*VRR[a2,c,m+1]) +CiO2zn*VRR[a1,c1,m+1]]
			       ,(* Taking down ell on ket side *)
			         If[q == 1, one = {1, 0, 0}; two={2,0,0}; QC=QCx; WQ=WQx;];
			         If[q == 2, one = {0, 1, 0}; two={0,2,0}; QC=QCy; WQ=WQy;];
			         If[q == 3, one = {0, 0, 1}; two={0,0,2}; QC=QCz; WQ=WQz;];
                                 Ci1=(c[[q]]-1)*r1x2E ; (* (c[i]-1)/(2 eta)   *)
			         AiO2zn=a[[q]]*HfxZpE; (* a[i]/(2(zeta+eta)) *)
                                 c1 = c - one;
                                 c2 = c - two;
                                 a1 = a - one;  
                                 Return[QC*VRR[a,c1,m]+WQ*VRR[a,c1,m+1]+Ci1*(VRR[a,c2,m]-ZxZpE*VRR[a,c2,m+1])+AiO2zn*VRR[a1,c1,m+1]]
                               ];
                             ];

(**************** THE HRR RELATIONS: THX HGP! ************************)

     HRR[a_List,b_List,c_List,d_List]:=Module[{AB,CD,pa,pb,pc,pd,MaxEll,a1,b1,c1,d1,one,two,adex,cdex},
                              pa=Position[a,Max[b]][[1, 1]];
                              pb=Position[b,Max[b]][[1, 1]];
                              pc=Position[c,Max[b]][[1, 1]];
                              pd=Position[d,Max[d]][[1, 1]];
                              MaxEll=Max[Join[b,d]];
			      (* Exit condition 1 *)
                              If[ a[[pa]] < 0 || b[[pb]] < 0 || c[[pc]] < 0 || d[[pd]] < 0 ,Return[0];]; 
			      (* Exit condition 2 *)
                              If[ b[[pb]]==0 && d[[pd]]==0,
                                 adex=LMNDex[a[[1]],a[[2]],a[[3]]];
                                 cdex=LMNDex[c[[1]],c[[2]],c[[3]]];
                                 BarExp=ToExpression[StringJoin["I",ToString[adex],"Bar",ToString[cdex]]];
                                 Return[BarExp]
                                ];
			      (* Recursion *)
			      If[b[[pb]]==MaxEll,(* Taking down ell on b *)
                                 If[pb==1, one = {1, 0, 0}; AB = ABx; ];
                                 If[pb==2, one = {0, 1, 0}; AB = ABy; ];
                                 If[pb==3, one = {0, 0, 1}; AB = ABz; ];
                                 a1 = a + one;
                                 b1 = b - one;                
                                 Return[ HRR[a1, b1, c, d] + AB  HRR[a, b1, c, d]]
                                ];
			      If[d[[pd]]==MaxEll,(* Taking down ell on d *)
                                 If[pd==1, one = {1, 0, 0}; CD = CDx; ];
                                 If[pd==2, one = {0, 1, 0}; CD = CDy; ];
                                 If[pd==3, one = {0, 0, 1}; CD = CDz; ];
                                 c1=c+one;
                                 d1=d-one;                
                                 Return[HRR[a,b,c1,d1]+CD HRR[a,b,c,d1]]
                                ];
                            ];

(* LOAD OPTIMIZING AND FORMATING ROUTINES FOR MMA AND SET ASSOCIATED OPTIONS *)

Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];
Get[StringJoin[MondoHome,"/MMA/Format.m"]];
Get[StringJoin[MondoHome,"/MMA/Optimize.m"]];

FF[x_] := ToString[FixedNumberForm[SetPrecision[N[x,32],32], 16, 2]];

SetOptions[FortranAssign,AssignOptimize->False,AssignMaxSize->350,AssignBreak->{300," & \n          "},AssignIndent->"      ",AssignTemporary->{W,Sequence}];
SetOptions[Optimize,OptimizeVariable->{V,Array},OptimizeTimes->True,OptimizePlus->True,OptimizeCoefficients->True,OptimizeFunction->False]; 
SetOptions[OpenWrite, PageWidth -> 200];

SetAttributes[o,NHoldAll];
SetAttributes[SSSS,NHoldAll];
SetAttributes[MBarN,NHoldAll];

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
						       MBarString=StringJoin["I",ToString[i],"Bar",ToString[k]];
                                                       oList=Append[oList,StringJoin["o(",ToString[Kount],")"]->MBarString];
                                                 ,{i,1,LEnd[BraEll]}];
                                                 ,{k,1,LEnd[KetEll]}];
                                                 Write[FileName,FortranAssign[o,IList,AssignReplace->oList]];
                                                 ];


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
                     WS[StringJoin["  AuxR0=Upq*",FStr[0]]];
            ];

         If[EllTot==1,
                     WS[StringJoin["  AuxR0=Upq*",FStr[0]]];
                     WS[StringJoin["  AuxR1=Upq*",FStr[1]]];
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
          WS[StringJoin["  AuxR",GSt,"=Upq*W",GSt,""]];
                           ,{L,0,EllTot}];
            ]; 

                     WS["ELSE"];
                     WS["  InvT=One/T"];
                     WS["  SqInvT=DSQRT(InvT)"];
                         Do[GSt=ToString[L];

           WS[StringJoin["  AuxR",GSt,"=",GammAss[L],"*Upq*SqInvT"]];

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
						 oList={" "->"","u"->"(","v"->",","w"->")","x1"->"CDx","y1"->"CDy","z1"->"CDz","x2"->"ABx","y2"->"ABy","z2"->"ABz","In"->"I"};
						 IList={};
						 Kount = 0;
                                                 Do[Do[Do[Do[
                                                 Do[Do[Do[Do[
                                                             Kount = Kount + 1;
                                                             a = {lx[i], my[i], nz[i]};
                                                             b = {lx[j], my[j], nz[j]};
                                                             c = {lx[k], my[k], nz[k]};
                                                             d = {lx[l], my[l], nz[l]};
                                                             IList=Append[IList,HRR[a,b,c,d]+In[OffSet]];
                                                             oList=Append[oList,StringJoin["o(",ToString[Kount],")"]-> 
                                                                                StringJoin["OffSet=(OA+",ToString[i-LBegin[il]],")*LDA",     \
                                                                                                 "+(OB+",ToString[j-LBegin[jl]],")*LDB",     \
                                                                                                 "+(OC+",ToString[k-LBegin[kl]],")*LDC",     \
                                                                                                 "+(OD+",ToString[l-LBegin[ll]],")*LDD \n", \
											   "      I(OffSet)"]];

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
           WriteString[Subroutine,StringJoin["   SUBROUTINE Int",ToString[IJKL],"(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, & \n", \
                                             "                              OA,LDA,OB,LDB,OC,LDC,OD,LDD,PBC,I) \n"]];

	   WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];

	   WS["USE DerivedTypes"];
	   WS["USE VScratch"];
	   WS["USE GlobalScalars"];
           (*WS["USE ONX2DataType"];*)
           WS["USE ShellPairStruct"];
	   If[LBra+LKet==1,
              WS["USE GammaF0"];
              WS["USE GammaF1"];,
              WS[StringJoin["USE GammaF",ToString[LBra+LKet]]]];
           WS["IMPLICIT REAL(DOUBLE) (A,I,W)"];
           WS["INTEGER        :: LBra,LKet"];
           WS["REAL(DOUBLE)   :: PrmBufB(10,LBra),PrmBufK(10,LKet)"];
	   WS["TYPE(SmallAtomInfo) :: ACInfo,BDInfo"];
           WS["TYPE(PBCInfo) :: PBC"];
	   LenBra=LEnd[LBra];
           LenKet=LEnd[LKet];

           WS[StringJoin["REAL(DOUBLE) :: I(*)"]];

           WS["REAL(DOUBLE)  :: Zeta,Eta,r1xZpE,HfxZpE,r1x2E,r1x2Z,ExZpE,ZxZpE,Omega,Up,Uq,Upq"];
           WS["REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz"];
           WS["REAL(DOUBLE)  :: QCx,QCy,QCz,PAx,PAy,PAz,PQx,PQy,PQz,WPx,WPy,WPz,WQx,WQy,WQz   "];
           WS["REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT,ABx,ABy,ABz,CDx,CDy,CDz"];

           WS["INTEGER       :: OffSet,OA,LDA,OB,LDB,OC,LDC,OD,LDD,J,K,L"];
           WS["REAL(DOUBLE)  :: FPQx,FPQy,FPQz"];


           Do[Do[
                 WS[StringJoin["I",ToString[i],"Bar",ToString[k],"=0.0d0"]];
            ,{i,1,LEnd[LBra]}],{k,1,LEnd[LKet]}];

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
           WS["   Qx =PrmBufK(2,J)"];
           WS["   Qy =PrmBufK(3,J)"];
           WS["   Qz =PrmBufK(4,J)"];
           WS["   Uq =PrmBufK(5,J)"];


           (* Add the conversion factor for SP shell *)
           If[IntegralClass[Classes[[kc]]]==2,
              Print["IntegralClass[Classes[[kc]]]=",IntegralClass[Classes[[kc]]] ];
              WS["   C1q=PrmBufK(6,J)"];
           ];
           If[IntegralClass[Classes[[lc]]]==2,
              Print["IntegralClass[Classes[[lc]]]=",IntegralClass[Classes[[lc]]] ];
              WS["   C2q=PrmBufK(7,J)"];
           ];
           If[IntegralClass[Classes[[kc]]]==2 && IntegralClass[Classes[[lc]]]==2,
              Print["IntegralClass[Classes[[kc]]]=",IntegralClass[Classes[[kc]]] ];
              Print["IntegralClass[Classes[[lc]]]=",IntegralClass[Classes[[lc]]] ];
              WS["   C3q=PrmBufK(8,J)"];
           ];



           WS["   QCx=Qx-Cx"];
           WS["   QCy=Qy-Cy"];
           WS["   QCz=Qz-Cz"];


           WS["   DO K=1,LBra ! K^2 VRR (M0| loop "];

           WS["      Zeta=PrmBufB(1,K)"];
           WS["      Px  =PrmBufB(2,K)"];
           WS["      Py  =PrmBufB(3,K)"];
           WS["      Pz  =PrmBufB(4,K)"];
           WS["      Up  =PrmBufB(5,K)"];

           (* Add the conversion factor for SP shell *)
           If[IntegralClass[Classes[[ic]]]==2,
              Print["IntegralClass[Classes[[ic]]]=",IntegralClass[Classes[[ic]]] ];
              WS["      C1p =PrmBufB(6,K)"];
           ];
           If[IntegralClass[Classes[[jc]]]==2,
              Print["IntegralClass[Classes[[jc]]]=",IntegralClass[Classes[[jc]]] ];
              WS["      C2p =PrmBufB(7,K)"];
           ];
           If[IntegralClass[Classes[[ic]]]==2 && IntegralClass[Classes[[jc]]]==2,
              Print["IntegralClass[Classes[[ic]]]=",IntegralClass[Classes[[ic]]] ];
              Print["IntegralClass[Classes[[jc]]]=",IntegralClass[Classes[[jc]]] ];
              WS["      C3p =PrmBufB(8,K)"];
           ];


           WS["      r1xZpE=One/(Zeta+Eta)"];
	   WS["      Upq=SQRT(r1xZpE)*Up*Uq"];					   
           WS["      HfxZpE=Half/(Zeta+Eta)"];
           WS["      r1x2E=Half/Eta"];
           WS["      r1x2Z=Half/Zeta"];
           WS["      ExZpE=Eta*r1xZpE"];
           WS["      ZxZpE=Zeta*r1xZpE"];
           WS["      Omega=Eta*Zeta*r1xZpE"];

           (*WS["      Wx=(Zeta*Px+Eta*Qx)*r1xZpE"];*)
           (*WS["      Wy=(Zeta*Py+Eta*Qy)*r1xZpE"];*)
           (*WS["      Wz=(Zeta*Pz+Eta*Qz)*r1xZpE"];*)

           WS["      PAx=Px-Ax"];
           WS["      PAy=Py-Ay"];
           WS["      PAz=Pz-Az"];

           WS["      PQx=Px-Qx"];
           WS["      PQy=Py-Qy"];
           WS["      PQz=Pz-Qz"];

           (* Should improve this part, scalar replacement... *)
           (*FPQx = PQx*PBC%InvBoxSh%D(1,1)+PQy*PBC%InvBoxSh%D(1,2)+PQz*PBC%InvBoxSh%D(1,3)
           FPQy = PQy*PBC%InvBoxSh%D(2,2)+PQz*PBC%InvBoxSh%D(2,3)
           FPQz = PQz*PBC%InvBoxSh%D(3,3)
           IF(PBC%AutoW%I(1)==1) FPQx = FPQx-ANINT(FPQx)
           IF(PBC%AutoW%I(2)==1) FPQy = FPQy-ANINT(FPQy)
           IF(PBC%AutoW%I(3)==1) FPQz = FPQz-ANINT(FPQz)
           PQx  = FPQx*PBC%BoxShape%D(1,1)+FPQy*PBC%BoxShape%D(1,2)+FPQz*PBC%BoxShape%D(1,3)
           PQy  = FPQy*PBC%BoxShape%D(2,2)+FPQz*PBC%BoxShape%D(2,3)
           PQz  = FPQz*PBC%BoxShape%D(3,3)*)

           WS["      INCLUDE 'ERIMIC.Inc'"];
(*         WS["! Need to be improve..."];
           WS["      FPQx = PQx*PBC%InvBoxSh%D(1,1)+PQy*PBC%InvBoxSh%D(1,2)+PQz*PBC%InvBoxSh%D(1,3)"];
           WS["      FPQy = PQy*PBC%InvBoxSh%D(2,2)+PQz*PBC%InvBoxSh%D(2,3)"];
           WS["      FPQz = PQz*PBC%InvBoxSh%D(3,3)"];
           WS["      IF(PBC%AutoW%I(1)==1) FPQx = FPQx-ANINT(ANINT(FPQx*1d9)*1d-9)"];
           WS["      IF(PBC%AutoW%I(2)==1) FPQy = FPQy-ANINT(ANINT(FPQy*1d9)*1d-9)"];
           WS["      IF(PBC%AutoW%I(3)==1) FPQz = FPQz-ANINT(ANINT(FPQz*1d9)*1d-9)"];
           WS["      PQx  = FPQx*PBC%BoxShape%D(1,1)+FPQy*PBC%BoxShape%D(1,2)+FPQz*PBC%BoxShape%D(1,3)"];
           WS["      PQy  = FPQy*PBC%BoxShape%D(2,2)+FPQz*PBC%BoxShape%D(2,3)"];
           WS["      PQz  = FPQz*PBC%BoxShape%D(3,3)"];
           WS["!"];*)


           (*WS["      WPx=Wx-Px"];*)
           (*WS["      WPy=Wy-Py"];*)
           (*WS["      WPz=Wz-Pz"];*)
           (*WS["      WQx=Wx-Qx"];*)
           (*WS["      WQy=Wy-Qy"];*)
           (*WS["      WQz=Wz-Qz"];*)

           WS["      WPx = -Eta*PQx*r1xZpE"];
           WS["      WPy = -Eta*PQy*r1xZpE"];
           WS["      WPz = -Eta*PQz*r1xZpE"];
           WS["      WQx = Zeta*PQx*r1xZpE"];
           WS["      WQy = Zeta*PQy*r1xZpE"];
           WS["      WQz = Zeta*PQz*r1xZpE"];

           WS["      T=Omega*(PQx*PQx+PQy*PQy+PQz*PQz)"];

           PunchGammas[Subroutine,LBra+LKet];
];
PunchVRRBack[Subroutine_,BKType_]:=Block[{WS},
      WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];
      WS["   ENDDO ! (M0| loop"];            					 
      WS["ENDDO ! |N0) loop"];            					 
];


MakeList="TwoEObjs= ";
RelsList="TwoERels= ";


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

(*            ijklType=1000000*IntegralClass[Classes[[ic]]] \
                      +10000*IntegralClass[Classes[[jc]]] \
                        +100*IntegralClass[Classes[[kc]]] \
	                    +IntegralClass[Classes[[lc]]];*)


          ijklFlag=1000000*IntegralClass[Classes[[ic]]] \
                  +10000*IntegralClass[Classes[[jc]]] \
	          +100*IntegralClass[Classes[[kc]]] \
	          +IntegralClass[Classes[[lc]]];


            ijklType=1000*IntegralClass[Classes[[ic]]] \
                      +100*IntegralClass[Classes[[jc]]] \
                        +10*IntegralClass[Classes[[kc]]] \
	                    +IntegralClass[Classes[[lc]]];


             WSI[StringJoin["CASE(",ToString[ijklFlag],")"]]
	     WSI[StringJoin["CALL Int",ToString[ijklType],"(ACAtmPair(CFAC)%SP%Cst(1,1),ACAtmPair(CFAC)%SP%L, & \n", 
                            "                       BDAtmPair(CFBD)%SP%Cst(1,1),BDAtmPair(CFBD)%SP%L, & \n",
                            "                       ACAtmPair(CFAC)%SP%AtmInfo,BDAtmPair(CFBD)%SP%AtmInfo, & \n",
                            "                       OffSet%A  ,1              , & \n",
                            "                       OffSet%C-1,NBFA           , & \n",
                            "                       OffSet%B-1,NBFA*NBFC      , & \n",
                            "                       OffSet%D-1,NBFA*NBFB*NBFC,GM%PBC,C(1)) \n"]];



Print["ijklType=",ijklType," i=",IntegralClass[Classes[[ic]]]," j=",IntegralClass[Classes[[jc]]]," k=",IntegralClass[Classes[[kc]]]," l=",IntegralClass[Classes[[lc]]]];

	   Subroutine=StringJoin["Int",ToString[ijklType],".F90"];

	   OpenWrite[Subroutine];
	   Print[" Openned ",Subroutine];

	   WS[String_]:=WriteString[Subroutine,"   ",String,"\n"];

	   WriteString[Subroutine,CommentLine]; 

	   MakeList=StringJoin[MakeList,StringJoin[" \\ \n Int",ToString[ijklType],".o"]];
           RelsList=StringJoin[RelsList,StringJoin[" \\ \n Int",ToString[ijklType],".x"]]; 

           BraEll=imax+jmax;
           KetEll=kmax+lmax;

           PunchFront[Subroutine,imax,jmax,kmax,lmax,ijklType]; 

           PunchVRRClass[Subroutine,BraEll,KetEll];

           PunchVRRBack[Subroutine,BKType];

           WS["   ! HRR "];
           PunchHRRClass[Subroutine,ic,jc,kc,lc]; 


           WS[StringJoin["END SUBROUTINE Int",ToString[ijklType]]];

           Close[Subroutine];
           Print[" Closed ",Subroutine];

	       ]; 

,{ic,1,LC}]
,{jc,1,LC}]
,{kc,1,LC}]
,{lc,1,LC}];

WSI["CASE DEFAULT"];
WSI["  STOP 'MISSED AN INTEGRAL' "];
WSI["END SELECT "];

Close[IncludeFile];

MakeList=StringJoin[MakeList," \n"];
RelsList=StringJoin[RelsList," \n"]; 

(**************** Print out the Makefile ************************)

Makefile="Makefile1"

OpenWrite[Makefile];

WriteString[Makefile,"include $(MONDO_COMPILER)\n"];
WriteString[Makefile,"include $(MONDO_HOME)/Includes/Suffixes\n"];
WriteString[Makefile,"include $(MONDO_HOME)/Includes/RemoveAll\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"CPPMISC =\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"EXTRA_INCLUDES=\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"SPObObjs=ShellPairStruct.o\n"];
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
WriteString[Makefile,"purge:clean \n","rm -f $(MONDO_LIB)/libTwoE.a\n","rm -f $(REMOVESRCE)\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"TwoE:$(SPObObjs) $(TwoEObjs)\n","$(AR) $(ARFLAGS) $(MONDO_LIB)/libTwoE.a $(?:.F90=.o)\n","$(RANLIB) $(MONDO_LIB)/libTwoE.a\n"];
WriteString[Makefile,"#\n"];
WriteString[Makefile,"CTwoE:\n","rm -f $(REMOVEMISC) $(REMOVEMODS)\n","rm -f \#*\n","rm -f *~\n","ln -s /dev/null core\n","ls -l\n"];

Close[Makefile];


