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
!    WRITE EXPLICIT CODE FOR CONVERSION OF HERMITE GAUSSIAN MULTIPOLE TO
!    SPHERICALS
!------------------------------------------------------------------------------
*)

(* GET THE MAX ANGULAR SYMMETRY TO BE USED *)

MondoHome=Environment["FREEON_HOME"];
If[MondoHome==$FAILED,
   Print["COULD NOT FIND $FREEON_HOME! CHECK YOUR .cshrc "];
   Abort[];
  ];
EllFile = StringJoin[MondoHome,"/Includes/Ell.m"];
Get[EllFile];

HGEll=7;

(* INDEXING OF ARRAY ELEMENTS *)

LHGTF[L_]:=(L+1)*(L+2)*(L+3)/6;                                    
LBegin[L_]:=(L*(L+1)*(L+2))/6+1;
LMNDex[L_,M_,N_]:=LBegin[L+M+N]+N*(2*(L+M+N)-N+3)/2+M
SPDex[L_,M_]:=Binomial[L+1,2]+M;
LSP[L_]:=L*(L+3)/2;   

(* THE IRREGULAR SPHERICAL HARMONICS *)
(*
p=Sqrt[px^2+py^2+pz^2];
cs=pz/p;
sn=Sqrt[1-cs^2];
qc[l_,m_]:= p^l Cos[ m ArcTan[py,px]] LegendreP[l,m,cs]/(l+m)! ;
qs[l_,m_]:= p^l Sin[ m ArcTan[py,px]] LegendreP[l,m,cs]/(l+m)! ;
*)
ClearAll[p,px,py,pz,cos,sin,e,y,q];

p=Sqrt[px^2+py^2+pz^2]; cos=pz/p;
qc[n_,m_]:= p^n Cos[ m ArcTan[py,px]] LegendreP[n,m,cos]/(n+m)! ;
qs[n_,m_]:= p^n Sin[ m ArcTan[py,px]] LegendreP[n,m,cos]/(n+m)! ;

(* TRANSFORMATION FUNCTIONS FROM HG TO SP *)

HGToSPC[L_,M_,l_,m_,n_]:=Module[{LMN,DC},
                                 DC=D[D[D[qc[L,M],{px,l}],{py,m}],{pz,n}];
                                 DC=Limit[DC,px->0]; 
                                 DC=Limit[DC,py->0]; 
                                 DC=Limit[DC,pz->0]; 
                                 LMN=LMNDex[l,m,n];
				DC*PiZ*HGCo[LMN,J]
                                ];

HGToSPS[L_,M_,l_,m_,n_]:=Module[{LMN,DS},
                                 DS=D[D[D[qs[L,M],{px,l}],{py,m}],{pz,n}];
                                 DS=Limit[DS,px->0]; 
                                 DS=Limit[DS,py->0]; 
                                 DS=Limit[DS,pz->0];                                     
                                 LMN=LMNDex[l,m,n];
				DS*PiZ*HGCo[LMN,J]
                                ];

(* GENERATE THE EXPLICIT TRANSFORMATIONS *)

Do[Do[
      LM=SPDex[L,M];
      SPC[LM]=0;
      SPS[LM]=0;
      Do[Do[Do[  
               SPC[LM]=SPC[LM]+HGToSPC[L,M,l,m,n];
               SPS[LM]=SPS[LM]+HGToSPS[L,M,l,m,n];
      ,{n,0,L-m-l}],{m,0,L-l}],{l,0,L}];
      SPC[LM]=Factor[Simplify[SPC[LM]]];
      SPS[LM]=Factor[Simplify[SPS[LM]]];
,{M,0,L}],{L,0,HGEll}];

(* LOAD OPTIMIZING AND FORMATING ROUTINES FOR MMA AND SET ASSOCIATED OPTIONS *)

Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];
Get[StringJoin[MondoHome,"/MMA/Format.m"]];
Get[StringJoin[MondoHome,"/MMA/Optimize.m"]];


SetOptions[Optimize,OptimizePower->Binary,OptimizeFunction->True,OptimizeVariable->{o,Array}];
SetOptions[FortranAssign,AssignOptimize->True,AssignIndent->"           ",AssignMaxSize->400,AssignPrecision->Infinity,AssignBreak->{132," & \n          "},AssignTemporary->{w,Array}];
SetOptions[OpenWrite, PageWidth -> 200];

(*
SetOptions[FortranAssign,AssignOptimize->False,AssignMaxSize->500,AssignBreak->{132," & \n          "},
                         AssignTemporary->{W,Array},AssignTemporary->{W,Array},AssignToArray->{Cp,Sp,HGCo}];

SetOptions[OpenWrite, PageWidth -> 200];
 *)

(* PUT THE TRANSFORMATIONS TO FILE *)

FileName="HGToSP_PoleNode.Inc";
Print[" Openned ",FileName];
OpenWrite[FileName];

Do[
   WriteString[FileName,"     CASE(",L,") \n"];
   WriteString[FileName,"        DO J=1,HG%NQ(",L,") \n"];
   WriteString[FileName,"           PiZ=(Pi/HG%Zeta(",L,")%D(J))**(ThreeHalves) \n"];
   Do[
      ls=ToString[l];
      Write[FileName,FortranAssign[c,SPC[l],AssignReplace->{"c"->StringJoin["CMTmp%D(J,",ls,")"],"HGCo"->StringJoin["HG%Coef(",ToString[L],")%D"]," "->""}]];
      Write[FileName,FortranAssign[s,SPS[l],AssignReplace->{"s"->StringJoin["SMTmp%D(J,",ls,")"],"HGCo"->StringJoin["HG%Coef(",ToString[L],")%D"]," "->""}]];
     ,{l,0,LSP[L]}];
   WriteString[FileName,"        ENDDO \n"];
  ,{L,0,HGEll}];
Close[FileName];          
Print[" Closed ",FileName];



FileName="HGToSP_Density.Inc";
Print[" Openned ",FileName];
OpenWrite[FileName];

Do[
   WriteString[FileName,"     CASE(",L,") \n"];
   WriteString[FileName,"        DO J=L,U \n"];
   WriteString[FileName,"           N=N+1 \n"];
   WriteString[FileName,"           IAdd=OffQ+J \n"];
   WriteString[FileName,"           JAdd=OffR+(J-1)*",ToString[LHGTF[L]]," \n"];
   WriteString[FileName,"           PTmp%D(:,N)=(/Rho%Qx%D(IAdd),Rho%Qy%D(IAdd),Rho%Qz%D(IAdd)/) \n"];
   Do[
      ls=ToString[l];
      Write[FileName,FortranAssign[c,SPC[l],AssignReplace->{"c"->StringJoin["CTmp%D(N,",ls,")"],"HGCo"->StringJoin["Rho%Co%D"],",J"->"+JAdd"," "->""}]];
      Write[FileName,FortranAssign[s,SPS[l],AssignReplace->{"s"->StringJoin["STmp%D(N,",ls,")"],"HGCo"->StringJoin["Rho%Co%D"],",J"->"+JAdd"," "->""}]];
     ,{l,0,LSP[L]}];
   WriteString[FileName,"        ENDDO \n"];
  ,{L,0,HGEll}];
Close[FileName];          
Print[" Closed ",FileName];


(*

           CASE(1) 
              
                 PiZ=(Pi/HG%Zeta(Ell)%D(J))**(ThreeHalves)
                 C(J,0)=PiZ*HG%Coef(Ell)%D(1,J)
                 S(J,0)=0
                 C(J,1)=PiZ*HG%Coef(Ell)%D(4,J)
                 S(J,1)=0
                 C(J,2)=-5.D-1*PiZ*HG%Coef(Ell)%D(3,J)
                 S(J,2)=-5.D-1*PiZ*HG%Coef(Ell)%D(2,J)
              ENDDO
           CASE(2) 

 *)
