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
!    WRITE INTERFACE FOR ERI'S.
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

(* Minimal 
   Classes = { {0,0},{1,1}} 
   Classes = { {0,0},{1,1}} 
*)
(* Maximal *)
 Classes = { {0,0},{1,1},{2,2},{3,3}}


CType[1]  = "S";
CType[2]  = "SP";
CType[3]  = "P";
CType[6]  = "D";
CType[10] = "F";

LC=Length[Classes];
MxEll = Max[Classes];


Print[" LC = ",LC];

(* GENERATE INTERFACE FOR ERIS *)

FileName="dERIInterfaceB.Inc";

SetOptions[OpenWrite, PageWidth -> 200];

WS[String_]:=WriteString[FileName,"   ",String,"\n"];

OpenWrite[FileName];
Print[" Openned ",FileName];

WS["SELECT CASE(IntType)"];

Do[Do[Do[Do[
          ijklType=1000000*IntegralClass[Classes[[ic]]] \
                  +10000*IntegralClass[Classes[[jc]]] \
	          +100*IntegralClass[Classes[[kc]]] \
	          +IntegralClass[Classes[[lc]]];


   imin = Classes[[ic, 1]]; imax = Classes[[ic, 2]];
   jmin = Classes[[jc, 1]]; jmax = Classes[[jc, 2]];
   kmin = Classes[[kc, 1]]; kmax = Classes[[kc, 2]];
   lmin = Classes[[lc, 1]]; lmax = Classes[[lc, 2]];



    WS[StringJoin["CASE(",ToString[ijklType],")"]];
(*    WS[StringJoin[" write(*,*) '",ToString[ijklType],"'"]]; *)

    If[IntegralClass[Classes[[ic]]]>=IntegralClass[Classes[[jc]]],
    itype=IntegralClass[Classes[[ic]]];
    jtype=IntegralClass[Classes[[jc]]];
    braswitch=0;
    ,
    jtype=IntegralClass[Classes[[ic]]];
    itype=IntegralClass[Classes[[jc]]];
    braswitch=1;
    ];


    If[IntegralClass[Classes[[kc]]]>=IntegralClass[Classes[[lc]]],
    ktype=IntegralClass[Classes[[kc]]];
    ltype=IntegralClass[Classes[[lc]]];
    ketswitch=0;
    ,
    ltype=IntegralClass[Classes[[kc]]];
    ktype=IntegralClass[Classes[[lc]]];
    ketswitch=1;
    ];


    If[itype*100+jtype>=ktype*100+ltype,
    braketswitch=0;
    ,
    itmp=itype;
    itype=ktype;
    ktype=itmp;
    itmp=jtype;
    jtype=ltype;
    ltype=itmp;
    braketswitch=1;
    ];


    Print["braketswitch ",braketswitch," braswitch ",braswitch," ketswitch ",ketswitch," ijkl ",itype,jtype,ktype,ltype];

          ijklType=1000000*itype \
                  + 10000*jtype \
	          +  100*ktype \
	          +     ltype;



                  ArgString1="";
                  ArgString2="";
                  ArgString3="";
                  ArgString4="";
                  ArgString5="";
                  ArgString6="";
                  ArgString7="";
                  ArgString8="";

                  If[braketswitch==0,
                      ArgString1="ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & \n";
                      ArgString2="                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & \n";
                      ArgString3="                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & \n ";
                      If[braswitch==0,
                          ArgString4="                 OffSet%A  ,1             ";
                          ArgString5="                 OffSet%C-1,NBFA*NBFB     ";
                          If[ketswitch==0,
                               ArgString6="                 OffSet%B-1,NBFA          ";
                               ArgString7="                 OffSet%D-1,NBFA*NBFB*NBFC";
                               ArgString8="                 1,4,7,10";
                          ,
                               ArgString6="                 OffSet%D-1,NBFA*NBFB*NBFC";
                               ArgString7="                 OffSet%B-1,NBFA          ";
                               ArgString8="                 1,4,10,7";
                          ];
                      ,
                          ArgString4="                 OffSet%C-1,NBFA*NBFB     ";
                          ArgString5="                 OffSet%A  ,1             ";
                          If[ketswitch==0,
                               ArgString6="                 OffSet%B-1,NBFA          ";
                               ArgString7="                 OffSet%D-1,NBFA*NBFB*NBFC";
                               ArgString8="                 4,1,7,10";
                          ,
                               ArgString6="                 OffSet%D-1,NBFA*NBFB*NBFC";
                               ArgString7="                 OffSet%B-1,NBFA          ";
                               ArgString8="                 4,1,10,7";
                          ];
                      ];
                  ,
                      ArgString1="BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & \n";
                      ArgString2="                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & \n";
                      ArgString3="                  BDAtmPair(iFBD)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & \n ";                  
                      If[ketswitch==0,
                         ArgString4="                 OffSet%B-1,NBFA          ";
                         ArgString5="                 OffSet%D-1,NBFA*NBFB*NBFC";
                         If[braswitch==0,
                            ArgString6="                 OffSet%A  ,1             ";
                            ArgString7="                 OffSet%C-1,NBFA*NBFB     ";
                            ArgString8="                 7,10,1,4";
                         ,
                            ArgString6="                 OffSet%C-1,NBFA*NBFB     ";
                            ArgString7="                 OffSet%A  ,1             ";
                            ArgString8="                 7,10,4,1";
                         ];
                      ,
                         ArgString4="                 OffSet%D-1,NBFA*NBFB*NBFC";
                         ArgString5="                 OffSet%B-1,NBFA          ";
                         If[braswitch==0,
                            ArgString6="                 OffSet%A  ,1             ";
                            ArgString7="                 OffSet%C-1,NBFA*NBFB     ";
                            ArgString8="                 10,7,1,4";
                         ,
                            ArgString6="                 OffSet%C-1,NBFA*NBFB     ";
                            ArgString7="                 OffSet%A  ,1             ";
                            ArgString8="                 10,7,4,1";
                         ];
                      ];
                  ];

                  ArgString=StringJoin[ArgString1,ArgString2,ArgString3, \
                                       ArgString4,", & \n ",ArgString5,", & \n ", \
                                       ArgString6,", & \n ",ArgString7,", & \n ", \
                                       ArgString8,",NIntBlk,GMc%PBC,C(1)"];



                  WS[StringJoin["  CALL dIntB",ToString[ijklType],"(",ArgString,")"]];

                  WS[StringJoin["  LocNInt=",ToString[12*(LEnd[imax]-LBegin[imin]+1)*(LEnd[jmax]-LBegin[jmin]+1)* \
                                                      (LEnd[kmax]-LBegin[kmin]+1)*(LEnd[lmax]-LBegin[lmin]+1)]]];

,{ic,1,LC}]
,{jc,1,LC}]
,{kc,1,LC}]
,{lc,1,LC}];

Print["We have ",LC*LC*LC*LC," integrals."];

WS["CASE DEFAULT"];
WS["   WRITE(*,*) 'We are in DERIInterface2.Inc'"];
WS["   WRITE(*,*) 'IntType=',IntType"];
WS["   STOP 'MISS AN INTEGRAL'"];
WS["END SELECT"];

Close[FileName];
Print[" Closed ",FileName];


(* GENERATE INTERFACE FOR THE LIST *)


FileName="DERIListInterface2.Inc";

SetOptions[OpenWrite, PageWidth -> 200];

WS[String_]:=WriteString[FileName,"   ",String,"\n"];

OpenWrite[FileName];
Print[" Openned ",FileName];

WS["SELECT CASE(IntType)"];

Do[Do[
   
   imin = Classes[[ic, 1]]; imax = Classes[[ic, 2]];
   jmin = Classes[[jc, 1]]; jmax = Classes[[jc, 2]];
   
   Print["imin =",imin," imax =",imax," nbr=",LEnd[imax]-LBegin[imin]+1];
   
   ijklType=100*IntegralClass[Classes[[ic]]]+IntegralClass[Classes[[jc]]];

   WS[StringJoin["CASE(",ToString[ijklType],")"]];
   
   If[IntegralClass[Classes[[ic]]]>=IntegralClass[Classes[[jc]]],
      itype=IntegralClass[Classes[[ic]]];
      jtype=IntegralClass[Classes[[jc]]];
      braswitch=0;
      ketswitch=0;
   ,
      jtype=IntegralClass[Classes[[ic]]];
      itype=IntegralClass[Classes[[jc]]];
      braswitch=1;
      ketswitch=1;
   ];

   Print["braketswitch ",braketswitch," braswitch ",braswitch," ketswitch ",ketswitch," ijkl ",itype,jtype,ktype,ltype];

   ijklType=1000000*itype+10000*jtype+100*itype+jtype;

   ArgString1="";
   ArgString2="";
   ArgString3="";
   ArgString4="";
   ArgString5="";
   ArgString6="";
   ArgString7="";
   
   ArgString1="ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & \n";
   ArgString2="                  ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & \n";
   ArgString3="                  ACAtmPair(iFAC)%SP%AtmInfo,ACAtmPair(iFAC)%SP%AtmInfo, & \n ";
   If[braswitch==0,

   Print["value=",(LEnd[imax]-LBegin[imin]+1)*(LEnd[jmax]-LBegin[jmin]+1)*(LEnd[imax]-LBegin[imin]+1)];

      ArgString4=StringJoin["                 1,1,0,",ToString[LEnd[imax]-LBegin[imin]+1]];
      ArgString5=StringJoin["0,",ToString[(LEnd[imax]-LBegin[imin]+1)*(LEnd[jmax]-LBegin[jmin]+1)],",",\
                            "0,",ToString[(LEnd[imax]-LBegin[imin]+1)^2*(LEnd[jmax]-LBegin[jmin]+1)]];
   ,
      ArgString4=StringJoin["                 1,1,0,",ToString[LEnd[jmax]-LBegin[jmin]+1]];
      ArgString5=StringJoin["0,",ToString[(LEnd[jmax]-LBegin[jmin]+1)*(LEnd[imax]-LBegin[imin]+1)],",",\
                            "0,",ToString[(LEnd[jmax]-LBegin[jmin]+1)^2*(LEnd[imax]-LBegin[imin]+1)]];
   ];


   ArgString=StringJoin[ArgString1,ArgString2,ArgString3, \
             ArgString4,",",ArgString5,",", \
             ToString[(LEnd[imax]-LBegin[imin]+1)^2*(LEnd[jmax]-LBegin[jmin]+1)^2], \
             ",1,4,7,10,GMc%PBC,C(1)"]; (* 1,4,7,10 -> THIS IS NOT RIGHT, BUT DOES NOT MATTER*)

   WS[StringJoin["  CALL dIntB",ToString[ijklType],"(",ArgString,")"]];

(*   WS[StringJoin["  LocNInt=",ToString[12*(LEnd[imax]-LBegin[imin]+1)^2*(LEnd[jmax]-LBegin[jmin]+1)^2]]]; *)

,{ic,1,LC}]
,{jc,1,LC}];

Print["We have ",LC*LC," integrals."];

WS["CASE DEFAULT"];
WS["   STOP 'MISS AN INTEGRAL'"];
WS["END SELECT"];

Close[FileName];
Print[" Closed ",FileName];




