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
!    GENERATE AND WRITE INTERPOLATION TABLES TO DISK
!    Author: Matt Challacombe
!------------------------------------------------------------------------------
*)
Do[
   Funk=FuncList[[m]];
   FileName = StringJoin[Funk,".Inc"];
   FileName2 = StringJoin[Funk,"_77.Inc"];
   Print["FileName = ",FileName];
   Print["FileName = ",FileName2];
   OpenWrite[FileName];
   OpenWrite[FileName2];
   WriteString[FileName2,StringJoin["      DATA ",FuncList[[m]],"Asymp   /",FF[Abs[2 n - 1]!!/(2 (2)^n) Sqrt[Pi]],"/ \n"]];

   ClearAll[b];
   time = SessionTime[];
   Do[ 
      x1 = 2 j Delta ;
      x2 = x1+(2) Delta ;
      Int = N[ Normal[ Series[ FunctionList[(x1+x2)/2+(x2-x1)x/2][[m]] ,{x,0,NTerms  +  30  }] ] , WP];
      Do[ a[i] = N[ (2/Pi) Int/.ChebyshevRules[i,x] , WP ] ,{i,0,NTerms}];
      PolyAnna = Expand[( a[0]/2  + Sum[ a[k] ChebyshevT[k,x] ,{k,1,NTerms}] )/.{x->(2*x-x1-x2)/(x2-x1)} ];
      Do[ b[i,j] = Coefficient[ PolyAnna , x^i ] , {i,NTerms} ];
      (* b[0,j] =  PolyAnna - Sum[ b[i,j] x^i ,{i,NTerms}] ; *)
      b[0,j] =  PolyAnna/.x->0;
     ,{j,0,NInterps}];

   (*
   Print["Error in Multipole Approx of F_",m," = ", Log[10,Abs[DTest[m-1,switch,30]]]];
   PlotMe = {};
   Do[ T =SetPrecision[ (i)*Delta*0.001, 16];
       j=Floor[T/(2 Delta)];             
       do[ b[k,j]=SetPrecision[b[k,j],16],{k,0,NTerms}];
       FF = Sum[b[k,j] T^k,{k,0,NTerms}] ;
       (*       Print[i," T = ",T," F = ",FF,F[m-1,T]]; *)
       xxx=Log[10,Abs[ FF-F[m-1,T] ]];
       PlotMe = Append[PlotMe,{T,xxx }]; 
       ,{i,0,1000}];
   ListPlot[PlotMe,PlotJoined->True];     
    *)

   steps=21;
   div  = Floor[NInterps/steps];
   Do[
      Do[
         begin = j steps;
         end   = (j+1)steps -1;
	 (* ----------------------------------------------------------------------------------------- *)
         DataString=StringJoin["      DATA (",FuncList[[m]],"_",ToString[i],
                               "(I",Funk,"),I",Funk,"=",ToString[begin],",",ToString[end],")/ & \n"];
         WriteString[FileName,DataString];
         Do[WriteString[FileName,StringJoin["      ",FF[b[i,k+0]],", ",FF[b[i,k+1]], ", ",FF[b[i,k+2]],", & \n"]];
            ,{k,begin,end-4,3}];
         WriteString[FileName,StringJoin["      ",FF[b[i,end-2]],", ",FF[b[i,end-1]], ", ",FF[b[i,end]],"/\n"]];
	 (* ----------------------------------------------------------------------------------------- *)
         DataString=StringJoin["      DATA (",FuncList[[m]],"_",ToString[i],
                               "(I",Funk,"),I",Funk,"=",ToString[begin],",",ToString[end],")/  \n"];
         WriteString[FileName2,DataString];
         Do[WriteString[FileName2,StringJoin["     >",FF[b[i,k+0]],", ",FF[b[i,k+1]], ", ",FF[b[i,k+2]],",  \n"]];
            ,{k,begin,end-4,3}];
         WriteString[FileName2,StringJoin[   "     >",FF[b[i,end-2]],", ",FF[b[i,end-1]], ", ",FF[b[i,end]],"/\n"]];

	 (* ----------------------------------------------------------------------------------------- *)
        ,{j,0,div-1}];
        begin = div steps;
        end   = NInterps;
        diff  = end-begin+1;
        (* ----------------------------------------------------------------------------------------- *)

        DataString=StringJoin["      DATA (",FuncList[[m]],"_",ToString[i],
                              "(I",Funk,"),I",Funk,"=",ToString[begin],",",ToString[end],")/ & \n"];
        WriteString[FileName,DataString];
        DataString=StringJoin["      DATA (",FuncList[[m]],"_",ToString[i],
                              "(I",Funk,"),I",Funk,"=",ToString[begin],",",ToString[end],")/  \n"];

        WriteString[FileName2,DataString];
       (* ----------------------------------------------------------------------------------------- *)

        If[diff>=4,
           Do[WriteString[FileName,StringJoin["      ",FF[b[i,k+0]],", ",FF[b[i,k+1]], ", ",FF[b[i,k+2]],", & \n"]];
              WriteString[FileName2,StringJoin["     >",FF[b[i,k+0]],", ",FF[b[i,k+1]], ", ",FF[b[i,k+2]],",  \n"]];
              ,{k,begin,end-3,3}];
           ];
         If[ Mod[diff,3]==1,
             WriteString[FileName,StringJoin["      ",FF[b[i,end]],"/ \n"]];
             WriteString[FileName2,StringJoin["     >",FF[b[i,end]],"/ \n"]];
            ,If[ Mod[diff,3]==2, 
                 WriteString[FileName,StringJoin["      ",FF[b[i,end-1]],", ",FF[b[i,end]],"/\n"]];
                 WriteString[FileName2,StringJoin["     >",FF[b[i,end-1]],", ",FF[b[i,end]],"/\n"]];
                ,If[ Mod[diff,3]==0, 
                      WriteString[FileName,StringJoin["      ",FF[b[i,end-2]],", ",FF[b[i,end-1]], ", ",FF[b[i,end]],"/\n"]];
                      WriteString[FileName2,StringJoin["     >",FF[b[i,end-2]],", ",FF[b[i,end-1]], ", ",FF[b[i,end]],"/\n"]];
                    ];
                ];
            ];
    ,{i,0,NTerms}];
   Print[FuncList[[m]],StringForm[" // Timing in hrs = ``",(SessionTime[]-time)/3600 ]];
   Close[FileName];
   Close[FileName2];
   ,{m,NFunctions}];
