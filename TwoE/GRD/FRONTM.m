
PunchFront[Subroutine_,ic_,jc_,kc_,lc_,IJKL_]:=Block[{WS,LBra,LKet,BKType,LenBra,LenKet},

	   imin = Classes[[ic, 1]]; imax = Classes[[ic, 2]];
	   jmin = Classes[[jc, 1]]; jmax = Classes[[jc, 2]];
	   kmin = Classes[[kc, 1]]; kmax = Classes[[kc, 2]];
           lmin = Classes[[lc, 1]]; lmax = Classes[[lc, 2]];

           LenI=LEnd[Classes[[ic,2]]]-LBegin[Classes[[ic,1]]]+1;
           LenJ=LEnd[Classes[[jc,2]]]-LBegin[Classes[[jc,1]]]+1;
           LenK=LEnd[Classes[[kc,2]]]-LBegin[Classes[[kc,1]]]+1;
           LenL=LEnd[Classes[[lc,2]]]-LBegin[Classes[[lc,1]]]+1;

           BEnd=LEnd[imax+jmax+1];
           KEnd=LEnd[kmax+lmax+1];

           BraMax=BEnd;
           If[ic==2&&jc==2,BraMax=BraMax+LenI+LenJ];
           If[ic==2&&jc!=2,BraMax=BraMax+1];
           If[jc==2&&ic!=2,BraMax=BraMax+LenI];

           KetMax=KEnd;
           If[kc==2&&lc==2,KetMax=KetMax+LenK+LenL];
           If[kc==2&&lc!=2,KetMax=KetMax+1];
           If[lc==2&&kc!=2,KetMax=KetMax+LenK];

           LBra=imax+jmax;
           LKet=kmax+lmax;
	   BKType=100*LBra+LKet;					   
	   LenBra=LEnd[LBra+1];
           LenKet=LEnd[LKet+1];



           WriteString[Subroutine,StringJoin["SUBROUTINE dIntB",ToString[IJKL],"(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, &\n ",
					     "OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NINT,PBC,GRADIENTS)\n "]];

	   WS[String_]:=WriteString[Subroutine,"      ",String,"\n"];

	   WS["USE DerivedTypes"];
	   WS["USE VScratchB"];
	   WS["USE GlobalScalars"];
           (*WS["USE ONX2DataType"];*)
           WS["USE ShellPairStruct"];
	   If[LBra+LKet+1==1,
              WS["USE GammaF0"];
              WS["USE GammaF1"];,
              WS[StringJoin["USE GammaF",ToString[LBra+LKet+1]]]];
           WS["IMPLICIT REAL(DOUBLE) (W)"]; 
           WS["INTEGER        :: LBra,LKet,NINT,CDOffSet"];
           WS["REAL(DOUBLE)   :: PrmBufB(10,LBra),PrmBufK(10,LKet)"];
	   WS["TYPE(SmallAtomInfo) :: ACInfo,BDInfo"];
           WS["TYPE(PBCInfo) :: PBC"];
         
           WS[StringJoin["REAL(DOUBLE)  :: GRADIENTS(NINT,12)"]];
           WS["REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz"];
           WS["REAL(DOUBLE)  :: Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz"];
           WS["REAL(DOUBLE)  :: PQx,PQy,PQz,FPQx,FPQy,FPQz"];
           WS["REAL(DOUBLE)  :: Zeta,Eta,Omega,Up,Uq,Upq"];
           WS["REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT"];
           WS["REAL(DOUBLE)  :: Alpha,Beta,Gamma"];
           WS[StringJoin["REAL(DOUBLE), DIMENSION(",ToString[BraMax],") :: HRRTmp "]];
           WS[StringJoin["REAL(DOUBLE), DIMENSION(",ToString[BraMax],",",ToString[KetMax],",",ToString[LEnd[lmax+1]],") :: HRR,HRRA,HRRB,HRRC "]];
           WS[StringJoin["REAL(DOUBLE)  :: VRR(",ToString[LenBra],",",ToString[LenKet],",0:",ToString[LBra+LKet+1],")"]];

           WS["INTEGER       :: OffSet,OA,LDA,GOA,OB,LDB,GOB,OC,LDC,GOC,OD,LDD,GOD,I,J,K,L"];
           WS["EXTERNAL InitDbl"];						     
           WS[StringJoin["CALL InitDbl(",ToString[BraMax],"*",ToString[KetMax],",HRR(1,1,1))"]];

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
           WS["   Qx=PrmBufK(2,J)"];
           WS["   Qy=PrmBufK(3,J)"];
           WS["   Qz=PrmBufK(4,J)"];
           WS["   Uq=PrmBufK(5,J)"];
           If[kc==2&&lc!=2,WS["   SpFnK=PrmBufK(6,J)"]];
           If[kc!=2&&lc==2,WS["   FnSpK=PrmBufK(6,J)"]];
           If[kc==2&&lc==2,WS["   SpSpK=PrmBufK(6,J)"];
                           WS["   FnSpK=PrmBufK(7,J)"];
                           WS["   SpFnK=PrmBufK(8,J)"]];
           WS["   Gamma =PrmBufK(9,J)"];
           WS["   QCx=Qx-Cx"];
           WS["   QCy=Qy-Cy"];
           WS["   QCz=Qz-Cz"];
           WS["   DO K=1,LBra ! K^2 VRR (M0| loop "];
           WS["      Zeta=PrmBufB(1,K)"];
           WS["      Px=PrmBufB(2,K)"];
           WS["      Py=PrmBufB(3,K)"];
           WS["      Pz=PrmBufB(4,K)"];
           WS["      Up=PrmBufB(5,K)"];
           If[ic==2&&jc!=2,WS["      SpFnB=PrmBufB(6,K)"]];
           If[ic!=2&&jc==2,WS["      FnSpB=PrmBufB(6,K)"]];
           If[ic==2&&jc==2,WS["      SpSpB=PrmBufB(6,K)"];
                           WS["      FnSpB=PrmBufB(7,K)"];
                           WS["      SpFnB=PrmBufB(8,K)"]];
           WS["      Alpha =PrmBufB(9,K)"];
           WS["      Beta  =PrmBufB(10,K)"];
           WS["      r1xZpE=One/(Zeta+Eta)"];
	   WS["      Upq=SQRT(r1xZpE)*Up*Uq"];					   
           WS["      HfxZpE=Half/(Zeta+Eta)"];
           WS["      r1x2E=Half/Eta"];
           WS["      r1x2Z=Half/Zeta"];
           WS["      ExZpE=Eta*r1xZpE"];
           WS["      ZxZpE=Zeta*r1xZpE"];
           WS["      Omega=Eta*Zeta*r1xZpE"];
           WS["      PAx=Px-Ax"];
           WS["      PAy=Py-Ay"];
           WS["      PAz=Pz-Az"];
           WS["      PQx=Px-Qx"];
           WS["      PQy=Py-Qy"];
           WS["      PQz=Pz-Qz"];
           WS["      ! Begin Minimum Image Convention "];
           WS["      FPQx = PQx*PBC%InvBoxSh%D(1,1)+PQy*PBC%InvBoxSh%D(1,2)+PQz*PBC%InvBoxSh%D(1,3)"];
           WS["      FPQy = PQy*PBC%InvBoxSh%D(2,2)+PQz*PBC%InvBoxSh%D(2,3)"];
           WS["      FPQz = PQz*PBC%InvBoxSh%D(3,3)"];
           WS["      IF(PBC%AutoW%I(1)==1)FPQx=FPQx-ANINT(FPQx-SIGN(1D-15,FPQx))"];
           WS["      IF(PBC%AutoW%I(2)==1)FPQy=FPQy-ANINT(FPQy-SIGN(1D-15,FPQy))"];
           WS["      IF(PBC%AutoW%I(3)==1)FPQz=FPQz-ANINT(FPQz-SIGN(1D-15,FPQz))"];
           WS["      PQx=FPQx*PBC%BoxShape%D(1,1)+FPQy*PBC%BoxShape%D(1,2)+FPQz*PBC%BoxShape%D(1,3)"];
           WS["      PQy=FPQy*PBC%BoxShape%D(2,2)+FPQz*PBC%BoxShape%D(2,3)"];
           WS["      PQz=FPQz*PBC%BoxShape%D(3,3)"];
           WS["      ! End MIC"];
           WS["      WPx = -Eta*PQx*r1xZpE"];
           WS["      WPy = -Eta*PQy*r1xZpE"];
           WS["      WPz = -Eta*PQz*r1xZpE"];
           WS["      WQx = Zeta*PQx*r1xZpE"];
           WS["      WQy = Zeta*PQy*r1xZpE"];
           WS["      WQz = Zeta*PQz*r1xZpE"];
           WS["      T=Omega*(PQx*PQx+PQy*PQy+PQz*PQz)"];
           PunchGammas[Subroutine,LBra+LKet+1];
];
