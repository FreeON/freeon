trace[m_List]:=Sum[ m[[i,i]],{i,Length[m]}]

ToAU=1.88972598857892;

it=2; n=2;
TwoDel=1/(2 10^(-3)*ToAU);

<<h2_p_out;
Pp=P[it];Kp=K[it];Jp=J[it];Sp=S;Vp=V;Tp=T;NukEp=NukE[it];
<<h2_p_ints;
Do[Do[Do[Do[ip[i,j,k,l]=Int[i,j,k,l],{i,n}],{j,n}],{k,n}],{l,n}];

<<h2_m_out;
Pm=P[it];Km=K[it];Jm=J[it];Sm=S;Vm=V;Tm=T;NukEm=NukE[it];
<<h2_m_ints;
Do[Do[Do[Do[im[i,j,k,l]=Int[i,j,k,l],{i,n}],{j,n}],{k,n}],{l,n}];

<<h2_0_out;
Po=P[it];Ko=K[it];Jo=J[it];So=S;Vo=V;Fo=F[it];To=T;NukEo=NukE[it];
<<h2_0_ints;
Do[Do[Do[Do[io[i,j,k,l]=Int[i,j,k,l],{i,n}],{j,n}],{k,n}],{l,n}];

Ep=N[trace[Pp.(2*Tp+Jp+Kp)],20]+NukEp;
Em=N[trace[Pm.(2*Tm+Jm+Km)],20]+NukEm;
Eo=N[trace[Po.(2*To+Jo+Ko)],20]+NukEo;

Print[" Eo = ",Eo];

dS=(Sp-Sm)*TwoDel;
dT=(Tp-Tm)*TwoDel;
dV=(Vp-Vm)*TwoDel;
dN=(NukEp-NukEm)*TwoDel;
dP=(Pp-Pm)*TwoDel;
dJ=(Jp-Jm)*TwoDel;
dK=(Kp-Km)*TwoDel;
dE=(Ep-Em)*TwoDel;
Do[Do[Do[Do[ dI[i,j,k,l]=(ip[i,j,k,l]-im[i,j,k,l])*TwoDel,{l,n}],{k,n}],{j,n}],{i,n}];
DJ=2*Table[Sum[ Po[[k,l]] dI[i,j,k,l] ,{l,n},{k,n}],{j,n},{i,n}];
DK=-Table[Sum[  Po[[k,l]] dI[i,k,j,l] ,{l,n},{k,n}],{j,n},{i,n}];

Print[" dE = ",dE," \n "];

dE1=N[trace[dP.(2*To+Jo+Ko)+Po.(2*dT+dJ+dK)],20]+dN;

Print[" dE1 = ",dE1," \n "];

dE2=N[trace[ -2 dS.Po.Fo.Po + Po.(2*dT+2 dV+DJ+DK)],20]+dN;

Print[" dE2 = ",dE2," \n "];

Print[" dS.W = ",trace[-2 dS.Po.Fo.Po]," \n"];

Print[" 2 dT.P = ",trace[dT.Po]," \n"];

Print[" (2 dJ+dK).P = ",trace[((2)dJ+dK).Po]]; 

Print[" (dK).P = ",trace[(dK).Po]]; 

Print[" dJ = ",MatrixForm[dJ]];

Print[" dK = ",MatrixForm[dK]];