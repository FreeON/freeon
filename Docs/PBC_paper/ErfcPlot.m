OutFile0 = "Plot.dat"
OpenWrite[OutFile0];

bet=0.1;

r1=0.001;
Do[
  a1 = 2.55-0.05*i;
  aa = 10^(a1);
  r1 = x /. FindRoot[Erfc[x]/x==(aa/Sqrt[bet]),{x,r1}, WorkingPrecision -> 20];
  r2 = Re[Sqrt[-Log[aa]]];
  Print[FortranForm[aa],"  ",FortranForm[r1]," ",FortranForm[r2]];
  WriteString[OutFile0,FortranForm[aa],"  ",FortranForm[r1],"  ",FortranForm[r2],"\n"];
,{i,1,251}];
 
Close[OutFile0];
