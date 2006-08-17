
SetOptions[FindRoot, WorkingPrecision -> 50, AccuracyGoal -> 20, MaxIterations -> 40]; 

Off[FindRoot::precw];
Unprotect[D, Derivative, Function];
ClearAll[F];
SetAttributes[F, Listable];
F[m_, T_] := Gamma[1/2 + m, 0, T]/(2*T^(1/2 + m)) /; T != 0
F[m_, T_] := 1/(2 m + 1) /; T == 0;
D[F[m_, x_], x_] := -F[m + 1, x];
Derivative[0, 1][F][m_, x_] := -F[m + 1, x];
Derivative[0, n_][F][m_, x_] := (-1)^n F[m + n, x];
Protect[D, Derivative, Function];

GammaInTheLimit[n_, T_] := Abs[2 n - 1]!!/(2 (2 T)^n) Sqrt[Pi/T];

DTest[M_,T_, WP_] := N[Abs[F[M, T] - GammaInTheLimit[M, T]], WP];

TInterval[M_,Err_] := x /. FindRoot[DTest[M,x, WP] == Err, {x, {1, 50}}]

ChebyshevRules[k_,x_] := { 
      a_Real x^m_ :> 
      a Pi*m!/(2^m ((m-k)/2)! (k+(m-k)/2)!) /;(m>=k&&Mod[m-k,2]==0),
      a_Real x :> 
      a Pi*1!/(2^1 ((1-k)/2)! (k+(1-k)/2)!) /;(1>=k&&Mod[1-k,2]==0),
      a_Real :>a Pi /;(k==0),
      a_Real x^m_ :> 0 /;Not[m>=k&&Mod[m-k,2]==0] ,
      a_Real x  :> 0 /;Not[1>=k&&Mod[1-k,2]==0] ,
      a_Real :>0 /;(k=!=0)
                       };

FF[x_]:=ToString[FixedNumberForm[x,16,2]];
