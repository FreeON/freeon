(*:Name: NumericalMath`Approximations` *)

(*:Title: Approximation of Functions *)

(*:Author: Jerry B. Keiper *)

(*:Summary:
This package provides tools for finding a rational approximation to
a differentiable function.  The approximation may be an interpolation
between specified abscissas or a minimax approximation over an interval.
The function may be specified explicitly or parametrically.
*)

(*:Context: NumericalMath`Approximations` *)

(*:Package Version: 2.0 *)

(* :Copyright: Copyright 1990-1994,  Wolfram Research, Inc.
*)

(*:History:
	Version 1.2 by Jerry B. Keiper, January 1989.
	Updated to 2.0 by Jerry B. Keiper, November 1990.
*)

(*:Keywords:
	functional approximation, Chebyshev approximation,
	rational approximation, minimax approximation
*)

(*Sources:
	Carl-Erik Froberg, Numerical Mathematics: Theory and Computer
		Applications, Benjamin/Cummings, 1985, pp. 250-266

	A. Ralston & P. Rabinowitz, A First Course in Numerical Analysis
		(2nd. ed.), McGraw-Hill, New York, 1978
*)

(*:Mathematica Version: 2.0 *)

(*:Limitation:
	Each of the tools will fail if you attempt to find an
	approximation to an even or odd function with numerator and
	denominator degrees that are impossible to achieve. For
	example, RationalInterpolation[Cos[x], {x, 3, 2}, {x, -1, 1},
	WorkingPrecision -> 30] or RationalInterpolation[Sin[x],
	{x, 3, 3}, {x, -1, 1}].  The denominator of such symmetric,
	rational interpolations should be even, and the numerator
	should be even or odd, as the function is even or odd.

	Each of the tools can be very sensitive to the form of the
	function being approximated.  For example, if the function
	has an irrational singularity (algebraic or transcendental),
	MiniMaxApproximation[ ] will be rather ineffective.  However,
	if you either subtract out or divide out such a singularity,
	MiniMaxApproximation[ ] will be much more effective.  Thus
	on the interval [1, 2], Integrate[1/Log[t], {t, x, 2}]
	should be approximated, as the approximation to (Log[x-1] +
	Integrate[1/Log[t], {t, x, 2}]) with Log[x-1] subtracted
	back off.

	Can occasionally get lost in its search for a minimax
	approximation.  With care and experience such problems can
	usually be dealt with by solving a sequence of approximations
	with parameters varying from "good" values to "bad" values.
	*)

(*:Discussion:
	This package provides the following tools for approximating
	differentiable functions:
	1. RationalInterpolation--finds a rational approximation to
	   a function that has no error at specified abscissas.  The
	   abscissas may be specified explicitly or implicitly.
	2. MiniMaxApproximation--finds a rational approximation to
	   a function.  The maximum of the relative error in the
	   approximation over the interval has the smallest possible value.
	3. GeneralRationalInterpolation--the same as RationalInterpolation 
	   except that the function may be specified parametrically.  This
	   simplifies the finding of approximations to inverse functions
	   that are difficult to specify explicitly (e.g., the inverse
	   of Erf[x]).
	4. GeneralMiniMaxApproximation--the same as MiniMaxApproximation,
	   except that the function may be specified parametrically, and the
	   error being minimized need not be the relative error.

RationalInterpolation and GeneralRationalInterpolation are for
interpolation.  MiniMaxApproximation and GeneralMiniMaxApproximation
are for approximating a function on an interval.  The "General" in
GeneralRationalInterpolation and GeneralMiniMaxApproximation refers
to the fact that the function need not be specified explicitly
and/or that the error being minimized is not required to be the
relative error.  Each of these functions has the option WorkingPrecision,
which specifies the number of decimal digits to use in calculating
the approximation.

MiniMaxApproximation uses RationalInterpolation to get an initial
approximation, which is then improved.  RationalInterpolation (in
the case of implicitly specified abscissas) and MiniMaxApproximation
both use the option Bias, which is a number between -1 and 1 and
which causes the abscissas to be chosen with a bias to the left or
right of the interval.  Bias -> 0 causes the abscissas to be chosen
symmetrically.

MiniMaxApproximation has several options in addition to WorkingPrecision and
Bias:
	Brake - {n1, n2} where both n1 and n2 are non-negative integers.
		MiniMaxApproximation tries to follow a "path" from the initial
		approximation provided by RationalInterpolation to the final
		answer.  This "path" is not well marked, and if the approx-
		imation changes too much from one iteration to the next,
		MiniMaxApproximation can lose its way.  Brake provides a way to
		restrict the changes.  n1 specifies how many of the changes are
		to be restricted and n2 specifies how much of a restriction is
		to be placed on the first change.  The restrictions decrease
		quadratically to no restriction at the n1+1-th iteration.
	MaxIterations - the maximum number of iterations allowed after the
		changes are no longer restricted by Brake.
	Derivatives - {func, D[func, x], D[func, x, 2]}.  This option may be
		left Automatic if Mathematica can find the derivatives itself.
		However, it may be more efficient to specify a function that
		returns the required list of values for any x in the interval.
	PrintFlag - if MiniMaxApproximation cannot be made to work (or the
		user would simply like to watch the convergence) setting
		PrintFlag -> True will cause certain information to be printed
		as it iterates.  More precisely, it prints the ordered pairs
		consisting of the abscissas of the extrema of the error and
		the value of the error at those extrema.
	PlotFlag - setting this to True causes a plot of the error to be
		made after each iteration.

It should be noted that in minimizing the maximum of the error,
MiniMaxApproximation works with the relative error, i.e., the maximum of
Abs[1 - approx[x]/func[x]] is minimized.  It will not normally be possible to
approximate a function with a zero in the interval, since the relative error
near the zero will approach infinity. More general forms of the error to be
minimized require the use of GeneralMiniMaxApproximation.

GeneralRationalInterpolation and GeneralMiniMaxApproximation differ from
RationalInterpolation and MiniMaxApproximation in that the function is
specified parametrically. For example, to approximate the inverse of the Gamma
function on some interval, we might specify the function as {Gamma[t], t}
or {Log[Gamma[t]], t}.  The latter is probably better, since the derivatives
of Log[Gamma[t]] are simpler to compute and the domain for which the
approximation is to be valid will be much smaller.  Another difference is
that the error that is minimized is 1 - approx[x[t]]/y[t], or as an 
alternative (y[t] - approx[x[t]])/g[t], where g is a third function
specified by the user.
*)

(*:Examples:
ri1 = RationalInterpolation[Exp[x],{x,2,2},{-1,-1/2,0,1/2,1}]
Plot[ri1 - Exp[x],{x,-2,2}]

ri2 = RationalInterpolation[Exp[x],{x,2,2},{x,-1,1}]
Plot[ri2 - Exp[x],{x,-2,2}]

(* Note the effect of changing the option Bias. *)
ri3 = RationalInterpolation[Exp[x],{x,2,2},{x,-1,1}, Bias -> .2]
Plot[ri3 - Exp[x],{x,-2,2}]

ri4 = RationalInterpolation[Exp[x],{x,2,2},{x,-1,1}, Bias -> -.7]
Plot[ri4 - Exp[x],{x,-2,2}]

mma1 = MiniMaxApproximation[Exp[x],{x,{-1,1},4,1}][[2,1]]
Plot[1 - mma1/Exp[x],{x,-1,1}]

mma2 = MiniMaxApproximation[Exp[x],{x,{-1,1},4,1},
	PlotFlag -> True, PrintFlag -> True]

mma3 = MiniMaxApproximation[Exp[x],{x,{-1,1},4,1},
	Bias -> -.09, Brake -> {0,0}, PrintFlag -> True]

(* If it cannot find the derivative of the function it fails. *)
mma4 = MiniMaxApproximation[Exp[Abs[x]],{x,{1,3},2,2}]

(* But we know how to find the derivative and we can tell it how. *)
mma5 = MiniMaxApproximation[Exp[Abs[x]],{x,{1,3},2,2},
	Derivatives -> Module[{exp = Exp[x]}, {exp, exp, exp}]]

gr1 = GeneralRationalInterpolation[{Exp[t],t},{t,2,3},x,{-1,-1/2,0,1/2,1,3/2}]
Plot[gr1 - Log[x],{x,Exp[-3/2],Exp[2]}]

gr2 = GeneralRationalInterpolation[{Exp[t],t},{t,2,3},x,{t,-1,3/2}]
Plot[gr2 - Log[x],{x,Exp[-3/2],Exp[2]}]

gr3 = GeneralRationalInterpolation[{Exp[t],t},{t,2,3},x,{t,-1,3/2},Bias->.4]
Plot[gr3 - Log[x],{x,Exp[-3/2],Exp[2]}]

gr4 = GeneralRationalInterpolation[{Exp[t],t},{t,2,3},x,{t,-1,3/2},Bias->-.3]
Plot[gr4 - Log[x],{x,Exp[-3/2],Exp[2]}]

(* If the weight function vanishes on the interval, there is a problem. *)
gmma1 = GeneralMiniMaxApproximation[{Exp[t],t},{t,{-1,1},3,0},x]

gmma2 = GeneralMiniMaxApproximation[{Exp[t],t},{t,{1,3},3,0},x][[2,1]]
Plot[1 - gmma2/Log[x],{x,Exp[1],Exp[3]}]

(* This is how we get it to minimize the maximum absolute error. *)
gmma3 = GeneralMiniMaxApproximation[{Exp[t],t,1},{t,{1,3},3,0},x][[2,1]]
Plot[gmma3 - Log[x],{x,Exp[1],Exp[3]}]

Finally we give a more involved use of MiniMaxApproximation.  We desire an
approximation to the complementary error function erfc[x] == Erf[x,Infinity]
on the interval 1 < x < Infinity.

We know that 

		Sqrt[Pi] x Exp[x^2] erfc[x] ~ 1

(cf. e.g., Abramowitz and Stegun, Handbook of Mathematical Functions, Dover,
eqn 7.1.23.)

If we can approximate the function x Exp[x^2] erfc[x] we can easily solve
for erfc[x].  (Note that the numerator and the denominator must be of the
same degree, since as x -> Infinity our function approaches neither 0 nor
Infinity.)  We will not be able to go all the way to Infinity; in fact we
can only go to a little more than 25 before MiniMaxApproximation complains
that it has a zero in the denominator of its approximation:

prec = 65;
sqrtpi = N[Sqrt[Pi],65];
f =  x  Erf[x,Infinity] Exp[x^2];
ffders = Module[{erf = Erf[x,Infinity], exp = Exp[x^2]},
	{x erf exp,erf exp (1+2x^2) - 2x/sqrtpi,
	-4(x^2+1)/sqrtpi + (4x^3+6x) erf exp}];
era = MiniMaxApproximation[f, {x,{1,30},3,3}, 
	Derivatives->ffders, WorkingPrecision->65, Brake->{40,30}, Bias->-.5];

MiniMaxApproximation::zeroden: The denominator of the approximation ... has a
	zero at 1.00456.

(* The interval is too big; we shorten it somewhat. *)

era = MiniMaxApproximation[f, {x,{1,25},3,3}, 
	Derivatives->ffders, WorkingPrecision->65, Brake->{40,30}, Bias->-.5];

(* Does find an approximation for the interval 1 < x < 25.  We now proceed to
stretch the interval to 1 < x < 70. *)

stretch[x1_] := (era = MiniMaxApproximation[f, era, {x,{1,x1},3,3},
			Derivatives->ffders, WorkingPrecision->25]);

stretch[27];
stretch[30];
stretch[33];
stretch[37];
stretch[42];
stretch[58];
stretch[55];
stretch[63];
stretch[70]

gives the approximation

x Exp[x^2] erfc[x] = (-0.034999362734701385 + 1.18540580001459306*x +
		1.07709037653126903*x^2 + 0.65975310784329749*x^3)/
	(1 + 2.6734238081308433*x + 1.9096936117306481*x^2 +
			1.1693737552523284*x^3)

The relative error on the interval 1 < x < 70 is less than .0000016.  On the
interval 70 < x < Infinity the relative error can get as large as .000007.
*)

BeginPackage["NumericalMath`Approximations`"]

RationalInterpolation::usage = 
"RationalInterpolation[func, {x, m, k}, {x1, x2, ..., xmk1}, (opts)] gives the
rational interpolant to func (a function of the variable x), where m and
k are the degrees of the numerator and denominator, respectively, and
{x1, x2, ..., xmk1} is a list of m+k+1 abscissas of the interpolation points.
An alternative form is \n
RationalInterpolation[func, {x, m, k}, {x, x0, x1}, (opts)], which
specifies the list of abscissas implicitly: the abscissas come from the
interval (x0,x1)."

Options[RationalInterpolation] = {WorkingPrecision -> Precision[N[1]],
				  Bias -> 0}

MiniMaxApproximation::usage = 
"MiniMaxApproximation[func, {x, {x0, x1}, m, k}, (opts)] finds the mini-max
approximation to func (a function of the variable x) on the interval (x0, x1),
where m and k are the degrees of the numerator and denominator, respectively.
The answer returned is {AbscissaList, {Approximation, MaxError}}, where
AbscissaList is a list of the abscissas where the maximum error occurs,
Approximation is the rational approximation desired, and MaxError is the
value of the mini-max error.\n
MiniMaxApproximation[f, approx, {x, {x0, x1}, m, k}, (opts)] is a form that
allows the user to start the iteration from a known approximation.  Here
approx must be in the form of an answer returned by MiniMaxApproximation."

Options[MiniMaxApproximation] = {
			Bias -> 0,
			Brake -> {5, 5},
			Derivatives -> Automatic,
			MaxIterations -> 20,
			WorkingPrecision -> Precision[N[1]],
			PrintFlag -> False,
			PlotFlag -> False
			}

GeneralRationalInterpolation::usage = 
"GeneralRationalInterpolation[{x, y}, {t, m, k}, xvar, {t1, t2, ..., tmk1}],
(opts)] gives the rational interpolant to the function y[x] on the interval
t0 < t < t1, where x and y are functions of t, which parametrically define y[x],
m and k are the degrees of the numerator and denominator, respectively, and xvar
is the variable to be used for x in the approximation, and {t1, t2, ..., tmk1}
is a list of m+k+1 abscissas (values of t) of the interpolation points.
An alternative form is \n
GeneralRationalInterpolation[{x, y}, {t, m, k}, xvar, {t, t0, t1}, (opts)],
which specifies the list of abscissas implicitly: the abscissas come from the
interval t0 < t < t1."

Options[GeneralRationalInterpolation] = {WorkingPrecision -> Precision[N[1]],
				  Bias -> 0}

GeneralMiniMaxApproximation::usage = 
"GeneralMiniMaxApproximation[{x, y, (g)}, {t, {t0, t1}, m, k}, xvar, (opts)]
finds the mini-max approximation to the function y[x] on the interval
t0 < t < t1 where x and y are functions of t, which parametrically define y[x],
m and k are the degrees of the numerator and denominator, respectively, and xvar
is the variable to be used for x in the approximation.  The optional g is also
a function of t, and the error to be minimized is (y[x[t]] - approx[x[t]])/g[t].
The answer returned is {AbscissaList, {Approximation, MaxError}}, where
AbscissaList is a list of the abscissas (values of t) where the maximum error
occurs. Approximation is the rational approximation desired, and MaxError is
the value of the mini-max error.\n
GeneralMiniMaxApproximation[{x, y, (g)}, approx, {t, {t0, t1}, m, k}, xvar,
(opts)] is a form that allows the user to start the iteration from a known
approximation.  Here approx must be in the form of an answer returned by
GeneralMiniMaxApproximation."

Options[GeneralMiniMaxApproximation] = {
			Bias -> 0,
			Brake -> {5, 5},
			Derivatives -> Automatic,
			MaxIterations -> 20,
			WorkingPrecision -> Precision[N[1]],
			PrintFlag -> False,
			PlotFlag -> False
			}

PrintFlag::usage = "PrintFlag is an option to MiniMaxApproximation and
GeneralMiniMaxApproximation that specifies whether data from the successive
iterates in the approximation algorithm are to be shown.  Its primary use
is to help understand what went wrong when the algorithm fails."

PlotFlag::usage = "PlotFlag is an option to MiniMaxApproximation and
GeneralMiniMaxApproximation that specifies whether plots of the successive
iterates in the approximation algorithm are to be drawn.  Its primary use
is to help understand what went wrong when the algorithm fails."

Bias::usage = "Bias is an option to RationalInterpolation,
GeneralRationalInterpolation, MiniMaxApproximation, and
GeneralMiniMaxApproximation. It is a number between -1 and 1
and causes the abscissas to be chosen with a bias to the left
or right of the interval.  Bias -> 0 causes the abscissas to
be chosen symmetrically."

Brake::usage = "Brake is an option to MiniMaxApproximation and
GeneralMiniMaxApproximation.  A valid choice for Brake is a list of two
nonnegative integers that control how the changes from one iteration to the
next are to be restricted.  The first integer indicates the number of
iterations that are to be restricted and the second indicates the magnitude
of the first restriction.  The restrictions decrease to 0 as the algorithm 
proceeds." 

Derivatives::usage = "Derivatives is an option to MiniMaxApproximation and
GeneralMiniMaxApproximation.  It can be Automatic or an expression that
evaluates to a list containing the function, the first derivative, and the
second derivative evaluated at the variable."
 
Begin["NumericalMath`Approximations`Private`"]

RationalInterpolation[func_, {x_, m_Integer, k_Integer}, xlist_, opts___] :=
    Module[{xinfo, bias, answer, biasOK = True,
	prec = WorkingPrecision /.{opts} /.Options[RationalInterpolation]},
	answer /;
	    (If[Head[xlist[[1]]] === Symbol,
		xinfo = {x, SetPrecision[N[{xlist[[2]], xlist[[3]]},
					prec], prec], m, k};
		bias = Bias /. {opts} /. Options[RationalInterpolation];
		bias = SetPrecision[N[bias, prec],prec];
		biasOK = (NumberQ[bias] && bias > -1 && bias < 1),
              (* else *)
		xinfo = {x, m, k};
		bias = SetPrecision[N[xlist,prec], prec]
	    ];
	    (biasOK &&
		(answer = RI[func, xinfo, prec, bias]; answer =!= Fail)))
    ];

MiniMaxApproximation[f_, {x_,{x0_, x1_}, m_Integer, k_Integer}, opts___] :=
    Module[{bias = Bias /. {opts} /. Options[MiniMaxApproximation],
        brake = Brake /. {opts} /. Options[MiniMaxApproximation],
        ffders = Derivatives /. {opts} /. Options[MiniMaxApproximation],
        maxit = MaxIterations /. {opts} /. Options[MiniMaxApproximation],
        prec = WorkingPrecision  /. {opts} /. Options[MiniMaxApproximation],
        sflag = PrintFlag /. {opts} /. Options[MiniMaxApproximation],
        pflag = PlotFlag /. {opts} /. Options[MiniMaxApproximation],
	answer},
	answer /; (If[SameQ[ffders,Automatic],
        	        ffders = D[f,x];
                	ffders = {f, ffders, D[ffders,x]}
                	];
        	bias = SetPrecision[N[bias,prec], prec];
        	answer = MMA[f, ffders, {x, {x0, x1, bias}, m, k},
			prec, maxit, brake, sflag, pflag];
		Fail =!= answer)
        ];

MiniMaxApproximation[f_, {xers_, erars_},
                        {x_,{x0_, x1_}, m_Integer, k_Integer}, opts___] :=
    Module[{brake = Brake /. {opts} /. Options[MiniMaxApproximation],
        ffders = Derivatives /. {opts} /. Options[MiniMaxApproximation],
        maxit = MaxIterations /. {opts} /. Options[MiniMaxApproximation],
        prec = WorkingPrecision  /. {opts} /. Options[MiniMaxApproximation],
        sflag = PrintFlag /. {opts} /. Options[MiniMaxApproximation],
        pflag = PlotFlag /. {opts} /. Options[MiniMaxApproximation],
	answer, mxers},
	answer /; (If[SameQ[ffders,Automatic],
                	ffders = D[f,x];
                	ffders = {f, ffders, D[ffders,x]}
                	];
		mxers = xers;
		mxers[[1]] = x0;
		mxers[[-1]] = x1;
        	answer = MMA[f, ffders, {mxers, erars}, {x, {x0, x1}, m, k},
			prec, maxit, brake, sflag, pflag];
		Fail =!= answer)
        ];

GeneralRationalInterpolation[flist_, {t_, m_Integer, k_Integer},
						xvar_, tlist_, opts___] :=
    Module[{tinfo, bias, answer,biasOK = True,
		prec = WorkingPrecision /. {opts} /.
				Options[GeneralRationalInterpolation]},
	answer /;
	    (If[SameQ[Head[tlist[[1]]],Symbol],
		tinfo = {t, SetPrecision[N[{tlist[[2]], tlist[[3]]},
					prec], prec], m, k};
		bias = Bias /. {opts} /. Options[GeneralRationalInterpolation];
		bias = SetPrecision[N[bias, prec], prec];
		biasOK = (NumberQ[bias] && bias > -1 && bias < 1),
	       (* else *)
		tinfo = {t, m, k};
		bias = SetPrecision[N[tlist, prec], prec];
	    ];
	    (biasOK &&
		(answer = GRI[flist,tinfo,xvar,prec,bias]; Fail =!= answer)))
	];

GeneralMiniMaxApproximation[flist_, {t_, {t0_, t1_}, m_Integer, k_Integer},
							xvar_, opts___] :=
    Module[{j,flen = Length[flist],
	bias = Bias /. {opts} /. Options[GeneralMiniMaxApproximation],
	brake = Brake /. {opts} /. Options[GeneralMiniMaxApproximation],
	ffders = Derivatives /. {opts} /. Options[GeneralMiniMaxApproximation],
	maxit = MaxIterations /. {opts} /.
					Options[GeneralMiniMaxApproximation],
	prec = WorkingPrecision  /. {opts} /.
					Options[GeneralMiniMaxApproximation],
	sflag = PrintFlag /. {opts} /. Options[GeneralMiniMaxApproximation],
	pflag = PlotFlag /. {opts} /. Options[GeneralMiniMaxApproximation],
	answer},
	answer /; (If[SameQ[ffders,Automatic],
			ffders = Array[1,{3,flen}];
			Do[ffders[[1,j]] = flist[[j]], {j,flen}];
			Do[ffders[[2,j]] = D[ffders[[1,j]], t], {j,flen}];
			Do[ffders[[3,j]] = D[ffders[[2,j]], t], {j,flen}];
			ffders = Transpose[ffders]
			];
		bias = SetPrecision[N[bias,prec], prec];
		answer = GMMA[flist, ffders, {t, {t0, t1, bias}, m, k},
			xvar, prec, maxit, brake, sflag, pflag];
		Fail =!= answer)
	];

GeneralMiniMaxApproximation[flist_, {ters_, erars_},
		{t_, {t0_, t1_}, m_Integer, k_Integer}, xvar_, opts___] :=
    Module[{j,flen = Length[flist],
	bias = Bias /. {opts} /. Options[GeneralMiniMaxApproximation],
	brake = Brake /. {opts} /. Options[GeneralMiniMaxApproximation],
	ffders = Derivatives /. {opts} /. Options[GeneralMiniMaxApproximation],
	maxit = MaxIterations /. {opts} /.
					Options[GeneralMiniMaxApproximation],
	prec = WorkingPrecision  /. {opts} /.
					Options[GeneralMiniMaxApproximation],
	sflag = PrintFlag /. {opts} /. Options[GeneralMiniMaxApproximation],
	pflag = PlotFlag /. {opts} /. Options[GeneralMiniMaxApproximation],
	answer, mters},
	answer /; (If[SameQ[ffders,Automatic],
			ffders = Array[1,{3,flen}];
			Do[ffders[[1,j]] = flist[[j]], {j,flen}];
			Do[ffders[[2,j]] = D[ffders[[1,j]], t], {j,flen}];
			Do[ffders[[3,j]] = D[ffders[[2,j]], t], {j,flen}];
			ffders = Transpose[ffders]
			];
		mters = ters;
		mters[[1]] = t0;
		mters[[-1]] = t1;
		answer = GMMA[flist, ffders, {mters, erars},
				{t, {t0, t1}, m, k},
				xvar, prec, maxit, brake, sflag, pflag];
		Fail =!= answer)
	];

RationalInterpolation::precloss = "Drastic loss of precision.  Try starting
	with higher precision."

RationalInterpolation::notnum = "Function is not numerical at the interpolation
	points: `1`."

pa[x_] := If[x == 0, Accuracy[x], Precision[x]];
precacc[l_List] := Min[pa /@ l]

RI[f_, xinfo_, prec_, bias_] :=
    Module[{i, mk1, xx, fx, mat, tempvec, x, x0, x1, m, k},
	If[(IntegerQ[prec] && (prec > 0)) =!= True, Return[Fail]];
	x = xinfo[[1]];
	If[SameQ[Length[xinfo],4],
		(x0 = xinfo[[2,1]];
		x1 = xinfo[[2,2]];
		m = xinfo[[3]];
		k = xinfo[[4]];
		mk1 = m+k+1;
		xx = Reverse[Table[Cos[N[Pi (i+1/2)/(mk1),prec]], {i,0,m+k}]];
		xx += bias (1 - xx xx);
		xx = xx*(x1-x0)/2 + (x1+x0)/2),
	    (* else *)
		(m = xinfo[[2]];
		k = xinfo[[3]];
		mk1 = m+k+1;
		xx = bias)
	];
	xx = SetPrecision[N[xx,prec],prec];
	fx = SetPrecision[N[f /. x->xx, prec], prec];
	If[!Apply[And,Map[NumberQ,fx]],
		Message[RationalInterpolation::notnum, fx];
		Return[Fail]
	];
	mat = Table[1,{i,mk1}];
	tempvec = mat;
	mat[[1]] = tempvec;
	Do[tempvec *= xx;
		If[i <= m,mat[[i+1]] = tempvec,,];
		If[i <= k,mat[[i+m+1]] = -tempvec*fx],{i,Max[m,k]}];
	xx = LinearSolve[Transpose[mat],fx];
	If[Head[xx] === LinearSolve, Return[Fail]];
	If[precacc[xx] < 5,
		Message[RationalInterpolation::precloss];
		Return[Fail]
	];
	xx = SetPrecision[xx, prec];
	(xx[[1]]+Sum[xx[[i+1]] x^i,{i,m}])/(1+Sum[xx[[i+m+1]] x^i,{i,k}])
    ];

MiniMaxApproximation::dervnotnum = "`1` is not numerical at `2` = `3` : `4`.";

MiniMaxApproximation::van = "Failed to locate the extrema in `1` iterations.
	The function `2` may be vanishing on the interval `3` or the
	WorkingPrecision may be insufficient to get convergence.";

LocateExtrema[ffders_, era_, xinfo_, prec_, xguess_, maxit_,sflag_] :=
    Module[{i, ii=1, mk1, xx, xd, fx, x, x0, x1, m, k, j=0,
		rnum, rnumd, rnumdd, rden, rdend, rdendd, dxmax, oldd,
		rnumx, rnum1, rnum2, rdenx, rden1, rden2, elist, r, repeat},
	If[Head[era] === List,
		r = era[[1]]; repeat = False,
		r = era; repeat = True,
	];
	x = xinfo[[1]];
	x0 = xinfo[[2,1]];
	x1 = xinfo[[2,2]];
	m = xinfo[[3]];
	k = xinfo[[4]];
	mk1 = m+k+1;
	If[SameQ[Head[xguess],List],xx=xguess,
		xx=Reverse[Table[Cos[N[Pi i/mk1,prec]],{i,0,mk1}]];
		xx += xguess (1 - xx xx);
		xx = xx*(x1-x0)/2 + (x1+x0)/2
	];
	If[Length[xx] != mk1+1, Return[Fail]];
	xx = SetPrecision[N[xx,prec],prec];
	xd = Take[xx,{2,mk1}];
	dxmax = (rnum1 = xd-Drop[xx,-2];
		 rnum2 = Drop[xx,2]-xd;
		 Map[Min,Transpose[{rnum1,rnum2}]]/6);
	oldd = Array[0&,{m+k}];
	rnum = Numerator[r];
	rnumd = D[rnum,x];
	rnumdd = D[rnumd,x];
	rden = Denominator[r];
	rdend = D[rden,x];
	rdendd = D[rdend,x];
	While[ii < 5+prec/10,
		If[j++ > 2 maxit,
			Message[MiniMaxApproximation::van, maxit,
					ffders[[1]], {x0,x1}];
			Return[Fail]
		];
		fx = SetPrecision[N[ffders /. x->xd,prec],prec];
		If[!Apply[And,Map[NumberQ,Flatten[fx]]],
			Message[MiniMaxApproximation::dervnotnum,
				ffders, x, xd, fx];
			Return[Fail]
		];
		rnumx = SetPrecision[N[rnum /. x->xd,prec],prec];
		rnum1 = SetPrecision[N[rnumd /. x->xd,prec],prec];
		rnum2 = SetPrecision[N[rnumdd /. x->xd,prec],prec];
		rdenx = SetPrecision[N[rden /. x->xd,prec],prec];
		rden1 = SetPrecision[N[rdend /. x->xd,prec],prec];
		rden2 = SetPrecision[N[rdendd /. x->xd,prec],prec];
		r = rnum1 rdenx - rnumx rden1;
		rden2 = (rdenx (rnum2 rdenx - rnumx rden2) - 
							2 r rden1)/(rdenx^3);
		rden1 = r/(rdenx^2);
		rdenx = rnumx/rdenx;
		elist = 1-rdenx/fx[[1]];
		rden1 = (rden1 fx[[1]] - rdenx fx[[2]]);
		rden2 = rden1/(rden2 fx[[1]] - rdenx fx[[3]]);
		rden1 = Sign[rden1] Sign[elist];
		Do[If[Sign[rden1[[i]]] != Sign[rden2[[i]]] ||
					Abs[rden2[[i]]] > dxmax[[i]],
			ii = 1;
			rden2[[i]] = dxmax[[i]] Sign[rden1[[i]]],,],
			{i,m+k}];
		Do[If[rden2[[i]]+oldd[[i]] == 0,
			rden2[[i]] /= 2;
			dxmax[[i]] /= 2],
			{i,m+k}];
		oldd = rden2;
		xd = xd - rden2;
		xd = SetPrecision[xd,prec];
		rden2 = N[Abs[rden2/xd]];
		If[sflag, Print[rden2]];
		If[repeat, If[Max[rden2] > .01, ii = 1, ii *= 2], ii = prec]
	];
	xd = Prepend[xd,xx[[1]]];
	xd = Append[xd,xx[[-1]]];
	xd
    ];

FindCoefficients[f_, era_, xinfo_, prec_, xlist_, brake_] :=
    Module[{mk2, x, m, k, i, aa, bb, suma, sumb, abe, ff, temp1, temp2, jac,
		r = If[SameQ[Head[era],List],era[[1]],era]},
	x = xinfo[[1]];
	m = xinfo[[3]];
	k = xinfo[[4]];
	mk2 = m+k+2;
	suma = Numerator[r] /. x->xlist;
	sumb = Denominator[r] /. x->xlist;
	abe = Join[CoefficientList[Numerator[r],x],
			Drop[CoefficientList[Denominator[r],x],1]];
	abe = Append[abe,If[SameQ[Head[era],List],era[[2]],0]];
	ff = SetPrecision[N[f /. x->xlist, prec],prec];
	If[!Apply[And,Map[NumberQ,ff]],
		Message[RationalInterpolation::notnum, ff];
		Return[Fail]];
	temp1 = -1/(ff sumb);
	temp2 = temp1;
	jac = Range[m+k+2];
	Do[jac[[i]] = temp2; temp2 *= xlist,{i,m+1}];
	temp1 *= -suma;
	temp2 = xlist temp1/sumb;
	Do[jac[[i]] = temp2; temp2 *= xlist,{i,m+2,mk2-1}];
	temp2 = Table[(-1)^i,{i,mk2}];
	jac[[mk2]] = temp2;
	temp2 = Table[1,{i,mk2}] - temp1 + abe[[mk2]] temp2;
	temp2 = SetPrecision[temp2,prec];
	temp2 = LinearSolve[Transpose[jac],temp2]/(brake+1);
	If[Head[temp2] === LinearSolve, Return[Fail]];
	abe = SetPrecision[abe-temp2,prec];
	{Sum[abe[[i]] x^(i-1),{i,m+1}]/
			(1+Sum[abe[[i]] x^(i-m-1),{i,m+2,mk2-1}]),
		abe[[mk2]]}
    ];


MiniMaxApproximation::zeroden = "The denominator of the approximation `1` has
	a zero at `2`.";

MMA[f_, ffders_, {x_,{x0_, x1_, bias_}, m_Integer, k_Integer},
				prec_, maxit_, brake_, sflag_, pflag_] :=
    Module[{i=0, root, xe, era, xinfo={x, {x0, x1}, m, k}},
	If[(IntegerQ[prec] && (prec > 0) && (bias < 1) && (bias > -1)
			&& IntegerQ[maxit] && (maxit > 0)) =!= True,
		Return[Fail]];
	era = RI[f, xinfo, prec, bias];
	If[era === Fail, Return[Fail]];
	If[k > 0,
		xe = NRoots[N[Denominator[era]]==0,x];
		While[i++ < k,
			root = If[k==1, xe[[2]], xe[[i,2]]];
			If[SameQ[Head[root],Real] && root<=x1 && root>=x0,
				Message[MiniMaxApproximation::zeroden,
					era, root];
				Return[Fail]
			]
		]
	];
	xe = LocateExtrema[ffders,era,xinfo,prec,bias,maxit,sflag];
	If[SameQ[xe, Fail], Return[Fail]];
	If[!Apply[And,Map[NumberQ,xe]],
		Message[MiniMaxApproximation::exnotnum, xe];
		Return[Fail]
	];
	MMA[f, ffders, {xe, era}, xinfo, prec, maxit, brake, sflag, pflag]
    ];

MiniMaxApproximation::exnotnum = "The extrema are not numerical: `1`."

MiniMaxApproximation::extalt =
"The extrema of the error do not alternate in sign.  It may be that 
MiniMaxApproximation has lost track of the extrema by going too fast.  If so
try increasing the values in the option Brake.  It may be that the
WorkingPrecision is insufficient.  Otherwise there is an extra extreme value of
the error, and MiniMaxApproximation cannot deal with this problem."

MiniMaxApproximation::extsort =
"MiniMaxApproximation has lost track of the extrema of the error; they are no
longer in sorted order.  Try increasing the values in the option Brake.";

MiniMaxApproximation::conv = "Warning: convergence was not complete.";

MMA[f_, ffders_, {xers_, erars_}, {x_,{x0_, x1_}, m_Integer, k_Integer},
				prec_, maxit_, brake_, sflag_, pflag_] :=
    Module[{xe=SetPrecision[N[xers,prec],prec],
		era=N[erars,prec], olderr=1, count=0, xinfo,
		brvar=brake[[1]], elist, notdone=True, temp, tmp, alpha},
	If[(IntegerQ[prec] && (prec > 0) && IntegerQ[maxit] &&
			(maxit > 0)) =!= True, Return[Fail]
	];
	xinfo = {x, SetPrecision[N[{x0, x1},prec],prec], m, k};
	alpha = If[brake[[2]] == 0, brvar = 0; 1, brake[[1]]^2/brake[[2]]];
	While[notdone,
		elist = SetPrecision[N[f /. x->xe, prec],prec];
		If[SameQ[List,Head[era]], tmp = era[[1]], tmp = era];
		elist = {xe, 1 - (tmp /. x->xe)/elist};
		temp = Sign[elist[[2]]];
		temp = Drop[temp,1]+Drop[temp,-1];
		If[Max[Abs[temp]]>0,
			Message[MiniMaxApproximation::extalt];
			Return[{xe, tmp, elist[[2]]}]
		];
		If[sflag, Print[N[Transpose[elist]]]];
		If[brvar != 0, count = 0; brvar--];
		era = FindCoefficients[f,era,xinfo,prec,xe,brvar^2/alpha];
		If[era === Fail, Return[Fail]];
		If[pflag,Plot[1-era[[1]]/f,{x,x0,x1}]];
		xe=LocateExtrema[ffders,era,xinfo,prec,xe,maxit,sflag];
		If[SameQ[xe, Fail], Return[Fail]];
		If[!Apply[And,Map[NumberQ,xe]],
			Message[MiniMaxApproximation::exnotnum, xe];
			Return[Fail]
		];
		If[!SortedQ[xe],
			Message[MiniMaxApproximation::extsort];
			Return[Fail],,];
		notdone = (brvar != 0) || (Abs[olderr-era[[2]]] >
				(Abs[olderr]+Abs[era[[2]]])/10000);
		If[notdone && (count == maxit),
			Message[MiniMaxApproximation::conv]
		];
		notdone = ((count++ < maxit) && notdone);
		olderr = era[[2]]
	];
	{xe, era}
    ];

GRI[f_, tinfo_, x_, prec_, bias_] :=
    Module[{i, mk1, tt, ft, mat, tempvec, t, t0, t1, m, k},
	If[(IntegerQ[prec] && (prec > 0)) =!= True, Return[Fail]];
	t = tinfo[[1]];
	If[SameQ[Length[tinfo],4],
		(t0 = tinfo[[2,1]];
		t1 = tinfo[[2,2]];
		m = tinfo[[3]];
		k = tinfo[[4]];
		mk1 = m+k+1;
		tt = Reverse[Table[Cos[N[Pi (i+1/2)/(mk1),prec]],{i,0,m+k}]];
		tt += bias (1 - tt tt);
		tt = tt*(t1-t0)/2 + (t1+t0)/2),
	    (* else *)
		(m = tinfo[[2]];
		k = tinfo[[3]];
		mk1 = m+k+1;
		tt = bias)
	];
	tt = SetPrecision[N[tt,prec],prec];
	ft = SetPrecision[N[f /. t->tt, prec], prec];
	If[!Apply[And,Map[NumberQ,Flatten[ft]]],
		Message[RationalInterpolation::notnum, ft];
		Return[Fail]
	];
	tt = ft[[1]];	(* tt is now a list of x-values *)
	ft = ft[[2]];	(* ft is now a list of y-values *)
	mat = Table[1,{i,mk1}];
	tempvec = mat;
	mat[[1]] = tempvec;
	Do[tempvec *= tt;
		If[i <= m,mat[[i+1]] = tempvec];
		If[i <= k,mat[[i+m+1]] = -tempvec*ft],
		{i,Max[m,k]}
	];
	tt = LinearSolve[Transpose[mat],ft];
	If[Head[tt] === LinearSolve, Return[Fail]];
	If[precacc[tt] < 5,
		Message[RationalInterpolation::precloss];
		Return[Fail]
	];
	tt = SetPrecision[tt, prec];
	(tt[[1]]+Sum[tt[[i+1]] x^i,{i,m}])/(1+Sum[tt[[i+m+1]] x^i,{i,k}])
    ];

GeneralMiniMaxApproximation::van = "Failed to locate the extrema in `1`
	iterations.  The weight function `2` may be vanishing on the
	interval `3`, or the WorkingPrecision may be insufficient to
	get convergence.";

GeneralLocateExtrema[ffders_,era_,tinfo_,x_,prec_,tguess_,maxit_,sflag_] :=
    Module[{i, ii = 1, mk1, tt, td, ft, gt, xd, t, t0, t1, m, k,
		rnum, rnumd, rnumdd, rden, rdend, rdendd, dtmax, oldd, j=0,
		rnumt, rnum1, rnum2, rdent, rden1, rden2, elist, r, repeat},
	If[SameQ[Head[era],List],
		r = era[[1]]; repeat = False,
		r = era; repeat = True,
	];
	t = tinfo[[1]];
	t0 = tinfo[[2,1]];
	t1 = tinfo[[2,2]];
	m = tinfo[[3]];
	k = tinfo[[4]];
	mk1 = m+k+1;
	If[SameQ[Head[tguess],List],tt=tguess,
		tt=Reverse[Table[Cos[N[Pi i/mk1,prec]],{i,0,mk1}]];
		tt += tguess (1 - tt tt);
		tt = tt*(t1-t0)/2 + (t1+t0)/2
	];
	If[Length[tt] != mk1+1,Return[Fail]];
	tt = SetPrecision[N[tt,prec],prec];
	td = Take[tt,{2,mk1}];
	dtmax = (rnum1 = td-Drop[tt,-2];
		 rnum2 = Drop[tt,2]-td;
		 Map[Min,Transpose[{rnum1,rnum2}]]/6);
	oldd = Array[0&,{m+k}];
	rnum = Numerator[r];
	rnumd = D[rnum,x];
	rnumdd = D[rnumd,x];
	rden = Denominator[r];
	rdend = D[rden,x];
	rdendd = D[rdend,x];
	While[ii < 5+prec/10,
		If[j++ > 2 maxit,
			Message[GeneralMiniMaxApproximation::van,
				maxit,ffders[[-1,1]], {t0,t1}];
			Return[Fail]
		];
		ft = SetPrecision[N[ffders /. t->td,prec],prec];
		If[!Apply[And,Map[NumberQ,Flatten[ft]]],
			Message[MiniMaxApproximation::dervnotnum,
				ffders, t, td, ft];
			Return[Fail]
		];
		xd = ft[[1]];	(* xd is a list of x val's & der's *)
		gt = ft[[-1]];	(* gt is a list of t val's & der's *)
		ft = ft[[2]];	(* ft is a list of y val's & der's *)
		rnumt = SetPrecision[N[rnum /. x->xd[[1]],prec],prec];
		rnum1 = SetPrecision[N[rnumd /. x->xd[[1]],prec],prec];
		rnum2 = SetPrecision[N[rnumdd /. x->xd[[1]],prec],prec];
		rdent = SetPrecision[N[rden /. x->xd[[1]],prec],prec];
		rden1 = SetPrecision[N[rdend /. x->xd[[1]],prec],prec];
		rden2 = SetPrecision[N[rdendd /. x->xd[[1]],prec],prec];
		r = rnum1 rdent - rnumt rden1;
		rden2 = (rdent(rnum2 rdent - rnumt rden2)-2 r rden1)/(rdent^3);
		rden1 = r/(rdent^2);
		(* now we convert derivatives wrt x to der's wrt t *)
		rden2 = rden2 (xd[[2]])^2 + rden1 xd[[3]];
		rden1 *= xd[[2]];
		rdent = rnumt/rdent;
		elist = (ft[[1]]-rdent)/gt[[1]];
		rden1 = (rden1 gt[[1]] - rdent gt[[2]]) -
				(ft[[2]] gt[[1]] - ft[[1]] gt[[2]]);
		rden2 = rden1/((rden2 gt[[1]] - rdent gt[[3]]) -
				(ft[[3]] gt[[1]] - ft[[1]] gt[[3]]));
		rden1 = Sign[rden1]  Sign[elist];
		Do[If[Sign[rden1[[i]]] != Sign[rden2[[i]]] ||
				Abs[rden2[[i]]] > dtmax[[i]],
			ii = 1;
			rden2[[i]] = dtmax[[i]]*Sign[rden1[[i]]]
			],
			{i,m+k}
		];
		Do[If[rden2[[i]]+oldd[[i]] == 0,
			rden2[[i]] /= 2;
			dtmax[[i]] /= 2],
			{i,m+k}
		];
		oldd = rden2;
		td = td - rden2;
		td = SetPrecision[td,prec];
		rden2 = N[Abs[rden2/td]];
		If[sflag, Print[rden2]];
		If[repeat, If[Max[rden2] > .01, ii = 1, ii *= 2], ii = prec]
	];
	td = Prepend[td,tt[[1]]];
	td = Append[td,tt[[-1]]];
	td
    ];

GeneralFindCoefficients[f_, era_, tinfo_, x_, prec_, tlist_, brake_] :=
    Module[{mk2, t, m, k, i, aa, bb, suma, sumb, abe, ff,
		xlist, temp1, temp2, jac, gt,
		r = If[SameQ[Head[era],List],era[[1]],era]},
	t = tinfo[[1]];
	m = tinfo[[3]];
	k = tinfo[[4]];
	mk2 = m+k+2;
	ff = SetPrecision[N[f /. t->tlist, prec],prec];
	If[!Apply[And,Map[NumberQ,Flatten[ff]]],
		Message[RationalInterpolation::notnum, ff];
		Return[Fail]
	];
	xlist = ff[[1]];   (* xlist are the abscissas of the extrema *)
	gt = ff[[-1]];     (* gt is the normalizing function *)
	ff = If[Length[ff] == 2, 1, ff[[2]]/gt];
				(* ff are the normalized ordinates *)
	abe = Join[CoefficientList[Numerator[r],x],
			Drop[CoefficientList[Denominator[r],x],1]];
	abe = Append[abe,If[SameQ[Head[era],List],era[[2]],0]];
	suma = Numerator[r] /. x->xlist;
	sumb = Denominator[r] /. x->xlist;
	temp1 = -1/(gt sumb);
	If[NumberQ[temp1], temp1 = Table[temp1,{i,m+k+2}]];
	temp2 = temp1;
	jac = Range[m+k+2];
	Do[jac[[i]] = temp2; temp2 *= xlist,{i,m+1}];
	temp1 *= -suma;
	temp2 = xlist temp1/sumb;
	Do[jac[[i]] = temp2; temp2 *= xlist,{i,m+2,mk2-1}];
	temp2 = Table[(-1)^i,{i,mk2}];
	jac[[mk2]] = temp2;
	temp2 = SetPrecision[ff-temp1 + abe[[mk2]] temp2,prec];
	temp2 = LinearSolve[Transpose[jac],temp2]/(brake+1);
	If[Head[temp2] === LinearSolve, Return[Fail]];
	abe = SetPrecision[abe-temp2,prec];
	{Sum[abe[[i]] x^(i-1),{i,m+1}]/
			(1+Sum[abe[[i]] x^(i-m-1),{i,m+2,mk2-1}]),
		abe[[mk2]]}
    ];

GMMA[f_, ffders_, {t_,{t0_, t1_, bias_}, m_Integer, k_Integer}, x_,
				prec_, maxit_, brake_, sflag_, pflag_] :=
    Module[{i=0, root, xe, te, era, tinfo = {t, {t0, t1}, m, k}},
	If[(IntegerQ[prec] && (prec > 0) && (bias < 1) && (bias > -1)
			&& IntegerQ[maxit] && (maxit > 0)) =!= True,
		Return[Fail]
	];
	era = GRI[f, tinfo, x, prec, bias];
	If[era === Fail, Return[Fail]];
	If[k > 0,
		te = NRoots[N[Denominator[era]]==0,x];
		xe = Sort[N[f /. t->{t0,t1}][[1]]];
		While[i++ < k,
			root = If[k==1, te[[2]], te[[i,2]]];
			If[Head[root]===Real && root>=xe[[1]] && root<=xe[[2]],
				Message[MiniMaxApproximation::zeroden,
						era, root];
				Return[Fail]]
			]
		];
	te = GeneralLocateExtrema[ffders,era,tinfo,x,prec,bias,maxit,sflag];
	If[SameQ[te, Fail], Return[Fail]];
	If[!Apply[And,Map[NumberQ,xe]],
		Message[MiniMaxApproximation::exnotnum, xe];
		Return[Fail]
	];
	GMMA[f, ffders, {te, era}, tinfo, x, prec, maxit, brake, sflag, pflag]
    ];

GMMA[f_, ffders_, {ters_, erars_}, {t_,{t0_, t1_}, m_Integer, k_Integer}, x_,
				prec_, maxit_, brake_, sflag_, pflag_] :=
    Module[{te=SetPrecision[N[ters,prec],prec],
		era=N[erars,prec], olderr=1, count=0, tinfo,
		brvar=brake[[1]], elist, notdone=True, temp, tmp, alpha},
	If[(IntegerQ[prec] && (prec > 0) && IntegerQ[maxit] &&
			(maxit > 0)) =!= True, Return[Fail]];
	tinfo = {t, SetPrecision[N[{t0, t1},prec],prec], m, k};
	alpha = If[brake[[2]] == 0, brvar = 0; 1, brake[[1]]^2/brake[[2]]];
	While[notdone,
		elist = SetPrecision[N[f /. t->te,prec],prec];
		If[SameQ[List,Head[era]], tmp = era[[1]], tmp = era];
		elist[[2]] = (elist[[2]] - SetPrecision[
			N[tmp /. x->elist[[1]],prec],prec])/elist[[-1]];
		temp = Sign[elist[[2]]];
		temp = Drop[temp,1]+Drop[temp,-1];
		If[Max[Abs[temp]]>0,
			Message[MiniMaxApproximation::extalt];
			Return[{elist[[1]], tmp, elist[[2]]}]
		];
		If[sflag, Print[N[Sort[Transpose[Take[elist,{1,2}]]]]]];
		If[brvar !=0,count = 0;brvar--];
		era = GeneralFindCoefficients[ f, era, tinfo, x, prec,
						te, brvar^2/alpha];
		If[era === Fail, Return[Fail]];
		If[pflag,
		    ParametricPlot[Evaluate[
			Module[{ft = f},
			    {ft[[1]],(ft[[2]]-(era[[1]] /. x->ft[[1]]))/
							ft[[-1]]}
			]],
			{t,t0,t1}]
		];
		te = GeneralLocateExtrema[
				ffders, era, tinfo, x, prec, te, maxit, sflag];
		If[SameQ[te, Fail], Return[Fail]];
		If[!Apply[And,Map[NumberQ,xe]],
			Message[MiniMaxApproximation::exnotnum, xe];
			Return[Fail]
		];
		If[!SortedQ[te],
			Message[MiniMaxApproximation::extsort];
			Return[Fail]
		];
		notdone = (brvar != 0) || (Abs[olderr-era[[2]]] >
					(Abs[olderr]+Abs[era[[2]]])/10000);
		If[notdone && (count == maxit),
				Message[MiniMaxApproximation::conv]
		];
		notdone = ((count++ < maxit) && notdone);
		olderr = era[[2]]
	];
	{te, era}
    ];

End[]  (* NumericalMath`Approximations`Private` *)

Protect[RationalInterpolation, MiniMaxApproximation, PlotFlag, PrintFlag,
GeneralRationalInterpolation, GeneralMiniMaxApproximation];

EndPackage[] (* NumericalMath`Approximations` *)


