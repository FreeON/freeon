(* :Title: ReIm *)

(* :Authors: Roman Maeder and Martin Buchholz *)

(* :Summary: 	Simplifications with Re[] and Im[]
                                                                      
  The basic idea is to allow variables to be assumed real.            
  To declare a variable x real, use                                   
                                                                      
  	x/: Im[x] = 0                                                   
                                                                      
  There are also rules for Arg[], Abs[], and Conjugate[].             
  Sometimes more simplifications are possible if variables are        
  known to be positive or negative (e.g., with Log[x]).                
  Such declarations can be made with                                 
                                                                      
  	p /: Positive[p] = True                                         
                                                                      
  There can be problems with multivalued functions, e.g., with
                                                                      
  			n /: Negative[n] = True                                        
  we get 
  			Im[Log[-n]] --> 2Pi                                         
                                                                      
  Functions can also be declared real-valued, meaning that f[x]       
  is real if x is real. To declare a function to be real valued, use                        
  	RealValued[f, g,...]                                            
                                                                      
  Reading in this file makes such definitions for some of the         
  built-in functions, for example the trigonometric ones.             
                                                                      
  See also ComplexExpand[].
*)

(* :Context: Algebra`ReIm` *)

(* :Package Version: 1.2 *)

(* :Copyright: 1990 Wolfram Research, Inc.  *)

(* :History: 
	Revised by Roman Maeder, October, 1990.
	Version 1.1 by Roman Maeder and Martin Buchholz, August, 1989.
*)

(* :Keywords: *)

(* :Source: none. *)

(* :Warning: Re, Im, Abs, Conjugate, and Arg redefined. *)

(* :Mathematica Version: 2.0 *)

(* :Limitation: *)

(* :Discussion: *)


BeginPackage["Algebra`ReIm`"]

ReIm::usage = "ReIm.m defines simplification rules for Re[], Im[], etc.
	Variables can be declared real by var/: Im[var] = 0."

RealValued::usage = "RealValued[f] declares f to be a real-valued function
	(for real-valued arguments)."

Begin["`Private`"]

protected = Unprotect[Re, Im, Abs, Conjugate, Arg]

(* fundamental rules, Im[x]==0 serves as our test for 'reality' *)

Re[x_] := x  /; Im[x] == 0
Arg[x_] := 0 /; Positive[x]
Arg[x_] :=Pi /; Negative[x]
Conjugate[x_] :=  x /; Im[x] == 0
Conjugate[x_] := -x /; Re[x] == 0

(* there must not be a rule for Im[x] in terms of Re[x] !! *)

(* things known to be real *)

Im[Re[_]] := 0
Im[Im[_]] := 0
Im[Abs[_]] := 0
Im[Arg[_]] := 0
Im[x_?Positive] = 0
Im[x_?Negative] = 0

Im[x_ ^ y_] := 0 /; Positive[x] && Im[y] == 0

Im[Log[r_?Positive]] := 0

(* arithmetic *)

Re[x_Plus] := Re /@ x
Im[x_Plus] := Im /@ x

Re[x_ y_Plus] := Re[Expand[x y]]
Im[x_ y_Plus] := Im[Expand[x y]]

Re[x_ y_] := Re[x] Re[y] - Im[x] Im[y]
Im[x_ y_] := Re[x] Im[y] + Im[x] Re[y]

(* hidden products *)
Re[(x_ y_)^k_] := Re[x^k y^k]
Im[(x_ y_)^k_] := Im[x^k y^k]

(* nested powers *)
Re[(x_^y_)^k_] := Re[x^(y k)]
Im[(x_^y_)^k_] := Im[x^(y k)]

Re[ 1/x_ ] :=  Re[x] / (Re[x]^2 + Im[x]^2)
Im[ 1/x_ ] := -Im[x] / (Re[x]^2 + Im[x]^2)

Im[x_^2] := 2 Re[x] Im[x]

Re[ x_^n_Integer ] :=
	Block[{a, b},
		a = Round[n/2]; b = n-a;
		Re[x^a] Re[x^b] - Im[x^a] Im[x^b]
	]

Im[ x_^n_Integer ] :=
	Block[{a, b},
		a = Round[n/2]; b = n-a;
		Re[x^a] Im[x^b] + Im[x^a] Re[x^b]
	]

Re[x_Integer^n_Rational] := 0                /; IntegerQ[2n] && Negative[x]
Im[x_Integer^n_Rational] := 
	(-x)^n (-1)^((Numerator[n]-1)/2)     /; IntegerQ[2n] && Negative[x]

(* functions *)

Re[Log[r_?Negative]] := Log[-r]
Im[Log[r_?Negative]] := Pi
Re[Log[z_]] := Log[Abs[z]] /; Im[z] == 0
Re[Log[z_]] := (1/2) Log[Re[z]^2 + Im[z]^2]
Im[Log[z_]] := Arg[z]

Re[Log[a_ b_]] := Re[Log[a] + Log[b]]
Im[Log[a_ b_]] := Im[Log[a] + Log[b]]
Re[Log[a_^c_]] := Re[c Log[a]]
Im[Log[a_^c_]] := Im[c Log[a]]

Re[E^x_] := Cos[Im[x]] Exp[Re[x]]
Im[E^x_] := Sin[Im[x]] Exp[Re[x]]

Re[Sin[x_]] := Sin[Re[x]] Cosh[Im[x]]
Im[Sin[x_]] := Cos[Re[x]] Sinh[Im[x]]

Re[Cos[x_]] :=  Cos[Re[x]] Cosh[Im[x]]
Im[Cos[x_]] := -Sin[Re[x]] Sinh[Im[x]]

Re[Sinh[x_]] := Sinh[Re[x]] Cos[Im[x]]
Im[Sinh[x_]] := Cosh[Re[x]] Sin[Im[x]]

Re[Cosh[x_]] := Cosh[Re[x]] Cos[Im[x]]
Im[Cosh[x_]] := Sinh[Re[x]] Sin[Im[x]]

(* conjugates *)

Re[Conjugate[z_]] :=  Re[z]
Im[Conjugate[z_]] := -Im[z]

Conjugate[x_Plus]:= Conjugate /@ x
Conjugate[x_Times]:= Conjugate /@ x
Conjugate[x_^n_Integer]:= Conjugate[x]^n
Conjugate[Conjugate[x_]]:= x

(* real-valued rules *)

Attributes[RealValued] = {Listable, HoldAll}
Attributes[RealValuedQ] = {HoldFirst}

RealValued[f_Symbol] := (f/: RealValuedQ[f] = True; f)
RealValued[f__] := RealValued /@ {f}

Im[ (_?RealValuedQ)[_?(Im[#]==0&)...] ] := 0

(* define built-in function to be real-valued *)

DoRules[flist_] :=
	Block[{protected},
		protected = Unprotect[flist];
		RealValued[flist];
		Protect[Evaluate[protected]]
	]

DoRules[{Sin, Cos, Tan, ArcSin, ArcCos, ArcTan, ArcCot,
	 Sinh, Cosh, Tanh, ArcSinh, ArcCosh, ArcTanh,
	 Floor, Ceiling, Round, Sign, Factorial}]

Protect[Evaluate[protected]]

End[]

Protect[RealValued]

EndPackage[]
