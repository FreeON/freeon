
(* :Title: Compile.m -- Inlining of functions *)

(* :Context: Utilities`Compile` *)

(* :Author: Markus Lischka (mlischka@physik.tu-muenchen.de) *)

(* :Summary:
   The Compile function is extended to provide a way of "inlining" functions 
   in compiled code.
*)

(* :Copyright: (c) 1997 by Markus Lischka

   Permission is granted to distribute this file for any purpose except 
   for inclusion in commercial software or program collections. 
   This copyright notice must remain intact.
*)

(* :Package Version: 1.3 *)

(* :Mathematica Version: 3.0 or higher *)

(* :History:
   1.0 Initial version (June 1997)
   1.1 Added support for OwnValues (September 1997)
   1.2 Corrected documentation (February 2000)
   1.3 Expanded documentation; initial submission to MathSource (2003-04-06)
*)

(* :Keywords: Compile *)

(* :Sources:
   Roman E. Maeder. Programming in Mathematica, 3rd ed. Addison-Wesley, 1996.
*)

(* :Discussion:
   Substitutes Up-, Down-, N- and OwnValues. 
   Check the compiled code for suspicious opcodes to see whether 
   the substitution was successful. For example, opcode 31 (Mathematica 3.0) 
   or opcode 26 (Mathematica 4.0) are used for function calls to 
   uncompiled code.

   Note that straightforward compilation of Mathematica statements often
   results in poor performance. For example, in both examples given below
   compilation without inlining the function definition results in a 
   significant performance penalty compared to the uncompiled code.
   
   See also
   <http://support.wolfram.com/mathematica/kernel/Symbols/System/Compile.html>.

   Cf. some usage examples at the end of the package.
*)


(* set up the package context, including public imports *)

BeginPackage["Utilities`Compile`"]

(* usage messages for the exported functions and the context itself *)

InlineFunction::usage = "InlineFunction is an option for compile 
that specifies which functions should be inlined before compilation."

Begin["`Private`"]    (* begin the private context (implementation part) *)

(* unprotect any system functions for which definitions will be made *)

protected = Unprotect[ Compile ]

(* definition of auxiliary functions and local (static) variables *)

(* cf. Maeder, Programming in Mathematica, ch. 5.3.2 *)
SetAttributes[WrapHold, HoldAll];
WrapHold[expr_] := Map[Hold, Hold[expr], {2}][[1]]

SetAttributes[values, HoldAll];

values[subs_List] := Flatten[ReleaseHold[{
    Map[Hold[UpValues], WrapHold[subs]],
    Map[Hold[DownValues], WrapHold[subs]],
    Map[Hold[NValues], WrapHold[subs]],
    Map[Hold[OwnValues], WrapHold[subs]]}]]

values[subs_] := Flatten[ReleaseHold[{
    Hold[UpValues][Hold[subs]],
    Hold[DownValues][Hold[subs]],
    Hold[NValues][Hold[subs]],
    Hold[OwnValues][Hold[subs]]}]]

(* definitions for system functions *)

Compile[vars_, expr_, subexpr_, InlineFunction -> subs_] :=
    ReleaseHold[Hold[Compile[vars, expr, subexpr]] //.
        values[subs]]
  
Compile[vars_, expr_, InlineFunction -> subs_] :=
    ReleaseHold[Hold[Compile[vars, expr]] //.
        values[subs]]

Protect[ Evaluate[protected] ]     (* restore protection of system symbols *)

End[ ]         (* end the private context *)

EndPackage[ ]  (* end the package context *)


(* :Examples:

Needs["Utilities`Compile`"];

x0 = Random[];
n = 100000;

f[x_] := 4 x (1 - x);
Timing[x = x0; Do[x = f[x], {n}]; x]

fc1 = Compile[{x0}, Module[{x = x0}, Do[x = f[x], {n}]; x]];
fc1 //FullForm
Timing[fc1[x0]]

fc2 = Compile[{x0}, Module[{x = x0}, Do[x = f[x], {n}]; x],
    InlineFunction -> f];
fc2 //FullForm
Timing[fc2[x0]]

f[x_] := Sin[x];
g[x_] := Cos[x];
Timing[x = x0; Do[x = f[x] + g[x], {n}]; x]

fc1 := Compile[{x0}, Module[{x = x0}, Do[x = f[x] + g[x], {n}]; x]];
fc1 //FullForm
Timing[fc1[x0]]

fc2 := Compile[{x0}, Module[{x = x0}, Do[x = f[x] + g[x], {n}]; x],
    InlineFunction -> {f, g}];
fc2 //FullForm
Timing[fc2[x0]]

*)
