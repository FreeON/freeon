(* :Title: Fixed Length Number Formatting - A Hack *)

(* :Context: FixedNumberForm`*)

(* :Author: John M. Novak *)

(* :Summary: 
	This packed introduces the hack FixedNumberForm.  It allows
	numbers to be formatted with a fixed length mantissa and
	exponent.  Numbers are also expressed in the "e" form.
	It is intended for output to files for import into
	applications that only understand this format...
*)

(* :Package Version: 1.0 *)

(* :Mathematica Version: 2.0 *)

(* :History:
	V.1.0, May 92 by John M. Novak
*)

(* :Keywords:
	fixed length numbers, exponent format
*)

(* :Warnings:
	All numbers will be converted to the precision
	specified (or machine precision, whichever is greater.)
	This also includes symbolic numbers (e.g., Pi, E) and
	integers.
*)

(* :Discussion:
	This is more of a hack than anything.  It's principle use
	should be limited to formatting numbers for export (where
	dealing with the conversion of integers and symbolic
	numbers isn't a problem.)  It must be emphasized that this
	conversion does modify the actual expression, not just
	the formatting (to do this just in the formatting would
	require that I be a good deal more clever... :-)
	
	The basic syntax is
	FixedNumberForm[expr, mantissalength, exponentlength]
	where all numbers will be changed to mantissalength
	precision, and will be formatted as <mantissa>e<exponent>.
*)

BeginPackage["FixedNumberForm`"]

FixedNumberForm::usage =
"FixedNumberForm[expr, mantissalength, exponentlength] converts
numbers in expr to mantissalength precision, and formats them
in a fixed length format.";

Begin["`Private`"]

padmantissa[num_, length_] :=
	Module[{str = num, zeros = length - (StringLength[num] - 1)},
		If[StringTake[num,1] === "-", zeros = zeros + 1,
			str = "+"<>str];
		str <> Table["0",{zeros}]
	]

padexponent[exp_, length_] :=
	Module[{str = exp, sign, zeros},
		If[StringTake[str,1] === "-",
			sign = "-"; str = StringDrop[str,1],
			sign = "+"];
		zeros = length - StringLength[str];
		If[zeros < 0,
			Message[FixedNumberForm::truncexp];
			str = StringTake[str, length];
			zeros = 0
		];
		sign <> Table["0",{zeros}] <> str
	]

FixedNumberForm::truncexp =
"Warning: length of exponent was too great for field,
truncating...";

FixedNumberForm[expr_, digits_Integer:6, expdigits_Integer:2] :=
	ScientificForm[N[expr,digits], digits,
		NumberFormat ->
			((SequenceForm[
				padmantissa[#1,digits],
				"D",padexponent[#3,expdigits]])&
			)]

End[]

EndPackage[]



