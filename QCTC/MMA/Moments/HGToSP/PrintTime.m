(*********************************************************************

        Adapted from
        Roman E. Maeder: Programming in Mathematica,
        Second Edition, Addison-Wesley, 1991.

 *********************************************************************)


PrintTime::usage = "PrintTime[expr] prints the time it takes
	to evaluate expr and returns the result of the evaluation."

Begin["`Private`"]

SetAttributes[PrintTime, HoldAll]

PrintTime[expr_] :=
	Module[{timing},
		timing = Timing[expr];
		Print[ timing[[1]] ];
		timing[[2]]
	]

End[]

Null
