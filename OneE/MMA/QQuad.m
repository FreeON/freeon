(*GET THE MONDO INFO*)MondoHome=Environment["MONDO_HOME"];
If[MondoHome\[Equal]$FAILED,
    Print["COULD NOT FIND $MONDO_HOME! CHECK YOUR .cshrc "];
    Abort[];];

(*LOAD FIXEDNUBERFORM*)

Get[StringJoin[MondoHome,"/MMA/FixedNumberForm.m"]];

<<NumericalMath`GaussianQuadrature`;

FF[x_]:=ToString[FixedNumberForm[x,16,2]];

FileName = "QQuad.Inc";
Print["FileName = ", FileName];
OpenWrite[FileName];

Do[ 

pw=GaussianQuadratureWeights[32,0,infinity,32];

len = Length[pw];
steps = 18;

(* WRITE OUT THE POINTS ARRAY *)

div = Floor[len/steps];
Do[
    begin = j steps+1;
    end = (j + 1)steps;
    DataString=StringJoin["      DATA(Points(IPts,",ToString[infinity],"),IPts=", ToString[begin], ",",ToString[end], ")/ & \n"];
    WriteString[FileName, DataString];
    Do[
       WriteString[FileName,StringJoin["      ", FF[pw[[k, 1]]], ", ", FF[pw[[k+1, 1]]],", ", FF[pw[[k+2, 1]]], ", & \n"]];
      ,{k,begin,end-3,3}];
    WriteString[FileName,StringJoin["      ", FF[pw[[end-2, 1]]], ", ", FF[pw[[end-1, 1]]],", ", FF[pw[[end, 1]]], "/\n"]];
  ,{j, 0, div - 1}];


begin = div steps+1;
end = len;
diff=end-begin+1;
DataString=StringJoin["      DATA(Points(IPts,",ToString[infinity],"),IPts=",ToString[begin], ",", ToString[end], ")/ & \n"];
WriteString[FileName, DataString];
If[diff >= 4,Do[
                WriteString[FileName,StringJoin["      ", FF[ pw[[k, 1]] ], ", ", FF[ pw[[k+1, 1]] ], ", ", FF[ pw[[k+2, 1]] ], ", & \n"]];
               ,{k,begin,end-3,3}];
  ];
If[Mod[diff, 3] == 1,WriteString[FileName, StringJoin["      ", FF[pw[[end, 1]] ], "/ \n"]]];
If[Mod[diff, 3] == 2,WriteString[FileName, StringJoin["      ", FF[pw[[end - 1, 1]]], ", ", FF[pw[[end, 1]]],"/\n"]]]; 
If[Mod[diff, 3] == 0,WriteString[FileName, StringJoin["      ", FF[pw[[end - 2, 1]]], ", ", FF[ pw[[end - 1, 1]]],", ", FF[pw[[end, 1]]], "/\n"]]];

(* WRITE OUT THE WEIGHTS ARRAY *)

div = Floor[len/steps];
Do[
    begin = j steps+1;
    end = (j + 1)steps;
    DataString=StringJoin["      DATA(Weights(IWts,",ToString[infinity],"),IWts=", ToString[begin], ",",ToString[end], ")/ & \n"];
    WriteString[FileName, DataString];
    Do[
       WriteString[FileName,StringJoin["      ", FF[pw[[k, 2]]], ", ", FF[pw[[k+1, 2]]],", ", FF[pw[[k+2, 2]]], ", & \n"]];
      ,{k,begin,end-3,3}];
     WriteString[FileName,StringJoin["      ", FF[pw[[end-2, 2]]], ", ", FF[pw[[end-1,2]]],", ", FF[pw[[end,2]]], "/\n"]];
  ,{j, 0, div - 1}];
begin = div steps+1;
end = len;
diff=end-begin+1;
DataString=StringJoin["      DATA(Weights(IWts,",ToString[infinity],"),IWts=",ToString[begin], ",", ToString[end], ")/ & \n"];
WriteString[FileName, DataString];
If[diff >= 4,Do[
                WriteString[FileName,StringJoin["      ", FF[ pw[[k, 2]] ], ", ", FF[ pw[[k+1, 2]] ], ", ", FF[ pw[[k+2, 2]] ], ", & \n"]];
               ,{k,begin,end-3,3}];
  ];
If[Mod[diff, 3] == 1,WriteString[FileName, StringJoin["      ", FF[pw[[end, 2]] ], "/ \n"]]];
If[Mod[diff, 3] == 2,WriteString[FileName, StringJoin["      ", FF[pw[[end - 1, 2]]], ", ", FF[pw[[end, 2]]],"/\n"]]]; 
If[Mod[diff, 3] == 0,WriteString[FileName, StringJoin["      ", FF[pw[[end - 2, 2]]], ", ", FF[ pw[[end - 1, 2]]],", ", FF[pw[[end, 2]]], "/\n"]]];

,{infinity,1,10}]

Close[FileName];

