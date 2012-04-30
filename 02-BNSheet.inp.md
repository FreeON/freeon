---
layout: default
title: 02-BNSheet.inp
---

    <BeginOptions>

    MPIInvoke="mpirun"
    MPIProcessors=8
    MPIProcFlag="-np"
    MPISpatialProc=8
    MPIMachineFlag="--nper"
    MPIMachineFile="2"

    Charge=0
    Multiplicity=1

    Grad=(Optimize,PrimInt,NoBackTr,BiSect,NoGDIIS)
    DebugAll=(MaxDebug,MinGeOp)

    Guess=SuperPos
    SCFMethod=(RH)
    BasisSets=(STO-2G-SPLIT)
    ModelChem=(B3LYP)
    Accuracy =(good)
    SCFConvergence=(DIIS)

    PBC=(T,T,F)
    PFFMaxEll=14
    Periodic=(AtomCoord)

    <EndOptions>

    <BeginPeriodic>
     7.260   7.260   1.000  90.000  90.000 120.000
    <EndPeriodic>

    <BeginGeometry>
       B          0.60500000          0.34929691          0.00000000
       N          0.59693333          1.76045644          0.00000000
       B         -0.60500000          2.44507839          0.00000000
       N         -0.61306667          3.85623792          0.00000000
       B         -1.81500000          4.54085987          0.00000000
       N         -1.82306667          5.95201940          0.00000000
       B          3.02500000          0.34929691          0.00000000
       N          3.01693333          1.76045644          0.00000000
       B          1.81500000          2.44507839          0.00000000
       N          1.80693333          3.85623792          0.00000000
       B          0.60500000          4.54085987          0.00000000
       N          0.59693333          5.95201940          0.00000000
       B          5.44500000          0.34929691          0.00000000
       N          5.43693333          1.76045644          0.00000000
       B          4.23500000          2.44507839          0.00000000
       N          4.22693333          3.85623792          0.00000000
       B          3.02500000          4.54085987          0.00000000
       N          3.01693333          5.95201940          0.00000000
    <EndGeometry>
