---
layout: default
title: Restricted-N2O.inp
---


    <BeginOptions>

    Charge=0
    Multiplicity=1

    Guess=Superpos
    OutPut=XYZ
    DebugAll=(MaxDebug)

    #Grad=(Optimize,PrimInt,BiSect,NoGDIIS)

    SCFMethod=(RH)
    BasisSets=(6-31G*)
    ModelChem=(B3LYP)
    Accuracy=(VeryTight)
    SCFConvergence=(DIIS)

    Geometry=InAngstrom

    <EndOptions>

    <BeginGeometry>
    N \t0.0000 \t0.0000 \t-1.2070
    N \t0.0000 \t0.0000 \t-0.0727
    O \t0.0000 \t0.0000 \t1.1198
    <EndGeometry>
