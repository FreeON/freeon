---
layout: default
title: Restricted-C.inp
---


    # NWChem: -37.846280 hartree

    <BeginOptions>

    Charge=0
    Multiplicity=3

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
    C 0.0 0.0 0.0
    <EndGeometry>
