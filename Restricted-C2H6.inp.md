---
layout: default
title: Restricted-C2H6.inp
---


    # NWChem with XFine setting: -79.830421 hartrees
    # which also agrees perfectly with the other code.

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
    C 0.0 0.0 0.7654
    C 0.0 0.0 -0.7654
    H 0.0 1.021 1.1645
    H -0.8842 -0.5105 1.1645
    H 0.8842 -0.5105 1.1645
    H 0.0 -1.021 -1.1645
    H -0.8842 0.5105 -1.1645
    H 0.8842 0.5105 -1.1645
    <EndGeometry>
