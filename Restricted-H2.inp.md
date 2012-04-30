---
layout: default
title: Restricted-H2.inp
---


    # NWChem: -1.175482 hartree

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
    H 0.078794 0.0 0.0
    H 0.821206 0.0 0.0
    <EndGeometry>
