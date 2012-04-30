---
layout: default
title: Restricted-CH3NH2.inp
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
    C \t0.0517 \t0.7038 \t0.0000
    N \t0.0517 \t-0.7609 \t0.0000
    H \t-0.9430 \t1.1848 \t0.0000
    H \t0.5941 \t1.0629 \t0.8815
    H \t0.5941 \t1.0629 \t-0.8815
    H \t-0.4584 \t-1.1037 \t-0.8123
    H \t-0.4584 \t-1.1037 \t0.8123
    <EndGeometry>
