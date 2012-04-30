---
layout: default
title: Restricted-CH3OH.inp
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
    C \t-0.0469 \t0.6608 \t0.0000
    O \t-0.0469 \t-0.7577 \t0.0000
    H \t-1.0944 \t0.9749 \t0.0000
    H \t0.4371 \t1.0864 \t0.8933
    H \t0.4371 \t1.0864 \t-0.8933
    H \t0.8762 \t-1.0514 \t0.0000
    <EndGeometry>
