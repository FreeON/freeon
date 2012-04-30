---
layout: default
title: Unrestricted-Be.inp
---



    GAMESS UPW91/3-21G  -0.9915542467
    FreeON UPW91/3-21G   -.9915540919968984D+00 (VeryTight)
    FreeON UPW91/3-21G   -.9915501733931975D+00 (Tight)

    <BeginOptions>

    Charge=0
    Multiplicity=1

    Guess=Superpos
    OutPut=XYZ
    DebugAll=(MaxDebug)

    #Grad=(Optimize,PrimInt,BiSect,NoGDIIS)

    SCFMethod=(RH,RH)
    BasisSets=(STO-3G,3-21G)
    SpinModel=(U,U)
    ModelChem=(HF,PW91xc)
    Accuracy=(Good,VeryTight)
    SCFConvergence=(DIIS,DIIS)

    Geometry=InAngstrom

    <EndOptions>

    <BeginGeometry>
    H 0.0 0.0 0.0
    H 5.0 0.0 0.0
    <EndGeometry>
