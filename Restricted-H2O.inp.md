---
layout: default
title: Restricted-H2O.inp
---


    <BeginOptions>

    Charge=0
    Multiplicity=1

    Pressure=0.0

    Guess=Superpos
    OutPut=XYZ
    #Grad=(Optimize,PrimInt,NoBackTr,NoGDIIS,BiSect,HBondOnly,NoFragmConnect,NonCovBend,NonCovTors)
    DebugAll=(MinDebug,CheckSums)

    SCFMethod=(TC2,TC2)
    BasisSets=(STO-2G,6-31G**-SPLIT)
    ModelChem=(HF,HF)
    Accuracy=(Loose,Good)
    SCFConvergence=(DIIS,DIIS)

    Geometry=InAngstrom

    <EndOptions>

    <BeginGeometry>
     O     0.0000000000000000     0.0000000000000000     0.0000000000000000
     H     0.9393538135187500    -0.1240976319899556     0.0000000000000000
     H    -0.1240979353667913     0.9393535662638797     0.0000000000000000
    <EndGeometry>
