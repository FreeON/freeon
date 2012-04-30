---
layout: default
title: Restricted-H2O-FD Smearing.inp
---


    <BeginOptions>

    Charge=0
    Multiplicity=1

    Pressure=0.0

    Guess=Superpos
    OutPut=XYZ
    DebugAll=(MaxDebug,CheckSums)

    SCFMethod=(RH)
    BasisSets=(STO-2G-SPLIT)
    ModelChem=(HF)
    Accuracy=(Good)
    SCFConvergence=(DIIS)
    Smearing=(Fermi-Dirac)
    SmearingTemperature=10000.0D0

    Geometry=InAngstrom

    <EndOptions>

    <BeginGeometry>
     O     0.0000000000000000     0.0000000000000000     0.0000000000000000
     H     0.9393538135187500    -0.1240976319899556     0.0000000000000000
     H    -0.1240979353667913     0.9393535662638797     0.0000000000000000
    <EndGeometry>
