---
layout: default
title: 15-LJ-2-CG.inp
---


    # Ar dimer with Lennard-Jones potential.

    <BeginOptions>
    Charge=0

    Multiplicity=1
    Guess=SuperPos
    Grad=(Optimize,Cartesian,CG)
    DebugAll=(MaxDebug,CheckSums,MaxGeOp)
    BasisSets=(STO-2G-SPLIT)
    SCFMethod=(RH)
    SCFConvergence=(DIIS)
    ModelChem=(HF)
    Accuracy=(verytight)
    Geometry=(InAngstroms)

    UseLennardJones=yes

    <EndOptions>

    <BeginGeometry>
    Ar 0.0 0.0 0.0
    Ar 2.0 0.0 0.0
    <EndGeometry>
