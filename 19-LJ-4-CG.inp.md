---
layout: default
title: 19-LJ-4-CG.inp
---

    # Ar cluster with Lennard-Jones potential.

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
    Ar 0.0 2.0 0.0
    Ar 0.0 0.0 2.0
    <EndGeometry>
