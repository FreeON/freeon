---
layout: default
title: 04-CH2Wire.inp
---

    <BeginTitle>
    Polyacetilene geometry optimization with fixed lattice vectors
    in STO-3G basis
    <EndTitle>


    <BeginOptions>
    MPIInvoke="mpiun"
    MPIProcFlag="-np"
    MPIProcessors=8
    MPISpatialProc=8
    MPIMachineFlag="--nper"
    MPIMachineFile="2"

    Charge=0
    Multiplicity=1

    Grad=(Optimize,PrimInt,BiSect,NoGDIIS)
    DebugAll=(MinDebug,CheckSums)
    Guess=SuperPos
    SCFMethod=(TC2)
    BasisSets=(STO-3G-SPLIT)
    ModelChem=(HF)
    Accuracy =(good)
    SCFConvergence=(DIIS)

    PBC=(T,F,F)
    PFFMaxEll=16
    Periodic=(AtomCoord)

    <EndOptions>

    <BeginPeriodic>
     12.000  1.0 1.0     90.0  90.0  90.0
    <EndPeriodic>

    <BeginGeometry>
     C    0.500   0.500   0.000
     H    0.500   1.300   0.800
     H    0.500   1.300  -0.800
     C    1.500  -0.500   0.000
     H    1.500  -1.300   0.800
     H    1.500  -1.300  -0.800
     C    2.500   0.500   0.000
     H    2.500   1.300   0.800
     H    2.500   1.300  -0.800
     C    3.500  -0.500   0.000
     H    3.500  -1.300   0.800
     H    3.500  -1.300  -0.800
     C    4.500   0.500   0.000
     H    4.500   1.300   0.800
     H    4.500   1.300  -0.800
     C    5.500  -0.500   0.000
     H    5.500  -1.300   0.800
     H    5.500  -1.300  -0.800
     C    6.500   0.500   0.000
     H    6.500   1.300   0.800
     H    6.500   1.300  -0.800
     C    7.500  -0.500   0.000
     H    7.500  -1.300   0.800
     H    7.500  -1.300  -0.800
     C    8.500   0.500   0.000
     H    8.500   1.300   0.800
     H    8.500   1.300  -0.800
     C    9.500  -0.500   0.000
     H    9.500  -1.300   0.800
     H    9.500  -1.300  -0.800
     C   10.500   0.500   0.000
     H   10.500   1.300   0.800
     H   10.500   1.300  -0.800
     C   11.500  -0.500   0.000
     H   11.500  -1.300   0.800
     H   11.500  -1.300  -0.800
    <EndGeometry>
