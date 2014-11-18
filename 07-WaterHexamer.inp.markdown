---
layout: default
title: 07-WaterHexamer.inp
---

    <BeginOptions>

    MPIInvoke="mpirun"
    MPIProcessors=4
    MPIProcFlag="-np"
    MPISpatialProc=4
    MPIMachineFlag="--nper"
    MPIMachineFile="2"

    Charge=0
    Multiplicity=1
    Guess=SuperPos

    OutPut=XYZ
    Grad=(Optimize,PrimInt,NoGDIIS,BiSect,HBondOnly,NoFragmConnect,NonCovBend,NonCovTors)

    DebugAll=(MinDebug,CheckSums)
    SCFMethod=(RH)
    BasisSets=(STO-3G-SPLIT)
    ModelChem=(HF)
    Accuracy =(Tight)
    SCFConvergence=(DIIS)

    <EndOptions>

    <BeginGeometry>
    O 15.166760 -3.618977 1.331220
    H 14.451769 -3.008255 1.582111
    H 14.933572 -4.497920 1.656471
    O -1.257823 1.704853 0.609162
    H -1.524830 2.098583 -0.232179
    H -0.880205 2.395651 1.169851
    O -1.639863 -0.218623 1.158248
    H -0.742152 -0.552159 1.225089
    H -1.639156 0.732185 1.290772
    O -0.777513 1.350434 -1.026719
    H -0.777513 1.350434 -0.066719
    H 0.127423 1.350434 -1.347174
    O -1.712020 3.116422 -1.091396
    H -1.712020 3.116422 -0.131396
    H -0.807084 3.116422 -1.411851
    O -1.633074 -2.742938 4.053256
    H -2.104468 -3.234595 4.729762
    H -0.688528 -2.878932 4.157846
    <EndGeometry>
