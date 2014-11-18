---
layout: default
title: 03-Caffeine.inp
---

    # vim: syntax=conf

    <BeginTitle>
                  ** Caffine **
    <EndTitle>

    <BeginOptions>

    Charge=0
    Multiplicity=1

    #DebugAll=(MaxDebug,CheckSums,MatDebug)
    DebugAll=(MinDebug)

    Grad=(Optimize,PrimInt,NoGDIIS,NoBackTr,DiagHess)

    Guess=SuperPos
    #Guess=Restart
    #HDFFile=/tmp/03-Caffeine_3826.hdf

    SCFMethod=(RH,RH)
    BasisSets=(STO-2G,6-31G**)
    ModelChem=(HF,BLYPxc)
    Accuracy=(Loose,Loose)
    SCFConvergence=(DIIS,DIIS)
    #DIISDamp=0.9
    #DIISDelay=2
    #DIISFirstSCF=2

    MPIInvoke="mpirun"
    MPIProcessors=4
    MPIProcFlag="-np"
    MPISpatialProc=4
    MPIMachineFlag="-machinefile"
    MPIMachineFile="$PBS_NODEFILE"

    <EndOptions>

    <BeginGeometry>
    C    0.00000   0.00000   0.00000
    C   -2.07229   0.51445   1.16368
    C   -2.31021  -0.77836   1.45023
    C   -1.39670  -1.80726   1.02520
    C   -3.83762   0.49274   2.30949
    C    0.70055  -2.32673  -0.17010
    C   -0.66824   2.36013   0.16415
    C   -4.22682  -1.95759   2.73238
    N   -0.32017  -1.35294   0.32996
    N   -0.90029   0.90864   0.44601
    N   -3.51598  -0.79506   2.16611
    N   -2.95349   1.30784   1.72590
    O    0.95316   0.29230  -0.66741
    O   -1.54357  -3.04162   1.22602
    H   -4.60867   0.79590   2.77433
    H    0.90850  -2.94850   0.51740
    H    0.35323  -2.79107  -0.92259
    H    1.48777  -1.85866  -0.42249
    H   -0.51555   2.82082   0.98082
    H    0.09007   2.45586  -0.40003
    H   -1.43198   2.72498  -0.26724
    H   -5.00451  -1.66093   3.19031
    H   -4.47780  -2.54754   2.03133
    H   -3.65576  -2.41270   3.34005
    <EndGeometry>
