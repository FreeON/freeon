---
layout: default
title: Restricted-C2H6.inp
---

![This example will show you how to run a calculation on [Ethane](http://en.wikipedia.org/Ethane) (C<sub>2</sub>H<sub>6</sub> using [FreeON](http://freeon.org) ](C2h6_eldens.png "This example will show you how to run a calculation on Ethane (C2H6 using FreeON ")

RHF calculation on Ethane
-------------------------

This page shows how to make an RHF calculation on [Ethane](http://en.wikipedia.org/Ethane) (C<sub>2</sub>H<sub>6</sub> using [FreeON](http://freeon.org)

### The complete input file

This file is available from the **Validate/SinglePoint/Restricted** subdirectory of the FreeON source code distribution.


    # NWChem with XFine setting: -79.830421 hartrees
    # which also agrees perfectly with the other code.

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
    C 0.0 0.0 0.7654
    C 0.0 0.0 -0.7654
    H 0.0 1.021 1.1645
    H -0.8842 -0.5105 1.1645
    H 0.8842 -0.5105 1.1645
    H 0.0 -1.021 -1.1645
    H -0.8842 0.5105 -1.1645
    H 0.8842 0.5105 -1.1645
    <EndGeometry>
