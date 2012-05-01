---
layout: default
title: Restricted-C2H4.inp
---

![This example will show you how to run a calculation on [Ethylene](http://en.wikipedia.org/wiki/Ethylene) (C<sub>2</sub>H<sub>4</sub>) using [FreeON](http://freeon.org)](c2h4_ESP_colored_VdW_surface.png "This example will show you how to run a calculation on Ethylene (C2H4) using FreeON")

Restricted Hartree-Fock calculation on Ethylene
-----------------------------------------------

This example will show you how to run a calculation on [Ethylene](http://en.wikipedia.org/wiki/Ethylene) (C<sub>2</sub>H<sub>4</sub>). Sample results obtained with [NWChem](http://www.nwchem-sw.org) are provided as reference, as well as commented instructions to compute a structure optimization.

### The complete input file

You can find this file in the **Validate/SinglePoint/Restricted** subdirectory of FreeON source code distribution.

    @ NWChem: -78.587459 hartree

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
    C 0.0 0.0 0.6656
    C 0.0 0.0 -0.6656
    H 0.0 0.9236 1.2397
    H 0.0 -0.9236 1.2397
    H 0.0 -0.9236 -1.2397
    H 0.0 0.9236 -1.2397
    <EndGeometry>
