---
layout: default
title: Restricted-H.inp
---

![This example performs an RHF calculation on the [Hydrogen](http://en.wikipedia.org/wiki/Hydrogen) atom (H) using [FreeON](http://freeon.org) ](H_1s.png "fig:This example performs an RHF calculation on the Hydrogen atom (H) using FreeON ") ![Hydrogen's 2s orbital ](H_2s.png "fig:Hydrogen's 2s orbital ") ![Hydrogen's 2px orbital ](H_2px.png "fig:Hydrogen's 2px orbital ") ![Hydrogen's 2px 2py and 2pz orbitals superposed ](H_2pxyz.png "fig:Hydrogen's 2px 2py and 2pz orbitals superposed ") ![Hydrogen's 3s orbital ](H_3s.png "fig:Hydrogen's 3s orbital ") ![Hydrogen's 3dxy orbital ](H_3dxy.png "fig:Hydrogen's 3dxy orbital ")

Hydrogen
--------

This example performs an RHF calculation on the [Hydrogen](http://en.wikipedia.org/wiki/Hydrogen) atom (H) using [FreeON](http://freeon.org). Reference final energy obtained with NWchem is provided as a reference.

### The complete input file

You can find this example configuration file under directory **Validate/SinglePoint/Restricted** of the FreeON source code distribution.


    # NWChem: -0.500273 hartree

    <BeginOptions>

    Charge=0
    Multiplicity=2

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
    H 0.0 0.0 0.0
    <EndGeometry>
