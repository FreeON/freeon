---
layout: default
title: Restricted-H2O-FD Smearing.inp
---

![This example drives a calculation on water using Fermi-Dirac smearing to simulate the canonical ensemble for electrons at 10000K using [FreeON](http://freeon.org) ](H2O_ESP.png "This example drives a calculation on water using Fermi-Dirac smearing to simulate the canonical ensemble for electrons at 10000K using FreeON ")

Water
-----

Simulation of the canonical ensemble for the electrons at nonzero temperature (Fermi-Dirac occupation numbers) using [FreeON](http://freeon.org). This calculation uses Fermi-Dirac occupation numbers to simulate a temperature of 10000 degrees Kelvin.

### The complete input file

You can find this file in the **Validate/SinglePoint/Restricted** subdirectory of FreeON source code distribution.

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
