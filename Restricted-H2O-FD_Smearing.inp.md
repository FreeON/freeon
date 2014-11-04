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

### Options section

Input to FreeON is divided in sections, each labeled by tags enclosed in less-than and greater-than signs. Every section in the file deals with a specific type of information: the options section allows you to enter your choices regarding the calculation and subsequent input:

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

You can get more detailed information in the manual and, possibly, in other examples.

    Charge=0
    Multiplicity=1

In this case, the file starts by defining the electronic state of the system, which has a net charge of zero, with an even number of electrons and all of them doubly filling the lowest level orbitals with no spin-uncoupled electrons.

    Pressure=0.0

Now, we tell FreeOn the pressure used in the Optimizer, MD or MC calculations, which in this case, is going to be none. This is one of the "miscellaneous" options in FreeON, and can be given in AU or GPa (Giga pascals).

    Guess=Superpos

Here FreeON is told to generate the initial guess from which the SCF calculation is to be started using the SuperPos (Superposition of Atomic Densities, SAD) method.

    OutPut=XYZ
    DebugAll=(MaxDebug)

These two lines help expand the output generated. The first one (Output=XYZ) requests that FreeON generate an XYZ formatted file with the geometry at the end of the calculation. The second (DebugAll=(MaxDebug) asks FreeON to output as much additional information as possible during the calculation. Normally, you will use the first line, but not the second.

    SCFMethod=(RH)

We are now specifying the kind of calculation: **RH** stands for [Restricted Hartree-Fock](http://en.wikipedia.org/wiki/Hartree–Fock_method). Notice that we are not specifying which kind of restricted Hartree-Fock calculation we want, the program will select the appropriate calculation for us.

    BasisSets=(STO-2G-SPLIT)

This is the [basis set](http://en.wikipedia.org/wiki/Basis_set_(chemistry)) to use. In this case we use a minimal basis set to speed up computation.

    ModelChem=(HF)

This option helps fine tune the kind of calculation, We are not further specifying it here, but could if needed (e. g. to specify the type of exchange-correlation functional for a DFT calculation).

    Accuracy=(good)

Here we indicate that we want to use convergence criteria that can be considered "good" yet not overly precise.

    SCFConvergence=(DIIS)

Now, we state the method to use for convergence calculations, which in this case is the [DIIS](http://en.wikipedia.org/wiki/DIIS) (direct inversion of iterated subspace) method.

    Smearing=(Fermi-Dirac)
    SmearingTemperature=10000.0D0

These two lines define the technique to use for smearing the electrons in order to account for non-ground temperature, in this case th Fermi-Dirac approach, and the temperature to simulate, which is set to 10000 degrees Kelvin.

    Geometry=InAngstrom

The last line in the options section prepares FreeON to interpret correctly the upcoming data: the atomic coordinates of the system geometry, which is here defined to give values in Angström units.

### The geometry section

Enclosed between the tags **<BeginGeometry>** and **<EndGeometry>** comes the definition of the system to be modeled: its format is very simple, it consists of the atom name followed by its XYZ coordinates, one atom per line. This format is similar enough to the well known XYZ format that it is trivial to convert from it (just remove the first two lines of an XYZ file leaving only the coordinates).

    <BeginGeometry>
     O     0.0000000000000000     0.0000000000000000     0.0000000000000000
     H     0.9393538135187500    -0.1240976319899556     0.0000000000000000
     H    -0.1240979353667913     0.9393535662638797     0.0000000000000000
    <EndGeometry>

[Jrvalverde](User:Jrvalverde "wikilink")
