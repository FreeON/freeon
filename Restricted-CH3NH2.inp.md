---
layout: default
title: Restricted-CH3NH2.inp
---

![This example shows how to run an RHF calculation on [Methylamine](http://en.wikipedia.org/wiki/Methylamine) (CH<sub>3</sub>NH<sub>2</sub>) using [FreeON](http://freeon.org/) ](CH3NH2.png "This example shows how to run an RHF calculation on Methylamine (CH3NH2) using FreeON ")

Restricted Hartree-Fock calculation on Methylamine
--------------------------------------------------

This example will show you how to run a calculation on [Methylamine](http://en.wikipedia.org/wiki/Methylamine) (CH<sub>3</sub>NH<sub>2</sub>). Commented instructions to compute a structure optimization are included.

### The complete input file

You can find this file in the **Validate/SinglePoint/Restricted** subdirectory of FreeON source code distribution.

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
    C \t0.0517 \t0.7038 \t0.0000
    N \t0.0517 \t-0.7609 \t0.0000
    H \t-0.9430 \t1.1848 \t0.0000
    H \t0.5941 \t1.0629 \t0.8815
    H \t0.5941 \t1.0629 \t-0.8815
    H \t-0.4584 \t-1.1037 \t-0.8123
    H \t-0.4584 \t-1.1037 \t0.8123
    <EndGeometry>

### Options section

Next comes the options section. Input to FreeON is divided in sections, each labeled by tags enclosed in less-than and greater-than signs. Every section in the file deals with a specific type of information: the options section allows you to enter your choices regarding the calculation and subsequent input:

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

You can get more detailed information in the manual and, possibly, in other examples.

    Charge=0
    Multiplicity=1

In this case, the file starts by defining the electronic state of the system, which has a net charge of zero, with an even number of electrons and all of them doubly filling the lowest level orbitals with no spin-uncoupled electrons.

    Guess=Superpos

Here FreeON is told to generate the initial guess from which the SCF calculation is to be started using the SuperPos (Superposition of Atomic Densities, SAD) method.

    OutPut=XYZ
    DebugAll=(MaxDebug)

These two lines help expand the output generated. The first one (Output=XYZ) requests that FreeON generate an XYZ formatted file with the geometry at the end of the calculation. The second (DebugAll=(MaxDebug) asks FreeON to output as much additional information as possible during the calculation. Normally, you will use the first line, but not the second.

    #Grad=(Optimize,PrimInt,BiSect,NoGDIIS)

This is a line that has been invalidated by using the \# at the beginning. Without it, it would ask FreeON to optimize the structure. You can reactivate it by removing the \#.

    SCFMethod=(RH)
    BasisSets=(6-31G*)
    ModelChem=(B3LYP)
    Accuracy=(VeryTight)
    SCFConvergence=(DIIS)

These lines define the computation to carry out: the SCF method will use the Restricted Hartree-Fock approach, with the 6-31G\* basis set, the B3LYP DFT functional, and it is to be iteratively applied until a very tight convergence is achieved using the Pulay solver DIIS (direct inversion of iterated subspace).

    Geometry=InAngstrom

This last option tells FreeOn that the geometry provided in the Geometry section uses Angströms as the measure unit,

### The geometry section

Enclosed between the tags **<BeginGeometry>** and **<EndGeometry>** comes the definition of the system to be modeled: its format is very simple, it consists of the atom name followed by its XYZ coordinates, one atom per line. This format is similar enough to the well known XYZ format that it is trivial to convert from it (just remove the first two lines of an XYZ file leaving only the coordinates).

    <BeginGeometry>
    C \t0.0517 \t0.7038 \t0.0000
    N \t0.0517 \t-0.7609 \t0.0000
    H \t-0.9430 \t1.1848 \t0.0000
    H \t0.5941 \t1.0629 \t0.8815
    H \t0.5941 \t1.0629 \t-0.8815
    H \t-0.4584 \t-1.1037 \t-0.8123
    H \t-0.4584 \t-1.1037 \t0.8123
    <EndGeometry>

[Jrvalverde](User:Jrvalverde "wikilink")