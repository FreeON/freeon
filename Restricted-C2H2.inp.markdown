---
layout: default
title: Restricted-C2H2.inp
---

![This example will show you how to calculate the energy of [acetylene](http://en.wikipedia.org/wiki/Acetylene) C<sub>2</sub>H<sub>2</sub> with DFT using [FreeON](http://freeon.org)](C2H2_EDS.png "This example will show you how to calculate the energy of acetylene C2H2 with DFT using FreeON")

Calculation on acetylene using DFT
----------------------------------

This configuration file runs a calculation on acetylene (C<sub>2</sub>H<sub>2</sub>) using DFT and the 6-31G\* basis set. The results obtained by NWchem are provided for reference so you can contrast your results.

Inside the code you will find instructions to run an optimization of the structure commented out. These will be ignored by the program.

### The complete input file

You can find the file below in the directory **Validate/SinglePoint/Restricted** of the FreeON source code distribution.

    @NWChem: -77.325646 hartree

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
    C 0.0 0.0 0.6025
    C 0.0 0.0 -0.6025
    H 0.0 0.0 1.6691
    H 0.0 0.0 -1.6691
    <EndGeometry>

### Introductory comments

As this is a validation example, the final energy obtained using [NWchem](http://www.nwchem-sw.org) to run a similar calculation is provided first. You can use this energy to verify that FreeON produces similar results (give or take an epsilon due to carried-on rounding errors) Note that any text outside one of the configuration sections (the text between *<Begin...>* and *<End...>* tags) is ignored by the program, so that you can use it to include comments on your files.

### Options section

In this section we describe the electronic properties of our system and the kind of calculation we want to perform. Please refer to the manual for a detailed description of all the options (or to the source code for the latest information).

The options section is enclosed by the **<BeginOptions>** and **<EndOptions>** tags.

    Charge=0
    Multiplicity=1

The file begins by defining the charge and multiplicty of the system, which, as should be obvious, is in the ground state, with no unpaired electrons and all charges equilibrated: the optimal situation for a Restricted Hartree-Fock (RHF) calculation.

    Guess=SuperPos

Next, we define how to guess the starting configuration. It is important to try to get this right as the closer it is to the final solution, the faster our calculation will converge (and the less time we will spend waiting for our results). Here we select the SuperPos (Superposition of Atomic Densities, SAD) method to generate the starting configuration.

    OutPut=XYZ
    DebugAll=(MaxDebug)

The following options ask FreeON to produce special output: first we ask to get a printout of the final conformation (the optimized structure) in XYZ format, and then we ask FreeON to log out additional information so we can verify the details of the calculation. You will normally **not** use this option (the **DebugAll**) unless you have reason to suspect the results (say because the calculation does not converge or gives a result you are uncomfortable with).

    #Grad=(Optimize,PrimInt,BiSect,NoGDIIS)

Any "option" that is not valid will be ignored by FreeON. This is a very handsome feature as long as you are careful with your work.

**Grad** tells FreeON what kind of gradients you want to compute and how. But as it has been written **\#Grad** it is now an invalid option and is ignored. This comes very handy: it means you can enable or disable any line in the configuration file by "invalidating" the corresponding keyword.

If you remove the \# from the line, then it will become valid again, and FreeON will use it, calculating an optimized structure of acetylene.

Note that you can do anything you like as long as you "invalidate" the line, but it is better to do something that is very obvious (see this and other examples).

The down side of this is that if you involuntarily "invalidate" a line (say by misspelling the keyword), then it will be ignored as well. You calculation might not be what you expected. This means that you should always check carefully your input file for misspelling errors, specially if the results you get depart wildly from your expectations.

    SCFMethod=(RH)
    BasisSets=(6-31G*)
    ModelChem=(B3LYP)
    Accuracy=(VeryTight)
    SCFConvergence=(DIIS)

The next block explains the type of calculation to perform (restricted Hartree-Fock, using the 6-31G\* basis and the B3LYP DFT functional, aim for a tight convergence using the DIIS converger).

    Geometry=InAngstrom

The final option we use tells FreeON that the geometry of the molecule (which goes in its own section) is to be specified in Angström.

### Molecular Geometry

The system geometry is specified in a separate section, delimited by the **<BeginGeometry>** and **<EndGeometry>** tags.

    <BeginGeometry>
    C 0.0 0.0 0.6025
    C 0.0 0.0 -0.6025
    H 0.0 0.0 1.6691
    H 0.0 0.0 -1.6691
    <EndGeometry>

As you can see it is very simple: all you need to provide is the atom names and their X, Y and Z coordinates. This is a nice feature as well: it means that you can take any structure from any database, convert it to XYZ format (e. g. with OpenBabel) and then take the coordinates off this XYZ file (which reduces to deleting the first two lines of the XYZ file).

[Jrvalverde](User:Jrvalverde "wikilink")
