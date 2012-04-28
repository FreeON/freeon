---
layout: default
title: CH4.inp
---

![This page will show you how to run a calculation on methane (CH<sub>4</sub>](CH4.png "This page will show you how to run a calculation on methane (CH4")

Calculation on methane CH<sub>4</sub>
-------------------------------------

This page will show you how to write a sample calculation for methane ground state.

### The complete file

You can find this file in the examples subdirectory of the [freeON](http://freeon.org) distribution.

    # NWChem: -40.518383 hartree

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
    C 0.0 0.0 0.0
    H 0.6312 0.6312 0.6312
    H -0.6312 -0.6312 0.6312
    H -0.6312 0.6312 -0.6312
    H 0.6312 -0.6312 -0.6312
    <EndGeometry>

### Comments

The first line in this file is a comment. FreeON will ignore any lines starting with a hash '\#' mark. This example starts with a comment where the energy calculated using [NWchem](http://www.nwchem-sw.org) is provided for your reference so that you can compare it with your own results using FreeON. Incidentally, the same calculation using [MPQC](http://www.mpqc.org) and the 6-311G\*\* basis set gives an energy of -105561.604 kJ/mol (-40.20629324 Hartree).

### Options

Calculations in FreeON are defined by specifying a number of options. This is done in a special section of the file, tagged between an opening *<BeginOptions>* tag and a closing *<EndOptions>* tag. Let us consider the options chosen for this calculation: first notice that they are stated in no predefined input order, you can arrange them as you wish.

    Charge=0

This is the total charge of the system, which in this case is zero.

    Multiplicity=1

This indicates FreeON the [multiplicity](http://en.wikipedia.org/wiki/Multiplicity_(chemistry)) of the system. Since we do not have any unpaired electrons, the multiplicity of our system is 1.

    Guess=Superpos

This telss FreeON that we want to start the calculation from an initial configuration computed using the SuperPos ([Superposition of Atomic Densities, SAD](http://igitur-archive.library.uu.nl/chem/2007-0302-200920/pdf18.pdf)) method.

    OutPut=XYZ

This option request from FreeON to produce an output file with the coordinates of the system in [XYZ format](http://en.wikipedia.org/wiki/XYZ_file_format). This is useful as it will allow us to visualize the system and use the coordinates for further computations. We can also convert these coordinates to other formats using [OpenBabel](http://en.wikipedia.org/wiki/OpenBabel).

    DebugAll=(MaxDebug)

This line tells FreeON to include in the output additional information. You will normally not use this. It is intended mainly for the developers of the software so they can track a calculation and spot where is any problem in the unlikely case that it went wrong, so they can fix the code. It may occasionally, though very unlikely, be useful for you to verify a calculation, but as mentioned, you will normally not use this option.

    #Grad=(Optimize,PrimInt,BiSect,NoGDIIS)

This line should have told FreeON which methods to use to compute analytical gradients during the calculations. However, since it starts with a hash mark (\#) it will be ignored and taken as a comment. Note that this is a very useful feaure: you can remove options from consideration by commenting them, so that during your normal work, you can try different approximations and keep the options you have tried. This is nice for instance to try various alternatives and save the previous choices while keeping them hidden from FreeON, e. g. so you can remember what you tried. It also implies another interesting feature: you can include comments anywhere in the FreeON file, and use those for documenting your work (e. g. adding explanatory comments).

    SCFMethod=(RH)

This line tells FreeON that we will use the Restricted Hartree-Fock method to solve the Self Consistent Field equations. Note that we can be more specific (say, indicating whether we want to perform an Open- or a Closed-shell calculation) but we need not to, and if we don't, the program will figure the correct option for us.

If you look at other examples, you will notice that you can give several values separated by commas. This is allowed so you can specify different computation methods for different phases of the calculation (e. g. an initial refinement at lower theory levels followed by a subsequent fine structure calculation). In this case there is only one value, meaning that all the calculation will be carried out using the same method in a single phase.

    BasisSets=(6-31G*)

Here is where you state the [basis set](http://en.wikipedia.org/wiki/Basis_set_(chemistry)) to be used for the calculation. In a multi-stage calculation you would state the different basis sets to use at each step separated by commas.

    ModelChem=(B3LYP)

We indicate here the model to use for the calculation. We use the [B3LYP](http://en.wikipedia.org/wiki/Hybrid_functional) method, which tells FreeON to run a [DFT](http://en.wikipedia.org/wiki/Density_functional_theory) calculation using the Becke, three-parameter, Lee-Yang-Parr (B3LYP) exchange-correlation functional.

    Accuracy=(VeryTight)

Next we indicate the precision requested from the calculation, in this case we want to use very tight convergence criteria, that is, the calculation will run until the error in between iterations is very small.

    SCFConvergence=(DIIS)

When iterating in search for the solution we need a method to achieve convergence. The [DIIS (direct inversion in the iterative subspace or direct inversion of the iterative subspace)](http://en.wikipedia.org/wiki/DIIS), also known as Pulay mixing, method is to be used in this calculation.

    Geometry=InAngstrom

Finally, our last option tells FreeON that we are going to provide the geometry using Ångströms as the measure unit.

### Molecular geometry

The last section in this input file defines the initial system geometry we will use to run the calculation,

    <BeginGeometry>
    C 0.0 0.0 0.0
    H 0.6312 0.6312 0.6312
    H -0.6312 -0.6312 0.6312
    H -0.6312 0.6312 -0.6312
    H 0.6312 -0.6312 -0.6312
    <EndGeometry>

The geometry is provided by stating each atom name followed by its X, Y and Z coordinates. The basic format is even simpler than the well known [XYZ format](http://en.wikipedia.org/wiki/XYZ_file_format), as we do not need to provide the number of atoms or the comment line. It is worth noting for one reason: if you have your molecule in XYZ format (as produced by OpenBabel), all you need to do is strip the first two lines (number of atoms and comment) to get a valid geometry for inclusion in this section of FreeON.
