---
layout: default
title: H2O.inp
---

![This page will show you how to write a configuration file for FreeON to optimize water structure](H2O_ESP.png "This page will show you how to write a configuration file for FreeON to optimize water structure")

Quick optimization of water structure
-------------------------------------

In this example we will see how can we run a quick calculation on water structure. The calculation will use a simple basis set and will be limited to 12 SCF cycles.

### The complete file

You can find this **examples** in the examples subdirectory of the FreeON source code distribution.


     MONDO:
    Frequencies= 2171.0366 4142.8791 4392.9252
                   Q               ln(Q)
      Elec     1.00000E+00       0.000000
      Trans    3.00540E+06      14.915920
      Rot      9.59594E+01       4.563926
      Vib      1.00003E+00       0.000028
      Tot      2.88404E+08      19.479874
                  E         H         G         CV        CP        S
               KJ/Mol    KJ/Mol    KJ/Mol    J/MolK    J/MolK    J/MolK
      Elec      0.000     0.000     0.000     0.000     0.000     0.000
      Trans     3.718     6.197   -36.976    12.472    20.786   144.805
      Rot       3.718     3.718   -11.314    12.472    12.472    50.419
      Vib      64.042    64.042    64.041     0.026     0.026     0.003
      Tot      71.479    73.958    15.751    24.969    33.284   195.227


     GAMESS:
     FREQUENCY:      2170.05     4139.84     4390.92
                   Q               LN Q
     ELEC.     1.00000E+00       0.000000
     TRANS.    3.00431E+06      14.915558
     ROT.      9.59495E+01       4.563822
     VIB.      1.00003E+00       0.000028
     TOT.      2.88270E+08      19.479409
                  E         H         G         CV        CP        S
               KJ/MOL    KJ/MOL    KJ/MOL   J/MOL-K   J/MOL-K   J/MOL-K
     ELEC.      0.000     0.000     0.000     0.000     0.000     0.000
     TRANS.     3.718     6.197   -36.975    12.472    20.786   144.800
     ROT.       3.718     3.718   -11.313    12.472    12.472    50.417
     VIB.      64.006    64.006    64.005     0.026     0.026     0.003
     TOTAL     71.443    73.922    15.717    24.969    33.283   195.220

    <BeginOptions>

    Charge=0
    Multiplicity=1

    G!rad=(NumFreq)
    Grad=(Optimize,PrimInt)

    DebugAll=(MaxDebug,CheckSums)

    MaxSCF=12

    Guess=SuperPos
    SCFMethod=(RH)
    BasisSets=(STO-3G-SPLIT)
    ModelChem=(HF)
    Accuracy=(verytight)
    SCFConvergence=(DIIS)

    <EndOptions>

    <BeginGeometry>
     o   0.95309668               0.00000000               0.63925771
     h   0.00278150               0.00000000               0.91468038
     h   0.84735582               0.00000000              -0.34449809
    <EndGeometry>

### Introductory comments

As you can see, the file is split in several sections. There are two parts clearly delimited by opening (Begin...) and closing (End...) tags, and an initial section that is not delimited by any tag at all. This first section, everything you write outside an appropriate tag, is ignored by FreeON. Hence you can use it to provide comments about the file, what it is about, what it does, etc...

It is advisable that you enter here as much information about the system being modeled and the calculations being run, as possible. Later, when you (or someone else) needs to refer to this file (say, to remember how things were done, or what you did) it will prove most useful to save work and headaches.


     MONDO:
    Frequencies= 2171.0366 4142.8791 4392.9252
                   Q               ln(Q)
      Elec     1.00000E+00       0.000000
      Trans    3.00540E+06      14.915920
      Rot      9.59594E+01       4.563926
      Vib      1.00003E+00       0.000028
      Tot      2.88404E+08      19.479874
                  E         H         G         CV        CP        S
               KJ/Mol    KJ/Mol    KJ/Mol    J/MolK    J/MolK    J/MolK
      Elec      0.000     0.000     0.000     0.000     0.000     0.000
      Trans     3.718     6.197   -36.976    12.472    20.786   144.805
      Rot       3.718     3.718   -11.314    12.472    12.472    50.419
      Vib      64.042    64.042    64.041     0.026     0.026     0.003
      Tot      71.479    73.958    15.751    24.969    33.284   195.227


     GAMESS:
     FREQUENCY:      2170.05     4139.84     4390.92
                   Q               LN Q
     ELEC.     1.00000E+00       0.000000
     TRANS.    3.00431E+06      14.915558
     ROT.      9.59495E+01       4.563822
     VIB.      1.00003E+00       0.000028
     TOT.      2.88270E+08      19.479409
                  E         H         G         CV        CP        S
               KJ/MOL    KJ/MOL    KJ/MOL   J/MOL-K   J/MOL-K   J/MOL-K
     ELEC.      0.000     0.000     0.000     0.000     0.000     0.000
     TRANS.     3.718     6.197   -36.975    12.472    20.786   144.800
     ROT.       3.718     3.718   -11.313    12.472    12.472    50.417
     VIB.      64.006    64.006    64.005     0.026     0.026     0.003
     TOTAL     71.443    73.922    15.717    24.969    33.283   195.220

This file can be used to check the accuracy of FreeON calculations. The developers of the package included therefore the results of running a similar calculation using two other packages: MondoSCF (an earlier version of FreeON) and [GAMESS-US](http://www.msg.ameslab.gov/gamess/). This way, you can verify how well the results you obtained after running this calculation performed in comparion.

### Calculation options

The first actual section containing instructions for FreeON, is the *Options* section. It is enclosed in between the *<BeginOptions>* and *<EndOptions>* tags, and contains the information needed to specify the calculations to carry out.

    Charge=0
    Multiplicity=1

These two lines define the [multiplicity](http://en.wikipedia.org/wiki/Multiplicity_(chemistry)) and charge of the system. We are interested in the ground state of the water molecule, so the number of electrons equilibrates the nuclear charges to a grand total charge of zero for the system. This gives an even number of electrons, all of them paired in the corresponding molecular orbitals, hence there are no spin-uncoupled electrons in the system, leading to a multiplicity of 1.

    G!rad=(NumFreq)
    Grad=(Optimize,PrimInt)

These two lines are interesting in themselves. In both cases we deal with the same option, but we use a trick here for our convenience. FreeON will ignore any option it does not understand. This allows us to use the same file as a template for different calculations by using a simple device: if we write out all the alternatives (two in this case), then we can have all at once in a single file, and then enable or disable the one we want to run by making it *invalid* somehow. In this case, the author wanted to keep the instructions for a frequency calculation and an optimization run in the same file. As we are now interested only in the optimization, we can invalidate the Frequency calculation inserting an special character (a bang ! here) to make it an invalid key word.

There are many ways to include comments in a FreeON configuration file (see other examples), and this is but one more.

    DebugAll=(MaxDebug,CheckSums)

This line tells FreeON to include in the output additional information. You will normally not use this. It is intended mainly for the developers of the software so they can track a calculation and spot where is any problem in the unlikely case that it went wrong, so they can fix the code. It may occasionally, though very unlikely, be useful for you to verify a calculation, but as mentioned, you will normally not use this option.

    MaxSCF=12

With this, we are telling FreeON to limit the number of SCF cycles to 12. This may come handy in complex calculations, where the system has trouble to converge, or where convergence is inherently difficult because the system oscillates between two similar energy configurations (e. g. in aromatic or resonance structures), or when we want to limit the time spend in calculation.

    Guess=SuperPos

Here we tell FreeON how to start the calculation, We could start our search for the optimal coefficients from any random value, but it is better to start from something that is closer to the solution, such as a prior quick Huckel calculation, as it will converge faster to the solution. In this case we are asking FreeON to use the SuperPos ([Superposition of Atomic Densities, SAD](http://igitur-archive.library.uu.nl/chem/2007-0302-200920/pdf18.pdf)) method to generate the starting configuration.

    SCFMethod=(RH)

We are now specifying the kind of calculation: **RH** stands for [Restricted Hartree-Fock](http://en.wikipedia.org/wiki/Hartree–Fock_method). Notice that we are not specifying which kind of restricted Hartree-Fock calculation we want, the program will select the appropriate calculation for us.

    BasisSets=(STO-3G-SPLIT)

This is the [basis set](http://en.wikipedia.org/wiki/Basis_set_(chemistry)) to use.

    ModelChem=(HF)

This option helps fine tune the kind of calculation, We are not further specifying it here, but could if needed (e. g. to specify the type of exchange-correlation functional for a DFT calculation).

    Accuracy=(verytight)

Here we indicate that we want to use very tight convergence criteria, or in other words to run the calculation to a high precision.

    SCFConvergence=(DIIS)

And finally we state the method to use for convergence calculations, which in this case is the [DIIS](http://en.wikipedia.org/wiki/DIIS) (direct inversion of iterated subspace) method.

### Molecule geometry

The final section in this example is the actual molecule geometry:

    <BeginGeometry>
     o   0.95309668               0.00000000               0.63925771
     h   0.00278150               0.00000000               0.91468038
     h   0.84735582               0.00000000              -0.34449809
    <EndGeometry>

The geometry is provided by stating each atom name followed by its X Y and Z coordinates. This is a very simple format. It is even simpler than the well known [XYZ format](http://en.wikipedia.org/wiki/XYZ_file_format), as we do not need to provide the number of atoms or the comment line. It is worth noting for one reason: if you have your molecule in XYZ format (as produced by [OpenBabel](http://en.wikipedia.org/wiki/OpenBabel)), all you need to do is strip the first two lines (number of atoms and comment) to get a valid geometry for inclusion in this section of the FreeON configuration file.

Note as well that FreeON does not care if the atom name is written in capitals or not. It will take and understand it all the same.

--[Jrvalverde](User:Jrvalverde "wikilink") 08:23, 28 April 2012 (MDT)
