---
layout: default
title: ICQCB Available software
---

Available software
------------------

### But now, man, like, uh, really, in practice, how do we go about this stuff?

Normally, you will rely on a graphical user interface to build your molecule or you system, and most often you will start by running a Molecular Mechanics or a Molecular Dynamics simulation to get an *educated* guess of the initial conformation of the system.

As an alternative you can get a starting structure from a database when it is known (or can readily be predicted), like **PDB**, **Zinc**, **Corina**, and many others.

Once we have our system set up (with its initial and final states if it is a reaction), the rest is a lot easier: there are many programs you can use to carry out the calculations, and all you need to do is specify what you want to do, letting the computer struggle with the problem. Often, you will be able to launch the calculation from a graphical user interface and to choose among several alternative programs (like with *Gabedit*, *Gdis* or *WebMO*).

The list of available programs is growing very quickly. The increased interest in computer simulations, partly driven by availability of linear scaling codes and of more capable computers has lead many groups to develop new algorithms and implement them often in their own program. Among the more popular and the ones that have been around for longer we could highlight (and keep in mind this is not exhaustive):

-   **FreeON** is distributed under the GNU license (i. e. you can get the source code in Fortran an C and modify it) and provides access to a large variety of calculations, both *ab-initio* and DFT, using linear scaling methods that make it most efficient. It can be applied to ground state calculations, to minimization problems, time-dependent calculations, and Quantum Molecular Dynamics using any of the methods available.

-   **GAMESS-US** is a program that is distributed for free and with source code, including a wide number of methods, bith *ab-initio* and semi-empirical, and in recent versions, also DFT. It is widely popular and can run on clusters. It has limited support for Molecular Dynamics and cannot perform mixed QM/MM calculations using popular biological oriented force fields, although there is an extension (*SIMOMM*) to integrate it with *Tinker*. A spin-off derivative is **Firefly/PC-Gamess**, which started as an optimized version for Intel architecture computers, but has now become fully independent although it maintains a mostly compatible input language.

-   **Gaussian** is doubtless the most popular program used in Chemistry, but it has a problem that is becoming a major hindrance: it is commercial, and by the way, terribly expensive. There are also reports of a discriminatory treatment of some relevant scientists that have used it to compare results or that work in competing codes. With the heavy calculations currently carried out today, and an very expensive license paid by CPU (and modern computers and clusters can have a lot of them), it is not the best option except for deep pockets aristocratic groups. On the bonus side, it can do almost anything you want and has many user contributed extensions.

-   **MOPAC** is doubtles the king of semi-empirical methods, at least in Biocomputing, to date. MOPAC uses parameter sets adapted for biochemical systems, can compute reaction paths, and is available in two main flavours: traditional MOPAC (up to version 7) was free software and distributed with source code. Modern MOPAC (versions 2000 and later) includes support for linear scaling calculations (MOZYME) and is distributed in executable form, being available for free for academic use.

-   **MPQC** the Massively Parallel Quantum Chemistry program, is distributed under GPL and supports a large number of *ab initio* calculations. It is written in C++ (whereas most other programs are written in various versions of Fortran), and has been designed to scale in multi-CPU computers and to be easy to integrate in other programs.

-   **NWchem** used to be free for academic use after registration, but is now distributed under a liberal license and can readily be got together with its source code without significant limitations. It can be integrated in a separate graphical user interface, **ECCE** and includes a wide range of methods, including mixed QM/MM and QM/MD models with biomolecular force fields. It can be run in clusters and grids.

-   **Orca** is another popular, rather complete, program that can carry out a wide variety of *ab initio* calculations. It is distributed for free to academics in executable form.

-   **SIESTA** or Spanish Initiative for for Electronic Simulations with Thousands of Atoms, was the first linear scaling DFT code to become widely popular. It can perform DFT calculations on the ground state, optimizations and Molecular Dynamics calculations. It is distributed for free, with source code, but requires registration.

All of these programs work in a similar way: since calculations may potentially take a very long time to complete, you do not want to tie your computer screen to the program for weeks or more. Instead, you will normally want to run them in a remote computer without any kind of permanent link to your desktop. The idea is that you will connect to said remote server, specify the calculation, disconnect, and then later, once it is finished recover the results for inspection. This rules out graphical interfaces for performing advanced QM calculations (they are still OK for small problems, but may be inconvenient for serious, large work).

So, in practice, what you will normally do is write down on a plain text file all you want to do, together with a description of the system to model. Then you will run the program telling it the name of this file, and forget about it until all the work is done. Usually, our main problem will be that there is no common format for this file, it changing from one program to the other, and so, we will need to learn how to specify our ideas to each different program. Some GUI environments can generate these files for us using a menu and dialog driven interface, but these are necessarily limited (they would be too complex if they were to include all options) and so some knowledge of how these files are written is always useful to fine tune our calculations.

Let us have a look at a few example input files to see how they look and the usual differences between them:

#### Sample FreeON input file

This input file will run a Quantum Molecular Dynamics simulation of four water molecules using FreeON. The calculation will use a minimal basis set, periodic boundary conditions and simulate 2000 steps.

As you can see, the coordinates of composing atoms are stated in a separate section as the atom symbols followed by their coordinates. And a simulation needs not be limited to a single molecule: it can encompass complex molecular systems as well.

    <BeginOptions>

    Charge=0
    Multiplicity=1

    Guess=SuperPos
    SCFConvergence=(ODA,DIIS)
    SCFMethod=(TRS4,TRS4)
    BasisSets=(STO-2G-SPLIT,STO-3G-SPLIT) 
    ModelChem=(HF,HF)
    Accuracy =(Loose,Good)

    Grad=(MolecularDynamics)
    InitTemp=100
    MDMethod=Verlet
    DeltaTime=0.5
    MaxMDStep=2000
    DMPOrder=2
    MinSCF=5
    MaxSCF=12

    PBC=(T,T,T)
    Periodic=(AtomCoords)
    <EndOptions>

    <BeginPeriodic>
     11.7124  11.7124  11.7124   90.00000   90.00000   90.00000 
    <EndPeriodic>

    <BeginGeometry>
      O       -3.41700        5.21400       -3.21100
      H       -4.04800        5.25600       -2.43700
      H       -3.10700        4.27100       -3.33600
      O        4.80800       -5.92700       -4.36200
      H        4.16600       -5.29200       -4.79200
      H        5.50800       -6.19200       -5.02500
      O        0.41400        1.88600        5.68900
      H        1.34900        1.80700        6.03400
      H        0.32900        1.36900        4.83700
      O       -5.03300       -6.23700       -5.37800
      H       -4.36700       -6.68900       -4.78600
      H       -4.65400       -5.37100       -5.70400
    <EndGeometry>

#### Sample GAMESS-US input file

This is a sample configuration file to minimize the energy of glycerin using RHF and PM3 with GAMESS-US

Like the former file, there options and coordinates are separate. Options are distributed in sections (with a different delimiter convention).

    $CONTRL SCFTYP=RHF RUNTYP=ENERGY ICHARG=0 MULT=1 COORD=ZMTMPC $END
     $BASIS GBASIS=PM3 $END
     $DATA
    Glycerin-pm3
    C1 1
    C 0.0000000 0 0.0000000 0 0.0000000 0 0 0 0
    C 1.5270663 1 0.0000000 0 0.0000000 0 1 0 0
    C 1.5210722 1 112.22985 1 0.0000000 0 2 1 0
    O 1.3986812 1 108.51506 1 -170.4171 1 3 2 1
    H 0.9471724 1 109.41726 1 -169.5354 1 4 3 2
    O 1.4000063 1 111.33784 1 -170.1437 1 1 2 3
    H 0.9486114 1 109.37466 1 -43.44499 1 6 1 2
    O 1.4019702 1 111.06335 1 -45.28770 1 2 1 6
    H 0.9464520 1 109.83009 1 -22.03894 1 8 2 1
    H 1.0873470 1 110.19312 1 66.070572 1 1 2 3
    H 1.0867309 1 108.39288 1 71.396516 1 2 1 6
    H 1.0905230 1 109.65714 1 -49.62054 1 3 2 1
    H 1.0871035 1 109.10518 1 68.624779 1 3 2 1
    H 1.0869613 1 110.97920 1 -52.58142 1 1 2 3
     $END 

#### Sample MOPAC input file

This is a sample configuration file to optimize the calculation of a transition state (TS) for Diels-Adler reaction with MOPAC using the AM1 parameter set.

This is a simple format: the first line contains all options (may be extended using specially tagged, continuation lines), the second line is a comment and then come the atoms and coordinates for the system.

    AM1 T=3600 TS PRECISE AIGOUT NOINTER
     starting point for ts optimization
     C 0.00000000 0 0.0000000 0 0.0000000 0 0 0 0 -0.1780
     C 1.40208483 1 0.0000000 0 0.0000000 0 1 0 0 -0.1609
     C 1.42259991 1 109.1668993 1 0.0000000 0 2 1 0 -0.0875
     C 2.15000000 1 95.2886004 1 72.3355385 1 3 2 1 -0.2117
     C 1.39345295 1 102.7414566 1 -68.9908735 1 4 3 2 -0.2003
     C 2.05000000 1 102.9376995 1 0.5014035 1 5 4 3 -0.1265
     C 1.52065200 1 106.6468717 1 -17.8627955 1 3 2 1 -0.1270
     H 1.08931724 1 125.8188251 1 -174.9635358 1 1 2 3 0.1381
     H 1.08951451 1 125.7796415 1 174.6084539 1 2 1 3 0.1376
     C 1.46768153 1 126.0338633 1 -174.4468082 1 3 2 1 -0.1710
     H 1.09908711 1 91.8502269 1 169.2505342 1 4 3 2 0.1030
     H 1.09782732 1 94.6228112 1 53.7146024 1 4 3 2 0.1057
     H 1.10096177 1 119.6664730 1 -102.0221978 1 5 4 3 0.1005
     H 1.09996024 1 119.9161149 1 105.9591396 1 5 4 3 0.1036
     H 1.08959020 1 125.3591475 1 172.7286299 1 6 1 2 0.1294
     H 1.11277814 1 113.6664178 1 149.3557267 1 7 3 2 0.1002
     H 1.12081000 1 110.3418915 1 -88.5033910 1 7 3 2 0.1045
     H 1.11840183 1 110.7890571 1 -174.2195130 1 10 3 2 0.0781
     H 1.11969567 1 110.6708390 1 65.7090700 1 10 3 2 0.0797
     H 1.11846345 1 110.5804688 1 -53.9527720 1 10 3 2 0.0825

### Glue programs

In addition to the programs that carry out the actual calculations, there are many others that although not able to do the calculations themselves, provide a nice user interface to facilitate creating the input files and running and controlling the calculations, or even to combine different methods. A few well known programs are:

-   **Gabedit** is a GPL'ed graphical user interface to build models and analyze cacluation results, It does not carry out the calculations itself, but relies on other programs to compute them. Gabedit can create input files and read output results for a growing variety of programs, which can be run locally or remotely on a server, all through the GUI. It is a good starting point for the most common calculations.

-   **Gdis** is another GPL'ed GUI which is easy to use and can drive calculations with other programs generating the input files and reading the results.

-   **Ghemical** is a graphical user interface that allows you to build models, perform calculations and visualize calculation results. It integrates tighly other codes like MOPAC or MPQC and is getting to be able to do mixed QM/MM calculations. It is easy to use and GPL, but requires you to do all the work through the GUI tying your screen to the calculations.

-   **Gromacs** is a GPL modeling system for Molecular Mechanics and Molecular Dynamics, likely the most popular one nowadays. In the latest versions it includes the ability to carry out mixed calculations using MOPAC, MPQC or CPMD. It is easy to use, but it is command-line and input-file driven, like CQC codes because it is designed to be used for very large (and long lasting) calculations on remote servers.

-   **Molden** is mainly a visualization tool, but has the ability to generate input files for some popular CQC codes. It is very popular and freely available in source an executable form.

-   **TRITON** is a program combining MODELLER, AUTODOCK and MOPAC. Modeller enables you to create homology models of protein structures (e. g. design *in silico* mutants and predict how the structure is affected), then Autodock will enable you to place ligands in the protein pocket, and finally, MOPAC can be used through a driving program to simulate the reaction process in the active site. The work environment is fully graphic, intuitive, easy and well documented, but its source code is not available.

There are many more programs, and the list of GUI environments and visualization systems is endless, but instead of being comprehensive, we think it is better to simply run a search on the Internet to locate the one you like better.

### Data preparation: Babel

At this point you may be wondering: how do I prepare my data if I do not have a nifty GUI program to do it for me? Well, as you have seen, expressing the calculation options is usually not difficult (a simple matter of reading the manual). However, there is an element that is not trivial, where you can easily make mistakes, is terribly boring and would take an inordinate amount of time for anything but the simpler models.

Or course we refer to entering the coordinates of atoms in the system. Each program will require them in a different format and changing them, atom by atom and coordinate by coordinate easily becomes a real pain (nay, torture!). Even if you have the data at hand (the most common case) it will most likely not be in the format you need (usually, in Biology, it will be in PDB format or mol2) and you will need to convert,

Enter **openbabel**, a GPL tool that is the workhorse of all Computational Chemistry and Biology. Babel can convert between a huge number of different formats (and even often fix the structure for you). Using Babel you can forget about file formats for specifying coordinates in most cases (there are still a few uncovered exceptions, but the most popular ones are there).

### Finally, the actual workflow, at last!

We finally reach what you were looking for: **how do we actually, really, in practice, carry out a calculation**. And in practice all it reduces to using available programs.

To actually run a calculation all we need is a description of the system structure and a way to specify the options.

The first step will usually be to obtain a system description, the molecular and atomic coordinates, normally gathered from a database like PDB or Zinc.

Usually, whenever there is one available, we will start from a *graphical user interface* program (ECCE, TRITON, WebMO, Gabedit, etc...) to obtain our initial input file with an easy menu-driven interface. Often we will be allowed to review the generated file and tweak it by hand adding, deleting or changing options if we want.

If we cannot find a suitable GUI, or if we are familiar with the program, we can directly use the programs from the command line (FreeON, ErgoSCF, NWchem, MPQC, MOPAC, GAMESS-US, etc...). The easiest way will be to take the system description file (the coordinates) and modify it with any text editor to add the calculation options we need.

The full process then, comes down to this:

-   Design or obtain the model system (the molecule/s)
-   Convert the file to the desired format using OpenBabel
-   Add to the file the desired options, either directly through a text editor or indirectly through a GUI
-   Run the chosen program
-   Wait...
-   Wait...
-   Wait...
-   And use a visualization program to have some fun after all this time (and quickly verify the results) although most often we will be interested in an specific numeric value.

Some popular visualization programs often used in Chemistry are:

-   Ghemical
-   GDIS
-   Gabedit
-   gOpenMol
-   Molden
-   OrbVis
-   QMView
-   VisMol

