---
layout: default
title: Restricted-HeH+-inp
---

![This example performs an RHF calculation on the [Helium hydride ion](http://en.wikipedia.org/wiki/Helium_hydride_ion) (HeH<sup>+\>) using [FreeON](http://freeon.org) ](HeH+_MO.png "This example performs an RHF calculation on the Helium hydride ion (HeH+>) using FreeON ")

Helium Hydride Ion
------------------

This example performs an RHF calculation on the [Helium hydride ion](http://en.wikipedia.org/wiki/Helium_hydride_ion) (HeH<sup>+</sup>) using [FreeON](http://freeon.org). Referece results calculated with MondoSCF and GAMESS-US as well as the GAMESS-US configuration file used to run the equivalent calculation are provided as comments.

### The complete input file

You can find this example configuration file under directory **Validate/SinglePoint/Restricted** of the FreeON source code distribution.


    Mondo                         GAMESS VERSION = 24 MAR 2007 (R3)

    <SCF>   = 1.36686714         ITER EX DEM     TOTAL ENERGY        E CHANGE  DENSITY CHANGE    DIIS ERROR
    <SCF>   = -2.78031168          1  0  0       -2.7803116840    -2.7803116840   0.294054433   0.000000000
    <SCF>   = -2.88688541          2  1  0       -2.8868854101    -0.1065737261   0.013459905   0.000000000
    <SCF>   = -2.88737780          3  2  0       -2.8873777957    -0.0004923856   0.002609729   0.000000000
    <SCF>   = -2.88738583          4  3  0       -2.8873858238    -0.0000080282   0.000617305   0.000000000
    <SCF>   = -2.88738612          5  0  0       -2.8873861149    -0.0000002910   0.000163981   0.000000000
    <SCF>   = -2.88738613          6  1  0       -2.8873861270    -0.0000000122   0.000000432   0.000000000
    <SCF>   = -2.88738613          7  2  0       -2.8873861270     0.0000000000   0.000000076   0.000000000
    <SCF>   = -2.88738613




    GAMESS input file:



     $CONTRL
      EXETYP=RUN
      SCFTYP=RHF
      RUNTYP=GRADIENT
      MPLEVL=0
      COORD=CART
      UNITS=BOHR
      ICHARG=1
      MULT=1
      NPRINT=9
     $END
     $SYSTEM
      TIMLIM=1
      MEMDDI=1
      PARALL=.TRUE.
     $END
     $BASIS
      GBASIS=N21
      NGAUSS=3
     $END
     $GUESS
      GUESS=HCORE
      PRTMO=.TRUE.
     $END
     $SCF
      DIIS=.TRUE.
      ETHRSH=10.0
     $END
     $DATA
    HeH+....
    C1      0
    HELIUM     2.0   0.0    0.0  0.0
    HYDROGEN   1.0   1.4632 0.0  0.0
     $END



    <BeginOptions>

    Charge=1
    Multiplicity=1

    DebugAll=(MaxDebug,MatDebug,CheckSums)

    Guess=Core
    SCFMethod=(RH)
    BasisSets=(3-21G-SPLIT)
    ModelChem=(HF)
    Accuracy=(verytight)
    SCFConvergence=(DIIS)
    DIISDimension=0
    Geometry=(InAU)

    <EndOptions>

    <BeginGeometry>
     He  0.0    0.0 0.0
     H   1.4632 0.0 0.0
    <EndGeometry>

### Prologue

Before actually entering in the meet of the configuration, the file starts with some free text. As it is not enclosed between any <Begin...> and <End...> tags, it is ignored by FreeON. This provides us a neat way to include comments about the calculation in the configuration file.

In this case, the content shows the results obtained after several iterations of the same calculation performed with MondoSCF (a previous version of FreeON) and GAMESS-US.

These sample results are followed by the input configuration file used to run a similar calculation with GAMESS-US. We will not enter in the details here, but if you are interested, they can be easily looked up in GAMESS-US documentation, which is available on the web.

### Options section

Next comes the options section. Input to FreeON is divided in sections, each labeled by tags enclosed in less-than and greater-than signs. Every section in the file deals with a specific type of information: the options section allows you to enter your choices regarding the calculation and subsequent input:

    <BeginOptions>

    Charge=1
    Multiplicity=1

    DebugAll=(MaxDebug,MatDebug,CheckSums)

    Guess=Core
    SCFMethod=(RH)
    BasisSets=(3-21G-SPLIT)
    ModelChem=(HF)
    Accuracy=(verytight)
    SCFConvergence=(DIIS)
    DIISDimension=0
    Geometry=(InAU)

    <EndOptions>

You can get more detailed information in the manual and, possibly, in other examples.

    Charge=1
    Multiplicity=1

In this case, the file starts by defining the electronic state of the system, which has a net charge of one. This would correspond to a bare proton, Since there are no singly occupied orbitals, the multiplicity is one.

    DebugAll=(MaxDebug,MatDebug,CheckSums)

This line requests additional output. You will not normally use it, unless you find yourself in trouble and want to know more details about the computation to find out the root of the problem and possibly how to fix it.

    Guess=Core
    SCFMethod=(RH)
    BasisSets=(3-21G-SPLIT)
    ModelChem=(HF)
    Accuracy=(verytight)
    SCFConvergence=(DIIS)
    DIISDimension=0
    Geometry=(InAU)

These lines define the kind of calculation to perform: **Guess=core** asks FreeON to diagonalize the core Hamiltonian in order to get the initial guess. The quality of the initial guess usually has a major impact on the speed at which the program can reach the solution in the iterative process. There are various other alternatives, but this is a good one for atomic calculations.

The following options request a Hartree-Fock family SCF calculation using the 3-21G-split basis set, using plain Hartree-Fock and iterating the SCF procedure until a very tight convergence is achieved (i. e. the energy change between one iteration and the next is very small), Convergence will use the Pulay solver (DIIS: direct inversion of iterated subspace) with a dimension of 0.

The last line prepares the next section and indicates that the geometry should be interpreted as stated in atomic units.

### The geometry section

Enclosed between the tags **<BeginGeometry>** and **<EndGeometry>** comes the definition of the system to be modeled: its format is very simple, it consists of the atom name followed by its XYZ coordinates, one atom per line. This format is similar enough to the well known XYZ format that it is trivial to convert from it (just remove the first two lines of an XYZ file leaving only the coordinates).

    <BeginGeometry>
     He  0.0    0.0 0.0
     H   1.4632 0.0 0.0
    <EndGeometry>
