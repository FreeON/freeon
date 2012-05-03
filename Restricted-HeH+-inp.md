---
layout: default
title: Restricted-HeH+-inp
---

![This example performs an RHF calculation on the [Helium hydride ion](http://en.wikipedia.org/wiki/Helium_hydride_ion) (HeH<sup>+\>) using [FreeON](http://freeon.org) ](HeH+_MO.png "This example performs an RHF calculation on the Helium hydride ion (HeH+>) using FreeON ")

Helium Hydride Ion
------------------

This example performs an RHF calculation on the [Helium hydride ion](http://en.wikipedia.org/wiki/Helium_hydride_ion) (HeH<sup>+\>) using [FreeON](http://freeon.org). Referece results calculated with MonoSCF and GAMESS-US as well as the GAMESS-US configuration file used to run the equivalent calculation are provided as comments.

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
