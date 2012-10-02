---
layout: default
title: Compiler Flags For Internal Use
---

GNU gcc Compiler
----------------

    ./configure

Lahey Compiler
--------------

    module load lahey-8.10b
    ./configure FC=lf95 F77=lf95 LIBS=-lblas\\ -llapack

NAG Compiler
------------

    module load nag-5.3
    ./configure FC=nagfor F77=nagfor LIBS=-lblas\\ -llapack

PGI Compiler
------------

    module load pgi-11.8
    ./configure CC=pgcc FC=pgf90 F77=pgf90 LIBS=-lblas\\ -llapack

Intel Compiler
--------------

    module load icc-2011_sp1.9.293
    ./configure CC=icc FC=ifort F77=ifort LIBS=-lblas\\ -llapack

Note, that running requires runtime libraries from the Intel compiler, if the binaries are not statically linked against those libraries. Simply sourcing the shell script before running is sufficient.
