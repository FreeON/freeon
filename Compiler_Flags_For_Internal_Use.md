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

NAG Compiler
------------

    module load nag-5.3

PGI Compiler
------------

    module load pgi-11.8

Intel Compiler
--------------

    module load icc-2011_sp1.9.293
    ./configure CC=icc FC=ifort F77=ifort

Note, that running requires runtime libraries from the Intel compiler, if the binaries are not statically linked against those libraries. Simply sourcing the shell script before running is sufficient.
