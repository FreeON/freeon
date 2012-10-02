---
layout: default
title: Compiler Flags For Internal Use
---

GNU gcc Compiler
----------------

    ./configure

Intel Compiler
--------------

    source /cluster/intel/composer_xe_2011_sp1.9.293/bin/compilervars.csh intel64
    ./configure CC=icc FC=ifort F77=ifort

Note, that running requires runtime libraries from the Intel compiler, if the binaries are not statically linked against those libraries. Simply sourcing the shell script before running is sufficient.
