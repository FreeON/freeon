---
layout: default
title: Developers Guide
---

Structure of the Program Suite
------------------------------

FreeON is split into a front-end and a back-end. The back-end involves numerically challenging tasks like matrix-multiplication and matrix-build operations. Back-end operations can always be run in stand-alone mode for tuning and debugging very large problems. In serial, the backend is spawned with calls to [`execvp()`](http://linux.die.net/man/3/execvp), and with parallel clones, through [`MPI_COMM_SPAWN()`](http://linux.die.net/man/3/mpi_comm_spawn).

This structure significantly simplifies the tuning and debugging of large calculations through access to the stand alone problem. Perhaps the greatest strength of this framework though, is that it enables rapid development and competition between algorithms.

### [APIs](APIs "wikilink")

Clones
------

NEB, parallel replica ..

The HDF5 Archive
----------------

Data is passed on the argument line to back-end subprograms, and structures are communicated via an archived HDF file.

Generalized *N*-Body Solvers
----------------------------

### Solver Stack

### Independent Solver Libraries

-   [SpAMMPACK](SpAMMPACK "wikilink")
-   Quantum Coulomb Tree Code

Task Parallelism
----------------

Code Base
---------

FreeON is hosted by the [Free Software Foundation](http://www.fsf.org/) at [<http://savannah.nongnu.org/projects/freeon>](http://savannah.nongnu.org/projects/freeon).

### Downloading and Building

We use the version control system [git](http://git-scm.com) to manage the FreeON source code. The master branch is hosted by the [Free Software Foundation](http://fsf.org) on [savanna.nongnu.org](http://savannah.nongnu.org/projects/freeon). If you would like to keep your local sources up to date with our development and/or are considering hacking the code, we recommend you use git yourself. A basic workflow for downloading and building the sources using git could look like the following:

    git clone http://git.savannah.gnu.org/r/freeon.git

which will copy the entire repository into \$PWD/freeon. This step has to be done only once.

### Validation

Bamboo waterfall, ...

### Working with Git

Simple commands. Building your own repository. Creating branches. Steal our code, please.

### Contributing

Patches, running "make check" etc...
