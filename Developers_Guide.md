---
layout: default
title: Developers Guide
---

Structure of the Program Suite
------------------------------

FreeON is split into a front-end and a back-end. The back-end involves numerically challenging tasks like matrix-multiplication and matrix-build operations. Back-end operations can always be run in stand-alone mode for tuning and debugging very large problems. In serial, the backend is spawned with calls to [`execvp()`](http://linux.die.net/man/3/execvp), and with parallel clones, through `MPI_COMM_SPAWN()`.

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

Git, ...

### Validation

Bamboo waterfall, ...

### Working with Git

Simple commands. Building your own repository. Creating branches. Steal our code, please.

### Contributing

Patches, running "make check" etc...
