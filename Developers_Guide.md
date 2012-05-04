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

We use the version control system [git](http://git-scm.com) to manage the FreeON source code. The master branch is hosted by the [Free Software Foundation](http://fsf.org) on [savanna.nongnu.org](http://savannah.nongnu.org/projects/freeon). If you would like to keep your local sources up to date with our development and/or are considering hacking the code, we recommend you use git yourself. A basic workflow of downloading and building the sources using git could look like the following:

    git clone http://git.savannah.gnu.org/r/freeon.git

which will copy the entire repository into \$PWD/freeon. This step has to be done only once. Now `cd` into freeon.

The command

    git branch

will list all available branches (after a clone operation, there will only be one, the master branch), with the current branch marked with a '\*'. It is good practice to create your own branch for development, so that the master branch stays in sync with the repository on savannah.

Let's assume for now that you want to build the master branch:

    autoreconf -fis

will recreate the necessary autoconf/automake/libtool scripts. This step requires that you have those tools installed on your system, and will fail otherwise. After this step you can configure the sources as described in the [Users Guide\#Downloading and Building the code](Users Guide#Downloading_and_Building_the_code "wikilink").

### Hacking the code

Now let's describe a basic workflow for hacking the code. First you create a new branch (which is local to your machine and is not shared with anyone unless you decide to) with

    git branch branchname

and switch to this branch with

    git checkout branchname

### Validation

Bamboo waterfall, ...

### Working with Git

The version control system (VCS) git is designed to be decentralized and is in this respect very different from more traditional VCS such as CVS and subversion. The fact that we use a server is simply a matter of convenience. Once the repository is cloned it is copied in its entirety and the copy could act as a server.

#### Cloning

Cloning our repository from savannah.

    git clone http://git.savannah.gnu.org/r/freeon.git

#### Creating a Branch

    git branch branchname

### Contributing

Patches, running "make check" etc...
