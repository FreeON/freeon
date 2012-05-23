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

#### Buildbot

We continuously build and test FreeON for verification, validation, and regression testing. You can see for yourself how FreeON is doing right now. Go to our [buildbot](http://www.freeon.org:8010) webpage.

#### Validating your local version

    make check

and

    make validate

### Working with Git

The version control system (VCS) git is designed to be decentralized and is in this respect very different from more traditional VCS such as CVS and subversion. The fact that we use a server is simply a matter of convenience. Once the repository is cloned it is copied in its entirety and the copy could act as a server.

#### Cloning

Cloning our repository from savannah.

    git clone http://git.savannah.gnu.org/r/freeon.git

#### Creating a Branch

    git branch branchname

### Contributing

If you decide to work with git, then we suggest you do the following. You first clone the repository as described above and then create your own topic branch. You can call this branch anything you like, it will only exist locally. Checkout that new branch and change whatever you feel like. You should periodically update your local branch with changes committed on the repository's master branch by running

    git pull --rebase master

### Eclipse

The [eclipse IDE](http://www.eclipse.org/) offers a powerful and convenient option for developing FreeON. Recommended plugins and settings can be found [here](Eclipse Settings "wikilink").

.adslot-overlay {position: absolute; font-family: arial, sans-serif; background-color: rgba(0,0,0,0.65); border: 2px solid rgba(0,0,0,0.65); color: white !important; margin: 0; z-index: 2147483647; text-decoration: none; box-sizing: border-box; text-align: left;}.adslot-overlay-iframed {top: 0; left: 0; right: 0; bottom: 0;}.slotname {position: absolute; top: 0; left: 0; right: 0; font-size: 13px; font-weight: bold; padding: 3px 0 3px 6px; vertical-align: middle; background-color: rgba(0,0,0,0.45); text-overflow: ellipsis; white-space: nowrap; overflow: hidden;}.slotname span {text-align: left; text-decoration: none; text-transform: capitalize;}.revenue {position: absolute; bottom: 0; left: 0; right: 0; font-size: 11px; padding: 3px 0 3px 6px; vertial-align: middle; text-align: left; background-color: rgba(0,0,0,0.45); font-weight: bold; text-overflow: ellipsis; overflow: hidden; white-space: nowrap;}.revenue .name {color: \#ccc;}.revenue .horizontal .metric {display: inline-block; padding-right: 1.5em;}.revenue .horizontal .name {padding-right: 0.5em;}.revenue .vertical .metric {display: block; line-height: 1.5em; margin-bottom: 0.5em;}.revenue .vertical .name, .revenue .vertical .value {display: block;}.revenue .square .metric, .revenue .button .metric {display: table-row;}.revenue .square .metric {line-height: 1.5em;}.revenue .square .name, .revenue .square .value, .revenue .button .value {display: table-cell;}.revenue .square .name {padding-right: 1.5em;}.revenue .button .name {display: block; margin-right: 0.5em; width: 1em; overflow: hidden; text-overflow: clip;}.revenue .button .name:first-letter {margin-right: 1.5em;}a.adslot-overlay:hover {border: 2px solid rgba(58,106,173,0.9);}a.adslot-overlay:hover .slotname {border-bottom: 1px solid rgba(81,132,210,0.9); background-color: rgba(58,106,173,0.9);}a.adslot-overlay:hover .revenue {border-top: 1px solid rgba(81,132,210,0.9); background-color: rgba(58,106,173,0.9);}div.adslot-overlay:hover {cursor: not-allowed; border: 2px solid rgba(64,64,64,0.9);}div.adslot-overlay:hover .slotname {border-bottom: 1px solid rgba(128,128,128,0.9); background-color: rgba(64,64,64,0.9);}div.adslot-overlay:hover .revenue {border-top: 1px solid rgba(128,128,128,0.9); background-color: rgba(64,64,64,0.9);}
