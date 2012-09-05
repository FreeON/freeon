---
layout: default
title: Developers Guide
---

Getting Started
---------------

FreeON is hosted by the [Free Software Foundation](http://www.fsf.org/) at [<http://savannah.nongnu.org/projects/freeon>](http://savannah.nongnu.org/projects/freeon) and by **[github](http://github.com) at [<http://github.com/FreeON/FreeON>](http://github.com/FreeON/FreeON)**.

Note that due to a technical issue the savannah repository only contains the commit history back to 2008, while the github repository contains the complete commit history back to 2000. Until savannah fixes this problem **we recommend using the github repository**.

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

### Working with git

The version control system (VCS) git is designed to be decentralized and is in this respect very different from more traditional VCS such as CVS and subversion. The fact that we use a server is simply a matter of convenience. Once the repository is cloned it is stand alone and can be used independently as a server. This feature allows extremely flexible development models involving networks and hierarchies of collaborators.

#### Cloning

Cloning our repository from savannah.

    git clone http://git.savannah.gnu.org/r/freeon.git

#### Creating local branches

Anyone is welcome to branch (fork) their own version of FreeON, and use it as they like within the GPL. To get started, we suggest the following. Clone the repository as described above and then create your own topic branch with

    git branch branchname

and switch to this branch with

    git checkout branchname

You can call this branch anything you like, it will only exist locally on your machine. Checkout the new branch and change whatever you feel like. Add your changes and new files to your local branch with

    git add filename

then commit them with

    git commit

Commit often and include descriptive commit messages to document your progress.

Also, you should periodically update your local branch with changes committed on the repository's master branch by running

    git pull --rebase

to make sure you don't get too far out of synch with the progress of others.

#### Contributing

In addition, we encourage the submission of new functionality for inclusion in the official branch, with the following stipulations:

1.  The code is understandable, commented and style conforming: (a) Fortran keywords are all capital. (b) Subroutines and Functions use mixed case naming. (c) Don't be afraid of long names, something like `i`, `ii`, `iii` is not acceptable.
2.  The submitted code is properly derived (GPL'able or GPL'd with attribution)
3.  The changes do not break validation
4.  Example files for regression testing are provided.

Once you are ready to contribute, update your branch a final time with master and then create patches with the following commands:

    git pull --rebase
    git format-patch master

Then, kindly send a brief description along with the attached patches in an email to <freeon-devel@nongnu.org>.

### Validation

#### Validating your local version

After building your local version, we suggest you test that it works properly in your local environment. You can do this with the commands:

    make check

and

    make validate

#### Buildbot

We continuously build and test FreeON for a variety of environments shown on the buildbot page [buildbot](http://www.freeon.org:8010) webpage. Note that we can rapidly test patches and work with authors to improve their submissions.

Structure of the Program Suite
------------------------------

FreeON is split into a front-end and a back-end. The back-end involves numerically challenging tasks like matrix-multiplication and matrix-build operations. Back-end operations can always be run in stand-alone mode for tuning and debugging. In serial, the backend is spawned with calls to [`execvp()`](http://linux.die.net/man/3/execvp), and with parallel clones, through [`MPI_COMM_SPAWN()`](http://linux.die.net/man/3/mpi_comm_spawn).

This structure significantly simplifies the tuning and debugging of large calculations through access to the stand alone problem. Perhaps the greatest strength of this framework though, is that it enables rapid development and competition between algorithms.

### API

We have started to document the [FreeON API](http://www.freeon.org/FreeON_API) in detail. This is still at an early stage and most APIs are missing, but we are rapidly adding information.

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

### Eclipse

The [eclipse IDE](http://www.eclipse.org/) offers a powerful and convenient option for developing FreeON. Recommended plugins and settings can be found [here](Eclipse Settings "wikilink").
