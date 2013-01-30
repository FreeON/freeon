---
layout: default
title: Main Page
---

\_\_NOTOC\_\_ \_\_NOTITLE\_\_

Introduction
------------

FreeON is an experimental, open source ([GPL](http://www.gnu.org/licenses/gpl.html)) suite of programs for linear scaling quantum chemistry. It is highly modular, and has been written from scratch for *N*-scaling SCF theory in Fortran95 and C. Platform independent I/O is supported with [HDF5](http://www.hdfgroup.org/HDF5/). FreeON should compile with most modern Linux distributions and OS X. FreeON performs Hartree-Fock, pure Density Functional, and hybrid HF/DFT calculations (e.g. B3LYP) in a Cartesian-Gaussian LCAO basis. All algorithms are *O(N)* or ''O(N *log* N)'' for non-metallic systems. Periodic boundary conditions in 1, 2 and 3 dimensions have been implemented through the Lorentz field (Γ-point), and an internal coordinate geometry optimizer allows full (atom+cell) relaxation using analytic derivatives. Effective core potentials for energies and forces have been implemented. Advanced features include *O(N)* static and dynamic response, as well as time reversible Born Oppenheimer Molecular Dynamics (MD).

Contents
--------

-   [History](History "wikilink")
-   [Vision](Vision "wikilink")
-   [Status](http://www.freeon.org:8010)
-   [Users](Users Guide "wikilink")
-   [Developers](Developers Guide "wikilink")
-   [Collaboratory](Collaboratory "wikilink")
-   [Solvers](Solvers "wikilink")
-   [Publications](Publications "wikilink")
-   [Citation](Citing FreeON "wikilink")
-   [Authors](Authors "wikilink")

Download
--------

We periodically release snapshots of the development tree in the [Savannah FreeON Download Area](http://savannah.nongnu.org/files/?group=freeon)

Here are a few recent ones:

-   [freeon-beta-2012-05-01.tar.bz2](http://download.savannah.gnu.org/releases/freeon/freeon-beta-2012-05-01.tar.bz2)
-   [freeon-beta-2012-05-03.tar.bz2](http://download.savannah.gnu.org/releases/freeon/freeon-beta-2012-05-03.tar.bz2)

Project Stats
-------------

<ohlohStats/>

Contact Us
----------

While you could e-mail one of the authors directly, we prefer if you used our mailing lists or our bug tracker instead. You will probably get a response quicker this way because your message will go to all of us as opposed to just one of us. We currently operate 2 mailing lists:

-   [freeon-users](http://lists.nongnu.org/mailman/listinfo/freeon-users), our public mailing list for users of FreeON.
-   [freeon-devel](http://lists.nongnu.org/mailman/listinfo/freeon-devel), our mailing list for developers of FreeON, and people who are interested in the development process of it.

Alternatively you can file a bug report or feature request with our issue tracker system:

-   [Issue Tracker](https://github.com/FreeON/freeon/issues)

Wikipedia
---------

We are also on wikipedia: [FreeON](http://en.wikipedia.org/wiki/FreeON)
