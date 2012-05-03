---
layout: default
title: Main Page
---

\_\_NOTOC\_\_

Introduction
------------

FreeON is an experimental, open source (GPL) suite of programs for linear scaling quantum chemistry. It is highly modular, and has been written from scratch for N-scaling SCF theory in Fortran95 and C. Platform independent I/O is supported with HDF5. FreeON should compile with most modern Linux distributions and OS X. FreeON performs Hartree-Fock, pure Density Functional, and hybrid HF/DFT calculations (e.g. B3LYP) in a Cartesian-Gaussian LCAO basis. All algorithms are *O(N)* or ''O(N *log* N)'' for non-metallic systems. Periodic boundary conditions in 1, 2 and 3 dimensions have been implemented through the Lorentz field (Γ-point), and an internal coordinate geometry optimizer allows full (atom+cell) relaxation using analytic derivatives. Effective core potentials for energies and forces have been implemented, but Effective Core Potential (ECP) lattice forces do not work yet. Advanced features include *O(N)* static and dynamic response, as well as time reversible Born Oppenheimer Molecular Dynamics (MD).

Contents
--------

-   [History](History "wikilink")
-   [Vision](Vision "wikilink")
-   [Users Guide](Users Guide "wikilink")
-   [Developers Guide](Developers Guide "wikilink")
-   [Collaboratory](Collaboratory "wikilink")
-   [Generalized Solvers](Generalized Solvers "wikilink")
-   [Publications](Publications "wikilink")
-   [Authors](Authors "wikilink")

Download
--------

We periodically produce new releases here: [Savannah FreeON Download Area](http://savannah.nongnu.org/files/?group=freeon)

Here are a few recent versions:

-   [freeon-beta-2012-05-01.tar.bz2](http://download.savannah.gnu.org/releases/freeon/freeon-beta-2012-05-01.tar.bz2)

Contact Us
----------

While you could e-mail one of the authors directly, we prefer if you used our mailing lists or our bug tracker instead. You will probably get a response quicker this way because your message will go to all of us as opposed to just one of us. We currently operate 2 mailing lists:

-   [freeon-users](http://lists.nongnu.org/mailman/listinfo/freeon-users), our public mailing list for users of FreeON.
-   [freeon-devel](http://lists.nongnu.org/mailman/listinfo/freeon-devel), our mailing list for developers of FreeON, and people who are interested in the development process of it.

Alternatively you can file a bug report or feature request with our bug tracker system:

-   [FreeON Bug Tracker](https://savannah.nongnu.org/bugs/?group=freeon)

When you connect to our bug tracker for the first time you might get an error related to the SSL certificate used by savannah.nongnu.org. Please follow the [instructions from savannah](http://savannah.nongnu.org/tls/tutorial/) to fix this error.

Citing FreeON
-------------

    @misc{FreeON,
      author = {Bock, Nicolas and Challacombe, Matt and Gan, Chee~Kwan
        and Henkelman, Graeme and Nemeth, Karoly and Niklasson, A.~M.~N.
        and Odell, Anders and Schwegler, Eric and Tymczak, C.~J.
        and Weber, Valery},
      title = {{\\sc FreeON}},
      year = 2011,
      note = {\\mbox{L}os Alamos National Laboratory (LA-CC 01-2; LA-CC-04-086),
        Copyright University of California.},
      url = {http://freeon.org/}
    }

Wikipedia
---------

We are also on wikipedia: [FreeON](http://en.wikipedia.org/wiki/FreeON)
