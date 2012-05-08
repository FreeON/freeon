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
-   [Users Guide](Users Guide "wikilink")
-   [Developers Guide](Developers Guide "wikilink")
-   [Collaboratory](Collaboratory "wikilink")
-   [Generalized Solvers](Generalized Solvers "wikilink")
-   [Publications](Publications "wikilink")
-   [Authors](Authors "wikilink")

Download
--------

We periodically release snapshots of the development tree in the [Savannah FreeON Download Area](http://savannah.nongnu.org/files/?group=freeon)

Here are a few recent ones:

-   [freeon-beta-2012-05-01.tar.bz2](http://download.savannah.gnu.org/releases/freeon/freeon-beta-2012-05-01.tar.bz2)
-   [freeon-beta-2012-05-03.tar.bz2](http://download.savannah.gnu.org/releases/freeon/freeon-beta-2012-05-03.tar.bz2)

SSL certificate
---------------

We force the login and the account creation page to go through SSL for security. Since we currently use a self-signed SSL certificate, your browser will likely complain that our certificate can not be trusted, but trust us, it can (we made it). You can verify that you are using the right certificate by comparing the certificat's fingerprint. It should be:

-   SHA-256 fingerprint: 44 4B 37 6E BA 1C BF A6 0E 1B 01 57 1C 82 9B 81 AB 7D 8B 5F EB CF D7 03 BD 46 6A 08 AC 70 8A 7F
-   SHA-1 fingerprint: 76 33 F9 C6 D0 67 3C 59 88 09 DF F6 7E A0 AA 04 5D 7B 4C F9

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

TeX Format

    @misc{FreeON,
      author = {Bock, Nicolas and Challacombe, Matt and Gan, Chee~Kwan
        and Henkelman, Graeme and Nemeth, Karoly and Niklasson, A.~M.~N.
        and Odell, Anders and Schwegler, Eric and Tymczak, C.~J.
        and Weber, Valery},
      title = {{\\sc FreeON}},
      year = 2012,
      note = {\\mbox{L}os Alamos National Laboratory (LA-CC 01-2; LA-CC-04-086),
        Copyright University of California.},
      url = {http://www.freeon.org/}
    }

EndNote/RIS

    TY  - GEN
    N1  - Los Alamos National Laboratory (LA-CC 01-2; LA-CC-04-086), Copyright University of California.
    AU  - Bock, Nicolas
    AU  - Challacombe, Matt
    AU  - Gan, Chee Kwan
    AU  - Henkelman, Graeme
    AU  - Nemeth, Karoly
    AU  - Niklasson, A. M. N.
    AU  - Odell, Anders
    AU  - Schwegler, Eric
    AU  - Tymczak, C. J.
    AU  - Weber, Valery
    T1  - FreeON
    PY  - 2012

Wikipedia
---------

We are also on wikipedia: [FreeON](http://en.wikipedia.org/wiki/FreeON)
