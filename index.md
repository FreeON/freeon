---
layout: default
title: FreeON
---

Introduction
------------

FreeON is an experimental, open source
([GPL](http://www.gnu.org/licenses/gpl.html)) suite of programs for linear
scaling quantum chemistry. It is highly modular, and has been written from
scratch for $$N$$-scaling SCF theory in Fortran and C. Platform independent
I/O is supported with [HDF5](http://www.hdfgroup.org/HDF5/). FreeON should
compile with most modern Linux distributions and OS X. FreeON performs
Hartree-Fock, pure Density Functional, and hybrid HF/DFT calculations (e.g.
B3LYP) in a Cartesian-Gaussian LCAO basis. All algorithms are $$O(N)$$ or
$$O(N \lg N)$$ for non-metallic systems. Periodic boundary conditions in 1, 2
and 3 dimensions have been implemented through the Lorentz field
($$\Gamma$$-point), and an internal coordinate geometry optimizer allows full
(atom+cell) relaxation using analytic derivatives. Effective core potentials
for energies and forces have been implemented. Advanced features include
$$O(N)$$ static and dynamic response, as well as time reversible Born
Oppenheimer Molecular Dynamics (MD).

Contents
--------

-   [History](History.html "wikilink")
-   [Vision](Vision.html "wikilink")
-   [Status](https://travis-ci.org/FreeON/freeon)
-   [Users](Users_Guide.html "wikilink")
-   [Developers](Developers_Guide.html "wikilink")
-   [Collaboratory](Collaboratory.html "wikilink")
-   [Solvers](Solvers.html "wikilink")
-   [Publications](Publications.html "wikilink")
-   [Citation](Citing_FreeON.html "wikilink")
-   [Authors](Authors.html "wikilink")

Download
--------

We periodically release snapshots of the development tree in the [Savannah
FreeON Download Area](http://savannah.nongnu.org/files/?group=freeon)

Here are a few recent ones:

-   [freeon-1.0.8.tar.bz2](http://download.savannah.gnu.org/releases/freeon/freeon-1.0.8.tar.bz2)

Project Stats
-------------

{% include openhub.html %}

Education
---------

[An introduction to Computational Quantum Chemistry for Biologists](An_introduction_to_Computational_Quantum_Chemistry_for_Biologists.html "wikilink")

Contact Us
----------

Please use our mailing lists or our bug tracker to ask questions or to report
bugs. The mailing lists can be found here:

-   [freeon-users](http://lists.nongnu.org/mailman/listinfo/freeon-users), our
    public mailing list for users of FreeON.
-   [freeon-devel](http://lists.nongnu.org/mailman/listinfo/freeon-devel), our
    mailing list for developers of FreeON, and people who are interested in
    the development process of it.

You are encouraged to file bug reports or feature requests with our issue
tracker system:

-   [Issue Tracker](https://github.com/FreeON/freeon/issues)

Wikipedia
---------

We are also on wikipedia: [FreeON](http://en.wikipedia.org/wiki/FreeON)
