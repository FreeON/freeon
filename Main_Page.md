---
layout: default
title: Main Page
---

Introduction
------------

FreeON is an experimental, open source (GPL) suite of programs for linear scaling quantum chemistry, formerly known as MondoSCF. It is highly modular, and has been written from scratch for N-scaling SCF theory in Fortran95 and C. Platform independent I/O is supported with HDF5. FreeON should compile with most modern Linux distributions. FreeON performs Hartree-Fock, pure Density Functional, and hybrid HF/DFT calculations (e.g. B3LYP) in a Cartesian-Gaussian LCAO basis. All algorithms are *O(N)* or ''O(N *log* N)'' for non-metallic systems. Periodic boundary conditions in 1, 2 and 3 dimensions have been implemented through the Lorentz field (Γ-point), and an internal coordinate geometry optimizer allows full (atom+cell) relaxation using analytic derivatives. Effective core potentials for energies and forces have been implemented, but Effective Core Potential (ECP) lattice forces do not work yet. Advanced features include *O(N)* static and dynamic response, as well as time reversible Born Oppenheimer Molecular Dynamics (MD).

Content of this wiki
--------------------

-   [API](API "wikilink") Documentation
-   [Input File Syntax](Input File Syntax "wikilink")
-   [Example Inputs](Example Inputs "wikilink")
-   [Vision Statement](Vision Statement "wikilink")

Features
--------

[coming soon]

Contact Us
----------

While you could e-mail one of the authors directly, we prefer if you used our mailing lists or our bug tracker instead. You will probably get a response quicker this way because your message will go to all of us as opposed to just one of us. We currently operate 2 mailing lists:

-   [freeon-users](http://lists.nongnu.org/mailman/listinfo/freeon-users), our public mailing list for users of FreeON.
-   [freeon-devel](http://lists.nongnu.org/mailman/listinfo/freeon-devel), our mailing list for developers of FreeON, and people who are interested in the development process of it.

Alternatively you can file a bug report or feature request with our bug tracker system:

-   [FreeON Bug Tracker](https://savannah.nongnu.org/bugs/?group=freeon)

When you connect to our bug tracker for the first time you might get an error related to the SSL certificate used by savannah.nongnu.org. Please follow the [instructions from savannah](http://savannah.nongnu.org/tls/tutorial/) to fix this error.

Developer Page
--------------

The FreeON source code is hosted by the [Free Software Foundation](http://www.fsf.org/) at

-   [<http://savannah.nongnu.org/projects/freeon>](http://savannah.nongnu.org/projects/freeon).

Related Electronic Structure Programs
-------------------------------------

Our group at Los Alamos National Laboratory is actively developing 3 independent electronic structure projects. You can find out more about the other 2 from their respective webpages:

-   [RSPt](http://www.rspt.net/) is an Open Source project and a code for band structure calculations. RSPt (Relativistic Spin Polarised (test)) is very robust and flexible and can be used to calculate the band structure and total energy for all elements, and combinations thereof, over a wide range of volumes and structures. The Full-Potential Linear Muffin-Tin Orbital (FP-LMTO) method allows for very small basis sets and fast calculations. Compared to LMTO-ASA there are no restrictions on the symmetry of the potential in FP-LMTO. RSPt allows for multiple energy sets (i.e. valence and semi-core states) with the same angular quantum numbers but different main quantum number.
-   [LATTE](http://savannah.nongnu.org/projects/latte) (Los Alamos Transferable Tight-binding for Energetics) enables quantum-based molecular dynamics simulations on many-core architectures of materials with mixed covalent and ionic bonding.

Authors
-------

The authors of FreeON, in alphabetical order:

|---|---|---|
|Nicolas Bock|\<[nbock@lanl.gov](http://mailto:nbock@lanl.gov)\>|Los Alamos National Laboratory, T-1|
|Matt Challacombe|\<[mchalla@lanl.gov](http://mailto:mchalla@lanl.gov)\>|Los Alamos National Laboratory, T-1|
|Chee Kwan Gan||Institute for High Performance Computing|
|Graeme Henkelman||University of Texas at Austin|
|Karoly Nemeth||Argonne National Laboratory|
|Anders Niklasson|\<[amn@lanl.gov](http://mailto:amn@lanl.gov)\>|Los Alamos National Laboratory, T-1|
|Anders Odell|\<[odell@mse.kth.se](http://mailto:odell@mse.kth.se)\>|KTH Stockholm|
|Eric Schwegler||Lawrence Livermore National Laboratory|
|C. J. Tymczak|\<[tymczakcj@tsu.edu](http://mailto:tymczakcj@tsu.edu)\>|Texas Southern University|
|Valery Weber||Universität Zürich|


