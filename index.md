---
layout: default
title: FreeON
---

FreeON
------

FreeON is an experimental, open source
([GPL](http://www.gnu.org/licenses/gpl.html)) suite of programs for linear
scaling quantum chemistry. It is highly modular, and has been written from
scratch for $$N$$-scaling SCF theory in Fortran95 and C. Platform independent I/O
is supported with HDF5. FreeON should compile with most modern Linux
distributions and OS X. FreeON performs Hartree-Fock, pure Density Functional,
and hybrid HF/DFT calculations (e.g. B3LYP) in a Cartesian-Gaussian LCAO
basis. All algorithms are $$\mathcal{O}(N)$$ or $$\mathcal{O}(N \lg N)$$ for
non-metallic systems.  Periodic boundary conditions in 1, 2 and 3 dimensions
have been implemented through the Lorentz field ($$\Gamma$$-point), and an
internal coordinate geometry optimizer allows full (atom+cell) relaxation
using analytic derivatives.  Effective core potentials for energies and forces
have been implemented.  Advanced features include $$\mathcal{O}(N)$$ static
and dynamic response, as well as time reversible Born Oppenheimer Molecular
Dynamics (MD).

Project Statistics
------------------

{% include openhub.html %}

Spammpack
---------

Please visit the [Spammpack Project](http://freeon.github.io/spammpack).

Authors
-------

The principal authors of FreeON are:

  - Matt Challacombe @mattchallacombe
  - Nicolas Bock @nicolasbock
