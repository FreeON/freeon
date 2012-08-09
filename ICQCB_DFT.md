---
layout: default
title: ICQCB DFT
---

Density Functional Theory
-------------------------

This is a model that has gained popularity in the last years because it is less demanding in computing capacity.

Hohenberg and Kohn established the basis for this method proving that molecular energy could be determined directly from the total electronic density of the molecule.

Next, Kohn and Sham defined the basic method: tota density can be computed as a linear combination of functions that are similar to the orbital wave functions we use in HF (*Kohn-Sham orbitals*). Once the density has been computed, the energy can be readily derived from it using a **functional**: a functional is a special kind of function that acts on other functions (not directly on the values of variables). For instance, we know that velocity is a function of two variables, time and space, velocity is a normal function, but we can imagine a function to compute derivatives: it does not make sense to compute the derivative of a variable's value (integral of "5"?), but it does to compute the derivative of a function: this would be a "functional".

The bottom line is that Khon and Sham. starting from the fact that energy depends on electronic density, developed all their calculations in terms of a final element, a "functional" that we need to determine, and that is a function of density (which is itself a function of many variables).

The problem lies in determining this functional, to determine how energy depends on density, and substitute it in Kohn-Sham calculations to obtain the values for energy. And it is a problem because we do not know it. We do not know how energy depends on density, the functional that is the only missing element in the equations.

If we knew this functional, then DFT would have a definitive advantage: it is simpler than HF as it depends only on orbital density, and therefore it is a problem that grows as N<sup>3</sup> on system size. If we knew this missing functional acting on density we would be able to perform our calculations a lot faster.

Until we can find this functional, we need to work with approximations to this unknown function. And, indeed, a number of approximations have been proposed, either based on previous *ab initio* calculations, or fitting experimental data. These approximations already include electronic correlation to some extent and, so, we do not need to compute it either:

-   **LDA (Local Density Approximation)** and **LSDA (Local Spin Density Approximation**) are two popular approaches that are based only on electronic density
-   **Gradient corrected methods** are based in the density gradient instead of in the density itself
-   **Hybrid methods** combine functionals of the two other methods with data from HF simulations

In the last years there have arisen a number of techniques that **scale linearly** with problem size: they basically neglect the influence of far electronic interactions, and so they grow linearly for large molecules. The most popular methods are likely *FMM (Fast Multipole Method)* and *CFMM (Continuous FMM)*.

The availability of linearly scaling methods for DFT has made them very attractive for large systems. However, we must remember that many traditional methods also impose distance limits when evaluating integrals to reduce computation, thus requiring less time too. As a rule, we can expect DFT to require 60-80% of the time spent by a sophisticated implementation using traditional methods.

In practice, our decision will depend on the basis set and the selected functional as the results may change considerably: DFT also requires basis sets (we will normally want at least 6-31G\* for precise calculations). The most commonly used method is currently the hybrid functional B3LYP.

The main problem with DFT is that we ignore the actual functional and need to resort to approximations that we cannot define rigorously, thus casting a shadow of doubt over the results.
