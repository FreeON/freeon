---
layout: default
title: ICQCB Basis Sets
---

Basis sets
----------

We compute molecular orbitals as a weighted sum of atomic orbitals:

Φ<sub>i</sub> = ∑<sub>k</sub> c<sub>ik</sub>Χ<sub>k</sub>

To solve this equation we take advantage of our extensive knowledge of atomic orbitals X<sub>k</sub>. To speed up calculations, we use functions that approximate the real value of each orbital wave function as the basis for our calculations. There are the **basis functions**.

Depending on how well these basis functions approximate the real wave functions, the quality of our calculations will conversely be better or worst. It is evident that, in principle, we want to use basis functions that are the better adjusted functions possible.

However, we have a problem: the calculations we wabt to carry out depend on the fourth power (or more) of the number of basis functions and we use an iterative process that will ve repeated until self consistency is accomplished, a process that might take an ideterminate amount of time, may be infinite. This means that we want to siplify the basis functions X<sub>k</sub> as much as possible, sacrificing some precision if needed, to speed up computation.

Indeed, it is like a game: as the basis functions get simpler, the result will be less precise but we will have it sooner. As they become more precise, the result will also become more accurate, but the calculation will proceed much slowlier (since the X<sub>k</sub> functions must be evaluated many times).

We can pre-compute these functions, but then we will need to increase our storage needs and, even so, we will still need to compute all the different linear combinations of them. We can also try to pre-compute all the more frequent combinations for our molecule, at an even higher cost in storage needs (we are talking combinatorial quantities here).

Thus, we need to find an equilibrium between time spent doing calculations, storage and precision, depending on the basis functions we chose.

There are many ways to simplify the functions X<sub>k</sub> to generate the basis functions can use. Each model needs to define all the basis functions for all the orbitals of all the atoms. This set of basis functions, computed according to some specific model is known as **basis set**.

### Slater Type Orbitals (STO)

Slater described a method to approximate analytically the numerical function that defined each atomic orbital. This method used orbitals that were similar to those of Hydrogen adding a series of rules to account for the type of orbital and the shielding produiced by the presence of other electrons.

To honor J. C. Slater, these functions are called Slater Type Orbitals (STO).

Slater used harmonic functions as the basis to approximate the atomic orbitals. This kind of functions raise issues when we deal with more than two atomic centres, and therefore are usually reserved for atoms and diatomic molecules or for semiempirical methods (which will be treated elsewhere).

### Gaussian Orbitals

An alternative approach consists in using Gauss functions. During the calculation we need to compute many integrals and combinations of them, and under these circumstances, the Gauss curve has tremendous advantages for speeding the calculation.

What happens is that while a Gauss curve can reasonably approach an s orbital, other orbitals do not lend themselves as well to be fitted by this kind of curve. What we do then is use a combination of gaussian curves (similar to what we do in a Fourier transform, when we substitute an arbitrary function by a sum of waves).

The speed up advantages obtained from defining orbitals as a sum of Gauss functions are more overwhelming when there are many d or higher orbitals (f, g, h...), speeding up the calculation enormously.

However, the Gauss curve depends on r<sub>2</sub> and this makes the orbitals defined as sums of gaussian curves (Gaussian Type Orbitals or GTO) inferior to STO in two aspects: en the nucleus, the Gauss curve has a slope of zero, but an STO incides sharply, so GTOs give worst results near the nucleus, and the other is that far from the nucleus the fall off very quickly. What this implies is that we need to include more Gaussian functions in the sum to achieve a better adjustment (as a rule, three times more than using STO). Even so, calculations proceed faster, making GTOs the preferred choice.

Once we have decided whether we will use STO (spherical harmonics) or GTP (sums of Gauss curves), the next step is deciding the number of functions to use and, if possible define a minimal basis set.

It seems evident that, as a minimum, we want to use one s function for H and He, two s functions for the first row elements, and so on. Each of these functions describing an orbita will be defined as a sum of Gauss curves (GTO) or a Sater type function (STO).

### n-Zeta sets

The next obvious enhancement is to double the number of all the basis functions. The term "zeta" stems from the exponen used in STOs, which is usually referred to as ζ. Thus, doubling the number of basis functions implies doubling the exponents, and we refer to this scheme as "double-zeta": in a "double-zeta" approach, we would use two s functions for H and He, and so on.

Using a **Double-Zeta (DZ)** basis set will increase precision, but duplicates the number of functions and increases the computation time. We can enhance it assuming that bonds affect specially the valence electrons, which require higher precision: if we use one function for internal electrons and two for valence electrons, we will obtain a set of basis functions that have been split by the valence electron, a **split valence basis**.

The next obvious step is to triple (**Triple-Zeta, TZ**), cuadruplicate (**Quadruple-Zeta, QZ**) or quintuplicate (**5Z**) the number of functions. Evidently, as we increase precision, computation time aslo grows, as a function usually of, at least, the fourth power of the number of basis functions, making these sets prohibitively expensive for large systems.

We can add in some functions to account for polarization effects to each of these basis sets, obtaining in this way **DZP** (Double Zeta plus Polarization), **TZ2P** (Triple Zeta with doubled polarization functions), etc...

### Contracted sets

We have already mentioned that inner electrons populate orbitals that require less precision, and contribute less to the chemical bond. During our variational calculation process they will change only slightly, becoming good candidates to speed up calculation by substituting them by a *fixed* combination of functions. This is to say, we compute them once and then consider them constant taking them out of our variational calculations (i. e. out of the iterative recalculation process).

When we use the **primitive** orbitals (PGTOs) to obtain a linear combination that will be kept fixed during all the calculations, we talk of *contraction*, and the resulting orbitals are called **contracted** (CGTOs). The degree of contraction is the number of primitive orbitals (PGTOs) used to compute the contracted orbitas (CGTOs).
