---
layout: default
title: ICQCB Ab initio methods
---

*Ab initio* methods
-------------------

### Hartree-Fock

The basic elements of this calculation are the use of the variational principle to adjust iteratively the value of the coefficients required by the orbital combinations. The method consists in using the self-consistent field theory to refine the values of the constants we need to compute the molecular orbitals using LCAO. The basic method is therefore referred to as SCF. As a rule, HF scales as N<sup>4</sup> where N is the number of basis functions (a calculation on a system twice as large requires 16 times longer).

When we apply the Hartree-Fock method, usually we use the same orbitals to define electrons with different spin, a method known as **Restricted Hartree-Fock (RHF)**. This implies that we use the same shape for both electrons (with different spin) in an orbital. As can be expected, this method works well as long as we have all orbitals fully occupied.

If we release this restriction (that is, we allow that each electron having a different spin uses an orbital with a different shape), we obtain the **Unrestricted Hartree-Fock (UHF)** method. The proble now is that paired electrons need not have the same spatial distribution, inducing an error known as *spin contamination*.

The alternative is to use the **Restricted Open Shell Hartree-Fock (ROHF)** approach: in this case we uncouple paired electrons considering their spin, but force them to share the same spatial orbital.

RHF is a good approximation for *closed shell* situations, i. e., when every orbital is fully occupied, and is faster than UHF, so we prefer it in these cases. However, when a system has orbitals which are not fully occupied (they have only one electron), the behaviour is different due to a polarization effect. In these *open shell* cases we can also use ROHF, but UHF is normally preferred.

On the other hand, UHF is less adequate for *singlet* situations (full orbital occupancy) since it may include contributions from states with a higher energy.

HF can account for up to 99% of the energy, but there is 1% that remains unaccounted for. The energy deficit in HF is known as *Electronic Correlation Energy*. This energy difference is due to the fact that electrons do not move independently.

We can visualize it imagining that we calculate the probability of each electron being in a given voxel: in principle, there is a probability that both are in the same voxel, but we know they have a strong repulstion, and so we do not expect both to be at the same place simultaneously. In other words, there are interactions among them that correlate their movement.

### Correlation

A way to increase accuracy is to account for these correlations in electron movement. Some of the methods that take it into account are those of **Moeller-Plesset (MPn)** (where N is a number indicating the level or precision), **General Valence Bond (GVB)**, **Multiconfigurational Self-Consistent Field (MCSCF)**, **Configuration Interaction (CI)** and **Coupled Cluster (CC)**.

Including correlation gives higher precision in our calculations. In general, it is not very important to consider it in biological problems, where we are more interested in cualitative aspects than in a high precision result. The exception are transition elements, where this effect becomes very important.

### Møller-Plesset Perturbation Theory

We can include correlation using a mathematical device known as perturbation theory. Basically, it consists in assuming that the solution to our problem does not differ too much from that or another problem that is slightly different but for which we know the solution, and then applying several correction levesl to that one in order to obtain our solution. Applied in Quantum Chemistry, this results in Many Body Perturbation Theory (MBPT).

In our case, the perturbation we use is doubly counting electronic repulstions (which certainly is not a minor correction, but works). Applying several correction levels we obtain various methods (MP0, MP1, MP2...) with increasing costs. MP2, for instance grows with the fifth power (N<sup>5</sup>) and accounts for 80-90% of the correlation energy. MP3 grows with N<sup>6</sup>, MP4 as N<sup>7</sup>, MP5 as N<sup>8</sup> to N<sup>10</sup> or worst...

A final consideration is that MP is not variational, and can yield energies below the real one.

### Configuration Interaction (CI)

=

Another approach to increase accuracy is to allow the interaction of various configurations. Formally this consists of adding a new linear combination to our process:

Ψ<sub>enhanced</sub> = aΨ<sub>0</sub> + bΨ<sub>1</sub> + ...

Where Ψ<sub>0</sub> is the initial wave function we find, and Ψ<sub>1</sub>, Ψ<sub>2</sub>, etc. are configurations for excited states with the same symmetry, and a, b, c... are coefficients chosen to increase precision.

What this implies is that the configurations with excited electron states affect the structure of all orbitals.

This method is known as Configuration Interaction. The problem here is knowing when to stop: as we increase the considered configurations, the results get better, but higher excited states also have a lower influence in the basal state. In principle, there is an infinite number of excited states, so although this method could potentially give exact results (if we consider all infinite possibilities), in practice we must restrict the number of configurations considered.

CI calculations grow very quickly as can be expected, being usually of the order of N<sup>8</sup> or worst. We normally refer to them by the number of excited configurations considered: CIS (single excitation), CISD (single and double), CISDT (..and triple), CISDTQ (..and quadruple), and so on.

### Multi-configuration self-consistent field (MCSCF)

We can consider this a variant of CI where we not only optimize the coefficients of the wave functions (determinants), but also the molecular orbitals used to build them. These methods can produce results that are more precise, as we can expect, but also at a singificantly higher cost.

It may be worth using this method if the HF calculations are not precise: to this end, we calculate both using HF and CIS, and if the result of HF is 0.9 or less than CI (i. e. if the contribution of electron correlation is bigger than 10%), then we can consider using MCSCF (or else, CISDTQ or greater).

The main problem is that this is no longer an automated method: we must decide which excited orbitals we need to include in the calculation. A popular approach is CASSCF (Complete Active Space SCF) where orbitals are divided into active and passive, treating completely only the active orbitals (which we pick by hand).

### Coupled Cluster

The mechanism is similar to CI: the final function is a combination of excited functions, but the selection is more complex. The basic idea is not to include all types of correction up to a given order, but one kind of corrections up to infinity. As in CI, we are forced to truncate the series and end up with different methods depending on the considered excitations: CCS, CCSD, CCSDT....

Precision is similar to CI, but it has the advantage of being a method *extensive* in size.

### Localized methods

As we have seen, generally the cost grows with N<sup>5</sup> to N<sup>8</sup> with respect to the size of the system. As most of the physical interactions are fundamentally local, this scale is scarcely reallistic. The reason is that canonical orbitals are delocalized on the molecule (all contribute). We can think therefore that using localized orbitals would be a nice starting point to carry out our calculations, considering that only a few orbitals participate in the main contribution and then ignoring all the rest. There are now MP2 and CC methods relying on localized molecular orbitals, but they are still far from being widely used. These methods are more difficult to develop, but has the potential to scale linearly with problem size and, therefore, deal with large scale systems.

We will learn more about linear scaling approaches elsewhere.

As a general rule, as we add more terms, accuracy increases, but on the other hand, computation times also grow exponentially. We may expect usually cost/precision relationships to be as follows:

HF \<\< MP2 \< CISD \< MP4(SDQ) \~ CCSD \< MP4 \< CCSD(T) \< CCSDT

### Other methods

There are many more *ab initio* methods (methods that calculate all the components of our system of equations). We are not going to describe all of them now. For the time being, we have more than enough to get started and, on the other hand, as we increase the theory level (precision), computation times also increase considerably. Given the current status of technology, we will normally not be concerned with more complex approaches.

In summary, we will normally use HF or methods incuding correlation interactions. With the later we can introduce additional terms to try to cover correlation energy. All in all, we should not forget that, although we can include many additional terms, we will still be starting from approximate basis functions to describe atomic orbitals, and that we are commonly ignoring relativistic effects.

Even so, we can cheat the game.
