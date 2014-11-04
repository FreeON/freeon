---
layout: default
title: ICQCB Molecular Orbitals
---

Molecular Orbitals
------------------

We can treat the molecular wave function much like we did with the atomic wave function. The main difference is that now we have several nuclei around which the electrons "orbit". The "orbitals" will be more complex and will span several nuclei, possibly the whole molecule, but just as we did before, we can obtain an approximation of the wave function as a linear combination of each of the molecular orbitals treated separately.

Ψ = Φ<sub>1</sub>Φ<sub>2</sub>...Φ<sub>n</sub>

Where we use Greek letter phi to indicate that this time the wave functions describe molecular orbitals.

As we did before for the Hydrogen atom, we could define each electron wave function as a combination of potential and kinetic energies, electrostatic and sping interactions of each orbital with respect to all other elements in the system (nuclei and electrons), and use that.

But for now we are interested in defining Ψ and defining each orbital in term of so many variables soon becomes impractical. It already was for atoms, wasn't it?

If we knew the molecular orbitals, it would be easy to compute Ψ, but this is not normally the case, and we need some method to compute the molecular orbitals. To do this, we assume that each molecular orbital can be expanded as a linear combination of its constituent atomic orbitals:

Φ<sub>i</sub> = ∑<sub>k</sub> c<sub>ik</sub>Χ<sub>k</sub>

In other words, we add all the constituent atomic orbitals, multiplying each of them by a specific constant. This method is known as **LCAO (Linear Combination of Atomic Orbitals)**.

Since we know well the atomic orbitals, our problem now resides in finding the appropriate values for the constants c<sub>ik</sub>. We use the *variational principle* to achieve this, which, states that as the wave function gets better (closer to the real solution), the energy will be lower.

Since we can consider that the interaction energy of two particles separated by an infinite distance is nil (0), as they get closer, the energy of the system will become more negative. Therefore, in our calculations, **the energy has negative values, and the bigger its absolute value, the lower the energy**.

In summary, we have a molecular wave function Ψ = Φ<sub>1</sub>Φ<sub>2</sub>...Φ<sub>n</sub> computed in terms of molecular orbitals. Each molecular orbital wave function in turn is defined as Φ<sub>i</sub> = ∑<sub>k</sub> c<sub>ik</sub>Χ<sub>k</sub> with the respective X<sub>i</sub> known. Akl we need is to solve the equations for the constants c<sub>ik</sub> of each atomic orbital for each molecular orbital.

Computing the constants c<sub>ik</sub> for each molecular orbital reduces to a relatively simple mathematical problem: a system of secular equations. As usual, we can write this system of equations as a matrix, and we now that the system will have a solution only if the matrix determinant is zero. The only thing we need to do then, is multiply the determinant terms to get a polynomial whose solution will yield the energies of the molecular orbitals. We can use these solutions then to find the constants (but we will not get into more detail for now).

Notice that in this process, we haven't considered the overlaps we talked about earlier and produced because we use approximate analytical functions to describe the atomic orbitals. If we want a realistic model, we need to include them (they are not difficult to compute) in our calculations.

To summarize up to this point, we have a problem where we want to:

-   Compute the molecular orbitals for all electrons
-   Include interactions among all particles
-   Include all possible electron exchanges (spin)
-   Consider overlaps between approximate atomic orbitals

We start from analytical approximations to the atomic orbitals and in principle we only need to find the constants c<sub>ik</sub>:

Ψ = Φ<sub>1</sub>Φ<sub>2</sub>...Φ<sub>n</sub>

where Ψ is a determinant where each

Φ<sub>i</sub> = ∑<sub>k</sub> c<sub>ik</sub>Χ<sub>k</sub>

and each

X<sub>k</sub> are substituted by predefined approximate functions.

We are missing only the c<sub>ik</sub> and then, working our way backwards, we will get a determinant where we simply substitute these constants to calculate it. This determinant can be solved using iterative methods and adjusting progressively the c<sub>ik</sub>
