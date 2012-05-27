---
layout: default
title: ICQCB Wavefunction
---

The wavefunction
----------------

All Quantum Mechanics methodology relies on Schrödinger's wave function:

<div class="center" style="width:auto; margin-left:auto; margin-right:auto;">
HΨ = EΨ

</div>
In this equation, H is an *operator* that we apply on Ψ, a *function*. and E is the energy of the system. It is worth noting that Ψ is nothing but a norma mathematic function, just like any other, and the fact we denote it with a Greek letter has no bearing or special meaning.

The only thing that makes Ψ special is that it describes a natural system, and therefore it is a continuous function defined in three dimensions (X, Y and Z). This has an important consequence: when we want to evaluate it, we need to do it in three dimensions, multiplying the space of values to compute. To get an idea: if it were a discrete function of X with values between 0 and N, we would have to evaluate N values. For a three-dimensional (X, Y, Z) function, with values between 0 and N, we need to evaluate N · N · N = N<sup>2</sup>. But in our case, we need to evaluate Ψ as a continuous function (using integrals) and, if that's not enough, it is defined over all the space.

In practice, the above is not the full truth: the above is usually implemented as a system of equations, where we have a function Ψ<sub>n</sub> for each energy level E<sub>n</sub>. Each of these functions will denote an energy level, i. e. a different excitation state of the system. However, this only implies that to solve the equation we need to use standard methods to solve systems of equations: more explicitly, we can (and usually do) define the system as a determinant to use the classical methods we learnt in the school.

This far, perhaps, the most strinking thing is that we can solve QM problems using pre-univeristy level mathematics.

Our problem lies in defining the function Ψ. To do this, we can simply use the standard formulas for potential and kinetic energy and for electrostatic interactions... again, basic formulas of pre-university physics, and the complexity lies in combining them and extending them for a system of many particles.

What happens then, is that as we keep adding elements to the system, calculations become more complex to the point where they can become unmanageable. However, once we have defined a function to describe the status of a system, we can use it to compute other observable properties using (for each of them) the appropriate operator:

<div class="center" style="width:auto; margin-left:auto; margin-right:auto;">
Observable = ∫ Ψ\* Operator Ψ dτ / ∫ Ψ\* Ψ dτ

</div>
Where Ψ\* is the complex conjugate of the function, and thus Ψ\* Ψ corresponds to the square of the function modulus, and where dτ represents (dx dy dz), i. e. we need to integrate over all space.

For Energy, whose operator is H, the function would become

<div class="center" style="width:auto; margin-left:auto; margin-right:auto;">
E = ∫ Ψ\* H Ψ dτ / ∫ Ψ\* Ψ dτ

</div>
What we are going to do is try to find ways to limit calculations (avoiding integration over all space) and to simplify the function (using others that are simpler and easier to evaluate instead of the fully explicit development).
