---
layout: default
title: ICQCB SCF Theory
---

Self Consistent Field Theory
----------------------------

Our next step to complete our calculations consists in finding a way to compute the missing coefficients as efficiently as possible.

Despite its esoteric name, what we are going to see is, at the bottom, a relatively simple calculation method. But let us not step in our toes.

We can use the variational principle to adjust our calculations. This is a method commonly used: the **Hartree-Fock method (HF)**. Basically, everything reduces to solving the equations for the unknowns (the c<sub>ik</sub> constants that define the relative weights of the constituent atomic orbitals of each molecular orbital). We can use them then to compute the atomic orbitals, and these in turn to solve the wave function for the whole system.

We want to find the best set of constants for our problem. The better we estimate them, the better the solution we get and (by the variational principle), the lower the energy ε of the model system. The estimate is biased because the atomic orbitals we use to start with are themselves approximations as well.

So, what we do is simple: we start with an "adequate" approximation for these constants and obtain a solution for the system that will give the values of energy ε. We now use these values in the secular equations to get a new value for the constants, which we hope will be closer to the real ones...

And so on: we take the new estimate and repeat the procedure to get new energies, etc.. until the system is stable, that is, until the c<sub>ik</sub> obtained in one iteration are equal to the ones obtained in the previous cycle (to some error limits). We then say that the solution is *self-consistent*.

In any case, once the function has converged, there is little more we can do, can we?

The results we get tend to agree with experiment to some extent, but often have a margin of error that can not be ignored. The reasons are basically two:

-   on one hand, the most obvious, we are starting from atomic orbitals expressed as analytical approximations of the real wave function
-   on the other hand, there are other factors that we are not taking into account in the calcuations.

Both factors incur in additional errors, and both can be enhanced. To the extent we use better approximations to the atomic orbitals, we will also get better results. And conversely, to the extent we include additional factors in our calculations, the result will also get better. And finally, we could try different approaches to obtain the results.

There are several factors that we have not considered and that normally do not bother us, but sometimes they are significant. Among them are the influence of excited states, the influence of external factors, nuclear movement and relativity theory.

We normally ignore nuclear movement: electrons move at much higher speed than nuclei and we can expect them to adapt to any nuclear movement that much faster as to allow us to neglect its impact. Neglecting nuclear movement allows us to dissociate it from orbital structure, simplifying calculations. This is known as **Born-Oppenheimer approximation**.

Relativistic effects are not so easy to deal with and we will consider them in a separate secion.
