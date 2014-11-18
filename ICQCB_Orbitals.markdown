---
layout: default
title: ICQCB Orbitals
---

![Hydrogen's 1s orbital ](images/H_1s.png "fig:Hydrogen's 1s orbital ")
![Hydrogen's 2s orbital ](images/H_2s.png "fig:Hydrogen's 2s orbital ")
![Hydrogen's 2px orbital ](images/H_2px.png "fig:Hydrogen's 2px orbital ")
![Hydrogen's 2px 2py and 2pz orbitals superposed ](images/H_2pxyz.png "fig:Hydrogen's 2px 2py and 2pz orbitals superposed ")
![Hydrogen's 3s orbital ](images/H_3s.png "fig:Hydrogen's 3s orbital ")
![Hydrogen's 3dxy orbital ](images/H_3dxy.png "fig:Hydrogen's 3dxy orbital ")

Atomic orbitals
---------------

We can solve exactly Schrödinger equation for the Hydrogen atom. When we do,
we observe that at least in its basal state (1s$${}^{1}$$ and in the first
excited state (2s$${}^{1}$$), the electronic density cloud has a spherical
distribution, which led to the use of the term *orbital* to describe these
states.

When we talk of an electronic cloud we do not mean that the electron is no
longer a particle/wave to become a diffuse cloud. What we refer to is the
probability of finding the electron in any position in space at a given time.
Let us remember that by Heisenberg's principle, we cannot know both its
location and momentum at the same time. Since, in practice, at a given time,
the electron can be in any point in space, we can consider its charge (over
time) as *diffused* over all the space. From this point of view, we consider
that the regions where it is more likely to find the electron have a bigger
*density of charge* than those where is it less likely found.

We are now in a better position to interpret $$\Psi$$, the wavefunction. But
first, let us take a different look.

### Statistics, amplitudes and probabilities

The interest of QM lies in its ability to predict the probability of different
events when considering large ensembles. Let us assume we take many measures
and always get the same outcome: that means the system is in a pure state
which will be associated to a set of describing quantum numbers.

Let us consider now some attribute of our system which can take various
values. We can take repeated measures and find how often we see each value in
these measures. For each possible value we will get a probability of getting
it. This probability may change if the system is in a different state. For
instance, you expect to get half heads and half tails on tossing a fair coin,
but only heads (or tails) if the coin is fake and have both sides the same.
The probability, of course will take a value between 0 and 1.

Simplifying a lot: We can denote the probability of each possible outcome by a
vector, and the values of each quantum number defining a state by another
vector. Since the outcomes will depend on the state, and we are talking
probabilities, we can define the *probability amplitude* (or simply,
*amplitude*) of finding an outcome as the vector product of both.

When we talk about an attribute that is a continuous function, we need to
define a *probability density*. And when we multiply this probability density
by the vector describing a state, what we get is the *probability density
amplitude*. This is what we call the wavefunction.

It is worth noting that this is an *amplitude*. To get the *probability*, we
need to square its modulus. Hence the common notation $$\Psi^{*} \cdot \Psi$$
to calculate the probability:

-  $$ \Psi(x, y, z, t) $$ = probability amplitude
-  $$ \Psi^{*} \cdot \Psi $$ = probability
-  $$ \int \Psi^{*} \Psi dt $$ = 1

The first line states the wave function in terms of space and time ($$t$$).

The last line simply summarizes that since we are calculating a probability,
the total sum of probabilities over all space (the integral) must be 1.

What, you say, are you telling me that QM is but a simple statistical problem?
Well, yes. You can think of QM as basic statistics and it will lead you to
deep understanding of chemical problems. But that is a longer story.

### Polielectronic atoms

Once we know the probability of finding an electron, it is straightforward to
imagine a simple way to compute the wave function of an element more complex
than Hydrogen. Deriving the full theoretical treatment is intractable, but we
can use our knowledge of basic statistics to derive a simple
**approximation**.

We can define the electronic charge distribution of polielectronic atoms
treating each electron separately. We know they fit in orbitals and that each
orbital can take up at most two electrons. So we can assign each electron its
own orbital with its associated wave function, $$X_{i}$$. Under this
approximation, we can describe the total wave function of the system $$\Psi$$,
combining the probabilities of finding each of the electrons in space as a
product of the wave functions of each of the electrons $$X_{i}$$.

$$ \Psi = Χ_{1} Χ_{2} \cdots Χ_{n} $$

### Spin and antisymmetry

Sadly for us, we know that things are not so simple: electrons have spin, and
we can only have two electrons in each orbital, each of them with different
spin. Using the standard notation, we add superindexes <sup>α</sup> and
<sup>β</sup> to indicate the two possible sping values (1/2 and -1/2). This
way we can describe Helium as 1s<sup>2</sup> and undertand it as an abbreviate
notation for 1s<sup>α</sup> and 1s<sup>β</sup>, or, written as a wave
function, Ψ = Χ<sup>α</sup> Χ<sup>β</sup>. It is also common to use a
simplified notation where the orbital name as such is assumed as associated to
an α spin electron and the orbital name covered by a macron as associated to
a β spin electron.

In order to include spin, we need to make things a little more complex in our
system of equations: if we stick to the former formula Ψ = Χ<sup>α</sup>
Χ<sup>β</sup>, we have a problem: we are considering both spin-electron wave
functions as equals, but Pauli's exclusion principle forces the function to
change sign when we exchange any pair of electrons (remember, it is an
amplitude, it won't matter when we square it to get a probability, but it does
affect the calculation).

We need to rewrite our system of equations to extend it in such a way that it
can accommodate any possible exchange of electrons.

We are not going to enter into details. It is enough to know that the result
of expanding the equations can be expressed as a determinant including all
possible exchanges. This is a special determinant, defined by diagonals, that
produces a multiplication of permutations, and which we need to normalize by
multiplying by 1/√n!. SInce the diagonals suffice to define the determinant,
and the normalization value is obvious, we usually write it in simplified form
as

$$ \Psi = X_{1} X_{2} \cdots X_{n} $$

and remember that we actually refer to a normalized determinant.

For instance, for Berillium (with 4 electrons and a structure 1s<sup>α</sup>,
1s<sup>β</sup>, 2s<sup>α</sup>, 2s<sup>β</sup>) we would have all the
following combinations of the four electrons in the four orbitals

    | 1sα(1) 1sα(2) 1sα(3) 1sα(4) |
    |                             |
    | 1sβ(1) 1sβ(2) 1sβ(3) 1sβ(4) |
    |                             |
    | 2sα(1) 2sα(3) 2sα(3) 2sα(4) |  1 / √(4!) = 1sα1 1sβ2 2sα3 2sβ4 = X1 X2 X3 X4 = Ψ
    |                             |
    | 2sβ(1) 2sβ(4) 2sβ(4) 2sβ(4) | 

### Defining the orbitals

Our next problem lies in defining each X<sub>i</sub>: the functions for each
spin-electron. Hartree defined several numerical functions that were rather
precise. These were then adapted into an analytical form that was easier to
compute by Slater. The last ones are commonly known as *Slater Type Orbitals
(STO)*

There are other approaches as we will see later, but for the time being, we
can remain with this concept: orbitals are nothing but distributions that we
can adjust mathematically.

### Adjusting the orbitals

When we have more than one electron, we may expect their behaviour to change:
the H electron only interacts wth its nucleus. But as soon as there is another
electron, they will also interact among themselves, repelling and damping the
nuclear attractive force. In Litium, for instance, the s<sup>1</sup> cloud
will slightly shield the nuclear attractive effect over the electron in
s<sup>2</sup>.

Atomic orbitals can be adjusted using the **variational principle**: the
better the quallity of the function we chose, the lower its energy will be.
Thanks to it, we now know rather well the orbitals in all the periodic system.

We will see the method later, but basically, it consists in starting with an
initial estimate and then adjusting progressively the calculations until we
can no longer decrease the energy.

### Orbital overlap

Finally, and before we jump into molecules, we need to understand a last
concept: atomic orbitals should be ortonormal.

$$ \int X_{i} X_{j} d\tau = 0 $$ if $$ i \neq j $$

However, since we normally use analytical functions to represent them, in
practice there is certain degree of overlap between the orbital wave
functions:

∫ ΧiΧj dτ = Sij

This is what we call overlap integral. We need to account for this property in
order to correct the calculations properly by substracting the overlap of each
orbital.

Note that we are referring to the fact that the wave functions should be
linearly independent (their product should be zero), we do not mean that
orbitals cannot overlap in space.
