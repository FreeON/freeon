---
layout: default
title: ICQCB Popular Basis Sets
---

Popular basis sets
------------------

### STO-nG

These are STOs that ise n primitive Gauss functions (without contraction). This is considered a minimal type basis set, where all exponents are obtained by curve fitting (not by the variational principle). In general, it has been shown that using more than 3 PGTOs does not improve significantly on the results, and so, the **STO-3G** set is widely used as a **minimal basis set**.

### k-nlmG

These are basis sets derived by *Pople et al.* of the kind of split valence sets. Here k states the number of PGTOs (and precision) used to approximate the core electrons, and nlm indicatre how the division has been done in the valence orbitals. When two numbers are used (nl) it means one valence division, if there are three numbers (nlm), two divisions (to adapt precision progressively). When there is a value after the G, it indicates the number of polarization functions added. We can see it better with some examples:

-   **3-21G** is a split basis set where the core orbitals are a contraction of 3 PGTOs, the inner part of the valence electrons uses 2 PGTOs and the external 1PGTO. This one would be similar to STO-3G but more preciseusing the double of functions in the valence region (and therefore more versatile)
-   **6-31G** uses a contraction of 6 PGTOs for the core, 3PGTOs for the inner valence electrons, and 1 for the external electrons. This is similar to 3-21G but more precise
-   **6-311G** is basis set divided in three parts: core orbitals are represented by a contraction of 6 functions, and the valence region is split in three parts with 3, 1 and 1 PGTOs each

We can add to each of them diffuse functions and/or polarization functions to increase precision. The **diffuse functions** are needed to deal with anions and are usually noted as + or ++ before the G. **Polarization functions** allow modeling the influence of an external electrical field and are noted after the G with \* (if it is one) or \*\* (if two) or designing their kind and number. Thus, a basis set 6-311++G(3d,2p) is a split basis set using 6 contracted GTOs for the core electrons, and dividing the valence shell in three levels, where the inner orbitals use a fixed combination of 3 GTOs, and two 2GTOs are left free to vary during the LCAO procedure. The ++ indicates that two extra diffuse functions have been added, and the (3d,2p) indicates that 3 functions have been added for d orbital polarization, and 2 for p orbitals.

### Other basis sets

Before going any forward, we should point out that there are many other basis sets, with higher precision levels. By extension too, they also have a higher computational cost. For this reason, at least by now, they are not commonly used in the Life and Health sciences, and their use is restricted to a great extent to theoretical chemistry, where high precision analysis of relatively small compounds is desired.

In Computational Biology, we have traditionally used semiempirical methods based on the use of valence bong electrons exclusively, with parameters adjusted to simplify calculations. Currently it is possible to use directly *ab initio* methods, such as RHF, UHF, MP, CI or others, using basis sets among the ones we mentioned above, and they are, indeed, becoming more popular each day.
