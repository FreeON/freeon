---
layout: default
title: Single point calculation examples
---

Single Point calculations
-------------------------

In subdirectory **SinglePoint**, under the **Validate** directory of FreeON source code distribution (./Validate/SinglePoint/.) you can find examples of various kinds of single point energy calculations.

A single point calculation solves the SCF equations to calculate the wavefunction of the system which can be used to derive a number of useful properties, most prominently, the energy of the system. It uses the geometry provided, without any initial geometry optimization.

The actual method used to solve the Schröedinger equations will depend on the characteristics of the system.

-   **Restricted**: this subdirectory contains examples for Restricted Hartree-Fock calculations:
    -   [CH2.inp](Restricted-C2H2.inp "wikilink") - Single point energy calculaton of C<sub>2</sub>H<sub>2</sub> using DFT (6-31G\* and B3LYP)
    -   [C2H4.inp](Restricted-C2H4.inp "wikilink") - Sinlge point energy calculaton of C<sub>2</sub>H<sub>4</sub> using DFT (6-31G\* and B3LYP)
    -   [C2H6.inp](Restricted-C2H6.inp "wikilink") - Single point energy calculaton of C<sub>2</sub>H<sub>6</sub> using DFT (6-31G\* and B3LYP)
    -   [CH3NH2.inp](Restricted-CH3NH2.inp "wikilink") - Single point energy calculaton of CH<sub>3</sub>NH<sub>2</sub> using DFT (6-31G\* and B3LYP)
    -   [CH3OH.inp](Restricted-CH3OH.inp "wikilink") - Single point energy calculaton of CH<sub>3</sub>OH using DFT (6-31G\* and B3LYP)
    -   [CH4.inp](Restricted-CH4.inp "wikilink") - Single point energy calculaton of CH<sub>4</sub> using DFT (6-31G\* and B3LYP)
    -   [C.inp](Restricted-C.inp "wikilink") - Single point calculation of C with two unpaired electrons using DFT (6-31G\* and B3LYP)
    -   [H2.inp](Restricted-H2.inp "wikilink") - Single point energy calculaton of H<sub>2</sub> using DFT (6-31G\* and B3LYP)
    -   [H2O-FD\_Smearing.inp](Restricted-H2O-FD_Smearing.inp "wikilink") - Single point energy calculation of H<sub>2</sub> with Fermi-Dirac smearing at a temperature of 10000 K using the split STO-2G basis set
    -   [H2O.inp](Restricted-H2O.inp "wikilink") - Two-step single point calculation on H<sub>2</sub>O, using an initial approximation to the wavefunction using STO-2G, followed by an accurate calculation using a split 6-31G\*\* basis set
    -   [HeH+-inp](Restricted-HeH+-inp "wikilink") - Single point calculation on HeH<sup>+</sup> using the 3-21G split basis set
    -   [H.inp](Restricted-H.inp "wikilink") - Single point calculation on H atom using DFT with the 6-31G\* basis set and B3LYP functional
    -   [N2.inp](Restricted-N2.inp "wikilink") - Single point calculation on N<sub>2</sub> using DFT with the 6-31G\* basis set and B3LYP functional
    -   [N2O.inp](Restricted-N2O.inp "wikilink") - Single point calculation on N<sub>2</sub>O using DFT with the 6-31G\* basis set and B3LYP functional
    -   [O2.inp](Restricted-O2.inp "wikilink") - Single point calculation on O<sub>2</sub> using DFT with the 6-31G\* basis set and B3LYP functional
    -   [Vancomycin.inp](Restricted-Vancomycin.inp "wikilink") - Single point calculation of Vancomycin, using TC2 with an STO-3G split basis set and the B3LYP/VWN5 functional
    -   [Vancomycin-quick.inp](Restricted-Vancomycin-quick.inp "wikilink") - Single point calculation of Vancomycin, using the STO-3G basis set

-   **Unrestricted**: this subdirectory contains examples for Unrestricted Hartree-Fock calculations:
    -   [Be.inp](Unrestricted-Be.inp "wikilink")
    -   [Furan.inp](Unrestricted-Furan.inp "wikilink")
    -   [H2.inp](Unrestricted-H2.inp "wikilink")


