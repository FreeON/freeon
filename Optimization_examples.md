---
layout: default
title: Optimization examples
---

In the Optimization subdirectory you can find a large number of examples of optimization calculations that can be carried out with FreeON.

Most of these examples include directives to use the parallel MPI implementation of FreeON if available. This means that if you compiled a version of FreeON with support for parallel computing using MPI, then the program will be able to use more resources in your computer to complete the calculation faster. If your version of FreeON does not support MPI parallel processing, these directives will be ignored and the calculation will proceed normally using only one CPU.

-   [01-Benzene.inp](01-Benzene.inp "wikilink") - reproduction of the optimization of benzenr published by G.A.Jeffrey, J.R.Ruble, R.K.McMullan, J.A.Pople (1987) Proc.R.Soc.London,Ser.A v414 p47, using DFT with STO3G and the BLYP exchange correlation functional
-   [02-BNSheet.inp](02-BNSheet.inp "wikilink") - Optimization of a BN sheet using DFT with split STO-2G and B3LYP
-   03-Caffeine.inp - Two step optimization of Caffeine, first using STO-2G, followed by DFT with 6-31G\*\* and BLYP
-   04-CH2Wire.inp - Polyacetilene geometry optimization with fixed lattice vectors in STO-3G basis
-   05-Ice.inp - Ice 1h-H2O-(P63/MMC) after Goto, A.;Hondoh, T.;Mae, S.(1990)
-   06-Sulfur.inp - Monoclinic gamma sulfur, after Gallacher, A.C.;Pinkerton, A.A. (1992) Phase Transition 38, 127-220, CSD entry 66517ICS
-   07-WaterHexamer.inp - Optimization of a water hexamer using STO-3G split basis set
-   08-Na2Cl2.inp - Optimization of Na<sub>2</sub>Cl<sub>2</sub>. Results can be comapred with reference in <http://cccbdb.nist.gov/>
-   09-TEMPO.inp - Optimization of 2,2,6,6-tetramethyl-piperidine-1-oxyl (TEMPO)
-   10-Malonaldahyde.inp - Transition state optimization from reactants and products of malonaldehyde.
-   11-CO2.inp - Solid CO\_2 C O Unit Cell, 83K, afetr J Phys Soc Japan, 53(3) 1176-1184 (1984)
-   12-BPA.inp - Calculation on Bisphenol A using STO-2G split basis, with commented example for a more accurate optimization
-   13-PhyticAcid.inp - Calculation on Phytic Acid using STO-3G split basis, with commented example for a more accurate optimization
-   14-LJ-2-SD.inp - Ar dimer with Lennard-Jones potential using Steepest Descent
-   15-LJ-2-CG.inp - Ar dimer with Lennard-Jones potential using Conjugate Gradients
-   16-LJ-3-SD.inp - Ar cluster (3 Ar) with Lennard-Jones potential using Steepest Descent
-   17-LJ-3-CG.inp - Ar cluster (3 Ar) with Lennard-Jones potential using Conjugate Gradients
-   18-LJ-4-SD.inp - Ar cluster (4 Ar) with Lennard-Jones potential using Steepest Descent
-   19-LJ-4-CG.inp - Ar cluster (4 Ar) with Lennard-Jones potential using Conjugate Gradients
-   20-LJ-5-SD.inp - Ar cluster (5 Ar) with Lennard-Jones potential using Steepest Descent
-   21-LJ-5-CG.inp - Ar cluster (5 Ar) with Lennard-Jones potential using Conjugate Gradients
-   22-LJ-100-SD.inp - 100 atom Lennard-Jones cluster using Steepest Descent. Structure taken from <http://physchem.ox.ac.uk/~doye/jon/structures/LJ/tables.150.html> and perturbed by ± 0.2 Å
-   23-LJ-100-CG.inp - 100 atom Lennard-Jones cluster using Conjugate Gradients. Structure taken from <http://physchem.ox.ac.uk/~doye/jon/structures/LJ/tables.150.html> and perturbed by ± 0.2 Å
-   24-Vancomycin-internal.inp - Vancomycin optimization using a split STO-2G basis and a maximum of 1000 optimization steps
-   25-Vancomycin-CG.inp - Vancomycin optimization using a split STO-2G basis and coarse tuned conjugate gradients
-   26-Vancomycin-CG.inp - Vancomycin optimization using a split STO-2G basis and medium tuned conjugate gradients
-   27-Vancomycin-CG.inp - Vancomycin optimization using a split STO-2G basis and fine tuned conjugate gradients
-   28-BPA-2.inp - two Bisphenol A optimization using STO-2G split with commented example of two step accurate calculation
-   29-BPA-3.inp - three Bisphenol A optimization using STO-2G split with commented example of two step accurate calculation

