---
layout: default
title: ICQCB CQC in practice
---

Computational Quantum Chemistry in Practice
-------------------------------------------

Although it may look very complex when seen from the point of view of the math involved and the sheer vastness of the options available, doing Computational Quantum Chemistry is not that difficult. All we need to do is keep in mind a few **basic concepts**:

-   we approximate atomic orbitals using a set of basis functions
-   we approximate molecular orbitals using a linear combination of atomic orbitals
-   we approximate the molecular wave function as a (determinantal) combination of molecular orbitals
-   we need to adjust parameters for each element in the calculation, usually using the Self Consistent Field method exploiting the variational principle
-   we can include electron correlation approximating the wave function as a combination of molecular wave functions for different excitation levels
-   traditional calculations depend on the fourth power of the number of basis functions
-   we can reduce the number of basis functions by concentrating on valence electrons
-   to reduce function complexity we can include empirical data
-   there are modern methods that can achieve linear scaling in the calculations

These points lead us to the full spectrum of options: when dealing with classical approaches, the most efficient are the **semiempirical methods** (they use valence electrons and substitute many calculations by empirical data or simplified approximations). They have traditionally been the most widely employed in Biology (MNDO, AM1, PM3, PM6).

More recently, DFT approaches have achieved similar efficiency and become very popular as well. **Density Functional Theory** is based on a different approach to electronic calculations: instead of reconstructing all orbitals, it concentrates on the overall electronic distribution. Without giving more detail, this methodology allows us to deal with bigger molecules and, because of that, it has gained momentum in the last decades, most specially its application in dynamic systems (*Time Dependent DFT, TD-DFT*) that simulate the dynamic aspects of electron behaviour (like excitation in answer to external stimuli such as radiation) and because of the early development of Quantum Molecular Dynamics methods based on the work of *Carr and Parrinello* (CPMD).

It is now feasible to carry out classical calculations using pure **ab initio** methods. For biological systems it is still common to use *STO-3G* or the more precise *6-311G* basis sets, using either *RHF* or *UHF* (depending on orbital occupancy). and, if we want to and our computing infrastructure does support it, we can also include correlation effects, like *CI* or *Møller-Plesset (MP)* at a higher computational cost.

**Linear scaling methods** use a number of novel shortcomings in the calculations to reduce computational cost, but software is becoming popular now, allowing us to carry out many of these calculations at a lower cost, a major achievement of interest in Molecular Biology. There are many tools available and choosing the best one for each given problem is still a difficult task.

A point worth noting is that although CPMD was an early development that has helped make DFT popular, in the last few years, similar Molecular Dynamics facilities have been added to most traditional methods as well, and so it is now becoming common to find studies using Semiempirical MD, or purely **ab-initio** MD.

<u>Regarding molecular analysis</u>, the results we will normally be interested in include finding the most stable conformation, the electronic structure (specially the frontier molecular orbitals HOMO and LUMO), the charge distribution, vibration frequencies, molecular flexibility and forces acting on different point of the molecule.

<u>Regarding functional chemical analysis</u>, we can also increase our understanding of reaction mechanisms by using saddle point and transition state calculations, as well as following the reaction paths. In so doing we will be able to spot relevant groups in drugs (helping us design better analogues or competitors) as well as amino acids involved and with a relevant role in the reaction (helping us design enhanced or inactive enzymes using site-directed mutagenesis).

Among many other things. Of course.

When dealing with macromolecules, it will often be sensible to reduce the analysis to the active site (taking the rest of the macromolecule as fixed) or to treat the active site usin QM and including the rest of the macromolecule represented with MM or MD.

We have already commented quickly that depending on the accuracy we want, we may need to use more or less computing power. We can provide gross estimates of the **computational cost** that can be expected for a given method, however, it is important to note as well that different authors are continuously working on enhancing them and increasing their efficiency, so that actual needs are reduced every year. And as a matter of fact, comparative reviews are being commonly published periodically (and you are advised to look for the most recent ones). However, we can, as a rule of thumb, expect something like the following table to hold for the default traditional methods:

|Method|Time|
|------|----|
|Molecular Mechanics|N\^2|
|Molecular Dynamics|N\^2|
|Hartree-Fock (HF)|N\^2 - N\^4 (depending on whether symmetry and cutoffs can be used)|
|MP2|N\^5|
|MP3, MP4|N\^6|
|MP5|N\^8|
|CC2|N\^5|
|CCSD|N\^6|
|CISD|N\^6|
|CASSF|A! (depends on orbital number)|
|Full CI|N!|
|Semi Empirical Methods|N\^2 (for small systems)|
|Semi Empirical Methods|N\^3 (for large systems)|
|DFT|N\^3|

Again, it is worth noting that improvements are continuously being made. And most importantly, that linear scaling methods are now becoming commonplace, significantly reducing computational needs.

Some commonly used examples of linear scaling programs you are likely to see are *MOPAC2009* for semi-empirical calculations, *SIESTA* for DFT, *ergoSCF* for ground state *ab-initio* calculations, and **FreeON** as an excellent general purpose, comprehensive tool.
