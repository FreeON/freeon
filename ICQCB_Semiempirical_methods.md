---
layout: default
title: ICQCB Semiempirical methods
---

![Aspirin EDM computed with PM3](Aspirin-elec_dens_map.png "Aspirin EDM computed with PM3")

Semiempirical methods
---------------------

Quantum Chemical calculations from first principles (*ab initio*) are purely theoretical, but endure some limitations that make them sometimes less desirable. On top of these, the number of integrals to evaluate is tremendous and grows exponentially in the number of basis functions (usually as N<sup>4</sup> or more) so that even dealing with a small molecule can involve millions of integrals.

There are approximations that can speed up computation taking advantage of the empirical knowledge we have gathered over time. Basically, what they do is substituting as many functions as possible by empirical valuies that have been measured or computed previously at higher precision levels, or by simple estimates derived from existing data.

All the methods commonly used rely on evaluating the valence bond electrons (VB), limiting the use of core electrons since their contribution to the chemical bond is minor. This approximation is specially interesting when the molecules to be studied contain a heavy atom (like metalloenzymes) as we can reduce significantly the number of functions to evaluate. The results are not too different from completely *ab initio* calculations, and we can obtain reductions in calculation times of about 50% or more.

Going back to the former paragraph, it is also worth noting that since we are substituting in the equation many properties by their empirically measured value, and the calculations are tuned and corrected to reproduce observed properties, and since these observed properties depend on the complete properties of the atoms, semiempirical methods can deal efficiently with transition metals -when they have been parametrized for them- without needing to include relativistic effects, as any needed correction is already included in the empirical parameters. Similarly, correlation, excitation and other effects are usually subsumed in the empirical corrections.

As we have mentioned. there are two approaches: one is substituting functions by empirical parameters, derived from experimental observations in compounds or bonds with similar properties. This implies that electron correlation effects are automatically considered, but also that, to the extend the parameters have been with obtained from experimental data that is more or less close to our problem, we will also get better or worst results. So, for example, one parameter set may have been derived to give good results for espectrometry calculations while another may have been parametrized for organic or inorganic chemistry...

The other approach is to use results derived from previous *ab initio* calculations. This approach is cleaner and does not introduce deviations to one kind of experiment or another, but on the other hand, it cannot yield better results than the *an initio* calculations themselves.

### Hückel

This was one of the first approaches derived. It considers only the PI valence electrons in a conjugated planar hydrocarbon. Its interest lies in its ability to produce some qualitative predictions, but it is seldomly used nowadays except for education (it can be computed by hand) or to obtain the initial guess for subsequent *ab initio* calculations.

Extended calculations can be performed modeling all valence orbitals using simple STO derived overlaps. The method is neither sofisticated nor precise, but it is still in use because it can be applied to all elements in the periodic table, specially in Inorganic Chemistry.

### Neglect of Differential Overlap

This is an approximation that neglects the differential overlap between atomic orbital basis functions centered around different nuclei. These overlaps tend to be small, but not so much as to allow us to ignore them without a penalty. Consequently, methods using this approach will usually add some parameters to adjust the calculation and make it comply with experiment.

There are various methods, from '''CNDO (Complete Neglect of Differential Overlap) ''' which is applied as such, to the more recents '''MNDO (Modified Neglect of Diatomic Overlap) ''', '''AM1 (Austin Model 1) ''', '''PM3 (Parameter Model 3) '''. Most of them are available in the public version of MOPAC (which is itself included in many commercial programs) and other *ab initio* packages. **PM6** has been developed recently and is currently avaiilable in a free -but closed and distribution-contolled- version of Mopac. The more recent methods tend to give generally satisfactory results and have been parametrized to include the most relevant atoms in Biochemistry (including several metals).

As these methods are much less costly in time and storage, and as they have been commonly parametrized for application in Biochemistry, **semiempirical methods have been the preferred methods to perform calculations in the Life and Health Sciences** (and often they have been the only available methods for many problems given their size). Things are changing now with the prevalence of linearly scaling methods, but that is quite a different story.

A large variety of semiempirical methods have been developed and parametrized for many different problems. In general we can expect that precision will be more or less

MNDO \< AM1 \< PM3

Today we normally use AM1 or PM3, or if we have the restricted version of MOPAC -which is free for academic use-, PM6. We can expect them to give good results in Organic Chemistry (for which they have been parametrized), the usually produce fairly good predictions of molecular geometries and energies, they can also predict vibrational mode and transition structures with less precision than *ab initio* methods and they usually give poor results for molecular Van der Waals forces and dispersion effects (but note that a special version of PM6 -PM6 DH2- has recently been proposed to deal specifically with these problems -and less well with everything else-).

As computing power has increased, it is now possible to use pure *ab initio* methods in many cases (mainly when analyzing drugs), allowing us to start considering the use of *ab initio* methods which are theoretically cleaner.
