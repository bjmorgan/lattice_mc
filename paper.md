---
title: 'lattice_mc: A Python Lattice-Gas Monte Carlo Module'  
tags:  
  - lattice-gas  
  - monte-carlo  
  - correlation-factors  
authors:  
 - name: Benjamin J. Morgan  
   email: b.j.morgan@bath.ac.uk  
   orcid: 0000-0002-3056-8233  
   affiliation: 1  
affiliations:  
 - name: Department of Chemistry, University of Bath, Bath, BA2 7AY, United Kingdom.  
   index: 1  
date: 13 April 2017  
bibliography: paper.bib
---

# Summary

`lattice_mc` is Python module for performing (kinetic) lattice-gas Monte Carlo (LGMC) simulations of ionic transport in solid electrolytes.

In solid electrolytes, ionic motion is typically effected by a series of discrete "jumps" where ions move between adjacent lattice sites \[[1](#Catlow_SolStatIonics1983)\]. For dilute mobile ions, ionic trajectories are random walks, and simple analytical expressions exits that relate macroscopic transport coefficients, i.e. diffusion coefficients and ionic conductivities, to the microscopic jump frequency of individual ions \[[2](#HowardAndLidiard_RepProgPhys1964),[3](#Stoneham_IonicSolids1989)\]. Practical solid electrolytes have high mobile ion concentrations, with significnat interparticle interactions producing deviations from the dilute limit random walk behaviour. In general, the *quantitative* effect of interparticle interactions cannot be determined analytically. As an alternative, numerical simulations, such as lattice-gas Monte Carlo methods, can be used to directly calculate these relationships. Lattice-gas Monte Carlo methods are particularly suited to studying how varying properties across broad classes of materials give quantitative differences in macroscopic ionic transport, and can be used to understand the different transport properties of materials with, for example, different crystal structures or mobile ion stoichiometries. 

`lattice_mc` has been written to allow materials scientists and solid-state chemists model how the microscopic physics of solid electrolytes (crystal structure, stoichiometry, interaction models) determine macroscopic transport behaviour (diffusion and ionic conduction), with the goal of understand the factors that make different materials more or less useful for specific applications (e.g. Li-ion batteries or fuel cells).

The code allows the programmatic construction of simple lattices (presently implemented are square, honeycomb, and cubic lattices). Lattices with arbitrary geometries can be constructed from a file format that defines the lattice sites and their connectivity, allowing models based on crystallographic data. The algorithms used and interaction models are described in more detail in Ref. \[[4](#Morgan_LLZO)\]. Calculated properties include tracer and &ldquo;jump&rdquo; diffusion coefficients; where the latter is proportional to the mobility (and hence the conductivity for charged particles) \[[5](#VanDerVenEtAl_AccChemRes2013)\]; and tracer (single particle) and collective correlation factors, *f* and *f*<sub>I</sub> \[6\]</a>. The simplest interaction model is for &ldquo;non-interacting&rdquo; particles, where the only restriction is volume exclusion (two particles cannot simultaneously occupy a single site) \[[7](#Kutner_PhysLett1981)\]. Additional interaction models include nearest-neighbour repulsion and on-site energies for inequivalent sites. Simulations are performed using an efficient rejection-free Monte Carlo scheme \[[8](#Voter_kMCmethod)\].

# Acknowledgements

BJM acknowledges support from the Royal Society (UF130329).

# References

## References
1. <a name="Catlow_SolStatIonics1983" />[C. R. A. Catlow, *Sol. Stat. Ionics* **8**, 89 (1983).](https://doi.org/10.1016/0167-2738%2883%2990069-3)
2. <a name="HowardAndLidiard1964" />[R. E. Howard and A. B. Lidiard, *Rep. Prog. Phys.* **27**, 161 (1964).](https://doi.org/10.1088/0034-4885/27/1/305)
3. <a name="Stoneham_IonicSolids1989" />[J. H. Harding, Defects and Transport in Ionic Solids, in *Ionic Solids at High Temperatures* ed. A. M. Stoneham, World Scientific (1989)](https://doi.org/10.1142/9789814503228_0003)
4. <a name="Morgan_LLZO" />B. J. Morgan, In Preparation.
5. <a name="VanDerVenEtAl_AccChemRes2013" />[A. Van der Ven *et al.* *Acc. Chem. Res.* **46**, 1216 (2013)](https://dx.doi.org/10.1021/ar200329r)
6. <a name="Murch_SolStationics1982" />[G. E. Murch *Sol. Stat. Ionics* **7**, 177 (1982)](https://dx.doi.org/10.1016/0167-2738%2882%2990050-9)
7. <a name="Kutner_PhysLett1981" />[R. Kutner *Phys. Lett.* **81A**, 239 (1981)](https://dx.doi.org/10.1016/0375-9601%2881%2990251-6)
8. <a name="Voter_kMCmethod" />[A. F. Voter, Introduction to the Kinetic Monte Carlo Method, in *Radiation Effects in Solids*, ed. K. E. Sicafus *et al.*, Springer( 2007)](https://doi.org/10.1007/978-1-4020-5295-8_1)
9. <a name="MorganAndMadden_JPhysCondensMatter2012" />[Morgan and Madden, *J. Phys. Condens. Matter* **24**, 275303 (2012)](http://www.iopscience.iop.org/article/10.1088/0953-8984/24/27/275303/)
10. <a name="MurchAndThorn_PhilMag1977" />[G. E. Murch & R. J. Thorn, *Phil. Mag.* **36** 529 (1977)](http://dx.doi.org/10.1080/14786437708239737)
