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

In solid electrolytes, ionic motion is typically effected by a series of discrete "jumps" where ions move between adjacent lattice sites [@Catlow_SolStatIonics1983]. For dilute mobile ions, ionic trajectories are random walks, and simple analytical expressions exits that relate macroscopic transport coefficients, i.e. diffusion coefficients and ionic conductivities, to the microscopic jump frequency of individual ions [@HowardAndLidiard_RepProgPhys1964;@Stoneham_IonicSolids1989]. Practical solid electrolytes have high mobile ion concentrations, with significnat interparticle interactions producing deviations from the dilute limit random walk behaviour. In general, the *quantitative* effect of interparticle interactions cannot be determined analytically. As an alternative, numerical simulations, such as lattice-gas Monte Carlo methods, can be used to directly calculate these relationships. Lattice-gas Monte Carlo methods are particularly suited to studying how varying properties across broad classes of materials give quantitative differences in macroscopic ionic transport, and can be used to understand the different transport properties of materials with, for example, different crystal structures or mobile ion stoichiometries. 

`lattice_mc` has been written to allow materials scientists and solid-state chemists model how the microscopic physics of solid electrolytes (crystal structure, stoichiometry, interaction models) determine macroscopic transport behaviour (diffusion and ionic conduction), with the goal of understand the factors that make different materials more or less useful for specific applications (e.g. Li-ion batteries or fuel cells).

The code allows the programmatic construction of simple lattices (presently implemented are square, honeycomb, and cubic lattices). Lattices with arbitrary geometries can be constructed from a file format that defines the lattice sites and their connectivity, allowing models based on crystallographic data. Calculated properties include tracer and &ldquo;jump&rdquo; diffusion coefficients; where the latter is proportional to the mobility (and hence the conductivity for charged particles) [@VanDerVenEtAl_AccChemRes2013]; and tracer (single particle) and collective correlation factors, $f$ and $f_I$ [@Murch_SolStatIonics1982] The simplest interaction model is for &ldquo;non-interacting&rdquo; particles, where the only restriction is volume exclusion (two particles cannot simultaneously occupy a single site) [@Kutner_PhysLett1981]. Additional interaction models include nearest-neighbour repulsion and on-site energies for inequivalent sites. Simulations are performed using an efficient rejection-free Monte Carlo scheme [@Voter_kMCmethod].

# Acknowledgements

BJM acknowledges support from the Royal Society (UF130329).

# References

