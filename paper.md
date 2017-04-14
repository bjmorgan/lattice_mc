---
title: 'lattice_mc: A Python Lattice-Gas Monte Carlo Module'
tags:
  - lattice-gas
  - monte-carlo
  - correlation-factors
authors:
 - name: Benjamin J. Morgan
   orcid: 0000-0002-3056-8233
   affiliation: 1
affiliations:
 - name: Department of Chemistry, University of Bath, Bath, BA2 7AY, United Kingdom.
   index: 1
date: 13 April 2017
bibliography: paper.bib
---

# Summary

`lattice_mc` is Python module for running (kinetic) lattice-gas Monte Carlo simulations. 
Simple lattices can be constructed programmatically (presently square, honeycomb, and cubic lattices).
Arbitrary lattices can be generated from files that define the lattice sites and their connectivity. The algorithms used and interaction models are described in <a href="#ref1">\[1\]</a>. Calculated properties include tracer and &ldquo;jump&rdquo; diffusion coefficients; where the latter is proportional to the mobility (and hence the conductivity for charged particles) <a href="#ref2">\[2\]</a>; and tracer (single particle) and collective correlation factors, f and f_I <a href="#ref3">\[3\]</a>. The simplest interaction model is for &ldquo;non-interacting&rdquo; particles, where the only restriction is volume exclusion (two particles cannot simultaneoulsy occupy a single site) <a href="#ref1">\[4\]</a>. Additional interaction models include nearest-neighbour repulsion and on-site energies for inequivalent sites.

# References

<span id='ref1'>[1] B. J. Morgan, In Preparation.</span>  
<span id='ref2'>[2] A. Van der Ven *et al.* [*Acc. Chem. Res.* **46**, 1216 (2013)](https://dx.doi.org/10.1021/ar200329r)</span>  
<span id='ref3'>[3] G. E. Murch [*Sol. Stat. Ionics* **7**, 177 (1982)](https://dx.doi.org/10.1016/0167-2738%2882%2990050-9)</span>  
<span id='ref4'>[4] R. Kutner [*Phys. Lett.* **81A**, 239 (1981)](https://dx.doi.org/10.1016/0375-9601%2881%2990251-6)
</span>  
<span id='ref5'>\[5\] Morgan and Madden, [*J. Phys. Condens. Matter* **24**, 275303 (2012)](http://www.iopscience.iop.org/article/10.1088/0953-8984/24/27/275303/)</span>.  
<span id='ref6'> \[6\] G. E. Murch & R. J. Thorn, [Phil. Mag. **36** 529 (1977)](http://dx.doi.org/10.1080/14786437708239737).
