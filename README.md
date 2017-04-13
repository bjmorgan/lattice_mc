# lattice_mc

`lattice_mc` is Python module for running (kinetic) lattice-gas Monte Carlo simulations. Simple lattices can be constructed programmatically (presently square, honeycomb, and cubic lattices). Arbitrary lattices can be generated from files that define the lattice sites and their connectivity. The algorithms used and interaction models are described in <a href="#ref1">\[1\]</a>. Calculated properties include tracer and &ldquo;jump&rdquo; diffusion coefficients; where the latter is proportional to the mobility (and hence the conductivity for charged particles) <a href="#ref2">\[2\]</a>; and tracer (single particle) and collective correlation factors, f and f_I <a href="#ref3">\[3\]</a>. The simplest interaction model is for &ldquo;non-interacting&rdquo; particles, where the only restriction is volume exclusion (two particles cannot simultaneously occupy a single site) <a href="#ref1">\[4\]</a>. Additional interaction models include nearest-neighbour repulsion and on-site energies for inequivalent sites.

## Installation

TODO

## Documentation

Full documentation and examples are contained in a [Jupyter notebook](http://jupyter-notebook.readthedocs.io/en/latest/#) at [examples/lattice_mc_example.ipynb](examples/lattice_mc_example.ipynb)

## Unit tests

```
python3 -m unittest discover
```

## References
[1] B. J. Morgan, In Preparation.  
[2] A. Van der Ven et al. Acc. Chem. Res. 46, 1216 (2013).  
[3] G. E. Murch Sol. Stat. Ionics 7, 177 (1982).  
[4] R. Kutner Phys. Lett. 81A, 239 (1981).  
[5] Morgan and Madden, J. Phys. Condens. Matter 24, 275303 (2012).  
[6] G. E. Murch & R. J. Thorn, Phil. Mag. 36 529 (1977).
