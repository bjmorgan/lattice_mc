# lattice_mc

`lattice_mc` is Python module for running (kinetic) lattice-gas Monte Carlo simulations. Simple lattices can be constructed programmatically (presently square, honeycomb, and cubic lattices). Arbitrary lattices can be generated from files that define the lattice sites and their connectivity. The algorithms used and interaction models are described in <a href="#ref1">\[1\]</a>. Calculated properties include tracer and &ldquo;jump&rdquo; diffusion coefficients; where the latter is proportional to the mobility (and hence the conductivity for charged particles) <a href="#ref2">\[2\]</a>; and tracer (single particle) and collective correlation factors, f and f_I <a href="#ref3">\[3\]</a>. The simplest interaction model is for &ldquo;non-interacting&rdquo; particles, where the only restriction is volume exclusion (two particles cannot simultaneously occupy a single site) <a href="#ref1">\[4\]</a>. Additional interaction models include nearest-neighbour repulsion and on-site energies for inequivalent sites.

## Installation

Download the latest release from [GitHub](https://github.com/bjmorgan/lattice_mc/releases)
```
https://github.com/bjmorgan/lattice_mc/archive/1.0.0.tar.gz
```
Then install
```
cd lattice_mc
python setup.py install
```

Or you can clone the latest development version:
```
git clone git@github.com:bjmorgan/lattice_mc.git
```
and install the same way.
```
cd lattice_mc
python setup.py install
```

Alternatively, you can install using `pip`, e.g.
```
pip3 install git+https://github.com/bjmorgan/lattice_mc.git
```

## Documentation

Full documentation and examples are contained in a [Jupyter notebook](http://jupyter-notebook.readthedocs.io/en/latest/#) at [examples/lattice_mc_example.ipynb](examples/lattice_mc_example.ipynb). This example notebook is also hosted on [GitHub](https://github.com/bjmorgan/lattice_mc/blob/master/examples/lattice_mc_examples.ipynb).

## Tests

```
python3 -m unittest discover
```

## References
<span id='ref1'>[1] B. J. Morgan, In Preparation.</span>  
<span id='ref2'>[2] A. Van der Ven *et al.* [*Acc. Chem. Res.* **46**, 1216 (2013)](https://dx.doi.org/10.1021/ar200329r)</span>  
<span id='ref3'>[3] G. E. Murch [*Sol. Stat. Ionics* **7**, 177 (1982)](https://dx.doi.org/10.1016/0167-2738%2882%2990050-9)</span>  
<span id='ref4'>[4] R. Kutner [*Phys. Lett.* **81A**, 239 (1981)](https://dx.doi.org/10.1016/0375-9601%2881%2990251-6)
</span>  
<span id='ref5'>\[5\] Morgan and Madden, [*J. Phys. Condens. Matter* **24**, 275303 (2012)](http://www.iopscience.iop.org/article/10.1088/0953-8984/24/27/275303/)</span>.  
<span id='ref6'> \[6\] G. E. Murch & R. J. Thorn, [Phil. Mag. **36** 529 (1977)](http://dx.doi.org/10.1080/14786437708239737).</span>
