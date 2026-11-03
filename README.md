# lattice_mc: A Python Lattice-Gas Monte Carlo Module

[![PyPI version](https://badge.fury.io/py/lattice-mc.svg)](https://badge.fury.io/py/lattice-mc)
[![DOI](https://zenodo.org/badge/75750912.svg)](https://zenodo.org/badge/latestdoi/75750912)
[![JOSS status](http://joss.theoj.org/papers/6940b7bb0d59be86b8823a10780caae0/status.svg)](http://joss.theoj.org/papers/6940b7bb0d59be86b8823a10780caae0)
[![Documentation Status](https://readthedocs.org/projects/lattice-mc/badge/?version=latest)](http://lattice-mc.readthedocs.io/en/latest/?badge=latest)


`lattice_mc` is a Python module for performing (kinetic) lattice-gas Monte Carlo (LGMC) simulations of ionic transport in solid electrolytes.

In solid electrolytes, ionic motion is typically effected by a series of discrete &ldquo;jumps&rdquo; where ions move between adjacent lattice sites \[[1](#Catlow_SolStatIonics1983)\]. For dilute mobile ions, ionic trajectories are random walks, and simple analytical expressions exist that relate macroscopic transport coefficients, i.e. diffusion coefficients and ionic conductivities, to the microscopic jump frequency of individual ions \[[2](#HowardAndLidiard_RepProgPhys1964),[3](#Stoneham_IonicSolids1989)\]. Practical solid electrolytes have high mobile ion concentrations, with significant interparticle interactions producing deviations from the dilute limit random walk behaviour. In general, the *quantitative* effect of interparticle interactions cannot be determined analytically. As an alternative, numerical simulations, such as lattice-gas Monte Carlo methods, can be used to directly calculate these relationships. Lattice-gas Monte Carlo methods are particularly suited to studying how varying properties across broad classes of materials give quantitative differences in macroscopic ionic transport, and can be used to understand the different transport properties of materials with, for example, different crystal structures or mobile ion stoichiometries.

`lattice_mc` has been written to allow materials scientists and solid-state chemists to model how the microscopic physics of solid electrolytes (crystal structure, stoichiometry, interaction models) determine macroscopic transport behaviour (diffusion and ionic conduction), with the goal of understanding the factors that make different materials more or less useful for specific applications (e.g. Li-ion batteries or fuel cells).

`lattice_mc` builds square, honeycomb and cubic lattices, and reads lattices of any geometry from a sites file that defines the sites and their connectivity. Atoms jump between neighbouring sites, with volume exclusion: no two atoms occupy the same site \[[4](#Kutner_PhysLett1981)\]. Site energies, a nearest-neighbour interaction, and coordination-number energies, which depend on the number of occupied neighbours with each label, can be added in any combination. Simulations use a rejection-free kinetic Monte Carlo scheme \[[5](#Voter_kMCmethod)\]; the algorithms are described in more detail in Ref. \[[6](#Morgan_LLZO)\]. A simulation reports the tracer and collective diffusion coefficients, the reduced ionic conductivity, the tracer and collective correlation factors, *f* and *f*<sub>I</sub> \[[7](#Murch_SolStatIonics1982)\], and the time-averaged number of atoms on the sites with each label.

## Installation

```
pip install lattice_mc
```

Or install the latest development version directly from GitHub:
```
pip install git+https://github.com/bjmorgan/lattice_mc.git
```

For development, clone the repository and install in editable mode with dev dependencies:
```
git clone git@github.com:bjmorgan/lattice_mc.git
cd lattice_mc
pip install -e ".[dev]"
```

## Usage

```python
import lattice_mc as lmc

lattice = lmc.cubic_lattice(4, 4, 4, spacing=1.0)

simulation = lmc.Simulation(
    lattice,
    params=lmc.SimulationParameters(temperature=298.0, rate_prefactor=1e13),
    n_atoms=16,
    nn_energy=0.1,
    lookup_table=True,
    seed=42,
)
simulation.run(n_jumps=1_000)  # equilibration
simulation.reset()
simulation.run(n_jumps=10_000)
results = simulation.results

print(f"D* = {results.tracer_diffusion_coefficient:.3e}")
print(f"f  = {results.tracer_correlation:.3f}")
```

## Documentation

Documentation, including a user guide and API reference, is at [lattice-mc.readthedocs.io](https://lattice-mc.readthedocs.io/).

## Tests

Tests can be run using
```
pytest
```

The code requires Python 3.12 or above.

## References
1. <a name="Catlow_SolStatIonics1983" />[C. R. A. Catlow, *Sol. Stat. Ionics* **8**, 89 (1983).](https://doi.org/10.1016/0167-2738%2883%2990069-3)
1. <a name="HowardAndLidiard_RepProgPhys1964" />[R. E. Howard and A. B. Lidiard, *Rep. Prog. Phys.* **27**, 161 (1964).](https://doi.org/10.1088/0034-4885/27/1/305)
1. <a name="Stoneham_IonicSolids1989" />[J. H. Harding, Defects and Transport in Ionic Solids, in *Ionic Solids at High Temperatures* ed. A. M. Stoneham, World Scientific (1989).](https://doi.org/10.1142/9789814503228_0003)
1. <a name="Kutner_PhysLett1981" />[R. Kutner *Phys. Lett.* **81A**, 239 (1981).](https://dx.doi.org/10.1016/0375-9601%2881%2990251-6)
1. <a name="Voter_kMCmethod" />[A. F. Voter, Introduction to the Kinetic Monte Carlo Method, in *Radiation Effects in Solids*, ed. K. E. Sickafus *et al.*, Springer (2007).](https://doi.org/10.1007/978-1-4020-5295-8_1)
1. <a name="Morgan_LLZO" />[B. J. Morgan, *R. Soc. Open Sci.* **4**, 170824 (2017).](http://dx.doi.org/10.1098/rsos.170824)
1. <a name="Murch_SolStatIonics1982" />[G. E. Murch *Sol. Stat. Ionics* **7**, 177 (1982).](https://dx.doi.org/10.1016/0167-2738%2882%2990050-9)
