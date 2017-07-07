lattice\_mc: A Python Lattice-Gas Monte Carlo Module
====================================================

|PyPI version| |Build Status| |DOI| |JOSS status|

``lattice_mc`` is Python module for performing (kinetic) lattice-gas
Monte Carlo (LGMC) simulations of ionic transport in solid electrolytes.

In solid electrolytes, ionic motion is typically effected by a series of
discrete “jumps” where ions move between adjacent lattice sites
[`1 <#Catlow_SolStatIonics1983>`__]. For dilute mobile ions, ionic
trajectories are random walks, and simple analytical expressions exits
that relate macroscopic transport coefficients, i.e. diffusion
coefficients and ionic conductivities, to the microscopic jump frequency
of individual ions
[`2 <#HowardAndLidiard_RepProgPhys1964>`__,\ `3 <#Stoneham_IonicSolids1989>`__].
Practical solid electrolytes have high mobile ion concentrations, with
significnat interparticle interactions producing deviations from the
dilute limit random walk behaviour. In general, the *quantitative*
effect of interparticle interactions cannot be determined analytically.
As an alternative, numerical simulations, such as lattice-gas Monte
Carlo methods, can be used to directly calculate these relationships.
Lattice-gas Monte Carlo methods are particularly suited to studying how
varying properties across broad classes of materials give quantitative
differences in macroscopic ionic transport, and can be used to
understand the different transport properties of materials with, for
example, different crystal structures or mobile ion stoichiometries.

``lattice_mc`` has been written to allow materials scientists and
solid-state chemists model how the microscopic physics of solid
electrolytes (crystal structure, stoichiometry, interaction models)
determine macroscopic transport behaviour (diffusion and ionic
conduction), with the goal of understand the factors that make different
materials more or less useful for specific applications (e.g. Li-ion
batteries or fuel cells).

The code allows the programmatic construction of simple lattices
(presently implemented are square, honeycomb, and cubic lattices).
Lattices with arbitrary geometries can be constructed from a file format
that defines the lattice sites and their connectivity, allowing models
based on crystallographic data. The algorithms used and interaction
models are described in more detail in Ref. [`4 <#Morgan_LLZO>`__].
Calculated properties include tracer and “jump” diffusion coefficients;
where the latter is proportional to the mobility (and hence the
conductivity for charged particles)
[`5 <#VanDerVenEtAl_AccChemRes2013>`__]; and tracer (single particle)
and collective correlation factors, *f* and *f*\ I
[`6 <#Murch_SolStatIonics1982>`__]. The simplest interaction model is
for “non-interacting” particles, where the only restriction is volume
exclusion (two particles cannot simultaneously occupy a single site)
[`7 <#Kutner_PhysLett1981>`__]. Additional interaction models include
nearest-neighbour repulsion and on-site energies for inequivalent sites.
Simulations are performed using an efficient rejection-free Monte Carlo
scheme [`8 <#Voter_kMCmethod>`__].

Installation
------------

::

    pip install lattice_mc

Or download the latest release from
`GitHub <https://github.com/bjmorgan/lattice_mc/releases>`__

::

    https://github.com/bjmorgan/lattice_mc/archive/1.0.0.tar.gz

Then install

::

    cd lattice_mc
    python setup.py install

Or you can clone the latest development version:

::

    git clone git@github.com:bjmorgan/lattice_mc.git

and install the same way.

::

    cd lattice_mc
    python setup.py install

Alternatively, you can install the latest build using ``pip``, direct
from GitHub, e.g.

::

    pip3 install git+https://github.com/bjmorgan/lattice_mc.git

Documentation
-------------

Full documentation and examples are contained in a `Jupyter
notebook <http://jupyter-notebook.readthedocs.io/en/latest/#>`__ at
`examples/lattice\_mc\_example.ipynb <examples/lattice_mc_example.ipynb>`__.
The example notebook is also hosted on
`GitHub <https://github.com/bjmorgan/lattice_mc/blob/master/examples/lattice_mc_examples.ipynb>`__.

Tests
-----

Automated testing of the latest build happens
`here <https://travis-ci.org/bjmorgan/lattice_mc>`__.

Manual tests can be run using

::

    python3 -m unittest discover

The code has been tested with Python versions 3.5 and above.

References
----------

1.  \ `C. R. A. Catlow, *Sol. Stat. Ionics* **8**, 89
    (1983). <https://doi.org/10.1016/0167-2738%2883%2990069-3>`__
2.  \ `R. E. Howard and A. B. Lidiard, *Rep. Prog. Phys.* **27**, 161
    (1964). <https://doi.org/10.1088/0034-4885/27/1/305>`__
3.  \ `J. H. Harding, Defects and Transport in Ionic Solids, in *Ionic
    Solids at High Temperatures* ed. A. M. Stoneham, World Scientific
    (1989). <https://doi.org/10.1142/9789814503228_0003>`__
4.  \ `B. J. Morgan, *arXiv*:
    1707.00491 <https://arxiv.org/abs/1707.00491>`__
5.  \ `A. Van der Ven *et al.* *Acc. Chem. Res.* **46**, 1216
    (2013). <https://dx.doi.org/10.1021/ar200329r>`__
6.  \ `G. E. Murch *Sol. Stat. Ionics* **7**, 177
    (1982). <https://dx.doi.org/10.1016/0167-2738%2882%2990050-9>`__
7.  \ `R. Kutner *Phys. Lett.* **81A**, 239
    (1981). <https://dx.doi.org/10.1016/0375-9601%2881%2990251-6>`__
8.  \ `A. F. Voter, Introduction to the Kinetic Monte Carlo Method, in
    *Radiation Effects in Solids*, ed. K. E. Sicafus *et al.*, Springer
    (2007). <https://doi.org/10.1007/978-1-4020-5295-8_1>`__
9.  \ `Morgan and Madden, *J. Phys. Condens. Matter* **24**, 275303
    (2012). <http://www.iopscience.iop.org/article/10.1088/0953-8984/24/27/275303/>`__
10. \ `G. E. Murch & R. J. Thorn, *Phil. Mag.* **36** 529
    (1977). <http://dx.doi.org/10.1080/14786437708239737>`__

.. |PyPI version| image:: https://badge.fury.io/py/lattice-mc.svg
   :target: https://badge.fury.io/py/lattice-mc
.. |Build Status| image:: https://travis-ci.org/bjmorgan/lattice_mc.svg?branch=master
   :target: https://travis-ci.org/bjmorgan/lattice_mc
.. |DOI| image:: https://zenodo.org/badge/75750912.svg
   :target: https://zenodo.org/badge/latestdoi/75750912
.. |JOSS status| image:: http://joss.theoj.org/papers/6940b7bb0d59be86b8823a10780caae0/status.svg
   :target: http://joss.theoj.org/papers/6940b7bb0d59be86b8823a10780caae0
