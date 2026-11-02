Introduction
============

``lattice_mc`` is a Python module for performing (kinetic) lattice-gas Monte Carlo
(LGMC) simulations of ionic transport in solid electrolytes.

Lattice-gas Monte Carlo
-----------------------

In solid electrolytes, ions typically move by a series of discrete jumps between
adjacent lattice sites [Catlow1983]_. For dilute mobile ions, ionic trajectories are
random walks, and simple analytical expressions relate macroscopic transport
coefficients, such as diffusion coefficients and ionic conductivities, to the jump
frequency of individual ions [HowardLidiard1964]_ [Harding1989]_. Practical solid
electrolytes have high mobile ion concentrations, and interactions between the
mobile ions give deviations from dilute random-walk behaviour. In general the size
of these deviations cannot be calculated analytically, but it can be calculated
directly from lattice-gas Monte Carlo simulations.

``lattice_mc`` models a lattice of sites, some of which are occupied by mobile atoms.
Each atom can jump to a neighbouring vacant site. The rate of each jump depends on
the change in energy it would cause, through an interaction model (see
:doc:`interaction_models`). Simulations use a rejection-free kinetic Monte Carlo
scheme [Voter2007]_. The algorithms and interaction models are described in more
detail in [Morgan2017]_.

What ``lattice_mc`` calculates
------------------------------

For atoms with displacements :math:`\mathbf{r}_i` after a simulation time :math:`t`:

Tracer diffusion coefficient, :math:`D^*`
    The mean squared displacement of individual atoms, from the Einstein relation

    .. math::

        \langle r_i^2 \rangle = 2 d D^* t.

Collective, or "jump", diffusion coefficient, :math:`D_\mathrm{J}`
    The squared displacement of the summed atom displacements,

    .. math::

        \left\langle \left( \sum_i \mathbf{r}_i \right)^2 \right\rangle = 2 d D_\mathrm{J} t.

    :math:`D_\mathrm{J}` is proportional to the mobility, and hence to the
    conductivity for charged particles [VanDerVen2013]_. No :math:`1/N`
    normalisation is applied; the per-atom value is also available.

Tracer correlation factor, :math:`f`, and collective correlation factor, :math:`f_\mathrm{I}`
    The ratios of the tracer and collective squared displacements to the sum of the
    squared lengths of the individual jumps [Murch1982]_. Both are 1 for uncorrelated
    random walks.

In the Einstein relations, :math:`d` is the dimensionality of the lattice, given by
:attr:`Lattice.dimensionality <lattice_mc.lattice.Lattice.dimensionality>`: the
number of non-zero cell lengths. Two-dimensional lattices, such as square and
honeycomb lattices, have a z cell length of 0, so :math:`d = 2`.

Installation
------------

.. code-block:: bash

    pip install lattice_mc

Or install the latest development version directly from GitHub:

.. code-block:: bash

    pip install git+https://github.com/bjmorgan/lattice_mc.git

``lattice_mc`` requires Python 3.11 or above.

References
----------

.. [Catlow1983] C. R. A. Catlow, *Sol. Stat. Ionics* **8**, 89 (1983).
   https://doi.org/10.1016/0167-2738(83)90069-3

.. [HowardLidiard1964] R. E. Howard and A. B. Lidiard, *Rep. Prog. Phys.* **27**, 161 (1964).
   https://doi.org/10.1088/0034-4885/27/1/305

.. [Harding1989] J. H. Harding, Defects and Transport in Ionic Solids, in *Ionic Solids
   at High Temperatures*, ed. A. M. Stoneham, World Scientific (1989).
   https://doi.org/10.1142/9789814503228_0003

.. [Morgan2017] B. J. Morgan, *R. Soc. Open Sci.* **4**, 170824 (2017).
   https://doi.org/10.1098/rsos.170824

.. [VanDerVen2013] A. Van der Ven *et al.*, *Acc. Chem. Res.* **46**, 1216 (2013).
   https://doi.org/10.1021/ar200329r

.. [Murch1982] G. E. Murch, *Sol. Stat. Ionics* **7**, 177 (1982).
   https://doi.org/10.1016/0167-2738(82)90050-9

.. [Voter2007] A. F. Voter, Introduction to the Kinetic Monte Carlo Method, in
   *Radiation Effects in Solids*, ed. K. E. Sickafus *et al.*, Springer (2007).
   https://doi.org/10.1007/978-1-4020-5295-8_1
