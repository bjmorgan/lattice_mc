Introduction
============

``lattice_mc`` is a Python module for kinetic lattice-gas Monte Carlo simulations of
ionic transport.

The model
---------

A simulation has a lattice of sites, some of which are occupied by atoms. Each atom
can jump to a neighbouring vacant site, at a rate set by the change in energy the jump
would cause under the chosen interaction model (see :doc:`interaction_models`). Jumps
are selected with a rejection-free kinetic Monte Carlo scheme [Voter2007]_. The
algorithms and interaction models are described in [Morgan2017]_.

What ``lattice_mc`` calculates
------------------------------

After a simulation of :math:`N` atoms for a time :math:`t`, where atom :math:`i` has
total displacement :math:`\mathbf{r}_i`, and :math:`\sum l^2` is the sum of the squared
lengths of every jump made:

Tracer diffusion coefficient, :math:`D^*`
    :attr:`~lattice_mc.simulation.Simulation.tracer_diffusion_coefficient`

    .. math::

        D^* = \frac{\sum_i r_i^2}{2 d N t}

Collective, or "jump", diffusion coefficient, :math:`D_\mathrm{J}`
    :attr:`~lattice_mc.simulation.Simulation.collective_diffusion_coefficient`

    .. math::

        D_\mathrm{J} = \frac{\left( \sum_i \mathbf{r}_i \right)^2}{2 d t}

    :attr:`~lattice_mc.simulation.Simulation.collective_diffusion_coefficient_per_atom`
    gives :math:`D_\mathrm{J}/N`.

Tracer correlation factor, :math:`f`
    :attr:`~lattice_mc.simulation.Simulation.tracer_correlation`

    .. math::

        f = \frac{\sum_i r_i^2}{\sum l^2}

Collective correlation factor, :math:`f_\mathrm{I}`
    :attr:`~lattice_mc.simulation.Simulation.collective_correlation`

    .. math::

        f_\mathrm{I} = \frac{\left( \sum_i \mathbf{r}_i \right)^2}{\sum l^2}

Diffusion coefficients are in units of (site spacing)\ :sup:`2` s\ :sup:`-1`.
:math:`d` is the dimensionality of the lattice, given by
:attr:`Lattice.dimensionality <lattice_mc.lattice.Lattice.dimensionality>`: the
number of non-zero cell lengths. Square and honeycomb lattices have a z cell length
of 0, so :math:`d = 2`.

Installation
------------

.. code-block:: bash

    pip install lattice_mc

Or install the latest development version directly from GitHub:

.. code-block:: bash

    pip install git+https://github.com/bjmorgan/lattice_mc.git

``lattice_mc`` requires Python 3.12 or above.

References
----------

.. [Morgan2017] B. J. Morgan, *R. Soc. Open Sci.* **4**, 170824 (2017).
   https://doi.org/10.1098/rsos.170824

.. [Voter2007] A. F. Voter, Introduction to the Kinetic Monte Carlo Method, in
   *Radiation Effects in Solids*, ed. K. E. Sickafus *et al.*, Springer (2007).
   https://doi.org/10.1007/978-1-4020-5295-8_1
