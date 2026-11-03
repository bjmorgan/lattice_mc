Introduction
============

``lattice_mc`` is a Python module for kinetic lattice-gas Monte Carlo simulations of
ionic transport.

The model
---------

A simulation has a lattice of sites, some of which are occupied by atoms. Each atom
can jump to a neighbouring vacant site, at a rate set by the change in energy the jump
would cause under the chosen interaction model (see :doc:`interaction_models`). Jumps
are selected with a rejection-free kinetic Monte Carlo scheme [Voter2007]_: each step
chooses one possible jump with probability proportional to its rate, and advances the
time by :math:`-\ln u / \sum_j k_j`, where :math:`k_j` are the rates of all possible
jumps and :math:`u` is drawn uniformly from :math:`(0, 1]`. The algorithms and
interaction models are described in [Morgan2017]_.

What ``lattice_mc`` calculates
------------------------------

:attr:`Simulation.results <lattice_mc.simulation.Simulation.results>` returns a
:class:`~lattice_mc.results.Results`, computed from every jump made since the
simulation was created or last reset. For :math:`N` atoms on a lattice of
:math:`N_\mathrm{sites}` sites, after a time :math:`t`, where atom :math:`i` has
total displacement :math:`\mathbf{r}_i`, and :math:`\sum l^2` is the sum of the
squared lengths of every jump made:

Tracer diffusion coefficient, :math:`D^*`
    :attr:`~lattice_mc.results.Results.tracer_diffusion_coefficient`

    .. math::

        D^* = \frac{\sum_i r_i^2}{2 d N t}

Collective, or "jump", diffusion coefficient, :math:`D_\mathrm{J}`
    :attr:`~lattice_mc.results.Results.collective_diffusion_coefficient`

    .. math::

        D_\mathrm{J} = \frac{\left( \sum_i \mathbf{r}_i \right)^2}{2 d N t}

Reduced ionic conductivity
    :attr:`~lattice_mc.results.Results.reduced_ionic_conductivity`

    .. math::

        \frac{N D_\mathrm{J}}{N_\mathrm{sites}}

Tracer correlation factor, :math:`f`
    :attr:`~lattice_mc.results.Results.tracer_correlation`

    .. math::

        f = \frac{\sum_i r_i^2}{\sum l^2}

Collective correlation factor, :math:`f_\mathrm{I}`
    :attr:`~lattice_mc.results.Results.collective_correlation`

    .. math::

        f_\mathrm{I} = \frac{\left( \sum_i \mathbf{r}_i \right)^2}{\sum l^2}

Average site occupations
    :attr:`~lattice_mc.results.Results.average_site_occupations`

    For each site label :math:`L`, the time-averaged number of atoms on sites with
    that label,

    .. math::

        \bar{n}_L = \frac{1}{t} \sum_{s \in L} t_s,

    where :math:`t_s` is the total time for which site :math:`s` has been occupied.

:attr:`~lattice_mc.results.Results.time` is :math:`t`, and
:attr:`~lattice_mc.results.Results.number_of_jumps` is the number of jumps made.

Diffusion coefficients and the reduced ionic conductivity are in units of (site
spacing)\ :sup:`2` s\ :sup:`-1`. :math:`d` is the dimensionality of the lattice, given
by :attr:`Lattice.dimensionality <lattice_mc.lattice.Lattice.dimensionality>`: the
number of non-zero cell lengths. Square and honeycomb lattices have a z cell length
of 0, so :math:`d = 2`. Jump vectors use the minimum-image convention, and have zero
components along axes with a cell length of 0.

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
