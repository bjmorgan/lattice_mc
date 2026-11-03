Interaction models
==================

The rate of each jump is :math:`\nu_0 P`, where :math:`\nu_0` is the rate prefactor
and :math:`P` is the Metropolis probability for the change in energy,
:math:`\Delta E`, that the jump would cause:

.. math::

    P = \begin{cases}
        1 & \Delta E \le 0 \\
        \exp(-\Delta E / kT) & \Delta E > 0
    \end{cases}

The interaction model defines :math:`\Delta E`. The energy terms below can be
combined, and all energies are in eV. Site energies are set on the lattice; the
nearest-neighbour and coordination-number energies are given when the
:class:`~lattice_mc.simulation.Simulation` is created.

The examples on this page share this set-up: a honeycomb lattice with alternating
``'A'`` and ``'B'`` sites, so that every ``'A'`` site has three ``'B'`` neighbours,
and vice versa.

.. testcode::

    import lattice_mc as lmc

    params = lmc.SimulationParameters(temperature=298.0, rate_prefactor=1e13)
    lattice = lmc.honeycomb_lattice(2, 2, 1.0, alternating_sites=True)

    def new_simulation(lattice, **interactions):
        return lmc.Simulation(lattice, params, n_atoms=4, seed=1, **interactions)

Non-interacting atoms
---------------------

By default, the only interaction is volume exclusion: two atoms cannot occupy the
same site. Unless the sites have occupation energies (see `Site energies`_), every
allowed jump has :math:`\Delta E = 0`.

.. testcode::

    simulation = new_simulation(lattice)
    simulation.run(n_jumps=200)
    print(simulation.results.number_of_jumps)

.. testoutput::

    200

Nearest-neighbour interactions
------------------------------

``nn_energy`` sets an interaction energy for each pair of occupied neighbouring
sites. A positive energy is a repulsion. For a jump from site :math:`i` to site
:math:`j`, with :math:`n_i` and :math:`n_j` occupied neighbours before the jump, this
adds :math:`E_\mathrm{nn} (n_j - n_i - 1)` to :math:`\Delta E`.

.. testcode::

    simulation = new_simulation(lattice, nn_energy=0.1, lookup_table=True)
    simulation.run(n_jumps=200)
    print(simulation.results.number_of_jumps)

.. testoutput::

    200

``lookup_table=True`` is optional. It precalculates the probability of every possible
jump, keyed by site label and numbers of occupied neighbours, when the simulation is
created. Without a table, probabilities are calculated as they are needed. The table
includes the site energies and the nearest-neighbour energy.

Because the table is keyed by site label, every site with a given label must have the
same energy; otherwise creating the simulation with ``lookup_table=True`` raises a
``ValueError``. A lookup table supports site energies and nearest-neighbour
interactions only, and giving ``cn_energies`` with ``lookup_table=True`` also raises
a ``ValueError``.

Coordination-number energies
----------------------------

``cn_energies`` gives each occupied site an energy that depends on how many of its
neighbours are occupied; vacant sites contribute no energy. The energies are nested by
site label, then by the label of the neighbouring sites, then by the number of those
neighbours that are occupied. The energy of an occupied site is the sum of its
energies for each neighbouring label. The mapping must include every site label,
every label of a neighbouring site, and every number of occupied neighbours from 0 to
the number of those neighbours; otherwise creating the simulation raises a
``ValueError``.

.. testcode::

    cn = {0: 0.0, 1: 0.1, 2: 0.2, 3: 0.3}
    simulation = new_simulation(lattice, cn_energies={"A": {"B": cn}, "B": {"A": cn}})
    simulation.run(n_jumps=200)
    print(simulation.results.number_of_jumps)

.. testoutput::

    200

Site energies
-------------

:meth:`~lattice_mc.lattice.Lattice.set_site_energies` gives sites an occupation energy,
by label or by site number (see :doc:`building_lattices`). A jump from site :math:`i`
to site :math:`j` adds :math:`E_j - E_i` to :math:`\Delta E`. Here ``'B'`` sites have
an occupation energy 0.1 eV higher than ``'A'`` sites:

.. testcode::

    lattice_with_energies = lmc.honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
    lattice_with_energies.set_site_energies(by_label={"A": 0.0, "B": 0.1})
    simulation = new_simulation(lattice_with_energies)
    simulation.run(n_jumps=200)
    print(simulation.results.number_of_jumps)

.. testoutput::

    200

Site energies can also be given per site in a sites file (see
:doc:`building_lattices`). The simulation reads the site energies from its copy of the
lattice when it is created, so changing the energies of the lattice afterwards has no
effect on it.
