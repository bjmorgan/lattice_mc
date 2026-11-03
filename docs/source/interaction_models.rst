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
combined, and all energies are in eV.

The examples on this page share this set-up: a honeycomb lattice with alternating
``'A'`` and ``'B'`` sites, so that every ``'A'`` site has three ``'B'`` neighbours,
and vice versa.

.. testcode::

    import lattice_mc
    from lattice_mc import init_lattice

    params = lattice_mc.SimulationParameters(temperature=298.0, rate_prefactor=1e13)

    def new_simulation():
        simulation = lattice_mc.Simulation(params)
        simulation.lattice = init_lattice.honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
        simulation.set_number_of_atoms(4)
        simulation.set_number_of_jumps(200)
        return simulation

Non-interacting atoms
---------------------

By default, the only interaction is volume exclusion: two atoms cannot occupy the
same site. Unless the sites have occupation energies (see `Site energies`_), every
allowed jump has :math:`\Delta E = 0`.

.. testcode::

    simulation = new_simulation()
    simulation.run()
    print(simulation.has_run)

.. testoutput::

    True

Nearest-neighbour interactions
------------------------------

:meth:`~lattice_mc.simulation.Simulation.set_nn_energy` sets an interaction energy for
each pair of occupied neighbouring sites. A positive energy is a repulsion.

.. testcode::

    simulation = new_simulation()
    simulation.set_nn_energy(0.1)
    simulation.setup_lookup_table()
    simulation.run()
    print(simulation.has_run)

.. testoutput::

    True

:meth:`~lattice_mc.simulation.Simulation.setup_lookup_table` is optional. It
precalculates the probability of every possible jump, which speeds up simulations
with nearest-neighbour interactions. Without a table, probabilities are calculated
as they are needed. The table includes the site energies, whether set with
:meth:`~lattice_mc.simulation.Simulation.set_site_energies` or read from a sites file,
and the nearest-neighbour energy, as they are when the table is created, so set those
first. A lookup table supports
nearest-neighbour interactions only, so do not use one with coordination-number
energies.

Coordination-number energies
----------------------------

:meth:`~lattice_mc.simulation.Simulation.set_cn_energies` gives each occupied site
an energy that depends on how many of its neighbours are occupied; vacant sites
contribute no energy. The energies are nested by site label, then by the label of the
neighbouring sites, then by the number of those neighbours that are occupied. The
energy of an occupied site is the sum of its energies for each neighbouring label, and
must be given for every possible number of occupied neighbours.

.. testcode::

    simulation = new_simulation()
    cn = {0: 0.0, 1: 0.1, 2: 0.2, 3: 0.3}
    simulation.set_cn_energies({"A": {"B": cn}, "B": {"A": cn}})
    simulation.run()
    print(simulation.has_run)

.. testoutput::

    True

Lookup tables do not support coordination-number energies.

Site energies
-------------

:meth:`~lattice_mc.simulation.Simulation.set_site_energies` gives every site with a
given label an occupation energy. Here ``'B'`` sites have an occupation energy
0.1 eV higher than ``'A'`` sites:

.. testcode::

    simulation = new_simulation()
    simulation.set_site_energies({"A": 0.0, "B": 0.1})
    simulation.run()
    print(simulation.has_run)

.. testoutput::

    True

Site energies can also be given per site in a sites file (see
:doc:`building_lattices`). ``set_site_energies`` overwrites any energies read from a
sites file for the labels it is given.
