Migrating from 1.x
==================

This page maps each 1.x call that has been removed or renamed in 2.0 to its
replacement. ``simulation`` is a :class:`~lattice_mc.simulation.Simulation`,
``lattice`` a :class:`~lattice_mc.lattice.Lattice`, and ``results`` the
:class:`~lattice_mc.results.Results` from ``simulation.results``.

Setting up a simulation
-----------------------

A simulation is now created from its lattice, parameters, atoms and interaction model
in one call, in place of the 1.x setters:

.. testcode::

    import lattice_mc as lmc

    lattice = lmc.square_lattice(4, 4, 1.0)
    params = lmc.SimulationParameters(temperature=298.0, rate_prefactor=1e13)
    simulation = lmc.Simulation(
        lattice,
        params,
        n_atoms=4,
        nn_energy=0.1,
        lookup_table=True,
        seed=1,
    )
    simulation.run(n_jumps=100)  # equilibration
    simulation.reset()
    simulation.run(n_jumps=500)
    results = simulation.results
    print(results.number_of_jumps)

.. testoutput::

    500

.. list-table::
    :header-rows: 1
    :widths: 50 50

    * - 1.x
      - 2.0
    * - ``Simulation(params)``, then ``simulation.lattice = lattice``
      - ``Simulation(lattice, params, ...)``
    * - ``simulation.define_lattice_from_file(filename, cell_lengths)``
      - ``lattice_from_sites_file(filename, cell_lengths)``, passed to ``Simulation``
    * - ``simulation.set_number_of_atoms(n)``
      - ``Simulation(..., n_atoms=n)``, or ``occupied_sites=[...]`` to place the
        atoms on given sites
    * - ``simulation.set_number_of_atoms(n, selected_sites=labels)``
      - ``Simulation(..., n_atoms=n, selected_labels=labels)``
    * - ``simulation.set_nn_energy(energy)``
      - ``Simulation(..., nn_energy=energy)``
    * - ``simulation.set_cn_energies(cn_energies)``
      - ``Simulation(..., cn_energies=cn_energies)``
    * - ``simulation.set_site_energies(energies)``
      - ``lattice.set_site_energies(by_label=energies)``, before creating the
        simulation
    * - ``simulation.setup_lookup_table()``
      - ``Simulation(..., lookup_table=True)``
    * - ``simulation.is_initialised()``
      - Removed: invalid arguments raise ``ValueError`` when the simulation is created
    * - ``simulation.number_of_atoms``
      - ``len(simulation.occupied_sites)``
    * - ``simulation.atoms``, ``simulation.for_time``,
        ``simulation.number_of_equilibration_jumps``
      - Removed

``seed`` (an integer or a ``numpy.random.Generator``) controls the initial
placement of atoms and the choice and timing of jumps. The simulation takes its own
copy of the lattice when it is created.

Running
-------

.. list-table::
    :header-rows: 1
    :widths: 50 50

    * - 1.x
      - 2.0
    * - ``simulation.set_number_of_jumps(n)``, then ``simulation.run()``
      - ``simulation.run(n_jumps=n)``
    * - ``simulation.run(for_time=t)``
      - ``simulation.run(for_time=t)``; the arguments of ``run`` are keyword-only
    * - ``simulation.set_number_of_equilibration_jumps(m)``, then ``simulation.run()``
      - ``simulation.run(n_jumps=m)``, ``simulation.reset()``,
        ``simulation.run(n_jumps=n)``
    * - ``simulation.has_run``
      - ``simulation.results is not None``; ``results`` is ``None`` from creation
        or a reset until the next jump
    * - ``simulation.lattice.is_blocked()``
      - ``simulation.is_blocked()``

Repeated calls to ``run`` accumulate: ``run(n_jumps=300)`` then ``run(n_jumps=700)``
gives the same results as ``run(n_jumps=1000)`` with the same seed.

Results
-------

.. list-table::
    :header-rows: 1
    :widths: 50 50

    * - 1.x
      - 2.0
    * - ``simulation.tracer_diffusion_coefficient``
      - ``results.tracer_diffusion_coefficient``
    * - ``simulation.collective_diffusion_coefficient``
      - ``results.collective_diffusion_coefficient * len(simulation.occupied_sites)``
    * - ``simulation.collective_diffusion_coefficient_per_atom``
      - ``results.collective_diffusion_coefficient``
    * - ``simulation.tracer_correlation``
      - ``results.tracer_correlation``
    * - ``simulation.collective_correlation``
      - ``results.collective_correlation``
    * - ``simulation.average_site_occupations``
      - ``results.average_site_occupations``
    * - ``simulation.number_of_jumps``
      - ``results.number_of_jumps``
    * - ``simulation.lattice.time``
      - ``results.time``

``collective_diffusion_coefficient`` is now per atom,
:math:`\left( \sum_i \mathbf{r}_i \right)^2 / (2 d N t)`, and
``collective_diffusion_coefficient_per_atom`` is removed.
``reduced_ionic_conductivity`` is new (see :doc:`introduction`).

Jump vectors are now zero along axes with a cell length of 0. In 1.x, those components
were left unchanged.

Lattices, sites and clusters
----------------------------

.. list-table::
    :header-rows: 1
    :widths: 50 50

    * - 1.x
      - 2.0
    * - ``lattice.occupied_sites()``, ``lattice.occupied_site_numbers()``,
        ``site.is_occupied``
      - ``simulation.occupied_sites``
    * - ``lattice.set_site_energies(energies)``
      - ``lattice.set_site_energies(by_label=energies)``; per-site energies with
        ``by_site={number: energy}``. Both arguments are keyword-only.
    * - ``lattice.transmute_sites(old_label, new_label, n)``
      - ``lattice.transmute_sites(old_label, new_label, n, seed=seed)``; ``seed``
        (an integer or a ``numpy.random.Generator``) is optional and
        keyword-only
    * - ``site.p_neighbours``
      - ``lattice.neighbours_of(site)``
    * - ``Cluster(sites)``
      - ``Cluster(sites, lattice)``
    * - ``cluster.is_periodically_contiguous(cell_lengths)``
      - ``cluster.is_periodically_contiguous()``

:attr:`Lattice.site_labels <lattice_mc.lattice.Lattice.site_labels>` and
:attr:`Lattice.site_populations <lattice_mc.lattice.Lattice.site_populations>` are
read-only properties, computed from the sites when read.

:class:`~lattice_mc.lattice.Lattice` and :class:`~lattice_mc.lattice_site.Site` no
longer hold simulation state. These members are removed:

- ``Lattice``: ``params``, ``nn_energy``, ``cn_energies``, ``jump_lookup_table``,
  ``number_of_occupied_sites``, ``time``, ``site_lookup``, ``reset()``,
  ``initialise_site_lookup_table()``, ``vacant_sites()``, ``occupied_sites()``,
  ``vacant_site_numbers()``, ``occupied_site_numbers()``, ``potential_jumps()``,
  ``update()``, ``populate_sites()``, ``jump()``, ``update_site_occupation_times()``,
  ``site_occupation_statistics()``, ``set_nn_energy()``, ``set_cn_energies()`` and
  ``is_blocked()``.
- ``Site``: the ``cn_energies`` constructor argument, ``index``, ``p_neighbours``,
  ``occupation``, ``atom``, ``is_occupied``, ``time_occupied``,
  ``cn_occupation_energies``, ``nn_occupation()``, ``site_specific_nn_occupation()``,
  ``site_specific_neighbours()``, ``set_cn_occupation_energies()`` and
  ``cn_occupation_energy()``.

Modules
-------

The ``lattice_mc.atom``, ``lattice_mc.species``, ``lattice_mc.jump``,
``lattice_mc.transitions`` and ``lattice_mc.lookup_table`` modules are removed.
:class:`~lattice_mc.simulation.Simulation`,
:class:`~lattice_mc.simulation.SimulationParameters`,
:class:`~lattice_mc.results.Results`, :class:`~lattice_mc.lattice.Lattice`,
:class:`~lattice_mc.lattice_site.Site`, :class:`~lattice_mc.error.BlockedLatticeError`
and the lattice builders ``square_lattice``, ``honeycomb_lattice``, ``cubic_lattice``
and ``lattice_from_sites_file`` can be imported from the top-level package.
