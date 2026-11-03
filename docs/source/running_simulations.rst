Running a simulation
====================

A :class:`~lattice_mc.simulation.Simulation` combines a lattice, a set of mobile
atoms, an interaction model, and the physical parameters for the jumps.

Simulation parameters
---------------------

:class:`~lattice_mc.simulation.SimulationParameters` holds the temperature, in K, and
the jump rate prefactor, :math:`\nu_0`, in s\ :sup:`-1`. Both must be finite and
positive. The thermal energy :math:`kT`, in eV, is derived from the temperature.

.. testcode::

    import lattice_mc as lmc

    params = lmc.SimulationParameters(temperature=298.0, rate_prefactor=1e13)
    print(f"{params.kT:.4f}")

.. testoutput::

    0.0257

Setting up
----------

A simulation is created from a lattice, the simulation parameters, and the atoms.
``n_atoms`` places that number of atoms on randomly chosen sites;
:attr:`~lattice_mc.simulation.Simulation.occupied_sites` gives the numbers of the
occupied sites.

.. testcode::

    lattice = lmc.square_lattice(4, 4, 1.0)
    simulation = lmc.Simulation(lattice, params, n_atoms=4, seed=1)
    print(len(simulation.occupied_sites))

.. testoutput::

    4

``occupied_sites`` instead places one atom on each of the listed sites, by site number:

.. testcode::

    placed = lmc.Simulation(lattice, params, occupied_sites=[1, 2, 3])
    print(placed.occupied_sites)

.. testoutput::

    (1, 2, 3)

Give exactly one of ``n_atoms`` and ``occupied_sites``. With ``n_atoms``,
``selected_labels`` restricts the random placement to sites with the given label or
labels. The atoms can then jump to any neighbouring vacant site, whatever its label.

.. testcode::

    honeycomb = lmc.honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
    on_a = lmc.Simulation(honeycomb, params, n_atoms=4, selected_labels="A", seed=1)
    print({honeycomb.site_with_id(number).label for number in on_a.occupied_sites})

.. testoutput::

    {'A'}

Invalid combinations of arguments, such as more atoms than sites, raise a
``ValueError`` when the simulation is created. The interaction model is also given
when the simulation is created (see :doc:`interaction_models`).

The simulation takes its own copy of the lattice, so later changes to the lattice
passed in have no effect on it.
:attr:`Simulation.lattice <lattice_mc.simulation.Simulation.lattice>` returns a new
copy of the simulation's lattice on each access, so changes to it have no effect on the
simulation.

Random numbers
--------------

``seed``, an integer or a ``numpy.random.Generator``, is the source of all
randomness in the simulation: the initial placement of atoms, and the choice and
timing of jumps. Two simulations created with the same arguments and the same integer
seed, and run in the same way, give the same results. A generator is used directly,
not copied, so the simulation advances it in place.
Without a seed, the generator is unseeded.

.. testcode::

    first = lmc.Simulation(lattice, params, n_atoms=4, seed=7)
    second = lmc.Simulation(lattice, params, n_atoms=4, seed=7)
    first.run(n_jumps=100)
    second.run(n_jumps=100)
    print(first.results == second.results)

.. testoutput::

    True

Running
-------

Until the simulation has made a jump, its results are ``None``:

.. testcode::

    print(simulation.results)

.. testoutput::

    None

:meth:`~lattice_mc.simulation.Simulation.run` makes a set number of jumps:

.. testcode::

    simulation.run(n_jumps=500)
    print(simulation.results.number_of_jumps)

.. testoutput::

    500

Instead of a number of jumps, ``for_time`` runs until the simulation time, in seconds,
reaches a set value:

.. testcode::

    timed = lmc.Simulation(lattice, params, n_atoms=4, seed=1)
    timed.run(for_time=1e-10)
    print(timed.results.time >= 1e-10)

.. testoutput::

    True

Give exactly one of ``n_jumps`` and ``for_time``. Calling
:meth:`~lattice_mc.simulation.Simulation.run` again continues from the current state,
and the time, displacements and other counters keep adding up across runs. With the
same integer seed, ``run(n_jumps=300)`` then ``run(n_jumps=700)`` gives the same
results as ``run(n_jumps=1000)``. ``for_time`` is the total time since the simulation
was created or last reset, not the time for that run.

:meth:`~lattice_mc.simulation.Simulation.reset` sets the time, the atom displacements
and the other counters back to zero, leaving the atoms where they are. The results are
then ``None`` until the next jump. To exclude a set of equilibration jumps from the
results, run them, reset, and then run again:

.. testcode::

    equilibrated = lmc.Simulation(honeycomb, params, n_atoms=4, selected_labels="A", seed=1)
    equilibrated.run(n_jumps=50)
    equilibrated.reset()
    print(equilibrated.results)
    equilibrated.run(n_jumps=100)
    print(equilibrated.results.number_of_jumps)

.. testoutput::

    None
    100

Results
-------

:attr:`~lattice_mc.simulation.Simulation.results` returns a new
:class:`~lattice_mc.results.Results`, a frozen dataclass, each time it is read. Its
fields are:

- :attr:`~lattice_mc.results.Results.tracer_diffusion_coefficient`, :math:`D^*`
- :attr:`~lattice_mc.results.Results.collective_diffusion_coefficient`,
  :math:`D_\mathrm{J}`
- :attr:`~lattice_mc.results.Results.reduced_ionic_conductivity`
- :attr:`~lattice_mc.results.Results.tracer_correlation`, :math:`f`
- :attr:`~lattice_mc.results.Results.collective_correlation`,
  :math:`f_\mathrm{I}`
- :attr:`~lattice_mc.results.Results.average_site_occupations`: the
  time-averaged number of atoms on the sites with each label, as a read-only mapping.
- :attr:`~lattice_mc.results.Results.time`: the simulation time, in seconds.
- :attr:`~lattice_mc.results.Results.number_of_jumps`: the number of jumps made.

These are defined in :doc:`introduction`.

.. testcode::

    results = simulation.results
    print(f"D*  = {results.tracer_diffusion_coefficient:.3e}")
    print(f"D_J = {results.collective_diffusion_coefficient:.3e}")
    print(f"f   = {results.tracer_correlation:.3f}")
    print(f"f_I = {results.collective_correlation:.3f}")
    print(f"reduced ionic conductivity = {results.reduced_ionic_conductivity:.3e}")
    print(sorted(results.average_site_occupations))

.. testoutput::

    D*  = ...
    D_J = ...
    f   = ...
    f_I = ...
    reduced ionic conductivity = ...
    ['L']
