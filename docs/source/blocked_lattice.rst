Blocked lattices
================

It is possible to construct a simulation where there are no possible jumps; such a
simulation is classified as "blocked". For example, a simulation where every site is
occupied is blocked.

:meth:`Simulation.is_blocked <lattice_mc.simulation.Simulation.is_blocked>` checks
for this before running. It returns ``True`` if no occupied site has a vacant
neighbour:

.. testcode::

    import lattice_mc as lmc

    params = lmc.SimulationParameters(temperature=298.0, rate_prefactor=1e13)
    lattice = lmc.square_lattice(4, 4, 1.0)
    simulation = lmc.Simulation(lattice, params, n_atoms=16)
    print(simulation.is_blocked())

.. testoutput::

    True

If a blocked simulation is run, a :class:`~lattice_mc.error.BlockedLatticeError`
exception is raised:

.. testcode::

    simulation.run(n_jumps=10)

.. testoutput::

    Traceback (most recent call last):
      ...
    lattice_mc.error.BlockedLatticeError: No jumps are possible.

:meth:`~lattice_mc.simulation.Simulation.run` also raises
:class:`~lattice_mc.error.BlockedLatticeError` if every possible jump has a
probability of zero, which ``is_blocked`` does not check. Jumps made in a run before
the error is raised are kept.
