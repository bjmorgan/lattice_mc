Blocked lattices
================

It is possible to construct a simulation lattice where there are no possible jumps;
such a lattice is classified as "blocked". For example, a lattice where every site is
occupied is blocked.

If a simulation with a blocked lattice is run, a
:class:`~lattice_mc.error.BlockedLatticeError` exception is raised:

.. testcode::

    import lattice_mc
    from lattice_mc import init_lattice

    params = lattice_mc.SimulationParameters(temperature=298.0, rate_prefactor=1e13)
    simulation = lattice_mc.Simulation(params)
    simulation.lattice = init_lattice.square_lattice(4, 4, 1.0)
    simulation.set_number_of_atoms(16)
    simulation.set_number_of_jumps(10)
    simulation.run()

.. testoutput::

    Traceback (most recent call last):
      ...
    lattice_mc.error.BlockedLatticeError: No moves are possible in this lattice

To test for this, a :class:`~lattice_mc.lattice.Lattice` object can be queried using
the :meth:`~lattice_mc.lattice.Lattice.is_blocked` method:

.. testcode::

    print(simulation.lattice.is_blocked())

.. testoutput::

    True

:meth:`~lattice_mc.lattice.Lattice.is_blocked` needs the simulation parameters to be
attached to the lattice. :meth:`Simulation.run <lattice_mc.simulation.Simulation.run>`
and :meth:`Simulation.setup_lookup_table
<lattice_mc.simulation.Simulation.setup_lookup_table>` both attach them, so call
``is_blocked()`` after one of these. Before that, it raises a ``RuntimeError``.
