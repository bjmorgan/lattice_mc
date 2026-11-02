Running a simulation
====================

A :class:`~lattice_mc.simulation.Simulation` combines a lattice, a number of mobile
atoms, and the physical parameters for the jumps.

Simulation parameters
---------------------

:class:`~lattice_mc.simulation.SimulationParameters` holds the temperature, in K, and
the jump rate prefactor, :math:`\nu_0`, in s\ :sup:`-1`. Both must be positive. The
thermal energy :math:`kT`, in eV, is derived from the temperature.

.. testcode::

    import lattice_mc
    from lattice_mc import init_lattice

    params = lattice_mc.SimulationParameters(temperature=298.0, rate_prefactor=1e13)
    print(f"{params.kT:.4f}")

.. testoutput::

    0.0257

Setting up and running
----------------------

Assign a lattice to the simulation, then add atoms and set the number of jumps.
:meth:`~lattice_mc.simulation.Simulation.set_number_of_atoms` places the atoms on
randomly chosen sites.

.. testcode::

    simulation = lattice_mc.Simulation(params)
    simulation.lattice = init_lattice.square_lattice(4, 4, 1.0)
    simulation.set_number_of_atoms(4)
    simulation.set_number_of_jumps(500)
    print(simulation.has_run)

.. testoutput::

    False

Until the simulation has run, its results are ``None``:

.. testcode::

    print(simulation.tracer_diffusion_coefficient)

.. testoutput::

    None

:meth:`~lattice_mc.simulation.Simulation.run` performs the jumps:

.. testcode::

    simulation.run()
    print(simulation.has_run)

.. testoutput::

    True

Instead of a number of jumps, a simulation can run until a set simulation time, in
seconds, has passed:

.. testcode::

    timed = lattice_mc.Simulation(params)
    timed.lattice = init_lattice.square_lattice(4, 4, 1.0)
    timed.set_number_of_atoms(4)
    timed.run(for_time=1e-10)
    print(timed.lattice.time >= 1e-10)

.. testoutput::

    True

Atoms can be restricted to sites with particular labels, and a number of
equilibration jumps can be run first. Equilibration jumps let the atoms move away
from their random starting positions, and are not included in the results.

.. testcode::

    equilibrated = lattice_mc.Simulation(params)
    equilibrated.lattice = init_lattice.honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
    equilibrated.set_number_of_atoms(4, selected_sites=["A"])
    equilibrated.set_number_of_jumps(100)
    equilibrated.set_number_of_equilibration_jumps(50)
    equilibrated.run()

:meth:`~lattice_mc.simulation.Simulation.reset` sets the time and the atom
displacements back to zero, leaving the atoms where they are. Call
:meth:`~lattice_mc.simulation.Simulation.run` again before reading the results: after a
reset, the diffusion coefficients and correlation factors cannot be calculated until
time has passed. A second run continues from the current atom positions.

Results
-------

After a run, the simulation reports:

- :attr:`~lattice_mc.simulation.Simulation.tracer_diffusion_coefficient`, :math:`D^*`
- :attr:`~lattice_mc.simulation.Simulation.collective_diffusion_coefficient`,
  :math:`D_\mathrm{J}`, and
  :attr:`~lattice_mc.simulation.Simulation.collective_diffusion_coefficient_per_atom`
- :attr:`~lattice_mc.simulation.Simulation.tracer_correlation`, :math:`f`
- :attr:`~lattice_mc.simulation.Simulation.collective_correlation`,
  :math:`f_\mathrm{I}`
- :attr:`~lattice_mc.simulation.Simulation.average_site_occupations`: the
  time-averaged number of atoms on the sites with each label.

These quantities are defined in :doc:`introduction`.

.. testcode::

    print(f"D*  = {simulation.tracer_diffusion_coefficient:.3e}")
    print(f"D_J = {simulation.collective_diffusion_coefficient:.3e}")
    print(f"f   = {simulation.tracer_correlation:.3f}")
    print(f"f_I = {simulation.collective_correlation:.3f}")
    print(sorted(simulation.average_site_occupations))

.. testoutput::

    D*  = ...
    D_J = ...
    f   = ...
    f_I = ...
    ['L']

The results of a single run are noisy. Converged values come from averaging over
many independent simulations, and from longer runs on larger lattices than these
examples use.
