lattice_mc: A Python Lattice-Gas Monte Carlo Module
===================================================

``lattice_mc`` runs lattice-gas Monte Carlo simulations.

Capabilities
------------

- **Lattices:** square, honeycomb and cubic lattices built programmatically, or
  lattices of any geometry read from a sites file.
- **Interaction models:** volume exclusion, nearest-neighbour interactions,
  coordination-number dependent energies, and site energies, in any combination.
- **Results:** tracer and collective diffusion coefficients, the reduced ionic
  conductivity, tracer and collective correlation factors, and average site
  occupations.
- **Runs** of a set number of jumps or a set simulation time, continued across
  calls, with counters that can be reset, and a seed for reproducible runs.

Quickstart
----------

A simulation of 16 atoms on a 4 x 4 x 4 cubic lattice, with nearest-neighbour
repulsion:

.. testcode::

    import lattice_mc as lmc

    lattice = lmc.cubic_lattice(4, 4, 4, spacing=1.0)

    simulation = lmc.Simulation(
        lattice,
        params=lmc.SimulationParameters(temperature=298.0, rate_prefactor=1e13),
        n_atoms=16,
        nn_energy=0.1,
        lookup_table=True,
        seed=42,
    )
    simulation.run(n_jumps=1_000)  # equilibration
    simulation.reset()
    simulation.run(n_jumps=10_000)
    results = simulation.results

    print(f"D* = {results.tracer_diffusion_coefficient:.3e}")
    print(f"f  = {results.tracer_correlation:.3f}")

.. testoutput::

    D* = ...
    f  = ...

.. toctree::
   :hidden:

   introduction
   building_lattices
   running_simulations
   interaction_models
   blocked_lattice
   migrating
   api
