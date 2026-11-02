lattice_mc: A Python Lattice-Gas Monte Carlo Module
===================================================

``lattice_mc`` runs kinetic lattice-gas Monte Carlo simulations of atoms jumping
between sites on a lattice, and calculates their transport coefficients.

Capabilities
------------

- **Lattices:** square, honeycomb and cubic lattices built programmatically, or
  lattices of any geometry read from a sites file.
- **Interaction models:** volume exclusion, nearest-neighbour interactions,
  coordination-number dependent energies, and site energies, in any combination.
- **Results:** tracer and collective diffusion coefficients, tracer and collective
  correlation factors, and average site occupations.
- **Runs** of a set number of jumps or a set simulation time, with optional
  equilibration jumps.

Quickstart
----------

A simulation of 50 atoms on a 10 x 10 square lattice:

.. testcode::

    import lattice_mc
    from lattice_mc import init_lattice

    params = lattice_mc.SimulationParameters(temperature=298.0, rate_prefactor=1e13)
    simulation = lattice_mc.Simulation(params)
    simulation.lattice = init_lattice.square_lattice(10, 10, 1.0)
    simulation.set_number_of_atoms(50)
    simulation.set_number_of_jumps(1000)
    simulation.run()

    print(f"D* = {simulation.tracer_diffusion_coefficient:.3e}")
    print(f"f  = {simulation.tracer_correlation:.3f}")

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
   modules
