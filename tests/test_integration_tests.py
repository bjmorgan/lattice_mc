import numpy as np

import lattice_mc

PARAMS = lattice_mc.SimulationParameters(temperature=298.0, rate_prefactor=1e13)


def three_site_chain(positions):
    """A chain of three 'A' sites, 1-2-3, at the given positions, in a large cell."""
    neighbours = [[2], [1, 3], [2]]
    sites = [
        lattice_mc.Site(number, np.array(r), neighbours[number - 1], 0.0, "A")
        for number, r in enumerate(positions, start=1)
    ]
    return lattice_mc.Lattice(sites, cell_lengths=np.array([10.0, 10.0, 10.0]))


class TestIntegration:
    def test_simulation_runs_with_selected_labels(self):
        lattice = lattice_mc.cubic_lattice(4, 4, 4, 1.0)
        lattice.transmute_sites("L", "X", 32, seed=0)
        s = lattice_mc.Simulation(lattice, PARAMS, n_atoms=32, selected_labels="L", lookup_table=True, seed=0)
        labels = {s.lattice.site_with_id(number).label for number in s.occupied_sites}
        assert labels == {"L"}
        assert s.results is None
        s.run(n_jumps=10)
        assert s.results is not None

    def test_simulation_runs_with_cn_energies(self):
        lattice = three_site_chain([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        s = lattice_mc.Simulation(
            lattice, PARAMS, n_atoms=1, cn_energies={"A": {"A": {0: 0.0, 1: -0.1, 2: -0.3}}}, seed=0
        )
        s.run(n_jumps=10)
        assert s.results is not None

    def test_simulation_runs_with_variable_coordination_numbers(self):
        lattice = three_site_chain([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 1.0, 1.0]])
        s = lattice_mc.Simulation(lattice, PARAMS, n_atoms=1, lookup_table=True, seed=0)
        s.run(n_jumps=10)
        assert s.results is not None
