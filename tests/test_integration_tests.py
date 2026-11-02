import numpy as np

import lattice_mc

PARAMS = lattice_mc.SimulationParameters(temperature=298.0, rate_prefactor=1e13)


class TestIntegration:
    def test_simulation_runs_with_selected_sites(self):
        a, b, c = 4, 4, 4
        spacing = 1.0
        n_atoms = 32
        selected_sites = "L"
        number_of_jumps = 10
        s = lattice_mc.Simulation(PARAMS)
        s.lattice = lattice_mc.init_lattice.cubic_lattice(a, b, c, spacing)
        s.lattice.transmute_sites("L", "X", n_atoms)
        s.set_number_of_atoms(n_atoms, selected_sites=selected_sites)
        for site in s.lattice.occupied_sites():
            assert site.label == selected_sites
        s.setup_lookup_table()
        s.set_number_of_jumps(number_of_jumps)
        assert s.tracer_correlation is None
        assert s.collective_correlation is None
        assert s.tracer_diffusion_coefficient is None
        assert s.collective_diffusion_coefficient_per_atom is None
        assert s.average_site_occupations is None
        s.run()
        assert s.tracer_correlation is not None
        assert s.collective_correlation is not None
        assert s.tracer_diffusion_coefficient is not None
        assert s.collective_diffusion_coefficient_per_atom is not None
        assert s.average_site_occupations is not None

    def test_simulation_runs_with_cn_energies(self):
        s = lattice_mc.Simulation(PARAMS)
        site_data = [
            [1, np.array([0.0, 0.0, 0.0]), [2], 0.0, "A"],
            [2, np.array([1.0, 0.0, 0.0]), [1, 3], 0.0, "A"],
            [3, np.array([2.0, 0.0, 0.0]), [2], 0.0, "A"],
        ]
        sites = [lattice_mc.lattice_site.Site(*d) for d in site_data]
        s.lattice = lattice_mc.lattice.Lattice(sites, cell_lengths=np.array([10.0, 10.0, 10.0]))
        s.set_cn_energies({"A": {"A": {0: 0.0, 1: -0.1, 2: -0.3}}})
        s.set_number_of_atoms(1)
        s.set_number_of_jumps(10)
        s.run()
        assert s.tracer_correlation is not None

    def test_simulation_runs_with_variable_coordination_numbers(self):
        s = lattice_mc.Simulation(PARAMS)
        site_data = [
            [1, np.array([0.0, 0.0, 0.0]), [2], 0.0, "A"],
            [2, np.array([1.0, 0.0, 0.0]), [1, 3], 0.0, "A"],
            [3, np.array([2.0, 1.0, 1.0]), [2], 0.0, "A"],
        ]
        sites = [lattice_mc.lattice_site.Site(*d) for d in site_data]
        s.lattice = lattice_mc.lattice.Lattice(sites, cell_lengths=np.array([10.0, 10.0, 10.0]))
        s.set_number_of_atoms(1)
        s.setup_lookup_table()
        s.set_number_of_jumps(10)
        s.run()
