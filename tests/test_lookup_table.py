import math
from unittest.mock import Mock, call, patch

import pytest

from lattice_mc.init_lattice import lattice_from_sites_file
from lattice_mc.lattice import Lattice
from lattice_mc.lookup_table import LookupTable, metropolis
from lattice_mc.simulation import SimulationParameters

PARAMS = SimulationParameters(temperature=298.0, rate_prefactor=1e13)


class TestLookupTableSupportFunctions:
    def test_metropolis_returns_one_for_negative_delta_E(self):
        assert metropolis(-0.02, PARAMS.kT) == 1.0

    def test_metropolis_returns_boltzmann_factor_for_positive_delta_E(self):
        expected = math.exp(-0.02 / PARAMS.kT)
        assert metropolis(+0.02, PARAMS.kT) == pytest.approx(expected, abs=5e-8)


class TestLookupTable:
    """Tests for LookupTable class"""

    @pytest.fixture(autouse=True)
    @patch("lattice_mc.lookup_table.LookupTable.generate_nearest_neighbour_lookup_table")
    def set_up(self, mock_generate_nearest_neighbour_lookup_table):
        self.lattice = Mock(spec=Lattice)
        self.lattice.params = PARAMS
        self.lattice.site_energies_by_label = Mock(return_value="foo")
        self.lattice.nn_energy = "bar"
        self.lattice.cn_energies = "baz"
        self.lattice.connected_site_pairs = Mock(return_value="qux")
        self.lattice.max_site_coordination_numbers = Mock(return_value={"A": 2})
        self.lattice.site_specific_coordination_numbers = Mock(return_value="quux")
        hamiltonian = "nearest-neighbour"
        self.table = LookupTable(self.lattice, hamiltonian)

    @patch("lattice_mc.lookup_table.LookupTable.generate_nearest_neighbour_lookup_table")
    def test_nearest_neighbour_lookup_table_is_initialised(self, mock_generate_nearest_neighbour_lookup_table):
        hamiltonian = "nearest-neighbour"
        lookup_table = LookupTable(self.lattice, hamiltonian)
        assert lookup_table.site_energies == "foo"
        assert lookup_table.nn_energy == "bar"
        assert lookup_table.cn_energy == "baz"
        assert lookup_table.connected_site_pairs == "qux"
        assert lookup_table.site_specific_coordination_per_site == "quux"
        assert mock_generate_nearest_neighbour_lookup_table.call_count == 1

    def test_lookup_table_init_with_invalid_hamiltonian_keywork_raises_ValueError(self):
        hamiltonian = "foo"
        with pytest.raises(ValueError):
            LookupTable(self.lattice, hamiltonian)

    @patch("lattice_mc.lookup_table.metropolis")
    def test_relative_probability_with_site_energies(self, mock_metropolis):
        self.table.site_energies = {"A": 1.0, "B": 2.0}
        self.table.nn_energy = None
        self.table.relative_probability("A", "B", 3, 1)
        mock_metropolis.assert_called_with(1.0, PARAMS.kT)

    @patch("lattice_mc.lookup_table.metropolis")
    def test_relative_probability_with_coordination_numbers(self, mock_metropolis):
        self.table.site_energies = {"A": 0.0, "B": 0.0}
        self.table.nn_energy = 1.0
        self.table.relative_probability("A", "B", 3, 1)
        mock_metropolis.assert_called_with(-3.0, PARAMS.kT)

    @patch("lattice_mc.lookup_table.metropolis")
    def test_relative_probability_with_site_energies_and_coordination_numbers(self, mock_metropolis):
        self.table.site_energies = {"A": 1.0, "B": 2.0}
        self.table.nn_energy = 1.0
        self.table.relative_probability("A", "B", 3, 1)
        mock_metropolis.assert_called_with(-2.0, PARAMS.kT)

    def test_generate_nearest_neighbour_lookup_table(self):
        self.table.connected_site_pairs = {"A": ["B"]}
        self.table.max_coordination_per_site = {"A": 2, "B": 3}
        self.table.relative_probability = Mock(side_effect=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        self.table.generate_nearest_neighbour_lookup_table()
        assert self.table.jump_probability == {"A": {"B": {0: {1: 1.0, 2: 2.0, 3: 3.0}, 1: {1: 4.0, 2: 5.0, 3: 6.0}}}}
        assert self.table.relative_probability.mock_calls == [
            call("A", "B", 0, 1),
            call("A", "B", 0, 2),
            call("A", "B", 0, 3),
            call("A", "B", 1, 1),
            call("A", "B", 1, 2),
            call("A", "B", 1, 3),
        ]


SITES_WITH_ENERGIES = """2

site: 1
centre: 0.0 0.0 0.0
neighbours: 2
label: A

site: 2
centre: 1.0 0.0 0.0
neighbours: 1
label: B
energy: 0.1
"""


def test_lookup_table_includes_site_energies_from_a_sites_file(tmp_path):
    path = tmp_path / "sites.dat"
    path.write_text(SITES_WITH_ENERGIES)
    lattice = lattice_from_sites_file(str(path), cell_lengths=[10.0, 10.0, 10.0])
    lattice.params = PARAMS
    table = LookupTable(lattice, "nearest-neighbour")
    assert table.relative_probability("A", "B", 0, 1) == pytest.approx(metropolis(0.1, PARAMS.kT))
