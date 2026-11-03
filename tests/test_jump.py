import math
from unittest.mock import Mock, patch

import numpy as np
import pytest

from lattice_mc.jump import Jump
from lattice_mc.lattice_site import Site
from lattice_mc.lookup_table import LookupTable
from lattice_mc.simulation import SimulationParameters

PARAMS = SimulationParameters(temperature=298.0, rate_prefactor=1e13)


class TestJump:
    """Test for Jump class"""

    @pytest.fixture(autouse=True)
    def set_up(self):
        self.mock_initial_site = Mock(spec=Site)
        self.mock_final_site = Mock(spec=Site)
        with patch("lattice_mc.jump.Jump.boltzmann_factor") as self.mock_boltzmann_factor:
            self.mock_boltzmann_factor.return_value = 0.13
            self.jump = Jump(self.mock_initial_site, self.mock_final_site, params=PARAMS)

    def test_jump_is_initialised_with_no_interactions(self):
        assert self.jump.initial_site == self.mock_initial_site
        assert self.jump.final_site == self.mock_final_site
        assert self.jump.nearest_neighbour_energy is None
        assert self.jump.coordination_number_energy is None
        assert self.jump.relative_probability == self.mock_boltzmann_factor.return_value

    def test_jump_is_initialised_with_nn_energy(self):
        nearest_neighbour_energy = "foo"
        with patch("lattice_mc.jump.Jump.boltzmann_factor") as mock_boltzmann_factor:
            mock_boltzmann_factor.return_value = 0.13
            jump = Jump(
                self.mock_initial_site,
                self.mock_final_site,
                nearest_neighbour_energy=nearest_neighbour_energy,
                params=PARAMS,
            )
        assert jump.initial_site == self.mock_initial_site
        assert jump.final_site == self.mock_final_site
        assert jump.nearest_neighbour_energy == nearest_neighbour_energy
        assert jump.coordination_number_energy is None
        assert jump.relative_probability == 0.13

    def test_jump_is_initialised_with_cn_energy(self):
        coordination_number_energy = "foo"
        with patch("lattice_mc.jump.Jump.boltzmann_factor") as mock_boltzmann_factor:
            mock_boltzmann_factor.return_value = 0.13
            jump = Jump(
                self.mock_initial_site,
                self.mock_final_site,
                coordination_number_energy=coordination_number_energy,
                params=PARAMS,
            )
        assert jump.initial_site == self.mock_initial_site
        assert jump.final_site == self.mock_final_site
        assert jump.nearest_neighbour_energy is None
        assert jump.coordination_number_energy == coordination_number_energy
        assert jump.relative_probability == 0.13

    def test_jump_is_initialised_with_lookup_table(self):
        mock_initial_site = Mock(spec=Site)
        mock_final_site = Mock(spec=Site)
        jump_lookup_table = "foo"
        with patch("lattice_mc.jump.Jump.relative_probability_from_lookup_table") as mock_p_from_lookup_table:
            mock_p_from_lookup_table.return_value = 0.73
            jump = Jump(mock_initial_site, mock_final_site, jump_lookup_table=jump_lookup_table, params=PARAMS)
        assert jump.initial_site == mock_initial_site
        assert jump.final_site == mock_final_site
        assert jump.nearest_neighbour_energy is None
        assert jump.coordination_number_energy is None
        assert jump.relative_probability == 0.73

    def test_rate(self):
        params = SimulationParameters(temperature=298.0, rate_prefactor=1e-3)
        with patch("lattice_mc.jump.Jump.boltzmann_factor") as mock_bf:
            mock_bf.return_value = 0.5
            j = Jump(self.mock_initial_site, self.mock_final_site, params=params)
        assert j.rate() == 5e-4

    def test_boltzmann_factor_lt_0(self):
        self.jump.delta_E = Mock(return_value=-2.5)
        assert self.jump.boltzmann_factor() == 1.0

    def test_boltzmann_factor_gt_0(self):
        self.jump.delta_E = Mock(return_value=0.2)
        expected = math.exp(-0.2 / PARAMS.kT)
        assert self.jump.boltzmann_factor() == pytest.approx(expected, abs=5e-8)

    def test_delta_E_non_interacting(self):
        initial_energy = 1.5
        final_energy = 2.1
        self.jump.initial_site.energy = initial_energy
        self.jump.final_site.energy = final_energy
        self.jump.nearest_neighbour_energy = None
        self.jump.coordination_number_energy = None
        delta_E = final_energy - initial_energy
        assert self.jump.delta_E() == delta_E

    def test_delta_E_nn_interactions(self):
        initial_energy = 1.5
        final_energy = 2.1
        self.jump.initial_site.energy = initial_energy
        self.jump.final_site.energy = final_energy
        self.jump.nearest_neighbour_energy = True
        self.jump.coordination_number_energy = None
        self.jump.nearest_neighbour_delta_E = Mock(return_value=0.3)
        delta_E = final_energy - initial_energy + self.jump.nearest_neighbour_delta_E.return_value
        assert self.jump.delta_E() == delta_E

    def test_delta_E_cn_dependent(self):
        initial_energy = 1.5
        final_energy = 2.1
        self.jump.initial_site.energy = initial_energy
        self.jump.final_site.energy = final_energy
        self.jump.nearest_neighbour_energy = None
        self.jump.coordination_number_energy = True
        self.jump.coordination_number_delta_E = Mock(return_value=0.5)
        delta_E = final_energy - initial_energy + self.jump.coordination_number_delta_E.return_value
        assert self.jump.delta_E() == delta_E

    def test_nearest_neighbour_delta_E(self):
        self.jump.initial_site.nn_occupation = Mock(return_value=3)
        self.jump.final_site.nn_occupation = Mock(return_value=1)
        self.jump.nearest_neighbour_energy = 0.7
        assert self.jump.nearest_neighbour_delta_E() == -3 * 0.7

    def test_dr(self):
        cell_lengths = np.array([10.0, 10.0, 10.0])
        self.jump.initial_site.r = np.array([3.0, 3.0, 3.0])
        self.jump.final_site.r = np.array([4.0, 5.0, 6.0])
        np.testing.assert_array_equal(self.jump.dr(cell_lengths), np.array([1.0, 2.0, 3.0]))

    def test_dr_minimum_image_convention(self):
        cell_lengths = np.array([10.0, 10.0, 10.0])
        self.jump.initial_site.r = np.array([9.0, 1.0, 3.0])
        self.jump.final_site.r = np.array([1.0, 9.0, 1.0])
        np.testing.assert_array_equal(self.jump.dr(cell_lengths), np.array([2.0, -2.0, -2.0]))

    def test_relative_probability_from_lookup_table(self):
        self.jump.initial_site.label = "A"
        self.jump.final_site.label = "B"
        self.jump.initial_site.nn_occupation = Mock(return_value=2)
        self.jump.final_site.nn_occupation = Mock(return_value=1)
        jump_lookup_table = Mock(spec=LookupTable)
        jump_lookup_table.jump_probability = {"A": {"B": {2: {1: 0.5}}}}
        assert self.jump.relative_probability_from_lookup_table(jump_lookup_table) == 0.5

    def test_relative_probability_getter(self):
        assert self.jump.relative_probability == self.jump._relative_probability

    def test_relative_probability_setter(self):
        self.jump.relative_probability = 0.5
        assert self.jump._relative_probability == 0.5


def total_cn_energy(sites):
    return sum(site.cn_occupation_energy() for site in sites if site.is_occupied)


def make_sites(neighbours, labels):
    """Build connected Site objects from a neighbour map and a label for each site."""
    numbers = {name: i for i, name in enumerate(neighbours, start=1)}
    sites = {
        name: Site(
            number=numbers[name],
            coordinates=np.zeros(3),
            neighbours=[numbers[n] for n in neighbours[name]],
            energy=0.0,
            label=labels[name],
        )
        for name in neighbours
    }
    for name, site in sites.items():
        site.p_neighbours = [sites[n] for n in neighbours[name]]
    return sites


def brute_force_cn_delta_E(initial_site, final_site, sites):
    """Total coordination-number energy of the occupied sites after the jump minus before."""
    before = total_cn_energy(sites)
    initial_site.is_occupied, final_site.is_occupied = False, True
    after = total_cn_energy(sites)
    initial_site.is_occupied, final_site.is_occupied = True, False
    return after - before


def test_coordination_number_delta_E_counts_a_site_next_to_both_ends_once():
    # Sites i and f share the neighbour s, which also neighbours t; u neighbours f:
    #
    #     i --- f --- u
    #      \   /
    #        s --- t
    #
    # i, s, t and u are occupied, and the atom on i jumps to f. s keeps two
    # occupied neighbours, so its energy does not change; counting it twice
    # would add E(1) + E(3) - 2 E(2) = 0.2 eV.
    sites = make_sites(
        {"i": ["f", "s"], "f": ["i", "s", "u"], "s": ["i", "f", "t"], "t": ["s"], "u": ["f"]},
        {name: "A" for name in "ifstu"},
    )
    cn = {"A": {0: 0.0, 1: 0.1, 2: 0.4, 3: 0.9}}
    for name, site in sites.items():
        site.set_cn_occupation_energies(cn)
        site.is_occupied = name in "istu"
    jump = Jump(sites["i"], sites["f"], coordination_number_energy={"A": cn}, params=PARAMS)
    expected = brute_force_cn_delta_E(sites["i"], sites["f"], sites.values())
    assert expected == pytest.approx(0.4)
    assert jump.coordination_number_delta_E() == pytest.approx(expected)


def test_coordination_number_delta_E_matches_brute_force_on_a_triangular_lattice():
    # Random occupations, labels and nonlinear energies on a periodic 4 x 4
    # triangular lattice, where every bond is part of a triangle.
    rng = np.random.default_rng(1)
    n = 4
    offsets = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
    neighbours = {(x, y): [((x + dx) % n, (y + dy) % n) for dx, dy in offsets] for x in range(n) for y in range(n)}
    for _ in range(5):
        labels = {key: str(rng.choice(["A", "B"])) for key in neighbours}
        sites = make_sites(neighbours, labels)
        cn = {a: {b: {c: float(rng.uniform(-0.5, 0.5)) for c in range(7)} for b in "AB"} for a in "AB"}
        for site in sites.values():
            site.set_cn_occupation_energies(cn[site.label])
            site.is_occupied = bool(rng.random() < 0.5)
        for initial_site in sites.values():
            for final_site in initial_site.p_neighbours:
                if initial_site.is_occupied and not final_site.is_occupied:
                    jump = Jump(initial_site, final_site, coordination_number_energy=cn, params=PARAMS)
                    expected = brute_force_cn_delta_E(initial_site, final_site, sites.values())
                    assert jump.coordination_number_delta_E() == pytest.approx(expected, abs=1e-12)
