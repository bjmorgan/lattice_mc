from unittest.mock import Mock, call, patch

import numpy as np
import pytest

from lattice_mc.atom import Atom
from lattice_mc.cluster import Cluster
from lattice_mc.error import BlockedLatticeError
from lattice_mc.jump import Jump
from lattice_mc.lattice import Lattice
from lattice_mc.lattice_site import Site
from lattice_mc.simulation import SimulationParameters
from lattice_mc.transitions import Transitions

PARAMS = SimulationParameters(temperature=298.0, rate_prefactor=1e13)


class TestLattice:
    """Tests for Lattice class"""

    @pytest.fixture(autouse=True)
    @patch("lattice_mc.lattice.Lattice.site_with_id")
    @patch("lattice_mc.lattice.Lattice.initialise_site_lookup_table")
    @patch("lattice_mc.lattice.Lattice.enforce_periodic_boundary_conditions")
    def set_up(self, pbc, islt, site_id):
        site_id.side_effect = (1, 2, 3, 4, 5, 6, 7, 8)
        self.site_id = site_id
        site_labels = ["A", "B", "A", "B", "C"]
        site_neighbours = [[2, 3], [1, 3], [1, 2], [5], [4]]
        self.mock_sites = [
            Mock(spec=Site, number=i, label=label, neighbours=n)
            for i, (label, n) in enumerate(zip(site_labels, site_neighbours), start=1)
        ]
        self.cell_lengths = np.array([7.0, 8.0, 9.0])
        self.lattice = Lattice(self.mock_sites, self.cell_lengths)

    def test_lattice_is_initialised(self):
        np.testing.assert_array_equal(self.lattice.cell_lengths, self.cell_lengths)
        assert self.lattice.sites == self.mock_sites
        assert self.lattice.number_of_sites == 5
        assert self.lattice.site_labels == {"A", "B", "C"}
        assert self.lattice.site_populations == {"A": 2, "B": 2, "C": 1}
        self.site_id.assert_has_calls([call(2), call(3), call(1), call(3), call(1), call(2), call(5), call(4)])
        assert self.mock_sites[0].p_neighbours == [1, 2]
        assert self.mock_sites[1].p_neighbours == [3, 4]
        assert self.mock_sites[2].p_neighbours == [5, 6]
        assert self.mock_sites[3].p_neighbours == [7]
        assert self.mock_sites[4].p_neighbours == [8]

    def test_dimensionality_is_three_for_three_non_zero_cell_lengths(self):
        assert self.lattice.dimensionality == 3

    def test_dimensionality_excludes_zero_cell_lengths(self):
        self.lattice.cell_lengths = np.array([7.0, 8.0, 0.0])
        assert self.lattice.dimensionality == 2

    def test_enforce_periodice_boundary_consitions(self):
        site_coordinates = [np.array([-1.0, 11.0, 3.0]), np.array([2.0, -3.0, 12.0])]
        mock_sites = [Mock(spec=Site, r=r) for r in site_coordinates]
        self.lattice.cell_lengths = np.array([10.0, 10.0, 10.0])
        self.lattice.sites = mock_sites
        self.lattice.enforce_periodic_boundary_conditions()
        np.testing.assert_array_equal(self.lattice.sites[0].r, np.array([9.0, 1.0, 3.0]))
        np.testing.assert_array_equal(self.lattice.sites[1].r, np.array([2.0, 7.0, 2.0]))

    def test_reset(self):
        self.lattice.reset()
        assert self.lattice.time == 0.0
        assert self.lattice.sites[0].time_occupied == 0.0

    def test_initialised_site_lookup_table(self):
        for i, site in enumerate(self.lattice.sites):
            site.number = i
        self.lattice.initialise_site_lookup_table()
        assert self.lattice.site_lookup[0] == self.mock_sites[0]

    def test_site_with_id(self):
        self.lattice.site_lookup = ["foo", "bar"]
        assert self.lattice.site_with_id(1) == "bar"

    def test_vacant_sites(self):
        occupied = [True, False, True, False, True]
        for o, site in zip(occupied, self.lattice.sites):
            site.is_occupied = o
        assert list(self.lattice.vacant_sites()) == self.mock_sites[1:5:2]

    def test_occupied_sites(self):
        occupied = [True, False, True, False, True]
        for o, site in zip(occupied, self.lattice.sites):
            site.is_occupied = o
        assert list(self.lattice.occupied_sites()) == self.mock_sites[0:5:2]

    def test_vacant_site_numbers(self):
        occupied = [True, False, True, False, True]
        numbers = [4, 5, 6, 7, 8]
        for n, o, site in zip(numbers, occupied, self.lattice.sites):
            site.is_occupied = o
            site.number = n
        assert self.lattice.vacant_site_numbers() == [5, 7]

    def test_occupied_site_numbers(self):
        occupied = [True, False, True, False, True]
        numbers = [4, 5, 6, 7, 8]
        for n, o, site in zip(numbers, occupied, self.lattice.sites):
            site.is_occupied = o
            site.number = n
        assert self.lattice.occupied_site_numbers() == [4, 6, 8]

    @patch("lattice_mc.jump.Jump")
    def test_potential_jumps_if_lattice_is_mostly_empty(self, mock_Jump):  # lattice is mostly vacant
        jumps = [Mock(spec=Jump), Mock(spec=Jump)]
        mock_Jump.side_effect = jumps
        self.lattice.number_of_occupied_sites = 1
        self.lattice_number_of_sites = 4
        self.lattice.params = PARAMS
        site = Mock(spec=Site)
        site.neighbours = [2, 3]
        occupied_sites = [site]
        unoccupied_sites = [Mock(spec=Site), Mock(spec=Site), Mock(spec=Site)]
        self.lattice.nn_energy = "A"
        self.lattice.cn_energies = "B"
        self.lattice.jump_lookup_table = "C"
        for s in unoccupied_sites:
            s.is_occupied = False
        with patch("lattice_mc.lattice.Lattice.occupied_sites") as mock_occupied_sites:
            with patch("lattice_mc.lattice.Lattice.site_with_id") as mock_site_with_id:
                mock_occupied_sites.return_value = occupied_sites
                mock_site_with_id.side_effect = unoccupied_sites[:2]
                potential_jumps = self.lattice.potential_jumps()
                assert potential_jumps == jumps
                assert mock_Jump.mock_calls[0][1] == (site, unoccupied_sites[0], "A", "B", "C")
                assert mock_Jump.mock_calls[0][2] == {"params": PARAMS}
                assert mock_Jump.mock_calls[1][1] == (site, unoccupied_sites[1], "A", "B", "C")
                assert mock_Jump.mock_calls[1][2] == {"params": PARAMS}
                mock_site_with_id.assert_has_calls([call(2), call(3)])

    @patch("lattice_mc.jump.Jump")
    def test_potential_jumps_if_lattice_is_mostly_filled(self, mock_Jump):  # lattice is mostly occupied
        mock_Jump.side_effect = ["jump1", "jump2"]
        self.lattice.number_of_occupied_sites = 3
        self.lattice_number_of_sites = 4
        self.lattice.params = PARAMS
        site = Mock(spec=Site)
        site.neighbours = [2, 3]
        vacant_sites = [site]
        occupied_sites = [Mock(spec=Site), Mock(spec=Site), Mock(spec=Site)]
        for s in occupied_sites:
            s.is_occupied = True
        self.lattice.nn_energy = "A"
        self.lattice.cn_energies = "B"
        self.lattice.jump_lookup_table = "C"
        with patch("lattice_mc.lattice.Lattice.vacant_sites") as mock_vacant_sites:
            with patch("lattice_mc.lattice.Lattice.site_with_id") as mock_site_with_id:
                mock_vacant_sites.return_value = vacant_sites
                mock_site_with_id.side_effect = occupied_sites[:2]
                jumps = self.lattice.potential_jumps()
                assert jumps == ["jump1", "jump2"]
                assert mock_Jump.mock_calls[0][1] == (occupied_sites[0], site, "A", "B", "C")
                assert mock_Jump.mock_calls[0][2] == {"params": PARAMS}
                assert mock_Jump.mock_calls[1][1] == (occupied_sites[1], site, "A", "B", "C")
                assert mock_Jump.mock_calls[1][2] == {"params": PARAMS}
                mock_site_with_id.assert_has_calls([call(2), call(3)])

    def test_update(self):
        atom = Mock(spec=Atom)  # the jumping atom
        atom.number = 28
        atom.number_of_hops = 4
        atom.dr = np.array([1.0, 2.0, 3.0])
        atom.summed_dr2 = 1.3
        jump = Mock(spec=Jump)  # the jump selected for the update
        jump.initial_site = Mock(spec=Site)  # atom jumps from this site
        jump.final_site = Mock(spec=Site)  # atom jumps to this site
        jump.initial_site.atom = atom
        jump.dr = Mock(return_value=np.array([2.0, 3.0, 4.0]))
        self.lattice.cell_lengths = np.array([1.0, 2.0, 3.0])
        self.lattice.update(jump)
        jump.dr.assert_called_with(self.lattice.cell_lengths)
        assert jump.final_site.occupation == atom.number
        assert jump.initial_site.occupation == 0
        assert jump.initial_site.atom is None
        assert jump.final_site.atom == atom
        assert jump.final_site.is_occupied is True
        assert jump.initial_site.is_occupied is False
        assert atom.site == jump.final_site
        assert atom.number_of_hops == 5
        np.testing.assert_array_equal(atom.dr, np.array([3.0, 5.0, 7.0]))
        assert atom.summed_dr2 == 30.3

    def test_populate_sites_raises_ValueError_with_too_many_atoms(self):
        with pytest.raises(ValueError):
            self.lattice.populate_sites(number_of_atoms=10)

    def test_populate_sites(self):
        number_of_atoms = 2
        with patch("random.sample") as mock_random_sample:
            with patch("lattice_mc.atom.Atom") as mock_Atom:
                mock_sites = [Mock(spec=Site), Mock(spec=Site)]
                mock_random_sample.return_value = mock_sites
                mock_Atom.return_value = [Mock(spec=Atom), Mock(spec=Atom)]
                self.lattice.populate_sites(2)
                assert mock_Atom.mock_calls[0][2]["initial_site"] == mock_sites[0]
                assert mock_Atom.mock_calls[1][2]["initial_site"] == mock_sites[1]
                assert self.lattice.number_of_occupied_sites == number_of_atoms

    def test_populate_sites_with_selected_sites(self):
        with patch("random.sample") as mock_random_sample:
            with patch("lattice_mc.atom.Atom") as mock_Atom:
                mock_sites = [Mock(spec=Site), Mock(spec=Site)]
                mock_random_sample.return_value = mock_sites
                mock_Atom.return_value = [Mock(spec=Atom), Mock(spec=Atom)]
                self.lattice.populate_sites(2, selected_sites=["A"])
                for site in mock_random_sample.call_args[0][0]:
                    assert site.label == "A"

    @patch("lattice_mc.transitions.Transitions")
    def test_jump(self, mock_Transitions):
        potential_jumps = [Mock(spec=Jump), Mock(spec=Jump)]
        self.lattice.potential_jumps = Mock(return_value=potential_jumps)
        self.lattice.params = PARAMS
        selected_jump = Mock(spec=Jump)
        mock_transitions = Mock(spec=Transitions)
        mock_transitions.random = Mock(return_value=selected_jump)
        mock_transitions.time_to_jump = Mock(return_value=5.0)
        mock_Transitions.return_value = mock_transitions
        self.lattice.time = 2.0
        self.lattice.update_site_occupation_times = Mock()
        self.lattice.update = Mock()
        self.lattice.jump()
        mock_Transitions.assert_called_with(potential_jumps, params=PARAMS)
        self.lattice.update.assert_called_with(selected_jump)
        self.lattice.update_site_occupation_times.assert_called_with(5.0)
        assert self.lattice.time == 2.0 + 5.0

    def test_potential_jumps_raises_RuntimeError_if_params_not_set(self):
        self.lattice.number_of_occupied_sites = 1
        self.lattice.params = None
        with pytest.raises(RuntimeError):
            self.lattice.potential_jumps()

    def test_jump_raises_BlockedLatticeError_if_no_possible_jumps(self):
        potential_jumps = []
        self.lattice.potential_jumps = Mock(return_value=potential_jumps)
        with pytest.raises(BlockedLatticeError):
            self.lattice.jump()

    def test_update_site_occupation_times(self):
        occupied_sites = [Mock(spec=Site), Mock(spec=Site)]
        occupied_sites[0].time_occupied = 2.0
        occupied_sites[1].time_occupied = 3.2
        self.lattice.occupied_sites = Mock(return_value=occupied_sites)
        self.lattice.update_site_occupation_times(4.2)
        assert occupied_sites[0].time_occupied == 2.0 + 4.2
        assert occupied_sites[1].time_occupied == 3.2 + 4.2

    def test_site_occupation_statistics(self):
        self.lattice.site_labels = ["A", "B"]
        sites = [Mock(spec=Site), Mock(spec=Site)]
        sites[0].time_occupied = 2.0
        sites[1].time_occupied = 3.2
        sites[0].label = "A"
        sites[1].label = "B"
        self.lattice.sites = sites
        self.lattice.time = 4.0
        occupation_stats = self.lattice.site_occupation_statistics()
        assert occupation_stats == {"A": 0.5, "B": 0.8}

    def test_site_energies_by_label(self):
        self.lattice.sites = [
            Site(number=i, coordinates=np.zeros(3), neighbours=[], energy=energy, label=label)
            for i, (label, energy) in enumerate([("A", 0.2), ("B", 0.4), ("A", 0.2)], start=1)
        ]
        assert self.lattice.site_energies_by_label() == {"A": 0.2, "B": 0.4}

    def test_site_energies_by_label_raises_ValueError_for_conflicting_energies(self):
        self.lattice.sites = [
            Site(number=i, coordinates=np.zeros(3), neighbours=[], energy=energy, label="A")
            for i, energy in enumerate([0.2, 0.3], start=1)
        ]
        with pytest.raises(ValueError, match="label A"):
            self.lattice.site_energies_by_label()

    def test_set_nn_energy(self):
        self.lattice.set_nn_energy("foo")
        assert self.lattice.nn_energy == "foo"

    def test_set_cn_energies(self):
        cn_energies = {"A": 0.5, "B": 0.3}
        sites = [Mock(spec=Site), Mock(spec=Site)]
        sites[0].label = "A"
        sites[1].label = "B"
        sites[0].set_cn_occupation_energies = Mock()
        sites[1].set_cn_occupation_energies = Mock()
        self.lattice.sites = sites
        self.lattice.set_cn_energies(cn_energies)
        sites[0].set_cn_occupation_energies.assert_called_with(0.5)
        sites[1].set_cn_occupation_energies.assert_called_with(0.3)
        assert self.lattice.cn_energies == cn_energies

    def test_site_coordination_numbers(self):
        sites = [Mock(spec=Site), Mock(spec=Site)]
        sites[0].label = "A"
        sites[1].label = "B"
        sites[0].neighbours = [1, 2, 3, 4]
        sites[1].neighbours = [1, 2, 3, 4, 5, 6]
        self.lattice.sites = sites
        self.lattice.site_labels = ("A", "B")
        assert self.lattice.site_coordination_numbers() == {"A": {4}, "B": {6}}

    def test_site_coordination_numbers_for_varied_coordination_numbers(self):
        sites = [Mock(spec=Site), Mock(spec=Site)]
        sites[0].label = "A"
        sites[1].label = "A"
        sites[0].neighbours = [1, 2, 3, 4]
        sites[1].neighbours = [1, 2]
        self.lattice.sites = sites
        self.lattice.site_labels = {"A"}
        assert self.lattice.site_coordination_numbers() == {"A": {2, 4}}

    def test_max_site_coordination_numbers(self):
        with patch("lattice_mc.lattice.Lattice.site_coordination_numbers") as mock_site_coordination_numbers:
            mock_site_coordination_numbers.return_value = {"A": {4}, "B": {2, 4}}
            assert self.lattice.max_site_coordination_numbers() == {"A": 4, "B": 4}

    def test_site_specific_coordination_numbers(self):
        sites = [Mock(spec=Site), Mock(spec=Site)]
        sites[0].label = "A"
        sites[1].label = "B"
        sites[0].site_specific_neighbours = Mock(return_value="foo")
        sites[1].site_specific_neighbours = Mock(return_value="bar")
        self.lattice.sites = sites
        assert self.lattice.site_specific_coordination_numbers() == {"A": "foo", "B": "bar"}

    def test_connected_site_pairs(self):
        sites = [Mock(spec=Site), Mock(spec=Site)]
        sites[0].label = "A"
        sites[1].label = "B"
        sites[0].p_neighbours = [sites[1]]
        sites[1].p_neighbours = [sites[0]]
        self.lattice.sites = sites
        assert self.lattice.connected_site_pairs() == {"A": ["B"], "B": ["A"]}

    def test_connected_site_pairs_combines_all_connections(self):
        sites = [Mock(spec=Site), Mock(spec=Site), Mock(spec=Site), Mock(spec=Site)]
        for s, label in zip(sites, ["A", "A", "B", "C"]):
            s.label = label
        sites[0].p_neighbours = [sites[1], sites[2]]
        sites[1].p_neighbours = [sites[0], sites[3]]
        sites[2].p_neighbours = [sites[0]]
        sites[3].p_neighbours = [sites[1]]
        self.lattice.sites = sites
        assert self.lattice.connected_site_pairs()["A"] == ["A", "B", "C"]

    @patch("random.sample")
    def test_transmute_sites_selects_from_old_site_label(self, mock_random_sample):
        sites = [Mock(spec=Site), Mock(spec=Site)]
        mock_random_sample.return_value = [sites[1]]
        sites[0].label = "A"
        sites[1].label = "B"
        self.lattice.sites = sites
        self.lattice.site_labels = set([s.label for s in self.lattice.sites])
        self.lattice.select_sites = Mock(return_value=[sites[1]])
        self.lattice.transmute_sites("B", "C", 1)
        assert set([s.label for s in self.lattice.sites]) == set(["A", "C"])
        assert self.lattice.site_labels == set(["A", "C"])
        mock_random_sample.assert_called_with([sites[1]], 1)

    @patch("lattice_mc.cluster.Cluster")
    def test_connected_sites(self, mock_cluster):
        sites = [Mock(spec=Site), Mock(spec=Site)]
        sites[0].label = "A"
        sites[1].label = "A"
        self.lattice.site_labels = set([s.label for s in sites])
        initial_clusters = [Mock(spec=Cluster), Mock(spec=Cluster)]
        initial_clusters[0].sites = sites[0]
        initial_clusters[1].sites = sites[1]
        initial_clusters[0].neighbours = sites[0]
        initial_clusters[1].neighbours = sites[1]
        merged_cluster = Mock(spec=Cluster)
        merged_cluster.neighbours = ()
        mock_cluster.side_effect = initial_clusters
        self.lattice.sites = sites
        initial_clusters[0].is_neighbouring.return_value = True
        initial_clusters[0].merge.return_value = merged_cluster
        assert self.lattice.connected_sites() == [merged_cluster]

    @patch("lattice_mc.cluster.Cluster")
    def test_connected_sites_with_site_labels(self, mock_cluster):
        sites = [Mock(spec=Site), Mock(spec=Site)]
        sites[0].label = "A"
        sites[1].label = "B"
        self.lattice.site_labels = set([s.label for s in sites])
        initial_clusters = [Mock(spec=Cluster)]
        initial_clusters[0].sites = sites[0]
        initial_clusters[0].neighbours = ()
        mock_cluster.side_effect = initial_clusters
        self.lattice.sites = sites
        self.lattice.select_sites = Mock(return_value=[sites[0]])
        initial_clusters[0].is_neighbouring.return_value = False
        site_labels = "B"
        assert self.lattice.connected_sites(site_labels) == initial_clusters
        self.lattice.select_sites.assert_called_with(site_labels)

    def test_select_sites(self):
        sites = [Mock(spec=Site), Mock(spec=Site), Mock(spec=Site)]
        sites[0].label = "A"
        sites[1].label = "AA"
        sites[2].label = "B"
        self.lattice.site_labels = set([s.label for s in sites])
        self.lattice.sites = sites
        assert self.lattice.select_sites("A") == [sites[0]]
        assert self.lattice.select_sites("AA") == [sites[1]]
        assert self.lattice.select_sites(["A", "B"]) == [sites[0], sites[2]]

    def test_detached_sites(self):
        sites = [Mock(spec=Site), Mock(spec=Site)]
        clusters = [Mock(spec=Cluster), Mock(spec=Cluster)]
        clusters[0].is_periodically_contiguous.return_value = (True, True, True)
        clusters[1].is_periodically_contiguous.return_value = (False, False, False)
        clusters[0].sites = set([sites[0]])
        clusters[1].sites = set([sites[1]])
        sites[0].label = "A"
        sites[1].label = "B"
        self.lattice.site_labels = set([s.label for s in sites])
        self.lattice.connected_sites = Mock(return_value=clusters)
        self.lattice.sites = sites
        assert self.lattice.detached_sites() == [sites[1]]

    @patch("lattice_mc.lattice.Lattice.site_with_id")
    @patch("lattice_mc.lattice.Lattice.initialise_site_lookup_table")
    @patch("lattice_mc.lattice.Lattice.enforce_periodic_boundary_conditions")
    @patch("lattice_mc.lattice.Lattice.potential_jumps")
    def test_is_blocked_returns_true(self, p_jumps, pbc, islt, site_id):
        p_jumps.return_value = []
        site_id.side_effect = (1, 2)
        site_labels = ["A", "A"]
        site_neighbours = [[], []]
        mock_sites = [
            Mock(spec=Site, number=i, label=label, neighbours=n)
            for i, (label, n) in enumerate(zip(site_labels, site_neighbours), start=1)
        ]
        cell_lengths = np.array([7.0, 8.0, 9.0])
        lattice = Lattice(mock_sites, cell_lengths)
        lattice.sites[0].is_occupied = True
        lattice.sites[1].is_occupied = True
        lattice.number_of_occupied_sites = 1
        assert lattice.is_blocked() is True

    @patch("lattice_mc.lattice.Lattice.site_with_id")
    @patch("lattice_mc.lattice.Lattice.initialise_site_lookup_table")
    @patch("lattice_mc.lattice.Lattice.enforce_periodic_boundary_conditions")
    @patch("lattice_mc.lattice.Lattice.potential_jumps")
    def test_is_blocked_returns_false(self, p_jumps, pbc, islt, site_id):
        p_jumps.return_value = [Mock(spec=Jump)]
        site_id.side_effect = (1, 2)
        site_labels = ["A", "A"]
        site_neighbours = [[2], [1]]
        mock_sites = [
            Mock(spec=Site, number=i, label=label, neighbours=n)
            for i, (label, n) in enumerate(zip(site_labels, site_neighbours), start=1)
        ]
        cell_lengths = np.array([7.0, 8.0, 9.0])
        lattice = Lattice(mock_sites, cell_lengths)
        lattice.sites[0].is_occupied = True
        lattice.sites[1].is_occupied = True
        lattice.number_of_occupied_sites = 1
        assert lattice.is_blocked() is False


def make_lattice(numbers, labels=None, neighbours=None, energies=None):
    """Build a Lattice from real Sites with the given site numbers."""
    n = len(numbers)
    labels = labels or ["A"] * n
    neighbours = neighbours or [[] for _ in numbers]
    energies = energies or [0.0] * n
    sites = [
        Site(number=number, coordinates=np.zeros(3), neighbours=nn, energy=e, label=label)
        for number, nn, e, label in zip(numbers, neighbours, energies, labels)
    ]
    return Lattice(sites, np.array([10.0, 10.0, 10.0]))


class TestLatticeIndexing:
    """Tests for site indices and neighbour indices"""

    def test_site_index_uses_position_not_number(self):
        lattice = make_lattice([-2, 5, 3])
        assert [lattice.site_index(number) for number in (-2, 5, 3)] == [0, 1, 2]

    def test_site_index_raises_KeyError_for_unknown_number(self):
        lattice = make_lattice([-2, 5, 3])
        with pytest.raises(KeyError):
            lattice.site_index(4)

    def test_neighbour_indices(self):
        lattice = make_lattice([-2, 5, 3], neighbours=[[5, 3], [-2], [-2, 5]])
        assert lattice.neighbour_indices == ((1, 2), (0,), (0, 1))


class TestLatticeSiteEnergies:
    """Tests for setting and reading per-site energies"""

    @pytest.fixture
    def lattice(self):
        return make_lattice([-2, 5, 3], labels=["A", "B", "A"])

    def test_set_site_energies_by_label(self, lattice):
        lattice.set_site_energies(by_label={"A": 0.2, "B": 0.4})
        assert [site.energy for site in lattice.sites] == [0.2, 0.4, 0.2]

    def test_set_site_energies_by_site(self, lattice):
        lattice.set_site_energies(by_site={5: 0.7, -2: 0.1})
        assert [site.energy for site in lattice.sites] == [0.1, 0.7, 0.0]

    def test_set_site_energies_by_site_overrides_by_label(self, lattice):
        lattice.set_site_energies(by_label={"A": 0.2, "B": 0.4}, by_site={3: 0.9})
        assert [site.energy for site in lattice.sites] == [0.2, 0.4, 0.9]

    def test_set_site_energies_requires_an_argument(self, lattice):
        with pytest.raises(ValueError):
            lattice.set_site_energies()

    def test_set_site_energies_unknown_label(self, lattice):
        with pytest.raises(ValueError, match="Unknown site labels"):
            lattice.set_site_energies(by_label={"C": 0.2})

    def test_set_site_energies_unknown_site(self, lattice):
        with pytest.raises(ValueError, match="Unknown site numbers"):
            lattice.set_site_energies(by_site={4: 0.2})

    def test_set_site_energies_unknown_key_leaves_energies_unchanged(self, lattice):
        with pytest.raises(ValueError):
            lattice.set_site_energies(by_label={"A": 0.2}, by_site={4: 0.5})
        assert [site.energy for site in lattice.sites] == [0.0, 0.0, 0.0]

    def test_set_site_energies_unknown_label_with_valid_site_leaves_energies_unchanged(self, lattice):
        with pytest.raises(ValueError):
            lattice.set_site_energies(by_label={"C": 0.2}, by_site={5: 0.5})
        assert [site.energy for site in lattice.sites] == [0.0, 0.0, 0.0]

    def test_site_energies_array(self):
        lattice = make_lattice([-2, 5, 3], energies=[0.1, 0.2, 0.3])
        energies = lattice.site_energies()
        assert energies.dtype == np.float64
        np.testing.assert_array_equal(energies, [0.1, 0.2, 0.3])
