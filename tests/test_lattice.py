import numpy as np
import pytest

from lattice_mc.lattice import Lattice
from lattice_mc.lattice_site import Site


def make_lattice(numbers, labels=None, neighbours=None, energies=None, cell_lengths=(10.0, 10.0, 10.0)):
    """Build a Lattice from Sites with the given site numbers, all at the origin."""
    n = len(numbers)
    labels = labels or ["A"] * n
    neighbours = neighbours or [[] for _ in numbers]
    energies = energies or [0.0] * n
    sites = [
        Site(number=number, coordinates=np.zeros(3), neighbours=nn, energy=e, label=label)
        for number, nn, e, label in zip(numbers, neighbours, energies, labels)
    ]
    return Lattice(sites, np.array(cell_lengths))


class TestLattice:
    """Tests for the Lattice geometry and labels"""

    @pytest.fixture
    def lattice(self):
        return make_lattice(
            [1, 2, 3, 4, 5], labels=["A", "B", "A", "B", "C"], neighbours=[[2, 3], [1, 3], [1, 2], [5], [4]]
        )

    def test_lattice_is_initialised(self, lattice):
        np.testing.assert_array_equal(lattice.cell_lengths, [10.0, 10.0, 10.0])
        assert [site.number for site in lattice.sites] == [1, 2, 3, 4, 5]
        assert lattice.number_of_sites == 5
        assert lattice.site_labels == {"A", "B", "C"}
        assert lattice.site_populations == {"A": 2, "B": 2, "C": 1}

    def test_dimensionality_is_three_for_three_non_zero_cell_lengths(self, lattice):
        assert lattice.dimensionality == 3

    def test_dimensionality_excludes_zero_cell_lengths(self):
        assert make_lattice([1], cell_lengths=(7.0, 8.0, 0.0)).dimensionality == 2

    def test_sites_are_wrapped_into_the_cell(self):
        sites = [Site(1, np.array([-1.0, 11.0, 3.0]), [], 0.0, "A"), Site(2, np.array([2.0, -3.0, 12.0]), [], 0.0, "A")]
        lattice = Lattice(sites, np.array([10.0, 10.0, 10.0]))
        np.testing.assert_array_equal(lattice.sites[0].r, [9.0, 1.0, 3.0])
        np.testing.assert_array_equal(lattice.sites[1].r, [2.0, 7.0, 2.0])

    def test_site_with_id(self):
        lattice = make_lattice([-2, 5, 3])
        assert lattice.site_with_id(5) is lattice.sites[1]

    def test_site_with_id_raises_KeyError_for_unknown_number(self):
        with pytest.raises(KeyError):
            make_lattice([-2, 5, 3]).site_with_id(4)

    def test_neighbours_of(self):
        lattice = make_lattice([-2, 5, 3], neighbours=[[3, 5], [-2], [-2, 5]])
        assert lattice.neighbours_of(lattice.sites[0]) == (lattice.sites[2], lattice.sites[1])

    def test_site_coordination_numbers(self):
        lattice = make_lattice([1, 2, 3], labels=["A", "A", "B"], neighbours=[[2, 3], [1], [1]])
        assert lattice.site_coordination_numbers() == {"A": {1, 2}, "B": {1}}

    def test_max_site_coordination_numbers(self):
        lattice = make_lattice([1, 2, 3], labels=["A", "A", "B"], neighbours=[[2, 3], [1], [1]])
        assert lattice.max_site_coordination_numbers() == {"A": 2, "B": 1}

    def test_site_specific_coordination_numbers(self, lattice):
        assert lattice.site_specific_coordination_numbers() == {"A": {"A": 1, "B": 1}, "B": {"C": 1}, "C": {"B": 1}}

    def test_connected_site_pairs(self, lattice):
        assert lattice.connected_site_pairs() == {"A": ["B", "A"], "B": ["A", "C"], "C": ["B"]}

    def test_transmute_sites(self, lattice):
        lattice.transmute_sites("B", "D", 1)
        assert sorted(site.label for site in lattice.sites) == ["A", "A", "B", "C", "D"]
        assert lattice.site_labels == {"A", "B", "C", "D"}
        assert lattice.site_populations == {"A": 2, "B": 1, "C": 1, "D": 1}

    def test_transmute_sites_with_the_same_seed_relabels_the_same_sites(self):
        lattices = [make_lattice(list(range(1, 21))) for _ in range(2)]
        for lattice in lattices:
            lattice.transmute_sites("A", "B", 10, seed=3)
        assert [s.label for s in lattices[0].sites] == [s.label for s in lattices[1].sites]

    def test_site_labels_and_populations_follow_changes_to_site_labels(self, lattice):
        lattice.sites[4].label = "A"
        assert lattice.site_labels == {"A", "B"}
        assert lattice.site_populations == {"A": 3, "B": 2}

    def test_select_sites(self):
        lattice = make_lattice([1, 2, 3], labels=["A", "AA", "B"])
        assert lattice.select_sites("A") == [lattice.sites[0]]
        assert lattice.select_sites("AA") == [lattice.sites[1]]
        assert lattice.select_sites(["A", "B"]) == [lattice.sites[0], lattice.sites[2]]

    def test_select_sites_raises_ValueError_for_other_types(self, lattice):
        with pytest.raises(ValueError):
            lattice.select_sites(("A",))

    def test_site_energies_by_label(self):
        lattice = make_lattice([1, 2, 3], labels=["A", "B", "A"], energies=[0.2, 0.4, 0.2])
        assert lattice.site_energies_by_label() == {"A": 0.2, "B": 0.4}

    def test_site_energies_by_label_raises_ValueError_for_conflicting_energies(self):
        lattice = make_lattice([1, 2], energies=[0.2, 0.3])
        with pytest.raises(ValueError, match="label A"):
            lattice.site_energies_by_label()


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
        with pytest.raises(ValueError, match="Unknown site labels"):
            lattice.set_site_energies(by_label={"C": 0.2}, by_site={3: 0.5})
        assert [site.energy for site in lattice.sites] == [0.0, 0.0, 0.0]

    def test_site_energies_array(self):
        lattice = make_lattice([-2, 5, 3], energies=[0.1, 0.2, 0.3])
        energies = lattice.site_energies()
        assert energies.dtype == np.float64
        np.testing.assert_array_equal(energies, [0.1, 0.2, 0.3])
