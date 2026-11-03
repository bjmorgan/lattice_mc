"""Tests for lattice_mc._energies."""

import math

import numpy as np
import numpy.typing as npt
import pytest

from lattice_mc._energies import CNEnergies, EnergyModel, metropolis
from lattice_mc.lattice import Lattice
from lattice_mc.lattice_site import Site


def make_lattice(
    neighbours: dict[str, list[str]],
    labels: dict[str, str] | None = None,
    energies: dict[str, float] | None = None,
) -> Lattice:
    """Build a Lattice from real Sites, from a map of site name to neighbour names.

    Site indices follow the order of `neighbours`.
    """
    numbers = {name: i for i, name in enumerate(neighbours, start=1)}
    sites = [
        Site(
            number=numbers[name],
            coordinates=np.zeros(3),
            neighbours=[numbers[n] for n in neighbours[name]],
            energy=(energies or {}).get(name, 0.0),
            label=(labels or {}).get(name, "A"),
        )
        for name in neighbours
    ]
    return Lattice(sites, np.array([10.0, 10.0, 10.0]))


def occupation(lattice: Lattice, occupied: list[int]) -> npt.NDArray[np.bool_]:
    """A boolean occupation array with the given site indices occupied."""
    result = np.zeros(lattice.number_of_sites, dtype=np.bool_)
    result[occupied] = True
    return result


def total_cn_energy(lattice: Lattice, occupied: npt.NDArray[np.bool_], cn: CNEnergies) -> float:
    """Total coordination-number energy, summed over the occupied sites."""
    total = 0.0
    for i, site in enumerate(lattice.sites):
        if not occupied[i]:
            continue
        neighbour_labels = [lattice.sites[j].label for j in lattice.neighbour_indices[i]]
        for label in set(neighbour_labels):
            count = sum(1 for j in lattice.neighbour_indices[i] if occupied[j] and lattice.sites[j].label == label)
            total += cn[site.label][label][count]
    return total


def brute_force_cn_delta_E(
    lattice: Lattice, occupied: npt.NDArray[np.bool_], initial: int, final: int, cn: CNEnergies
) -> float:
    """Total coordination-number energy after the jump minus before."""
    after = occupied.copy()
    after[initial], after[final] = False, True
    return total_cn_energy(lattice, after, cn) - total_cn_energy(lattice, occupied, cn)


class TestMetropolis:
    def test_negative_delta_E_is_always_accepted(self) -> None:
        assert metropolis(-0.3, 0.025) == 1.0

    def test_zero_delta_E_is_always_accepted(self) -> None:
        assert metropolis(0.0, 0.025) == 1.0

    def test_positive_delta_E_gives_boltzmann_factor(self) -> None:
        assert metropolis(0.05, 0.025) == pytest.approx(math.exp(-2.0))


class TestSiteEnergies:
    def test_delta_E_is_the_site_energy_difference(self) -> None:
        lattice = make_lattice({"i": ["f"], "f": ["i"]}, energies={"i": 1.5, "f": 2.1})
        model = EnergyModel(lattice)
        assert model.delta_E(occupation(lattice, [0]), 0, 1) == pytest.approx(0.6)


class TestNearestNeighbour:
    def test_delta_E_counts_neighbours_of_each_end(self) -> None:
        # The initial site has three occupied neighbours; the final site's only occupied
        # neighbour is the jumping atom, which is not counted after the jump.
        lattice = make_lattice({"i": ["f", "a", "b", "c"], "f": ["i"], "a": ["i"], "b": ["i"], "c": ["i"]})
        model = EnergyModel(lattice, nn_energy=0.7)
        assert model.delta_E(occupation(lattice, [0, 2, 3, 4]), 0, 1) == pytest.approx(-3 * 0.7)

    def test_delta_E_does_not_modify_the_occupation(self) -> None:
        lattice = make_lattice({"i": ["f"], "f": ["i"]})
        occupied = occupation(lattice, [0])
        EnergyModel(lattice, nn_energy=0.7).delta_E(occupied, 0, 1)
        assert occupied.tolist() == [True, False]

    def test_delta_E_adds_the_site_energy_difference(self) -> None:
        # The final site has two occupied neighbours, one of which is the jumping atom.
        lattice = make_lattice({"i": ["f"], "f": ["i", "x"], "x": ["f"]}, energies={"i": 0.1, "f": 0.4})
        model = EnergyModel(lattice, nn_energy=0.7)
        assert model.delta_E(occupation(lattice, [0, 2]), 0, 1) == pytest.approx(0.3 + 0.7)


class TestCoordinationNumber:
    def test_site_next_to_both_ends_is_counted_once(self) -> None:
        # Sites i and f share the neighbour s, which also neighbours t; u neighbours f:
        #
        #     i --- f --- u
        #      \   /
        #        s --- t
        #
        # i, s, t and u are occupied, and the atom on i jumps to f. s keeps two
        # occupied neighbours, so its energy does not change; counting it twice
        # would add E(1) + E(3) - 2 E(2) = 0.2 eV.
        lattice = make_lattice({"i": ["f", "s"], "f": ["i", "s", "u"], "s": ["i", "f", "t"], "t": ["s"], "u": ["f"]})
        cn = {"A": {"A": {0: 0.0, 1: 0.1, 2: 0.4, 3: 0.9}}}
        occupied = occupation(lattice, [0, 2, 3, 4])
        expected = brute_force_cn_delta_E(lattice, occupied, 0, 1, cn)
        assert expected == pytest.approx(0.4)
        assert EnergyModel(lattice, cn_energies=cn).delta_E(occupied, 0, 1) == pytest.approx(expected)

    def test_matches_brute_force_on_a_triangular_lattice(self) -> None:
        # Random occupations, labels and nonlinear energies on a periodic 4 x 4
        # triangular lattice, where every bond is part of a triangle.
        rng = np.random.default_rng(1)
        n = 4
        offsets = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
        neighbours = {
            f"{x}{y}": [f"{(x + dx) % n}{(y + dy) % n}" for dx, dy in offsets] for x in range(n) for y in range(n)
        }
        checked = 0
        for _ in range(5):
            labels = {name: str(rng.choice(["A", "B"])) for name in neighbours}
            lattice = make_lattice(neighbours, labels)
            cn = {a: {b: {c: float(rng.uniform(-0.5, 0.5)) for c in range(7)} for b in "AB"} for a in "AB"}
            occupied = rng.random(lattice.number_of_sites) < 0.5
            model = EnergyModel(lattice, cn_energies=cn)
            for initial in np.flatnonzero(occupied):
                for final in lattice.neighbour_indices[initial]:
                    if not occupied[final]:
                        expected = brute_force_cn_delta_E(lattice, occupied, int(initial), final, cn)
                        assert model.delta_E(occupied, int(initial), final) == pytest.approx(expected, abs=1e-12)
                        checked += 1
        assert checked > 50

    def test_repeated_neighbours_match_brute_force(self) -> None:
        # In very small periodic lattices a site can appear more than once in a neighbour list.
        lattice = make_lattice({"a": ["b", "b"], "b": ["a", "a"]})
        cn = {"A": {"A": {0: 0.0, 1: 0.3, 2: 0.8}}}
        occupied = occupation(lattice, [0])
        expected = brute_force_cn_delta_E(lattice, occupied, 0, 1, cn)
        assert EnergyModel(lattice, cn_energies=cn).delta_E(occupied, 0, 1) == pytest.approx(expected)

    def test_all_three_contributions_are_summed(self) -> None:
        lattice = make_lattice(
            {"i": ["f", "s"], "f": ["i", "s", "t"], "s": ["i", "f"], "t": ["f"]}, energies={"i": 0.1, "f": 0.5}
        )
        cn = {"A": {"A": {0: 0.0, 1: 0.2, 2: 0.7, 3: 1.1}}}
        occupied = occupation(lattice, [0, 2, 3])
        model = EnergyModel(lattice, nn_energy=0.3, cn_energies=cn)
        nn_part = 0.3 * (3 - 1 - 1)
        cn_part = brute_force_cn_delta_E(lattice, occupied, 0, 1, cn)
        assert model.delta_E(occupied, 0, 1) == pytest.approx(0.4 + nn_part + cn_part)


class TestCoordinationNumberValidation:
    LATTICE = {"a": ["b", "c", "d"], "b": ["a"], "c": ["a"], "d": ["a"]}
    LABELS = {"a": "A", "b": "B", "c": "B", "d": "B"}

    def lattice(self) -> Lattice:
        return make_lattice(self.LATTICE, self.LABELS)

    def test_complete_energies_are_accepted(self) -> None:
        cn = {
            "A": {"B": {0: 0.0, 1: 0.1, 2: 0.2, 3: 0.3}},
            "B": {"A": {0: 0.0, 1: 0.1}},
        }
        EnergyModel(self.lattice(), cn_energies=cn)

    def test_extra_neighbour_labels_are_accepted(self) -> None:
        cn = {
            "A": {"B": {0: 0.0, 1: 0.1, 2: 0.2, 3: 0.3}, "C": {0: 0.0}},
            "B": {"A": {0: 0.0, 1: 0.1}},
        }
        EnergyModel(self.lattice(), cn_energies=cn)

    def test_missing_site_label_raises(self) -> None:
        cn = {"A": {"B": {0: 0.0, 1: 0.1, 2: 0.2, 3: 0.3}}}
        with pytest.raises(ValueError, match="site label 'B'"):
            EnergyModel(self.lattice(), cn_energies=cn)

    def test_missing_neighbour_label_raises(self) -> None:
        cn = {"A": {"B": {0: 0.0, 1: 0.1, 2: 0.2, 3: 0.3}}, "B": {}}
        with pytest.raises(ValueError, match="site label 'B'.*neighbour label 'A'"):
            EnergyModel(self.lattice(), cn_energies=cn)

    def test_missing_count_raises(self) -> None:
        cn = {"A": {"B": {0: 0.0, 1: 0.1, 2: 0.2}}, "B": {"A": {0: 0.0, 1: 0.1}}}
        with pytest.raises(ValueError, match="site label 'A'.*neighbour label 'B'.*count 3"):
            EnergyModel(self.lattice(), cn_energies=cn)

    def test_lattice_without_cn_energies_needs_no_checking(self) -> None:
        EnergyModel(self.lattice(), nn_energy=0.2)
