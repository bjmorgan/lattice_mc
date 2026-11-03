"""Tests for lattice_mc._lookup."""

from pathlib import Path

import numpy as np
import pytest

from lattice_mc._energies import EnergyModel, metropolis
from lattice_mc._lookup import LookupTable
from lattice_mc.init_lattice import honeycomb_lattice, lattice_from_sites_file
from lattice_mc.lattice import Lattice
from lattice_mc.lattice_site import Site

KT = 0.025

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


def lattice_from_text(tmp_path: Path, text: str) -> Lattice:
    path = tmp_path / "sites.dat"
    path.write_text(text)
    return lattice_from_sites_file(str(path), cell_lengths=[10.0, 10.0, 10.0])


def star_lattice() -> Lattice:
    """A central 'A' site, with three neighbouring 'B' sites that each have only the centre as a neighbour."""
    neighbours = {1: [2, 3, 4], 2: [1], 3: [1], 4: [1]}
    sites = [Site(number, np.zeros(3), neighbours[number], 0.0, "A" if number == 1 else "B") for number in neighbours]
    return Lattice(sites, np.array([10.0, 10.0, 10.0]))


def alternating_honeycomb() -> Lattice:
    """A honeycomb lattice with 'A' and 'B' sites of different energies."""
    lattice = honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
    lattice.set_site_energies(by_label={"A": 0.05, "B": -0.03})
    return lattice


class TestRelativeProbability:
    def test_matches_the_energy_model_for_every_jump_in_random_configurations(self) -> None:
        lattice = alternating_honeycomb()
        model = EnergyModel(lattice, nn_energy=0.04)
        table = LookupTable(model, KT)
        rng = np.random.default_rng(0)
        checked = 0
        for _ in range(10):
            occupied = rng.random(lattice.number_of_sites) < 0.5
            for initial in np.flatnonzero(occupied):
                for final in lattice.neighbour_indices[initial]:
                    if occupied[final]:
                        continue
                    expected = metropolis(model.delta_E(occupied, int(initial), final), KT)
                    assert table.relative_probability(occupied, int(initial), final) == pytest.approx(expected)
                    checked += 1
        assert checked > 0

    def test_without_nn_energy_depends_only_on_site_energies(self) -> None:
        lattice = honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
        lattice.set_site_energies(by_label={"A": 0.0, "B": 0.05})
        table = LookupTable(EnergyModel(lattice), KT)
        occupied = np.zeros(lattice.number_of_sites, dtype=np.bool_)
        initial = next(i for i, site in enumerate(lattice.sites) if site.label == "A")
        occupied[initial] = True
        final = lattice.neighbour_indices[initial][0]
        assert table.relative_probability(occupied, initial, final) == pytest.approx(metropolis(0.05, KT))

    def test_includes_site_energies_from_a_sites_file(self, tmp_path: Path) -> None:
        lattice = lattice_from_text(tmp_path, SITES_WITH_ENERGIES)
        table = LookupTable(EnergyModel(lattice), KT)
        occupied = np.array([True, False])
        assert table.relative_probability(occupied, 0, 1) == pytest.approx(metropolis(0.1, KT))

    def test_matches_the_energy_model_when_labels_have_different_coordination_numbers(self) -> None:
        lattice = star_lattice()
        lattice.set_site_energies(by_label={"A": 0.02, "B": -0.01})
        model = EnergyModel(lattice, nn_energy=0.05)
        table = LookupTable(model, KT)
        largest_counts = set()
        for configuration in range(2**lattice.number_of_sites):
            occupied = np.array([(configuration >> i) & 1 for i in range(lattice.number_of_sites)], dtype=np.bool_)
            for initial in np.flatnonzero(occupied):
                for final in lattice.neighbour_indices[initial]:
                    if occupied[final]:
                        continue
                    expected = metropolis(model.delta_E(occupied, int(initial), final), KT)
                    assert table.relative_probability(occupied, int(initial), final) == pytest.approx(expected)
                    # Record jumps from the 'A' site with 2 occupied neighbours, and onto it with 3.
                    largest_counts.add(
                        (int(initial) == 0 and occupied[1:].sum() == 2, final == 0 and occupied[1:].sum() == 3)
                    )
        assert (True, False) in largest_counts
        assert (False, True) in largest_counts

    def test_uses_the_site_energies_the_model_was_built_with(self) -> None:
        lattice = honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
        lattice.set_site_energies(by_label={"A": 0.0, "B": 0.05})
        model = EnergyModel(lattice, nn_energy=0.04)
        lattice.set_site_energies(by_label={"A": 0.2, "B": -0.2})
        table = LookupTable(model, KT)
        initial = next(i for i, site in enumerate(lattice.sites) if site.label == "A")
        occupied = np.zeros(lattice.number_of_sites, dtype=np.bool_)
        occupied[initial] = True
        final = lattice.neighbour_indices[initial][0]
        expected = metropolis(model.delta_E(occupied, initial, final), KT)
        assert table.relative_probability(occupied, initial, final) == pytest.approx(expected)


class TestConstruction:
    @pytest.mark.parametrize("kT", [0.0, -0.025])
    def test_raises_for_non_positive_kT(self, kT: float) -> None:
        lattice = honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
        with pytest.raises(ValueError, match="kT"):
            LookupTable(EnergyModel(lattice), kT)

    def test_raises_for_a_model_with_coordination_number_energies(self) -> None:
        lattice = honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
        cn = {0: 0.0, 1: 0.1, 2: 0.2, 3: 0.3}
        model = EnergyModel(lattice, cn_energies={"A": {"B": cn}, "B": {"A": cn}})
        with pytest.raises(ValueError, match="coordination"):
            LookupTable(model, KT)

    def test_raises_for_one_label_with_different_energies(self, tmp_path: Path) -> None:
        lattice = lattice_from_text(tmp_path, SITES_WITH_ENERGIES.replace("label: B", "label: A"))
        with pytest.raises(ValueError, match="label A"):
            LookupTable(EnergyModel(lattice), KT)
