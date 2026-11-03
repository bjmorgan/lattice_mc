"""Tests for lattice_mc._lookup."""

from pathlib import Path

import numpy as np
import pytest

from lattice_mc._energies import EnergyModel, metropolis
from lattice_mc._lookup import LookupTable
from lattice_mc.init_lattice import honeycomb_lattice, lattice_from_sites_file
from lattice_mc.lattice import Lattice

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


class TestRelativeProbability:
    def test_matches_the_energy_model_for_every_jump_in_random_configurations(self) -> None:
        lattice = honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
        lattice.set_site_energies(by_label={"A": 0.05, "B": -0.03})
        model = EnergyModel(lattice, nn_energy=0.04)
        table = LookupTable(lattice, model, KT)
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
        table = LookupTable(lattice, EnergyModel(lattice), KT)
        occupied = np.zeros(lattice.number_of_sites, dtype=np.bool_)
        initial = next(i for i, site in enumerate(lattice.sites) if site.label == "A")
        occupied[initial] = True
        final = lattice.neighbour_indices[initial][0]
        assert table.relative_probability(occupied, initial, final) == pytest.approx(metropolis(0.05, KT))

    def test_includes_site_energies_from_a_sites_file(self, tmp_path: Path) -> None:
        lattice = lattice_from_text(tmp_path, SITES_WITH_ENERGIES)
        table = LookupTable(lattice, EnergyModel(lattice), KT)
        occupied = np.array([True, False])
        assert table.relative_probability(occupied, 0, 1) == pytest.approx(metropolis(0.1, KT))


class TestConstruction:
    def test_raises_for_a_model_with_coordination_number_energies(self) -> None:
        lattice = honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
        cn = {0: 0.0, 1: 0.1, 2: 0.2, 3: 0.3}
        model = EnergyModel(lattice, cn_energies={"A": {"B": cn}, "B": {"A": cn}})
        with pytest.raises(ValueError, match="coordination"):
            LookupTable(lattice, model, KT)

    def test_raises_for_one_label_with_different_energies(self, tmp_path: Path) -> None:
        lattice = lattice_from_text(tmp_path, SITES_WITH_ENERGIES.replace("label: B", "label: A"))
        with pytest.raises(ValueError, match="label A"):
            LookupTable(lattice, EnergyModel(lattice), KT)
