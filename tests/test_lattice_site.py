from unittest.mock import Mock

import pytest

from lattice_mc.lattice_site import Site


class TestSite:
    """Test for Site class"""

    @pytest.fixture(autouse=True)
    def set_up(self):
        self.number = "foo1"
        self.coordinates = "foo2"
        self.neighbours = "foo3"
        self.energy = "foo4"
        self.label = "foo5"
        self.cn_energies = {"A": [1.0, 2.0, 3.0], "B": [2.0, 4.0, 6.0]}
        self.site = Site(
            self.number, self.coordinates, self.neighbours, self.energy, self.label, cn_energies=self.cn_energies
        )
        self.neighbouring_sites = [
            Mock(spec=Site, is_occupied=True, label="A"),
            Mock(spec=Site, is_occupied=False, label="A"),
            Mock(spec=Site, is_occupied=True, label="B"),
            Mock(spec=Site, is_occupied=False, label="B"),
        ]

    def test_site_is_initialised(self):
        assert self.site.number == self.number
        assert self.site.index == Site.index - 1
        assert self.site.r == self.coordinates
        assert self.site.neighbours == self.neighbours
        assert self.site.p_neighbours is None
        assert self.site.energy == self.energy
        assert self.site.occupation == 0
        assert self.site.atom is None
        assert self.site.is_occupied is False
        assert self.site.label == self.label
        assert self.site.time_occupied == 0.0
        assert self.site.cn_occupation_energies == self.cn_energies

    def test_nn_occupation(self):
        self.site.p_neighbours = self.neighbouring_sites
        assert self.site.nn_occupation() == 2

    def test_site_specific_nn_occupation(self):
        self.site.p_neighbours = self.neighbouring_sites
        assert self.site.site_specific_nn_occupation() == {"A": 1, "B": 1}

    def test_site_specific_neighbours(self):
        self.site.p_neighbours = self.neighbouring_sites
        assert self.site.site_specific_neighbours() == {"A": 2, "B": 2}

    def test_set_cn_occupation_energies(self):
        self.site.set_cn_occupation_energies("foo7")
        assert self.site.cn_occupation_energies == "foo7"

    def test_unchanged_cn_occupation_energy(self):
        self.site.p_neighbours = self.neighbouring_sites
        self.site_specific_nn_occupation = Mock(return_value={"A": 2, "B": 2})
        assert self.site.cn_occupation_energy() == 6.0

    def test_changed_cn_occupation_energy(self):
        delta_occupation = {"A": 1, "B": -1}
        self.site.p_neighbours = self.neighbouring_sites
        self.site_specific_nn_occupation = Mock(return_value={"A": 2, "B": 2})
        assert self.site.cn_occupation_energy(delta_occupation=delta_occupation) == 5.0
