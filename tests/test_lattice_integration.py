from unittest.mock import Mock, patch

import numpy as np
import pytest

from lattice_mc.init_lattice import honeycomb_lattice, square_lattice
from lattice_mc.lattice import Lattice
from lattice_mc.lattice_site import Site


class TestLatticeIntegration:
    """Tests for integration of Lattice class with the rest of the code"""

    @pytest.fixture(autouse=True)
    @patch("lattice_mc.lattice.Lattice.site_with_id")
    @patch("lattice_mc.lattice.Lattice.initialise_site_lookup_table")
    @patch("lattice_mc.lattice.Lattice.enforce_periodic_boundary_conditions")
    def set_up(self, pbc, islt, site_id):
        site_id.side_effect = (1, 2, 3, 4, 5, 6, 7, 8)
        self.site_id = site_id
        site_labels = ["A", "A", "A", "A", "A"]
        site_neighbours = [[2, 3], [1, 3], [1, 2], [5], [4]]
        self.mock_sites = [
            Mock(spec=Site, number=i, label=label, neighbours=n)
            for i, (label, n) in enumerate(zip(site_labels, site_neighbours), start=1)
        ]
        self.cell_lengths = np.array([7.0, 8.0, 9.0])
        self.lattice = Lattice(self.mock_sites, self.cell_lengths)

    def test_connected_sites(self):
        sites = [Mock(spec=Site), Mock(spec=Site)]
        sites[0].p_neighbours = [sites[1]]
        sites[1].p_neighbours = [sites[0]]
        for s in sites:
            s.label = "A"
        self.lattice.sites = sites
        connected_sites = self.lattice.connected_sites()
        assert len(connected_sites) == 1
        for s in sites:
            assert s in connected_sites[0].sites
        assert list(connected_sites[0].neighbours) == []

    def test_connected_sites_2(self):
        sites = [Mock(spec=Site), Mock(spec=Site), Mock(spec=Site), Mock(spec=Site)]
        sites[0].p_neighbours = [sites[1]]
        sites[1].p_neighbours = [sites[0]]
        sites[2].p_neighbours = [sites[3]]
        sites[3].p_neighbours = [sites[2]]
        for s in sites:
            s.label = "A"
        self.lattice.sites = sites
        connected_sites = self.lattice.connected_sites()
        assert len(connected_sites) == 2
        for c in connected_sites:
            if sites[0] in c.sites:
                assert sites[1] in c.sites
            else:
                assert sites[2] in c.sites
                assert sites[3] in c.sites

    def test_connected_sites_3(self):
        sites = [Mock(spec=Site), Mock(spec=Site), Mock(spec=Site), Mock(spec=Site), Mock(spec=Site)]
        sites[0].p_neighbours = [sites[1], sites[2]]
        sites[1].p_neighbours = [sites[0], sites[2]]
        sites[2].p_neighbours = [sites[0], sites[1]]
        sites[3].p_neighbours = [sites[4]]
        sites[4].p_neighbours = [sites[3]]
        self.lattice.sites = sites
        for s in sites:
            s.label = "A"
        connected_sites = self.lattice.connected_sites()
        assert len(connected_sites) == 2
        for c in connected_sites:
            if sites[0] in c.sites:
                assert sites[1] in c.sites
                assert sites[2] in c.sites
                assert c.size() == 3
            else:
                assert sites[3] in c.sites
                assert sites[4] in c.sites
                assert c.size() == 2


class TestDetachedSitesIntegration:
    """Tests for Lattice.detached_sites using real Site and Lattice objects."""

    def test_disconnected_cliques_are_all_detached(self):
        clique_coordinates = [
            np.array([0.0, 0.0, 1.0]),
            np.array([1.0, 1.0, 0.0]),
            np.array([0.0, 1.0, 1.0]),
            np.array([1.0, 0.0, 0.0]),
        ]
        sites = []
        for clique, offset in enumerate([np.zeros(3), np.array([5.0, 5.0, 5.0])]):
            numbers = [clique * 4 + i + 1 for i in range(4)]
            for number, r in zip(numbers, clique_coordinates):
                neighbours = [n for n in numbers if n != number]
                sites.append(Site(number=number, coordinates=r + offset, neighbours=neighbours, energy=0.0, label="A"))
        lattice = Lattice(sites, cell_lengths=np.array([10.0, 10.0, 10.0]))
        assert {s.number for s in lattice.detached_sites()} == {s.number for s in sites}

    def test_winding_ring_has_no_detached_sites(self):
        sites = [
            Site(number=1, coordinates=np.array([0.0, 0.0, 0.0]), neighbours=[2, 3], energy=0.0, label="A"),
            Site(number=2, coordinates=np.array([1.0, 0.0, 0.0]), neighbours=[1, 3], energy=0.0, label="A"),
            Site(number=3, coordinates=np.array([2.0, 0.0, 0.0]), neighbours=[1, 2], energy=0.0, label="A"),
        ]
        lattice = Lattice(sites, cell_lengths=np.array([3.0, 10.0, 10.0]))
        assert lattice.detached_sites() == []

    def test_cluster_closed_only_through_another_label_is_detached(self):
        # The ring 1-2-3 winds along x, but only through the B site 3, so the
        # A cluster {1, 2} does not wind.
        sites = [
            Site(number=1, coordinates=np.array([0.0, 0.0, 0.0]), neighbours=[2, 3], energy=0.0, label="A"),
            Site(number=2, coordinates=np.array([1.0, 0.0, 0.0]), neighbours=[1, 3], energy=0.0, label="A"),
            Site(number=3, coordinates=np.array([2.0, 0.0, 0.0]), neighbours=[1, 2], energy=0.0, label="B"),
        ]
        lattice = Lattice(sites, cell_lengths=np.array([3.0, 10.0, 10.0]))
        assert {s.number for s in lattice.detached_sites("A")} == {1, 2}

    @pytest.mark.parametrize(
        "lattice", [square_lattice(3, 3, 1.0), honeycomb_lattice(2, 2, 1.0)], ids=["square", "honeycomb"]
    )
    def test_fully_occupied_two_dimensional_lattice_winds_along_x_and_y(self, lattice):
        (cluster,) = lattice.connected_sites()
        assert cluster.is_periodically_contiguous(lattice.cell_lengths) == (True, True, False)
        assert lattice.detached_sites() == []

    def test_dimer_on_a_square_lattice_is_detached(self):
        lattice = square_lattice(4, 4, 1.0)
        relabel(lattice, {square_site(0, 0, 4), square_site(1, 0, 4)})
        assert {s.number for s in lattice.detached_sites("A")} == {square_site(0, 0, 4), square_site(1, 0, 4)}

    def test_block_straddling_the_boundary_is_detached(self):
        lattice = square_lattice(4, 4, 1.0)
        block = {square_site(3, 3, 4), square_site(0, 3, 4), square_site(3, 0, 4), square_site(0, 0, 4)}
        relabel(lattice, block)
        assert {s.number for s in lattice.detached_sites("A")} == block

    def test_diagonal_winding_is_reported_along_both_axes(self):
        # (0,0)-(1,0)-(1,1)-(2,1)-(2,2), then across the x boundary to (0,2) and
        # across the y boundary back to (0,0): one winding along x and y together.
        lattice = square_lattice(3, 3, 1.0)
        relabel(lattice, {square_site(x, y, 3) for x, y in [(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (0, 2)]})
        (cluster,) = lattice.connected_sites("A")
        assert cluster.is_periodically_contiguous(lattice.cell_lengths) == (True, True, False)


def square_site(x, y, a):
    """Number of the site at (x, y) in square_lattice(a, b, spacing)."""
    return x + a * y + 1


def relabel(lattice, a_site_numbers):
    """Label the given sites 'A' and every other site 'B'."""
    for site in lattice.sites:
        site.label = "A" if site.number in a_site_numbers else "B"
    lattice.site_labels = {site.label for site in lattice.sites}
