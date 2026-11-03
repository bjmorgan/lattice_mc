import numpy as np
import pytest

from lattice_mc.init_lattice import honeycomb_lattice, square_lattice
from lattice_mc.lattice import Lattice
from lattice_mc.lattice_site import Site


def lattice_with_neighbours(neighbours):
    """A lattice of 'A' sites at the origin, with the given zero-based neighbour indices."""
    sites = [
        Site(number=index + 1, coordinates=np.zeros(3), neighbours=[n + 1 for n in nn], energy=0.0, label="A")
        for index, nn in enumerate(neighbours)
    ]
    return Lattice(sites, cell_lengths=np.array([10.0, 10.0, 10.0]))


class TestConnectedSites:
    """Tests for Lattice.connected_sites using real Site and Lattice objects."""

    def test_connected_pair_is_one_cluster(self):
        lattice = lattice_with_neighbours([[1], [0]])
        (cluster,) = lattice.connected_sites()
        assert cluster.sites == set(lattice.sites)
        assert cluster.neighbours == set()

    def test_two_separate_pairs_are_two_clusters(self):
        lattice = lattice_with_neighbours([[1], [0], [3], [2]])
        clusters = lattice.connected_sites()
        assert sorted(sorted(s.number for s in c.sites) for c in clusters) == [[1, 2], [3, 4]]
        assert all(c.neighbours == set() for c in clusters)

    def test_triangle_and_pair_are_clusters_of_three_and_two(self):
        lattice = lattice_with_neighbours([[1, 2], [0, 2], [0, 1], [4], [3]])
        assert sorted(c.size() for c in lattice.connected_sites()) == [2, 3]

    def test_site_labels_restrict_the_search(self):
        lattice = lattice_with_neighbours([[1], [0, 2], [1]])
        lattice.sites[1].label = "B"
        clusters = lattice.connected_sites("A")
        assert sorted(sorted(s.number for s in c.sites) for c in clusters) == [[1], [3]]


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
        "build_lattice",
        [lambda: square_lattice(3, 3, 1.0), lambda: honeycomb_lattice(2, 2, 1.0)],
        ids=["square", "honeycomb"],
    )
    def test_fully_occupied_two_dimensional_lattice_winds_along_x_and_y(self, build_lattice):
        lattice = build_lattice()
        (cluster,) = lattice.connected_sites()
        assert cluster.is_periodically_contiguous() == (True, True, False)
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
        assert cluster.is_periodically_contiguous() == (True, True, False)


def square_site(x, y, a):
    """Number of the site at (x, y) in square_lattice(a, b, spacing)."""
    return x + a * y + 1


def relabel(lattice, a_site_numbers):
    """Label the given sites 'A' and every other site 'B'."""
    for site in lattice.sites:
        site.label = "A" if site.number in a_site_numbers else "B"
