import numpy as np
import pytest

from lattice_mc.cluster import Cluster
from lattice_mc.lattice import Lattice
from lattice_mc.lattice_site import Site


def make_lattice(positions, neighbours, cell_lengths=(10.0, 10.0, 10.0), labels=None):
    """Build a Lattice from site positions and zero-based neighbour indices."""
    labels = labels or ["A"] * len(positions)
    sites = [
        Site(index + 1, np.array(position, dtype=float), [n + 1 for n in neighbours[index]], 0.0, labels[index])
        for index, position in enumerate(positions)
    ]
    return Lattice(sites, np.array(cell_lengths))


def fully_connected(positions, cell_lengths=(10.0, 10.0, 10.0)):
    """A lattice in which every site neighbours every other site."""
    n = len(positions)
    return make_lattice(positions, [[j for j in range(n) if j != i] for i in range(n)], cell_lengths)


def chain(n):
    """An open chain of n sites along x, spaced 1 apart, in a large cell."""
    neighbours = [[j for j in (i - 1, i + 1) if 0 <= j < n] for i in range(n)]
    return make_lattice([[float(i), 0.0, 0.0] for i in range(n)], neighbours)


class TestCluster:
    """Test for Cluster class"""

    def test_cluster_is_initialised(self):
        lattice = make_lattice([[0, 0, 0]] * 4, [[1, 2, 3], [0, 2], [0, 1], [0]])
        sites = lattice.sites
        cluster = Cluster(sites[0:2], lattice)
        assert cluster.sites == set(sites[0:2])
        assert cluster.neighbours == set(sites[2:4])

    def test_cluster_merge(self):
        lattice = chain(4)
        sites = lattice.sites
        combined_cluster = Cluster([sites[1]], lattice).merge(Cluster([sites[2]], lattice))
        assert combined_cluster.sites == {sites[1], sites[2]}
        assert combined_cluster.neighbours == {sites[0], sites[3]}

    def test_cluster_merge_with_blocked_neighbour(self):
        lattice = chain(4)
        sites = lattice.sites
        cluster1 = Cluster([sites[1]], lattice)
        cluster1.neighbours = cluster1.neighbours - {sites[0]}
        combined_cluster = cluster1.merge(Cluster([sites[2]], lattice))
        assert combined_cluster.sites == {sites[1], sites[2]}
        assert combined_cluster.neighbours == {sites[3]}

    def test_cluster_is_neighbouring(self):
        lattice = chain(2)
        clusters = [Cluster([site], lattice) for site in lattice.sites]
        assert clusters[0].is_neighbouring(clusters[1]) is True
        assert clusters[1].is_neighbouring(clusters[0]) is True

    def test_cluster_is_not_neighbouring(self):
        lattice = chain(3)
        clusters = [Cluster([lattice.sites[0]], lattice), Cluster([lattice.sites[2]], lattice)]
        assert clusters[0].is_neighbouring(clusters[1]) is False
        assert clusters[1].is_neighbouring(clusters[0]) is False

    def test_cluster_size(self):
        lattice = chain(3)
        assert Cluster(lattice.sites[0:2], lattice).size() == 2

    @pytest.mark.parametrize("axis", [0, 1, 2])
    def test_cluster_is_periodically_contiguous_for_a_winding_ring(self, axis):
        positions = np.zeros((3, 3))
        positions[:, axis] = [0.0, 1.0, 2.0]
        cell_lengths = np.full(3, 10.0)
        cell_lengths[axis] = 3.0
        lattice = fully_connected(positions, cell_lengths)
        expected = [False, False, False]
        expected[axis] = True
        assert Cluster(lattice.sites, lattice).is_periodically_contiguous(cell_lengths) == tuple(expected)

    @pytest.mark.filterwarnings("error::RuntimeWarning")
    def test_cluster_ignores_an_axis_with_a_zero_cell_length(self):
        # Regression: a zero z cell length (as for two-dimensional lattices)
        # gave 0/0 in the bond wrap calculation, and the NaN cast to int is
        # platform dependent. Warnings are errors so this fails on every platform.
        cell_lengths = np.array([3.0, 10.0, 0.0])
        lattice = fully_connected([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], cell_lengths)
        assert Cluster(lattice.sites, lattice).is_periodically_contiguous(cell_lengths) == (True, False, False)

    def test_cluster_ignores_a_zero_length_axis_with_buckled_sites(self):
        # Sites that are not all at z = 0 give a non-zero z bond component, which
        # divided by a zero cell length is infinite rather than NaN.
        cell_lengths = np.array([3.0, 10.0, 0.0])
        lattice = fully_connected([[0.0, 0.0, 0.0], [1.0, 0.0, 0.3], [2.0, 0.0, 0.0]], cell_lengths)
        with np.errstate(divide="ignore", invalid="ignore"):
            assert Cluster(lattice.sites, lattice).is_periodically_contiguous(cell_lengths) == (True, False, False)

    def test_is_periodically_contiguous_raises_ValueError_for_an_empty_cluster(self):
        with pytest.raises(ValueError):
            Cluster([], chain(1)).is_periodically_contiguous(np.array([3.0, 3.0, 3.0]))

    def test_cluster_is_not_periodically_contiguous_for_a_fully_connected_cluster(self):
        # Regression: a 4-clique has neighbouring coordinate-extreme sites on
        # every axis, but does not wind round the cell.
        lattice = fully_connected([[0.0, 0.0, 1.0], [1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 0.0]])
        cluster = Cluster(lattice.sites, lattice)
        assert cluster.is_periodically_contiguous(lattice.cell_lengths) == (False, False, False)

    def test_cluster_is_not_periodically_contiguous_when_straddling_a_boundary(self):
        # Regression: a bond crossing the wrapped boundary does not imply
        # that the cluster winds round the cell.
        lattice = fully_connected([[0.5, 5.0, 5.0], [9.5, 5.0, 5.0]])
        cluster = Cluster(lattice.sites, lattice)
        assert cluster.is_periodically_contiguous(lattice.cell_lengths) == (False, False, False)

    def test_cluster_is_not_periodically_contiguous_for_an_open_chain(self):
        # An open chain spans the cell but does not reconnect through the
        # periodic boundary.
        cell_lengths = np.array([3.0, 10.0, 10.0])
        lattice = make_lattice([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], [[1], [0, 2], [1]], cell_lengths)
        assert Cluster(lattice.sites, lattice).is_periodically_contiguous(cell_lengths) == (False, False, False)

    def test_cluster_is_periodically_contiguous_ignores_neighbours_outside_the_cluster(self):
        # The ring winds along x, but only through the third site, which is not in the cluster.
        cell_lengths = np.array([3.0, 10.0, 10.0])
        lattice = fully_connected([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], cell_lengths)
        cluster = Cluster(lattice.sites[0:2], lattice)
        assert cluster.is_periodically_contiguous(cell_lengths) == (False, False, False)

    def test_cluster_is_not_periodically_contiguous_for_a_single_site(self):
        lattice = chain(1)
        cluster = Cluster(lattice.sites, lattice)
        assert cluster.is_periodically_contiguous(lattice.cell_lengths) == (False, False, False)

    def test_cluster_remove_sites_from_neighbours(self):
        lattice = make_lattice([[0, 0, 0]] * 3, [[1, 2], [0], [0]], labels=["A", "A", "B"])
        sites = lattice.sites
        cluster = Cluster([sites[0]], lattice)
        cluster.remove_sites_from_neighbours("A")
        assert cluster.neighbours == {sites[2]}
        cluster.remove_sites_from_neighbours(["B"])
        assert cluster.neighbours == set()
