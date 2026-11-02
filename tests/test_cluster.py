from unittest.mock import Mock

import numpy as np

from lattice_mc.cluster import Cluster
from lattice_mc.lattice_site import Site


class TestCluster:
    """Test for Cluster class"""

    def test_cluster_is_initialised(self):
        sites = [Mock(spec=Site), Mock(spec=Site), Mock(spec=Site), Mock(spec=Site)]
        cluster_sites = sites[0:2]
        neighbour_sites = sites[2:4]
        sites[0].p_neighbours = [sites[1], sites[2], sites[3]]
        sites[1].p_neighbours = [sites[0], sites[2]]
        cluster = Cluster(cluster_sites)
        assert cluster.sites == set(cluster_sites)
        assert cluster.neighbours == set(neighbour_sites)

    def test_cluster_merge(self):
        sites = [Mock(spec=Site), Mock(spec=Site), Mock(spec=Site), Mock(spec=Site)]
        sites[0].p_neighbours = [sites[1], sites[2]]
        sites[2].p_neighbours = [sites[0], sites[3]]
        cluster1 = Cluster([sites[0]])
        cluster2 = Cluster([sites[2]])
        combined_cluster = cluster1.merge(cluster2)
        assert combined_cluster.sites == set([sites[0], sites[2]])
        assert combined_cluster.neighbours == set([sites[1], sites[3]])

    def test_cluster_merge_with_blocked_neighbour(self):
        sites = [Mock(spec=Site), Mock(spec=Site), Mock(spec=Site), Mock(spec=Site)]
        sites[0].p_neighbours = [sites[1], sites[2]]
        sites[2].p_neighbours = [sites[0], sites[3]]
        cluster1 = Cluster([sites[0]])
        cluster2 = Cluster([sites[2]])
        cluster1.neighbours = cluster1.neighbours - set([sites[1]])
        combined_cluster = cluster1.merge(cluster2)
        assert combined_cluster.sites == set([sites[0], sites[2]])
        assert combined_cluster.neighbours == set([sites[3]])

    def test_cluster_is_neighbouring(self):
        sites = [Mock(spec=Site), Mock(spec=Site)]
        sites[0].p_neighbours = [sites[1]]
        sites[1].p_neighbours = [sites[0]]
        clusters = [Cluster([sites[0]]), Cluster([sites[1]])]
        assert clusters[0].is_neighbouring(clusters[1]) is True
        assert clusters[1].is_neighbouring(clusters[0]) is True

    def test_cluster_is_not_neighbouring(self):
        sites = [Mock(spec=Site), Mock(spec=Site)]
        sites[0].p_neighbours = [Mock(spec=Site)]
        sites[1].p_neighbours = [Mock(spec=Site)]
        clusters = [Cluster([sites[0]]), Cluster([sites[1]])]
        assert clusters[0].is_neighbouring(clusters[1]) is False
        assert clusters[1].is_neighbouring(clusters[0]) is False

    def test_cluster_size(self):
        sites = [Mock(spec=Site), Mock(spec=Site)]
        sites[0].p_neighbours = [Mock(spec=Site)]
        sites[1].p_neighbours = [Mock(spec=Site)]
        cluster = Cluster(sites)
        assert cluster.size() == 2

    def test_cluster_edges(self):
        sites = [Mock(spec=Site), Mock(spec=Site), Mock(spec=Site), Mock(spec=Site)]
        sites[0].r = np.array([0.0, 0.0, 1.0])
        sites[1].r = np.array([1.0, 1.0, 0.0])
        sites[2].r = np.array([0.0, 1.0, 1.0])
        sites[3].r = np.array([1.0, 0.0, 0.0])
        sites[0].p_neighbours = [sites[1], sites[2], sites[3]]
        sites[1].p_neighbours = [sites[0], sites[2], sites[3]]
        sites[2].p_neighbours = [sites[0], sites[1], sites[3]]
        sites[3].p_neighbours = [sites[0], sites[1], sites[2]]
        min_x = (sites[0], sites[2])
        max_x = (sites[1], sites[3])
        min_y = (sites[0], sites[3])
        max_y = (sites[1], sites[2])
        min_z = (sites[1], sites[3])
        max_z = (sites[0], sites[2])
        cluster = Cluster(sites)
        for s in cluster.sites_at_edges()[0]:
            assert s in max_x
        for s in cluster.sites_at_edges()[1]:
            assert s in min_x
        for s in cluster.sites_at_edges()[2]:
            assert s in max_y
        for s in cluster.sites_at_edges()[3]:
            assert s in min_y
        for s in cluster.sites_at_edges()[4]:
            assert s in max_z
        for s in cluster.sites_at_edges()[5]:
            assert s in min_z

    def test_cluster_is_periodically_contiguous(self):
        sites = [Mock(spec=Site), Mock(spec=Site), Mock(spec=Site)]
        sites[0].r = np.array([0.0, 0.0, 0.0])
        sites[1].r = np.array([1.0, 0.0, 0.0])
        sites[2].r = np.array([2.0, 0.0, 0.0])
        sites[0].p_neighbours = [sites[1], sites[2]]
        sites[1].p_neighbours = [sites[0], sites[2]]
        sites[2].p_neighbours = [sites[0], sites[1]]
        cluster = Cluster(sites)
        assert cluster.is_periodically_contiguous()[0] is True

    def test_cluster_is_not_periodically_contiguous(self):
        sites = [Mock(spec=Site), Mock(spec=Site), Mock(spec=Site)]
        sites[0].r = np.array([0.0, 0.0, 0.0])
        sites[1].r = np.array([1.0, 0.0, 0.0])
        sites[2].r = np.array([2.0, 0.0, 0.0])
        sites[0].p_neighbours = [sites[1]]
        sites[1].p_neighbours = [sites[0], sites[2]]
        sites[2].p_neighbours = [sites[1]]
        cluster = Cluster(sites)
        assert cluster.is_periodically_contiguous()[0] is False

    def test_cluster_remove_sites_from_neighbours(self):
        sites = [Mock(spec=Site), Mock(spec=Site), Mock(spec=Site)]
        sites[0].p_neighbours = [sites[1], sites[2]]
        sites[1].label = "A"
        sites[2].label = "B"
        cluster = Cluster([sites[0]])
        cluster.remove_sites_from_neighbours("A")
        assert cluster.neighbours == set([sites[2]])
        cluster.remove_sites_from_neighbours(["B"])
        assert cluster.neighbours == set()
