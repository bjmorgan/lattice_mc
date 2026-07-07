import unittest
from unittest.mock import Mock, patch

import numpy as np

from lattice_mc.lattice import Lattice
from lattice_mc.lattice_site import Site


class LatticeIntegrationTestCase(unittest.TestCase):
    """Tests for integration of Lattice class with the rest of the code"""

    @patch("lattice_mc.lattice.Lattice.site_with_id")
    @patch("lattice_mc.lattice.Lattice.initialise_site_lookup_table")
    @patch("lattice_mc.lattice.Lattice.enforce_periodic_boundary_conditions")
    def setUp(self, pbc, islt, site_id):
        site_id.side_effect = (1, 2, 3, 4, 5, 6, 7, 8)
        self.site_id = site_id
        site_labels = ["A", "A", "A", "A", "A"]
        site_neighbours = [[2, 3], [1, 3], [1, 2], [5], [4]]
        self.mock_sites = [Mock(spec=Site, label=label, neighbours=n) for label, n in zip(site_labels, site_neighbours)]
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
        self.assertEqual(len(connected_sites), 1)
        for s in sites:
            self.assertEqual(s in connected_sites[0].sites, True)
        self.assertEqual(list(connected_sites[0].neighbours), [])

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
        self.assertEqual(len(connected_sites), 2)
        for c in connected_sites:
            if sites[0] in c.sites:
                self.assertEqual(sites[1] in c.sites, True)
            else:
                self.assertEqual(sites[2] in c.sites, True)
                self.assertEqual(sites[3] in c.sites, True)

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
        self.assertEqual(len(connected_sites), 2)
        for c in connected_sites:
            if sites[0] in c.sites:
                self.assertEqual(sites[1] in c.sites, True)
                self.assertEqual(sites[2] in c.sites, True)
                self.assertEqual(c.size(), 3)
            else:
                self.assertEqual(sites[3] in c.sites, True)
                self.assertEqual(sites[4] in c.sites, True)
                self.assertEqual(c.size(), 2)


class DetachedSitesIntegrationTestCase(unittest.TestCase):
    """Regression tests for Lattice.detached_sites using real Site objects.

    The old coordinate-extremes heuristic in Cluster.is_periodically_contiguous
    reported fully connected 4-site clusters as periodically contiguous on all
    axes, so a lattice of disconnected cliques wrongly reported no detached
    sites.

    The winding-ring test guards the complementary direction: genuinely
    winding clusters must not be reported as detached.
    """

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
        self.assertEqual({s.number for s in lattice.detached_sites()}, {s.number for s in sites})

    def test_winding_ring_has_no_detached_sites(self):
        sites = [
            Site(number=1, coordinates=np.array([0.0, 0.0, 0.0]), neighbours=[2, 3], energy=0.0, label="A"),
            Site(number=2, coordinates=np.array([1.0, 0.0, 0.0]), neighbours=[1, 3], energy=0.0, label="A"),
            Site(number=3, coordinates=np.array([2.0, 0.0, 0.0]), neighbours=[1, 2], energy=0.0, label="A"),
        ]
        lattice = Lattice(sites, cell_lengths=np.array([3.0, 10.0, 10.0]))
        self.assertEqual(lattice.detached_sites(), [])


if __name__ == "__main__":
    unittest.main()
