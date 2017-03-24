import unittest
from lattice_mc.cluster import Cluster
from lattice_mc.lattice_site import Site
from unittest.mock import Mock, patch
import numpy as np

class ClusterTestCase( unittest.TestCase ):
    """Test for Cluster class"""

    def test_cluster_is_initialised( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ) ]
        cluster_sites = sites[0:2]
        neighbour_sites = sites[2:4]
        sites[0].p_neighbours = [ sites[1], sites[2], sites[3] ]
        sites[1].p_neighbours = [ sites[0], sites[2] ]
        cluster = Cluster( cluster_sites )
        self.assertEqual( cluster.sites, set( cluster_sites ) )
        self.assertEqual( cluster.neighbours, set( neighbour_sites ) )
    
    def test_cluster_merge( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].p_neighbours = [ sites[1], sites[2] ]
        sites[2].p_neighbours = [ sites[0], sites[3] ]
        cluster1 = Cluster( [ sites[0] ] )
        cluster2 = Cluster( [ sites[2] ] )
        combined_cluster = cluster1.merge( cluster2 )
        self.assertEqual( combined_cluster.sites, set( [ sites[0], sites[2] ] ) )
        self.assertEqual( combined_cluster.neighbours, set( [ sites[1], sites[3] ] ) )

    def test_cluster_is_neighbouring( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].p_neighbours = [ sites[1] ]
        sites[1].p_neighbours = [ sites[0] ]
        clusters = [ Cluster( [ sites[0] ] ), Cluster( [ sites[1] ] ) ]
        self.assertEqual( clusters[0].is_neighbouring( clusters[1] ), True )
        self.assertEqual( clusters[1].is_neighbouring( clusters[0] ), True )

    def test_cluster_is_not_neighbouring( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].p_neighbours = [ Mock( spec=Site ) ]
        sites[1].p_neighbours = [ Mock( spec=Site ) ]
        clusters = [ Cluster( [ sites[0] ] ), Cluster( [ sites[1] ] ) ]
        self.assertEqual( clusters[0].is_neighbouring( clusters[1] ), False )
        self.assertEqual( clusters[1].is_neighbouring( clusters[0] ), False )

    def test_cluster_size( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].p_neighbours = [ Mock( spec=Site ) ]
        sites[1].p_neighbours = [ Mock( spec=Site ) ]
        cluster = Cluster( sites )
        self.assertEqual( cluster.size(), 2 )

    def test_cluster_edges( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].r = np.array( [ 0.0, 0.0, 1.0 ] )
        sites[1].r = np.array( [ 1.0, 1.0, 0.0 ] )
        sites[2].r = np.array( [ 0.0, 1.0, 1.0 ] )
        sites[3].r = np.array( [ 1.0, 0.0, 0.0 ] )
        sites[0].p_neighbours = [ sites[1], sites[2], sites[3] ]
        sites[1].p_neighbours = [ sites[0], sites[2], sites[3] ]
        sites[2].p_neighbours = [ sites[0], sites[1], sites[3] ]
        sites[3].p_neighbours = [ sites[0], sites[1], sites[2] ]
        min_x = ( sites[0], sites[2] )
        max_x = ( sites[1], sites[3] )
        min_y = ( sites[0], sites[3] )
        max_y = ( sites[1], sites[2] )
        min_z = ( sites[1], sites[3] )
        max_z = ( sites[0], sites[2] )
        cluster = Cluster( sites )
        for s in cluster.sites_at_edges()[0]: 
            self.assertEqual( s in min_x, True )
        for s in cluster.sites_at_edges()[1]:
            self.assertEqual( s in max_x, True )
        for s in cluster.sites_at_edges()[2]:
            self.assertEqual( s in min_y, True )
        for s in cluster.sites_at_edges()[3]:
            self.assertEqual( s in max_y, True )
        for s in cluster.sites_at_edges()[4]:
            self.assertEqual( s in min_z, True )
        for s in cluster.sites_at_edges()[5]:
            self.assertEqual( s in max_z, True )

    def test_cluster_is_periodically_contiguous( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].r = np.array( [ 0.0, 0.0, 0.0 ] )
        sites[1].r = np.array( [ 1.0, 0.0, 0.0 ] )
        sites[2].r = np.array( [ 2.0, 0.0, 0.0 ] )
        sites[0].p_neighbours = [ sites[1], sites[2] ]
        sites[1].p_neighbours = [ sites[0], sites[2] ]
        sites[2].p_neighbours = [ sites[0], sites[1] ]
        cluster = Cluster( sites )
        self.assertEqual( cluster.is_periodically_contiguous()[0], True )

    def test_cluster_is_not_periodically_contiguous( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].r = np.array( [ 0.0, 0.0, 0.0 ] )
        sites[1].r = np.array( [ 1.0, 0.0, 0.0 ] )
        sites[2].r = np.array( [ 2.0, 0.0, 0.0 ] )
        sites[0].p_neighbours = [ sites[1] ]
        sites[1].p_neighbours = [ sites[0], sites[2] ]
        sites[2].p_neighbours = [ sites[1] ]
        cluster = Cluster( sites )
        self.assertEqual( cluster.is_periodically_contiguous()[0], False )

if __name__ == '__main__':
    unittest.main()
