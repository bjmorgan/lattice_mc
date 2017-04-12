import unittest
from lattice_mc.lattice import Lattice
from lattice_mc.lattice_site import Site
from unittest.mock import Mock, patch
import numpy as np

class LatticeIntegrationTestCase( unittest.TestCase ):
    """Tests for integration of Lattice class with the rest of the code"""

    @patch( 'lattice_mc.lattice.Lattice.site_with_id' )
    @patch( 'lattice_mc.lattice.Lattice.initialise_site_lookup_table' )
    @patch( 'lattice_mc.lattice.Lattice.enforce_periodic_boundary_conditions' )
    def setUp( self, pbc, islt, site_id ):
        site_id.side_effect = ( 1,2,3,4,5,6,7,8 )
        self.site_id = site_id
        site_labels = [ 'A', 'A', 'A', 'A', 'A' ]
        site_neighbours = [ [ 2, 3 ], [ 1, 3 ], [ 1, 2 ], [ 5 ], [ 4 ] ]
        self.mock_sites = [ Mock( spec=Site, label = l, neighbours = n ) for l, n in zip( site_labels, site_neighbours ) ]
        self.cell_lengths = np.array( [ 7.0, 8.0, 9.0 ] )
        self.lattice = Lattice( self.mock_sites, self.cell_lengths )

    def test_connected_sites( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].p_neighbours = [ sites[1] ]
        sites[1].p_neighbours = [ sites[0] ]
        for s in sites:
            s.label = 'A'
        self.lattice.sites = sites
        connected_sites = self.lattice.connected_sites() 
        self.assertEqual( len( connected_sites ), 1 ) 
        for s in sites:
            self.assertEqual( s in connected_sites[0].sites, True )
        self.assertEqual( list( connected_sites[0].neighbours ), [] )
  
    def test_connected_sites_2( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].p_neighbours = [ sites[1] ]
        sites[1].p_neighbours = [ sites[0] ]
        sites[2].p_neighbours = [ sites[3] ]
        sites[3].p_neighbours = [ sites[2] ]
        for s in sites:
            s.label = 'A'
        self.lattice.sites = sites
        connected_sites = self.lattice.connected_sites()
        self.assertEqual( len( connected_sites ), 2 )
        for c in connected_sites:
            if sites[0] in c.sites:
                self.assertEqual( sites[1] in c.sites, True )
            else:
                self.assertEqual( sites[2] in c.sites, True )
                self.assertEqual( sites[3] in c.sites, True )

    def test_connected_sites_3( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].p_neighbours = [ sites[1], sites[2] ]
        sites[1].p_neighbours = [ sites[0], sites[2] ]
        sites[2].p_neighbours = [ sites[0], sites[1] ]
        sites[3].p_neighbours = [ sites[4] ]
        sites[4].p_neighbours = [ sites[3] ]
        self.lattice.sites = sites
        for s in sites:
            s.label = 'A'
        connected_sites = self.lattice.connected_sites()
        self.assertEqual( len( connected_sites ), 2 )
        for c in connected_sites:
            if sites[0] in c.sites:
                self.assertEqual( sites[1] in c.sites, True )
                self.assertEqual( sites[2] in c.sites, True )
                self.assertEqual( c.size(), 3 )
            else:
                self.assertEqual( sites[3] in c.sites, True )
                self.assertEqual( sites[4] in c.sites, True )
                self.assertEqual( c.size(), 2 )

if __name__ == '__main__':
    unittest.main()
