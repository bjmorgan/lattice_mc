import unittest
from unittest.mock import Mock, MagicMock, patch, call, mock_open
from lattice_mc import init_lattice
from lattice_mc.lattice_site import Site
from lattice_mc.lattice import Lattice
import numpy as np

#TODO write integration tests for creating lattices

class InitLatticeTestCase( unittest.TestCase ):
    """Test for Specific Lattice initialisation routines"""

    @patch( 'lattice_mc.lattice_site.Site' )
    @patch( 'lattice_mc.lattice.Lattice' )
    def test_square_lattice( self, mock_lattice, mock_site ):
        mock_site.side_effect = range(2*3)
        mock_lattice.return_value = 'foo'
        a = 2
        b = 3
        spacing = 1.0
        lattice = init_lattice.square_lattice( a, b, spacing )
        expected_site_calls = [ [ 1, np.array([ 0.,  0.,  0.]), [2, 2, 5, 3], 0.0, 'L' ],
                                [ 2, np.array([ 1.,  0.,  0.]), [1, 1, 6, 4], 0.0, 'L' ],
                                [ 3, np.array([ 0.,  1.,  0.]), [4, 4, 1, 5], 0.0, 'L' ],
                                [ 4, np.array([ 1.,  1.,  0.]), [3, 3, 2, 6], 0.0, 'L' ],
                                [ 5, np.array([ 0.,  2.,  0.]), [6, 6, 3, 1], 0.0, 'L' ],
                                [ 6, np.array([ 1.,  2.,  0.]), [5, 5, 4, 2], 0.0, 'L' ] ]
        for c, e in zip( mock_site.mock_calls, expected_site_calls ):
            self.assertEqual( c[1][0], e[0] ) # site number
            np.testing.assert_array_equal( c[1][1], e[1] ) # site coordinates
            self.assertEqual( c[1][2], e[2] ) # neighbour lists
            self.assertEqual( c[1][3], e[3] ) # ?
            self.assertEqual( c[1][4], e[4] ) # site label
        self.assertEqual( mock_lattice.mock_calls[0][1][0], list(range(2*3)) )
        self.assertEqual( lattice, 'foo' )
        np.testing.assert_array_equal( mock_lattice.mock_calls[0][2]['cell_lengths'], np.array( [ 2.0, 3.0, 0.0 ] ) )

    @patch( 'lattice_mc.lattice_site.Site' )
    @patch( 'lattice_mc.lattice.Lattice' )
    def test_honeycomb_lattice( self, mock_lattice, mock_site ):
        a, b = 2, 1
        mock_site.side_effect = range( a*b*4 )
        mock_lattice.return_value = 'foo'
        spacing = 1.0
        lattice = init_lattice.honeycomb_lattice( a, b, spacing )
        expected_site_calls = [ [ 1, np.array([ 0., 0., 0. ]),          [2, 4, 6], 0.0, 'L' ],
                                [ 2, np.array([ 0.8660254, 0.5, 0. ]),  [1, 3, 5], 0.0, 'L' ],
                                [ 3, np.array([ 0.8660254, 1.5, 0. ]),  [2, 4, 8], 0.0, 'L' ],
                                [ 4, np.array([ 0., 2., 0. ]),          [1, 3, 7], 0.0, 'L' ],
                                [ 5, np.array([ 1.73205081, 0., 0. ]),  [2, 6, 8], 0.0, 'L' ],
                                [ 6, np.array([ 2.59807621, 0.5, 0. ]), [5, 1, 7], 0.0, 'L' ],
                                [ 7, np.array([ 2.59807621, 1.5, 0. ]), [8, 6, 4], 0.0, 'L' ],
                                [ 8, np.array([ 1.73205081, 2., 0. ]),  [3, 5, 7], 0.0, 'L' ] ]
        for c, e in zip( mock_site.mock_calls, expected_site_calls ):
            self.assertEqual( c[1][0], e[0] ) # site number
            np.testing.assert_array_almost_equal( c[1][1], e[1] ) # site coordinates
            self.assertEqual( c[1][2].sort(), e[2].sort() ) # neighbour lists
            self.assertEqual( c[1][3], e[3] ) # ?
            self.assertEqual( c[1][4], e[4] ) # site label
        self.assertEqual( mock_lattice.mock_calls[0][1][0], list(range(2*1*4)) )
        self.assertEqual( lattice, 'foo' )
        np.testing.assert_array_almost_equal( mock_lattice.mock_calls[0][2]['cell_lengths'], np.array( [ 3.464102,  3., 0. ] ) )

    @patch( 'lattice_mc.lattice_site.Site' )
    @patch( 'lattice_mc.lattice.Lattice' )
    def test_honeycomb_lattice_periodic_boundaries( self, mock_lattice, mock_site ):
        a, b = 3, 3
        mock_site.side_effect = range( a*b*4 )
        mock_lattice.return_value = 'foo'
        spacing = 1.0
        lattice = init_lattice.honeycomb_lattice( a, b, spacing )	
        expected_neighbours = { 1  : [ 2, 26, 12 ],
                                12 : [ 1, 35, 11 ],
                                35 : [ 36, 34, 12 ],
                                26 : [ 27, 25, 1 ] }
        for c in mock_site.mock_calls:
            site_number = c[1][0]
            if site_number in expected_neighbours:
                site_neighbours = c[1][2]
                self.assertEqual( site_neighbours.sort(), expected_neighbours[ site_number ].sort() )

    @patch( 'lattice_mc.lattice_site.Site' )
    @patch( 'lattice_mc.lattice.Lattice' )
    def test_honeycomb_lattice_alternating_sites( self, mock_lattice, mock_site ):
        a, b, = 3, 3
        mock_site.side_effect = range( a*b*4 )
        mock_lattice.return_value = 'foo'
        spacing = 1.0
        lattice = init_lattice.honeycomb_lattice( a, b, spacing, alternating_sites=True )
        for c in mock_site.mock_calls:
            site_label = c[1][4]
            site_neighbours = c[1][2]
            for n in site_neighbours:
                neighbour_label = mock_site.mock_calls[n-1][1][4]
                self.assertEqual( set( [ site_label, neighbour_label ] ), set( [ 'A', 'B' ] ) )
  
    @patch( 'lattice_mc.lattice_site.Site' )
    @patch( 'lattice_mc.lattice.Lattice' )
    def test_cubic_lattice( self, mock_lattice, mock_site ):
        a, b, c, = 2, 2, 2
        spacing = 1.0
        mock_site.side_effect = range(2*2*2)
        mock_lattice.return_value = 'bar'
        lattice = init_lattice.cubic_lattice( a, b, c, spacing )
        expected_site_calls = [ [ 1, np.array( [ 0., 0., 0. ] ), [ 2, 2, 3, 3, 5, 5 ], 0.0, 'L' ],
                                [ 2, np.array( [ 1., 0., 0. ] ), [ 1, 1, 4, 4, 6, 6 ], 0.0, 'L' ],
                                [ 3, np.array( [ 0., 1., 0. ] ), [ 4, 4, 1, 1, 7, 7 ], 0.0, 'L' ],
                                [ 4, np.array( [ 1., 1., 0. ] ), [ 3, 3, 2, 2, 8, 8 ], 0.0, 'L' ],
                                [ 5, np.array( [ 0., 0., 1. ] ), [ 6, 6, 7, 7, 1, 1 ], 0.0, 'L' ],
                                [ 6, np.array( [ 1., 0., 1. ] ), [ 5, 5, 8, 8, 2, 2 ], 0.0, 'L' ],
                                [ 7, np.array( [ 0., 1., 1. ] ), [ 8, 8, 5, 5, 3, 3 ], 0.0, 'L' ],
                                [ 8, np.array( [ 1., 1., 1. ] ), [ 7, 7, 6, 6, 4, 4 ], 0.0, 'L' ] ]
        for call, e in zip( mock_site.mock_calls, expected_site_calls ):
            self.assertEqual( call[1][0], e[0] ) # site number
            np.testing.assert_array_equal( call[1][1], e[1] ) # site coordinates
            self.assertEqual( call[1][2], e[2] ) # neighbour lists
            self.assertEqual( call[1][3], e[3] ) # ?
            self.assertEqual( call[1][4], e[4] ) # site label
        self.assertEqual( mock_lattice.mock_calls[0][1][0], list(range(2*2*2)) )
        self.assertEqual( lattice, 'bar' )
        np.testing.assert_array_equal( mock_lattice.mock_calls[0][2]['cell_lengths'], np.array( [ a, b, c, ] ) )

    @patch( 'lattice_mc.lattice_site.Site' )
    @patch( 'lattice_mc.lattice.Lattice' )
    def test_lattice_from_file( self, mock_lattice, mock_site ):
        cell_lengths = np.array( [ 2.0, 3.0, 4.0 ] )
        example_file = """1\n
                          site: 2
                          centre: 21.4669 -1.37 6.1334
                          vertices: 1 315 435 649
                          neighbours: 4 5 6
                          label: H"""
        mock_site.return_value = 'site'
        mock_lattice.return_value = 'lattice' 
        with patch( 'builtins.open', mock_open( read_data=example_file), create=True) as m:
            lattice = init_lattice.lattice_from_sites_file( 'filename', cell_lengths )
        site_calls = mock_site.mock_calls[0][1]
        self.assertEqual( site_calls[0], 2 ) # site number
        np.testing.assert_array_equal( site_calls[1], np.array( [ 21.4669, -1.37, 6.1334 ] ) ) # r
        self.assertEqual( site_calls[2], [ 4, 5, 6 ] ) # neighbour list
        self.assertEqual( site_calls[3], 0.0 ) # ?
        self.assertEqual( site_calls[4], 'H' ) # site label 
        self.assertEqual( mock_lattice.mock_calls[0][1][0], [ 'site' ] )
        np.testing.assert_array_equal( mock_lattice.mock_calls[0][2]['cell_lengths'], cell_lengths )
        self.assertEqual( lattice, 'lattice' )
 
    @patch( 'lattice_mc.lattice_site.Site' )
    @patch( 'lattice_mc.lattice.Lattice' )
    def test_lattice_from_file_new( self, mock_lattice, mock_site ):
        cell_lengths = np.array( [ 2.0, 3.0, 4.0 ] )
        example_file = """1\n
                          site: 2
                          centre: 21.4669 -1.37 6.1334
                          vertices: 1 315 435 649
                          neighbours: 4 5 6
                          label: H"""
        mock_site.return_value = 'site'
        mock_lattice.return_value = 'lattice'
        with patch( 'builtins.open', mock_open( read_data=example_file), create=True) as m:
            lattice = init_lattice.lattice_from_sites_file_new( 'filename', cell_lengths )
        site_calls = mock_site.mock_calls[0][1]
        self.assertEqual( site_calls[0], 2 ) # site number
        np.testing.assert_array_equal( site_calls[1], np.array( [ 21.4669, -1.37, 6.1334 ] ) ) # r
        self.assertEqual( site_calls[2], [ 4, 5, 6 ] ) # neighbour list
        self.assertEqual( site_calls[3], 0.0 ) # ?
        self.assertEqual( site_calls[4], 'H' ) # site label 
        self.assertEqual( mock_lattice.mock_calls[0][1][0], [ 'site' ] )
        np.testing.assert_array_equal( mock_lattice.mock_calls[0][2]['cell_lengths'], cell_lengths )
        self.assertEqual( lattice, 'lattice' ) 

    @patch( 'lattice_mc.lattice_site.Site' )
    @patch( 'lattice_mc.lattice.Lattice' )
    def test_lattice_from_file_with_energy( self, mock_lattice, mock_site ):
        cell_lengths = np.array( [ 2.0, 3.0, 4.0 ] )
        example_file = """1\n
                          site: 2
                          centre: 21.4669 -1.37 6.1334
                          vertices: 1 315 435 649
                          neighbours: 4 5 6
                          energy: -1.0
                          label: H"""
        mock_site.return_value = 'site'
        mock_lattice.return_value = 'lattice'
        with patch( 'builtins.open', mock_open( read_data=example_file), create=True) as m:
            lattice = init_lattice.lattice_from_sites_file_new( 'filename', cell_lengths )
        site_calls = mock_site.mock_calls[0][1]
        self.assertEqual( site_calls[0], 2 ) # site number
        np.testing.assert_array_equal( site_calls[1], np.array( [ 21.4669, -1.37, 6.1334 ] ) ) # r
        self.assertEqual( site_calls[2], [ 4, 5, 6 ] ) # neighbour list
        self.assertEqual( site_calls[3], -1.0 ) # site energy
        self.assertEqual( site_calls[4], 'H' ) # site label 
        self.assertEqual( mock_lattice.mock_calls[0][1][0], [ 'site' ] )
        np.testing.assert_array_equal( mock_lattice.mock_calls[0][2]['cell_lengths'], cell_lengths )
        self.assertEqual( lattice, 'lattice' )

if __name__ == '__main__':
    unittest.main()
