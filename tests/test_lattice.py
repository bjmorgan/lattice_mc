import unittest
from lattice_mc.lattice import Lattice
from lattice_mc.lattice_site import Site
from lattice_mc.atom import Atom
from lattice_mc.jump import Jump
from lattice_mc.transitions import Transitions
from unittest.mock import Mock, patch, call
import numpy as np

class LatticeTestCase( unittest.TestCase ):
    """Tests for Lattice class"""

    @patch( 'lattice_mc.lattice.Lattice.site_with_id' )
    @patch( 'lattice_mc.lattice.Lattice.initialise_site_lookup_table' )
    @patch( 'lattice_mc.lattice.Lattice.enforce_periodic_boundary_conditions' )
    def setUp( self, pbc, islt, site_id ):
        site_id.side_effect = ( 1,2,3,4,5,6,7,8 )
        self.site_id = site_id
        site_labels = [ 'A', 'B', 'A', 'B', 'C' ]
        site_neighbours = [ [ 2, 3 ], [ 1, 3 ], [ 1, 2 ], [ 5 ], [ 4 ] ]
        self.mock_sites = [ Mock( spec=Site, label = l, neighbours = n ) for l, n in zip( site_labels, site_neighbours ) ]
        self.cell_lengths = np.array( [ 7.0, 8.0, 9.0 ] )
        self.lattice = Lattice( self.mock_sites, self.cell_lengths )

    def test_lattice_is_initialised( self ):
        np.testing.assert_array_equal( self.lattice.cell_lengths, self.cell_lengths )
        self.assertEqual( self.lattice.sites, self.mock_sites )
        self.assertEqual( self.lattice.number_of_sites, 5 )
        self.assertEqual( self.lattice.site_labels, { 'A', 'B', 'C' } )
        self.assertEqual( self.lattice.site_populations, { 'A':2, 'B':2, 'C':1 } )
        self.site_id.assert_has_calls( [ call(2), call(3), call(1), call(3), call(1), call(2), call(5), call(4) ] )
        self.assertEqual( self.mock_sites[0].p_neighbours, [ 1, 2 ] )        
        self.assertEqual( self.mock_sites[1].p_neighbours, [ 3, 4 ] )        
        self.assertEqual( self.mock_sites[2].p_neighbours, [ 5, 6 ] )        
        self.assertEqual( self.mock_sites[3].p_neighbours, [ 7 ] )        
        self.assertEqual( self.mock_sites[4].p_neighbours, [ 8 ] )        

    def test_enforce_periodice_boundary_consitions( self ):
        site_coordinates = [ np.array( [ -1.0, 11.0, 3.0 ] ), np.array( [ 2.0, -3.0, 12.0 ] ) ]
        mock_sites = [ Mock( spec=Site, r = r ) for r in site_coordinates ]
        self.lattice.cell_lengths = np.array( [ 10.0, 10.0, 10.0 ] )
        self.lattice.sites = mock_sites
        self.lattice.enforce_periodic_boundary_conditions()
        np.testing.assert_array_equal( self.lattice.sites[0].r, np.array( [ 9.0, 1.0, 3.0 ] ) )
        np.testing.assert_array_equal( self.lattice.sites[1].r, np.array( [ 2.0, 7.0, 2.0 ] ) )

    def test_reset( self ):
        self.lattice.reset()
        self.assertEqual( self.lattice.time, 0.0 )
        self.assertEqual( self.lattice.sites[0].time_occupied, 0.0 )

    def test_initialised_site_lookup_table( self ):
        for i, site in enumerate( self.lattice.sites ):
            site.number = i
        self.lattice.initialise_site_lookup_table()
        self.assertEqual( self.lattice.site_lookup[0], self.mock_sites[0] )

    def test_site_with_id( self ):
        self.lattice.site_lookup = [ 'foo', 'bar' ]
        self.assertEqual( self.lattice.site_with_id( 1 ), 'bar' )

    def test_vacant_sites( self ):
        occupied = [ True, False, True, False, True ]
        for o, site in zip( occupied, self.lattice.sites ):
            site.is_occupied = o
        self.assertEqual( list( self.lattice.vacant_sites() ), self.mock_sites[ 1:5:2 ] )

    def test_occupied_sites( self ):
        occupied = [ True, False, True, False, True ]
        for o, site in zip( occupied, self.lattice.sites ):
            site.is_occupied = o
        self.assertEqual( list( self.lattice.occupied_sites() ), self.mock_sites[ 0:5:2 ] )

    def test_vacant_site_numbers( self ):
        occupied = [ True, False, True, False, True ] 
        numbers = [ 4, 5, 6, 7, 8 ]
        for n, o, site in zip( numbers, occupied, self.lattice.sites ):
            site.is_occupied = o
            site.number = n
        self.assertEqual( self.lattice.vacant_site_numbers(), [ 5, 7 ] )

    def test_occupied_site_numbers( self ):
        occupied = [ True, False, True, False, True ]
        numbers = [ 4, 5, 6, 7, 8 ]
        for n, o, site in zip( numbers, occupied, self.lattice.sites ):
            site.is_occupied = o
            site.number = n
        self.assertEqual( self.lattice.occupied_site_numbers(), [ 4, 6, 8 ] )

    @patch( 'lattice_mc.jump.Jump' )
    def test_potential_jumps_if_lattice_is_mostly_empty( self, mock_Jump ): # lattice is mostly vacant
        jumps = [ Mock( spec=Jump ), Mock( spec=Jump ) ]
        mock_Jump.side_effect = jumps
        self.lattice.number_of_occupied_sites = 1 
        self.lattice_number_of_sites = 4
        site = Mock( spec=Site )
        site.neighbours = [ 2, 3 ]
        occupied_sites = [ site ]
        unoccupied_sites = [ Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ) ]
        self.lattice.nn_energy = 'A'
        self.lattice.cn_energies = 'B'
        self.lattice.jump_lookup_table = 'C'
        for s in unoccupied_sites:
            s.is_occupied = False
        with patch( 'lattice_mc.lattice.Lattice.occupied_sites' ) as mock_occupied_sites:
            with patch( 'lattice_mc.lattice.Lattice.site_with_id' ) as mock_site_with_id:
                mock_occupied_sites.return_value = occupied_sites
                mock_site_with_id.side_effect = unoccupied_sites[:2]
                potential_jumps = self.lattice.potential_jumps()
                self.assertEqual( potential_jumps, jumps )
                self.assertEqual( mock_Jump.mock_calls[0][1], ( site, unoccupied_sites[0], 'A', 'B', 'C' ) )
                self.assertEqual( mock_Jump.mock_calls[1][1], ( site, unoccupied_sites[1], 'A', 'B', 'C' ) )
                mock_site_with_id.assert_has_calls( [ call(2), call(3) ] )

    @patch( 'lattice_mc.jump.Jump' )
    def test_potential_jumps_if_lattice_is_mostly_filled( self, mock_Jump ): # lattice is mostly occupied
        mock_Jump.side_effect = [ 'jump1', 'jump2' ]
        self.lattice.number_of_occupied_sites = 3
        self.lattice_number_of_sites = 4
        site = Mock( spec=Site )
        site.neighbours = [ 2, 3 ]
        vacant_sites = [ site ]
        occupied_sites = [ Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ) ]
        for s in occupied_sites:
            s.is_occupied = True
        self.lattice.nn_energy = 'A'
        self.lattice.cn_energies = 'B'
        self.lattice.jump_lookup_table = 'C'
        with patch( 'lattice_mc.lattice.Lattice.vacant_sites' ) as mock_vacant_sites:
            with patch( 'lattice_mc.lattice.Lattice.site_with_id' ) as mock_site_with_id:
                mock_vacant_sites.return_value = vacant_sites
                mock_site_with_id.side_effect = occupied_sites[:2]
                jumps = self.lattice.potential_jumps()
                self.assertEqual( jumps, [ 'jump1', 'jump2' ] )
                self.assertEqual( mock_Jump.mock_calls[0][1], ( occupied_sites[0], site, 'A', 'B', 'C' ) )
                self.assertEqual( mock_Jump.mock_calls[1][1], ( occupied_sites[1], site, 'A', 'B', 'C' ) )
                mock_site_with_id.assert_has_calls( [ call(2), call(3) ] )

    def test_update( self ):
        atom = Mock( spec=Atom ) # the jumping atom
        atom.number = 28
        atom.number_of_hops = 4
        atom.dr = np.array( [ 1.0, 2.0, 3.0 ] )
        atom.summed_dr2 = 1.3
        jump = Mock( spec=Jump ) # the jump selected for the update
        jump.initial_site = Mock( spec=Site ) # atom jumps from this site
        jump.final_site   = Mock( spec=Site ) # atom jumps to this site
        jump.initial_site.atom = atom
        jump.dr = Mock( return_value=np.array( [ 2.0, 3.0, 4.0 ] ) )
        self.lattice.cell_lengths = np.array( [ 1.0, 2.0, 3.0 ] )
        self.lattice.update( jump )
        jump.dr.assert_called_with( self.lattice.cell_lengths )
        self.assertEqual( jump.final_site.occupation, atom.number )
        self.assertEqual( jump.initial_site.occupation, 0 )
        self.assertEqual( jump.initial_site.atom, None )
        self.assertEqual( jump.final_site.atom, atom )
        self.assertEqual( jump.final_site.is_occupied, True )
        self.assertEqual( jump.initial_site.is_occupied, False )
        self.assertEqual( atom.site, jump.final_site )
        self.assertEqual( atom.number_of_hops, 5 )
        np.testing.assert_array_equal( atom.dr, np.array( [ 3.0, 5.0, 7.0 ] ) )
        self.assertEqual( atom.summed_dr2, 30.3 )

    def test_populate_sites_raises_ValueError_with_too_many_atoms( self ):
        with self.assertRaises( ValueError ):
            self.lattice.populate_sites( number_of_atoms=10 )

    def test_populate_sites( self ):
        number_of_atoms = 2
        with patch( 'random.sample' ) as mock_random_sample:
            with patch( 'lattice_mc.atom.Atom' ) as mock_Atom:
                mock_sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
                mock_random_sample.return_value = mock_sites
                mock_Atom.return_value = [ Mock( spec=Atom ), Mock( spec=Atom ) ]
                self.lattice.populate_sites( 2 )
                self.assertEqual( mock_Atom.mock_calls[0][2]['initial_site'], mock_sites[0] ) 
                self.assertEqual( mock_Atom.mock_calls[1][2]['initial_site'], mock_sites[1] ) 
                self.assertEqual( self.lattice.number_of_occupied_sites, number_of_atoms )           

    @patch( 'lattice_mc.transitions.Transitions' )
    def test_jump( self, mock_Transitions ):
        potential_jumps = [ Mock( spec=Jump ), Mock( spec=Jump ) ] 
        self.lattice.potential_jumps = Mock( return_value=potential_jumps )
        selected_jump = Mock( spec=Jump )
        mock_transitions = Mock( spec=Transitions )
        mock_transitions.random = Mock( return_value=selected_jump )
        mock_transitions.time_to_jump = Mock( return_value=5.0 )
        mock_Transitions.return_value = mock_transitions
        self.lattice.time = 2.0
        self.lattice.update_site_occupation_times = Mock()
        self.lattice.update = Mock()
        self.lattice.jump()
        mock_Transitions.assert_called_with( potential_jumps )
        self.lattice.update.assert_called_with( selected_jump )
        self.lattice.update_site_occupation_times.assert_called_with( 5.0 )
        self.assertEqual( self.lattice.time, 2.0 + 5.0 )

    def test_update_site_occupation_times( self ):
        occupied_sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        occupied_sites[0].time_occupied = 2.0
        occupied_sites[1].time_occupied = 3.2
        self.lattice.occupied_sites = Mock( return_value=occupied_sites )
        self.lattice.update_site_occupation_times( 4.2 )
        self.assertEqual( occupied_sites[0].time_occupied, 2.0 + 4.2 )
        self.assertEqual( occupied_sites[1].time_occupied, 3.2 + 4.2 )
    
    def test_site_occupation_statistics( self ):
        self.lattice.site_labels = [ 'A', 'B' ]
        sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].time_occupied = 2.0
        sites[1].time_occupied = 3.2
        sites[0].label = 'A'
        sites[1].label = 'B'
        self.lattice.sites = sites
        self.lattice.time = 4.0
        occupation_stats = self.lattice.site_occupation_statistics()
        self.assertEqual( occupation_stats, { 'A' : 0.5, 'B' : 0.8 } )
   
    def test_set_site_energies( self ):
        site_energies = { 'A' : 0.2, 'B' : 0.4 }
        sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].label = 'A'
        sites[1].label = 'B'
        self.lattice.sites = sites
        self.lattice.set_site_energies( site_energies )
        self.assertEqual( sites[0].energy, 0.2 )
        self.assertEqual( sites[1].energy, 0.4 )
    
    def test_set_nn_energy( self ):
        self.lattice.set_nn_energy( 'foo' )
        self.assertEqual( self.lattice.nn_energy, 'foo' )
  
    def test_set_cn_energies( self ):
        cn_energies = { 'A' : 0.5, 'B' : 0.3 }
        sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].label = 'A'
        sites[1].label = 'B'
        sites[0].set_cn_occupation_energies = Mock()
        sites[1].set_cn_occupation_energies = Mock()
        self.lattice.sites = sites
        self.lattice.set_cn_energies( cn_energies )
        sites[0].set_cn_occupation_energies.assert_called_with( 0.5 )
        sites[1].set_cn_occupation_energies.assert_called_with( 0.3 )
        self.assertEqual( self.lattice.cn_energies, cn_energies )

    def test_site_coordination_numbers( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].label = 'A'
        sites[1].label = 'B'
        sites[0].neighbours = [ 1,2,3,4 ]
        sites[1].neighbours = [ 1,2,3,4,5,6 ]
        self.lattice.sites = sites
        self.assertEqual( self.lattice.site_coordination_numbers(), { 'A' : 4, 'B' : 6 } )

    def test_mixed_site_coordination_numbers_raises_ValueError( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].label = 'A'
        sites[1].label = 'A'
        sites[0].neighbours = [ 1,2,3,4 ]
        sites[1].neighbours = [ 1,2,3,4,5,6 ]
        self.lattice.sites = sites
        with self.assertRaises( ValueError ):
            self.lattice.site_coordination_numbers()

    def test_site_specific_coordination_numbers( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].label = 'A'
        sites[1].label = 'B'
        sites[0].site_specific_neighbours = Mock( return_value='foo' )
        sites[1].site_specific_neighbours = Mock( return_value='bar' )
        self.lattice.sites = sites
        self.assertEqual( self.lattice.site_specific_coordination_numbers(), { 'A' : 'foo', 'B' : 'bar' } ) 

    def test_mixed_site_specific_coordination_numbers_raises_ValueError( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].label = 'A'
        sites[1].label = 'A'
        sites[0].site_specific_neighbours = Mock( return_value='foo' )
        sites[1].site_specific_neighbours = Mock( return_value='bar' )
        self.lattice.sites = sites
        with self.assertRaises( ValueError ):
            self.lattice.site_specific_coordination_numbers()

    def test_connected_site_pairs( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].label = 'A'
        sites[1].label = 'B'
        sites[0].p_neighbours = [ sites[1] ]
        sites[1].p_neighbours = [ sites[0] ]
        self.lattice.sites = sites
        self.assertEqual( self.lattice.connected_site_pairs(), { 'A' : [ 'B' ], 'B' : [ 'A' ] } )

    def test_conflicting_connected_site_pairs_raises_ValueError( self ):
        sites = [ Mock( spec=Site ), Mock( spec=Site ), Mock( spec=Site ) ]
        sites[0].label = 'A'
        sites[1].label = 'A'
        sites[2].label = 'B'
        sites[0].p_neighbours = [ sites[2] ]
        sites[1].p_neighbours = [ sites[2], sites[0] ]
        sites[2].p_neighbours = [ sites[1] ]
        self.lattice.sites = sites
        with self.assertRaises( ValueError ):
            self.lattice.connected_site_pairs()

if __name__ == '__main__':
    unittest.main()
