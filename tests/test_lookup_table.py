import unittest
from lattice_mc.lookup_table import LookupTable, metropolis
from lattice_mc.lattice import Lattice
from unittest.mock import Mock, patch, call
import numpy as np

class LookupTableSupportFunctionsTestCase( unittest.TestCase ):

    def test_metropolis_returns_one_for_negative_delta_E( self ):
        self.assertEqual( metropolis( -0.02 ), 1.0 )

    def test_metropolis_returns_boltzmann_factor_for_positive_delta_E( self ):
        self.assertAlmostEqual( metropolis( +0.02 ), 0.45894412999 )

class LookupTableTestCase( unittest.TestCase ):
    """Tests for LookupTable class"""

    @patch( 'lattice_mc.lookup_table.LookupTable.generate_nearest_neighbour_lookup_table' )
    def setUp( self, mock_generate_nearest_neighbour_lookup_table ):
        self.lattice = Mock( spec=Lattice )
        self.lattice.site_energies = 'foo'
        self.lattice.nn_energy = 'bar'
        self.lattice.cn_energies = 'baz'
        self.lattice.connected_site_pairs = Mock( return_value='qux' )
        self.lattice.site_specific_coordination_numbers = Mock( return_value='quux' )
        hamiltonian = 'nearest-neighbour'
        self.table = LookupTable( self.lattice, hamiltonian )
    
    @patch( 'lattice_mc.lookup_table.LookupTable.generate_nearest_neighbour_lookup_table' )
    def test_nearest_neighbour_lookup_table_is_initialised( self, mock_generate_nearest_neighbour_lookup_table ):
        hamiltonian = 'nearest-neighbour'
        lookup_table = LookupTable( self.lattice, hamiltonian )
        self.assertEqual( lookup_table.site_energies, 'foo' )
        self.assertEqual( lookup_table.nn_energy, 'bar' )
        self.assertEqual( lookup_table.cn_energy, 'baz' )
        self.assertEqual( lookup_table.connected_site_pairs, 'qux' )
        self.assertEqual( lookup_table.site_specific_coordination_per_site, 'quux' )
        mock_generate_nearest_neighbour_lookup_table.assert_called_once()
       
    def test_lookup_table_init_with_invalid_hamiltonian_keywork_raises_ValueError( self ):
        hamiltonian = 'foo'
        with self.assertRaises( ValueError ):
            LookupTable( self.lattice, hamiltonian )

    @patch( 'lattice_mc.lookup_table.metropolis' )
    def test_relative_probability_with_site_energies( self, mock_metropolis ):
        self.table.site_energies = { 'A' : 1.0, 'B' : 2.0 }   
        self.table.nn_energy = None
        self.table.relative_probability( 'A', 'B', 3, 1 )
        mock_metropolis.assert_called_with( 1.0 )

    @patch( 'lattice_mc.lookup_table.metropolis' )
    def test_relative_probability_with_coordination_numbers( self, mock_metropolis ):
        self.table.site_energies = None
        self.table.nn_energy = 1.0
        self.table.relative_probability( 'A', 'B', 3, 1 )
        mock_metropolis.assert_called_with( -3.0 )

    @patch( 'lattice_mc.lookup_table.metropolis' )
    def test_relative_probability_with_site_energies_and_coordination_numbers( self, mock_metropolis ):
        self.table.site_energies = { 'A' : 1.0, 'B' : 2.0 }
        self.table.nn_energy = 1.0
        self.table.relative_probability( 'A', 'B', 3, 1 )       
        mock_metropolis.assert_called_with( -2.0 )

    def test_generate_nearest_neighbour_lookup_table( self ):
        self.table.connected_site_pairs = { 'A' : [ 'B' ] }
        self.table.max_coordination_per_site = { 'A' : 2, 'B' : 3 }
        self.table.relative_probability = Mock( side_effect = [ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] )
        self.table.generate_nearest_neighbour_lookup_table()
        self.assertEqual( self.table.jump_probability, {'A': {'B': {0: {1: 1.0, 2: 2.0, 3: 3.0}, 1: {1: 4.0, 2: 5.0, 3: 6.0}}}} )
        self.assertEqual( self.table.relative_probability.mock_calls, [call('A', 'B', 0, 1), call('A', 'B', 0, 2), call('A', 'B', 0, 3), call('A', 'B', 1, 1), call('A', 'B', 1, 2), call('A', 'B', 1, 3)] )      

if __name__ == '__main__':
    unittest.main()

