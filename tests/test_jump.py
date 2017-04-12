import unittest
from lattice_mc import jump as jump_module
from lattice_mc.lattice_site import Site
from lattice_mc.lookup_table import LookupTable
from unittest.mock import Mock, patch
import numpy as np

Jump = jump_module.Jump

class JumpTestCase( unittest.TestCase ):
    """Test for Jump class"""

    def setUp( self ):
        self.mock_initial_site = Mock( spec=Site )
        self.mock_final_site = Mock( spec=Site )
        with patch( 'lattice_mc.jump.Jump.boltzmann_factor' ) as self.mock_boltzmann_factor:
            self.mock_boltzmann_factor.return_value = 0.13
            self.jump = Jump( self.mock_initial_site, self.mock_final_site )

    def test_jump_is_initialised_with_no_interactions( self ):
        self.assertEqual( self.jump.initial_site, self.mock_initial_site )
        self.assertEqual( self.jump.final_site, self.mock_final_site )
        self.assertEqual( self.jump.nearest_neighbour_energy, False )
        self.assertEqual( self.jump.coordination_number_energy, False )
        self.assertEqual( self.jump.relative_probability, self.mock_boltzmann_factor.return_value )

    def test_jump_is_initialised_with_nn_energy( self ):
        nearest_neighbour_energy = 'foo'
        with patch( 'lattice_mc.jump.Jump.boltzmann_factor' ) as mock_boltzmann_factor:
            mock_boltzmann_factor.return_value = 0.13
            jump = Jump( self.mock_initial_site, self.mock_final_site, nearest_neighbour_energy = nearest_neighbour_energy )
        self.assertEqual( jump.initial_site, self.mock_initial_site )
        self.assertEqual( jump.final_site, self.mock_final_site )
        self.assertEqual( jump.nearest_neighbour_energy, nearest_neighbour_energy )
        self.assertEqual( jump.coordination_number_energy, False )
        self.assertEqual( jump.relative_probability, 0.13 )

    def test_jump_is_initialised_with_cn_energy( self ):
        coordination_number_energy = 'foo'
        with patch( 'lattice_mc.jump.Jump.boltzmann_factor' ) as mock_boltzmann_factor:
            mock_boltzmann_factor.return_value = 0.13
            jump = Jump( self.mock_initial_site, self.mock_final_site, coordination_number_energy = coordination_number_energy )
        self.assertEqual( jump.initial_site, self.mock_initial_site )
        self.assertEqual( jump.final_site, self.mock_final_site )
        self.assertEqual( jump.nearest_neighbour_energy, False )
        self.assertEqual( jump.coordination_number_energy, coordination_number_energy )
        self.assertEqual( jump.relative_probability, 0.13 )

    def test_jump_is_initialised_with_lookup_table( self ):
        mock_initial_site = Mock( spec=Site )
        mock_final_site = Mock( spec=Site )
        jump_lookup_table = 'foo'
        with patch( 'lattice_mc.jump.Jump.relative_probability_from_lookup_table' ) as mock_p_from_lookup_table:
            mock_p_from_lookup_table.return_value = 0.73
            jump = Jump( mock_initial_site, mock_final_site, jump_lookup_table = jump_lookup_table )
        self.assertEqual( jump.initial_site, mock_initial_site )
        self.assertEqual( jump.final_site, mock_final_site )
        self.assertEqual( jump.nearest_neighbour_energy, False )
        self.assertEqual( jump.coordination_number_energy, False )
        self.assertEqual( jump.relative_probability, 0.73 )

    def test_rate( self ):
        self.jump.relative_probability = 0.5
        jump_module.rate_prefactor = 1e-3
        self.assertEqual( self.jump.rate(), 5e-4 )

    def test_boltzmann_factor_lt_0( self ):
        self.jump.delta_E = Mock( return_value=-2.5 )
        self.assertEqual( self.jump.boltzmann_factor(), 1.0 )
 
    def test_boltzmann_factor_gt_0( self ):
        self.jump.delta_E = Mock( return_value=0.2 )
        self.assertEqual( self.jump.boltzmann_factor(), 0.000414569521855495 )

    def test_delta_E_non_interacting( self ):
        initial_energy = 1.5
        final_energy = 2.1
        self.jump.initial_site.energy = initial_energy
        self.jump.final_site.energy = final_energy
        self.jump.nearest_neighbour_energy = False
        self.jump.coordination_number_energy = False
        delta_E = final_energy - initial_energy
        self.assertEqual( self.jump.delta_E(), delta_E )

    def test_delta_E_nn_interactions( self ):
        initial_energy = 1.5
        final_energy = 2.1
        self.jump.initial_site.energy = initial_energy
        self.jump.final_site.energy = final_energy
        self.jump.nearest_neighbour_energy = True
        self.jump.coordination_number_energy = False
        self.jump.nearest_neighbour_delta_E = Mock( return_value = 0.3 )
        delta_E = final_energy - initial_energy + self.jump.nearest_neighbour_delta_E.return_value
        self.assertEqual( self.jump.delta_E(), delta_E )

    def test_delta_E_cn_dependent( self ):
        initial_energy = 1.5
        final_energy = 2.1
        self.jump.initial_site.energy = initial_energy
        self.jump.final_site.energy = final_energy
        self.jump.nearest_neighbour_energy = False
        self.jump.coordination_number_energy = True
        self.jump.coordination_number_delta_E = Mock( return_value = 0.5 )
        delta_E = final_energy - initial_energy + self.jump.coordination_number_delta_E.return_value
        self.assertEqual( self.jump.delta_E(), delta_E )

    def test_nearest_neighbour_delta_E( self ):
        self.jump.initial_site.nn_occupation = Mock( return_value = 3 )
        self.jump.final_site.nn_occupation = Mock( return_value = 1 )
        self.jump.nearest_neighbour_energy = 0.7
        self.assertEqual( self.jump.nearest_neighbour_delta_E(), -3 * 0.7 )

    def test_coordination_number_delta_E( self ):
        self.jump.initial_site.cn_occupation_energy = Mock( return_value = 1.5 )
        self.jump.initial_site.label = 'A'
        self.jump.final_site.cn_occupation_energy = Mock( return_value = 0.2 )
        self.jump.final_site.label = 'A'
        occupied = [ True, False, True, True, False ]
        cn_occupation_energy = [ 1.0, 2.0, 3.0, 4.0, 5.0 ]
        mock_sites = [ Mock( spec=Site, is_occupied=o ) for o in occupied ]
        for m, e in zip( mock_sites, cn_occupation_energy ):
            m.cn_occupation_energy.return_value = e
        self.jump.initial_site.p_neighbours = mock_sites[:2]
        self.jump.final_site.p_neighbours = mock_sites[2:]
        delta_E = ( 0.2 + 1.0 + 3.0 + 4.0 ) - ( 1.5 + 1.0 + 3.0 + 4.0 )
        self.assertEqual( self.jump.coordination_number_delta_E(), delta_E )

    def test_dr( self ):
        cell_lengths = np.array( [ 10.0, 10.0, 10.0 ] )
        self.jump.initial_site.r = np.array( [ 3.0, 3.0, 3.0 ] )
        self.jump.final_site.r = np.array( [ 4.0, 5.0, 6.0 ] )
        np.testing.assert_array_equal( self.jump.dr( cell_lengths ), np.array( [ 1.0, 2.0, 3.0 ] ) )

    def test_dr_minimum_image_convention( self ):
        cell_lengths = np.array( [ 10.0, 10.0, 10.0 ] )
        self.jump.initial_site.r = np.array( [ 9.0, 1.0, 3.0 ] )
        self.jump.final_site.r = np.array( [ 1.0, 9.0, 1.0 ] )
        np.testing.assert_array_equal( self.jump.dr( cell_lengths ), np.array( [ 2.0, -2.0, -2.0 ] ) )

    def test_relative_probability_from_lookup_table( self ):
        self.jump.initial_site.label = 'A'
        self.jump.final_site.label = 'B'
        self.jump.initial_site.nn_occupation = Mock( return_value=2 )
        self.jump.final_site.nn_occupation = Mock( return_value=1 )
        jump_lookup_table = Mock( spec=LookupTable )
        jump_lookup_table.jump_probability = { 'A': { 'B': { 2: { 1: 0.5 } } } }
        self.assertEqual( self.jump.relative_probability_from_lookup_table( jump_lookup_table ), 0.5 )

    def test_relative_probability_getter( self ):
        self.assertEqual( self.jump.relative_probability, self.jump._relative_probability )

    def test_relative_probability_setter( self ):
        self.jump.relative_probability = 0.5
        self.assertEqual( self.jump._relative_probability, 0.5 )
 
if __name__ == '__main__':
    unittest.main()
