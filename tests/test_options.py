import unittest
from lattice_mc.options import Options
from unittest.mock import Mock

class OptionsTestCase( unittest.TestCase ):
    """Test for Options class"""

    def setUp( self ):
        self.options = Options()

    def test_options_is_initialised( self ):
        self.assertIs( self.options.site_energies, None )
        self.assertIs( self.options.nn_energy_scaling, None )
        self.assertEqual( self.options.number_of_equilibration_jumps, 0 )

    def test_set_number_of_atoms( self ):
        self.options.set_number_of_atoms( 36 )
        self.assertEqual( self.options.number_of_atoms, 36 )

    def test_set_nn_energy_scaling( self ):
        self.options.set_nn_energy_scaling( 2.0 )
        self.assertEqual( self.options.nn_energy_scaling, 2.0 )

    def test_cn_energies( self ):
        self.options.set_cn_energies( 10.0 )
        self.assertEqual( self.options.cn_energies, 10.0 )

    def test_set_site_energies( self ):
        self.options.set_site_energies( { 'O' : 0.0 } )
        self.assertEqual( self.options.site_energies, { 'O' : 0.0 } )

    def test_set_number_of_jumps( self ):
        self.options.set_number_of_jumps( 1000 )
        self.assertEqual( self.options.number_of_jumps, 1000 )

    def test_set_number_of_equilibration_jumps( self ):
        self.options.set_number_of_equilibration_jumps( 1000 )
        self.assertEqual( self.options.number_of_equilibration_jumps, 1000 )

    def test_read_lattice_from_file( self ):
        self.options.read_lattice_from_file( 'foo' )
        self.assertEqual( self.options.lattice_site_file, 'foo' )

    def test_set_lattice_cell_lengths( self ):
        self.options.set_lattice_cell_lengths( 1.0 )
        self.assertEqual( self.options.lattice_cell_lengths, 1.0 )

if __name__ == '__main__':
    unittest.main()
