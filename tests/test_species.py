import unittest
from lattice_mc.species import Species
from lattice_mc.atom import Atom
from unittest.mock import Mock
import numpy as np

class SpeciesTestCase( unittest.TestCase ):
    """Test for Species class"""

    def test_species_is_initialised( self ):
        mock = Mock( spec=Atom )
        atoms = [ mock ] * 3
        species = Species( atoms )
        self.assertEqual( species.atoms, atoms )

    def test_sites_occupied( self ):
        mock = Mock( spec=Atom )
        mock.site.number = 3
        atoms = [ mock ] * 3
        species = Species( atoms )
        self.assertEqual( species.sites_occupied(), [ 3, 3, 3 ] )

    def test_sum_dr_squared( self ):
        mock = Mock( spec=Atom )
        mock.dr_squared.return_value = 2.0
        atoms = [ mock ] * 3
        species = Species( atoms )
        self.assertEqual( species.sum_dr_squared(), 6.0 )

    def test_collective_dr_squared( self ):
        mock = Mock( spec=Atom )
        mock.dr = np.array( [ 2.0, 1.0, 1.0 ] )
        atoms = [ mock ] * 3
        species = Species( atoms )
        self.assertEqual( species.collective_dr_squared(), 54.0 ) # (2+2+2)**2 + (1+1+1)**2 + (1+1+1)**2

    def test_occupation( self ):
        mock = Mock( spec=Atom )
        mock.site.label = 'L'
        atoms = [ mock ] * 3
        species = Species( atoms )
        self.assertEqual( species.occupations( 'L' ), 3 )

    def test_summed_dr2( self ):
        mock = Mock( spec=Atom )
        mock.summed_dr2 = 5.0
        n_atoms = 3
        atoms = [ mock ] * n_atoms
        species = Species( atoms )
        self.assertEqual( species.summed_dr2(), mock.summed_dr2 * n_atoms )

    def test_tracer_correlation( self ):
        species = Species( atoms = [ 1, 2, 3 ] )
        species.sum_dr_squared = Mock( return_value=10.0 )
        species.summed_dr2 = Mock( return_value=5.0 )
        self.assertEqual( species.tracer_correlation(), 2.0 )

    def test_collective_correlation( self ):
        species = Species( atoms = [ 1, 2, 3 ] )
        species.collective_dr_squared = Mock( return_value=10.0 )
        species.summed_dr2 = Mock( return_value=5.0 )
        self.assertEqual( species.collective_correlation(), 2.0 )

if __name__ == '__main__':
    unittest.main()
