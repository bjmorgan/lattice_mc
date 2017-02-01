import unittest
from lattice_mc.transitions import Transitions
from lattice_mc.jump import Jump
from unittest.mock import Mock, patch
import numpy as np

class TranstitionsTestCase( unittest.TestCase ):
    """Test for Transitions class"""

    def setUp( self ):
        self.jumps = [ Mock( spec=Jump, relative_probability=0.25 ) for i in range(4) ] 
        self.transitions = Transitions( self.jumps )

    def test_transitions_is_initialised( self ):
        self.assertEqual( self.transitions.jumps, self.jumps )
        np.testing.assert_array_equal( self.transitions.p, np.array( [ 0.25, 0.25, 0.25, 0.25 ] ) )

    @patch( 'random.random' )
    def test_random( self, mock_random ):
        self.transitions.p = np.array( [ 0.1, 0.2, 0.3, 0.4 ] )
        mock_random.return_value = 0.05
        self.assertIs( self.transitions.random(), self.jumps[0] )
        mock_random.return_value = 0.15
        self.assertIs( self.transitions.random(), self.jumps[1] )
        mock_random.return_value = 0.60
        self.assertIs( self.transitions.random(), self.jumps[2] )
        mock_random.return_value = 0.75
        self.assertIs( self.transitions.random(), self.jumps[3] )

    @patch( 'random.random' )
    def test_time_to_jump( self, mock_random ):
        self.transitions.p = np.array( [ 0.1, 0.2, 0.3, 0.4 ] )
        mock_random.return_value = 0.15
        self.assertAlmostEqual( self.transitions.time_to_jump(), 1.89711998489e-13 )
 
if __name__ == '__main__':
    unittest.main()
