import numpy as np
import math
import random
from lattice_mc.global_vars import rate_prefactor

class Transitions:
    """
    Transitions class

    Contains methods that operate on sets of Jumps.
    """

    def __init__( self, jumps ):
        """
        Initialise a Transitions object.

        Args:
            jumps (List(Jump)): List of jumps to be contained in this Transitions object.

        Returns:
            None
        """
        self.jumps = jumps
        self.p = np.array( [ jump.relative_probability for jump in self.jumps ] )
        
    def cumulative_probabilities( self ):
        """
        Cumulative sum of the relative probabilities for all possible jumps.

        Args:
            None

        Returns:
            (np.array): Cumulative sum of relative jump probabilities.
        """
        partition_function = np.sum( self.p )
        return np.cumsum( self.p ) / partition_function

    def random( self ):
        """
        Select a jump at random with appropriate relative probabilities.

        Args:
            None

        Returns:
            (Jump): The randomly selected Jump.
        """
        j = np.searchsorted( self.cumulative_probabilities(), random.random() )
        return self.jumps[ j ]

    def time_to_jump( self ):
        """
        The timestep until the next jump.

        Args:
            None

        Returns:
            (Float): The timestep until the next jump.
        """ 
        k_tot = rate_prefactor * np.sum( self.p )
        return -( 1.0 / k_tot ) * math.log( random.random() )
