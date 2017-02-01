import numpy as np
import math
import random

prefactor = 1e13 # standard vibrational frequency

class Transitions:

    def __init__( self, jumps ):
        self.jumps = jumps
        self.p = np.array( [ jump.relative_probability for jump in self.jumps ] )
        
    def random( self ):
        partition_function = np.sum( self.p )
        cumulative_probabilities = np.cumsum( self.p ) / partition_function
        j = np.searchsorted( cumulative_probabilities, random.random() )
        return self.jumps[ j ]

    def time_to_jump( self ):
        k_tot = prefactor * np.sum( self.p )
        return -( 1.0 / k_tot ) * math.log( random.random() )
