import numpy as np
import math
import random
from lattice_mc.global_vars import rate_prefactor

class Transitions:

    def __init__( self, jumps ):
        self.jumps = jumps
        self.p = np.array( [ jump.relative_probability for jump in self.jumps ] )
        
    def cumulative_probabilities( self ):
        partition_function = np.sum( self.p )
        return np.cumsum( self.p ) / partition_function

    def random( self ):
        j = np.searchsorted( self.cumulative_probabilities(), random.random() )
        return self.jumps[ j ]

    def time_to_jump( self ):
        k_tot = rate_prefactor * np.sum( self.p )
        return -( 1.0 / k_tot ) * math.log( random.random() )
