import numpy as np
import random

class Transitions:

    def __init__( self, jumps ):
        self.jumps = jumps

    def random( self ):
        p = np.array( [ jump.relative_probability for jump in self.jumps ] )
        partition_function = np.sum( p )
        cumulative_probabilities = np.cumsum( p ) / partition_function
        j = np.searchsorted( cumulative_probabilities, random.random() )
        return self.jumps[ j ]
