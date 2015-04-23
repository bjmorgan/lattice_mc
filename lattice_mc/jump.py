class Jump:

    k_boltzmann = 8.6173324e-5

    def __init__( self, initial_site, final_site ):
        self.initial_site = initial_site
        self.final_site = final_site
        self.relative_probability = self.boltzmann_factor( 298.0 )

    def boltzmann_factor( self, temp ):
        return 1.0
        #return math.exp( -self.delta_E() / ( k_boltzmann * temp ) )

    def delta_E( self ):
        return self.final_site.energy - self.initial_site.energy

    def dr( self, cell_lengths ):
        half_cell_lengths = cell_lengths / 2.0
        this_dr = self.final_site.r - self.initial_site.r
        for i in range( 3 ):
            if this_dr[ i ] > half_cell_lengths[ i ]:
                this_dr[ i ] -= cell_lengths[ i ]
            if this_dr[ i ] < -half_cell_lengths[ i ]:
                this_dr[ i ] += cell_lengths[ i ]
        return this_dr
