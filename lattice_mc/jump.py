import math

k_boltzmann = 8.6173324e-5
temperature = 298.0
kT = k_boltzmann * temperature
rate_prefactor = 1e13 # standard vibrational frequency

class Jump:

    def __init__( self, initial_site, final_site, nearest_neighbour_energy = False, coordination_number_energy = False, jump_lookup_table = False ):
        self.initial_site = initial_site
        self.final_site = final_site
        self.nearest_neighbour_energy = nearest_neighbour_energy
        self.coordination_number_energy = coordination_number_energy
        if jump_lookup_table:
            self.relative_probability = self.relative_probability_from_lookup_table( jump_lookup_table )
        else:
            self.relative_probability = self.boltzmann_factor()
        
    def rate( self ):
        return rate_prefactor * self.relative_probability
    
    def boltzmann_factor( self ):
        if self.delta_E() <= 0.0:
            return 1.0
        else:
            return math.exp( -self.delta_E() / kT )

    def delta_E( self ):
        site_delta_E = self.final_site.energy - self.initial_site.energy
        if self.nearest_neighbour_energy:
            site_delta_E += self.nearest_neighbour_delta_E()
        if self.coordination_number_energy:
            site_delta_E += self.coordination_number_delta_E()
        return site_delta_E

    def nearest_neighbour_delta_E( self ):
        delta_nn = self.final_site.nn_occupation() - self.initial_site.nn_occupation() - 1 # -1 because the hopping ion is not counted in the final site occupation number
        return ( delta_nn * self.nearest_neighbour_energy )
    
    def coordination_number_delta_E( self ):
        initial_site_neighbours = [ s for s in self.initial_site.p_neighbours if s.is_occupied ] # excludes final site, since this is always uncocupied
        final_site_neighbours = [ s for s in self.final_site.p_neighbours if s.is_occupied and s is not self.initial_site ] # excludes initial site
        initial_cn_occupation_energy = ( self.initial_site.cn_occupation_energy() + 
            sum( [ site.cn_occupation_energy() for site in initial_site_neighbours ] ) +
            sum( [ site.cn_occupation_energy() for site in final_site_neighbours ] ) )
        final_cn_occupation_energy = ( self.final_site.cn_occupation_energy( delta_occupation = { self.initial_site.label : -1 } ) +
            sum( [ site.cn_occupation_energy( delta_occupation = { self.initial_site.label : -1 } ) for site in initial_site_neighbours ] ) +
            sum( [ site.cn_occupation_energy( delta_occupation = { self.final_site.label : +1 } ) for site in final_site_neighbours ] ) )
        return ( final_cn_occupation_energy - initial_cn_occupation_energy )

    def dr( self, cell_lengths ):
        half_cell_lengths = cell_lengths / 2.0
        this_dr = self.final_site.r - self.initial_site.r
        for i in range( 3 ):
            if this_dr[ i ] > half_cell_lengths[ i ]:
                this_dr[ i ] -= cell_lengths[ i ]
            if this_dr[ i ] < -half_cell_lengths[ i ]:
                this_dr[ i ] += cell_lengths[ i ]
        return this_dr

    def relative_probability_from_lookup_table( self, jump_lookup_table ):
        l1 = self.initial_site.label
        l2 = self.final_site.label
        c1 = self.initial_site.nn_occupation()
        c2 = self.final_site.nn_occupation()
        return jump_lookup_table.jump_probability[ l1 ][ l2 ][ c1 ][ c2 ]

    @property
    def relative_probability( self ):
        return self.relative_probability
