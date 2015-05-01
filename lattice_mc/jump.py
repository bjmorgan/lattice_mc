import math

k_boltzmann = 8.6173324e-5
temperature = 298.0
kT = k_boltzmann * temperature
rate_prefactor = 1e13 # standard vibrational frequency

class Jump:

    def __init__( self, initial_site, final_site, nearest_neighbour_energy = False, jump_lookup_table = False ):
        self.initial_site = initial_site
        self.final_site = final_site
        self.nearest_neighbour_energy = nearest_neighbour_energy
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
            delta_nn = self.final_site.nn_occupation() - self.initial_site.nn_occupation() -1 # -1 because the hopping ion is not counted in the final site occupation number
            site_delta_E += delta_nn * self.nearest_neighbour_energy
        # HARD CODED: penalise octahedral sites having both neighbours occupied
        #if self.final_site.label in [ 'T1', 'T2', 'T3' ]:
        #    for n_site in self.final_site.p_neighbours:
        #        if n_site is not self.initial_site:
        #            if n_site.is_occupied and n_site.nn_occupation() == 1:
        #                return site_delta_E + 10 * kT
        # END
        return site_delta_E

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
