import math
import sys
import itertools as it

k_boltzmann = 8.6173324e-5
temperature = 298.0
kT = k_boltzmann * temperature
rate_prefactor = 1e13 # standard vibrational frequency

def cn_combinations( cn ):
    return it.product(*( range( cn[i]+1) for i in sorted( cn ) ) )

class LookupTable: #TODO if nearest-neighbour and coordination number dependent look-up tables have different data structures, they should each subclass this general class: the different implementations for setting these up and accessing the jump probabilities can then be self-contained

    def __init__( self, lattice, hamiltonian ):
        expected_hamiltonian_values = [ 'nearest-neighbour', 'coordination_number' ]
        assert( hamiltonian in expected_hamiltonian_values )
        self.site_energies = lattice.site_energies
        self.nn_energy = lattice.nn_energy
        self.cn_energy = lattice.cn_energies
        self.connected_site_pairs = lattice.connected_site_pairs()
        self.max_coordination_per_site = lattice.site_coordination_numbers()
        self.site_specific_coordination_per_site = lattice.site_specific_coordination_numbers()
        if hamiltonian == 'nearest-neighbour':
            self.generate_nearest_neighbour_lookup_table()
        elif hamiltonian == 'coordination_number':
            self.generate_coordination_number_lookup_table()
        
    def relative_probability( self, l1, l2, c1, c2 ):
        if self.site_energies:
            site_delta_E = self.site_energies[ l2 ] - self.site_energies[ l1 ]
        else:
            site_delta_E = 0.0
        if self.nn_energy:
            delta_nn = c2 - c1 - 1 # -1 because the hopping ion is not counted in the final site occupation number
            site_delta_E += delta_nn * self.nn_energy
        if site_delta_E <= 0.0:
            return 1.0
        else:
            return math.exp( -site_delta_E / ( kT ) )

    def generate_nearest_neighbour_lookup_table( self ):
        self.jump_probability = {}
        for site_label_1 in self.connected_site_pairs:
            self.jump_probability[ site_label_1 ] = {}
            for site_label_2 in self.connected_site_pairs[ site_label_1 ]:
                self.jump_probability[ site_label_1 ][ site_label_2 ] = {}
                for coordination_1 in range( self.max_coordination_per_site[ site_label_1 ] ):
                    self.jump_probability[ site_label_1 ][ site_label_2 ][ coordination_1 ] = {}
                    for coordination_2 in range( 1, self.max_coordination_per_site[ site_label_2 ] + 1 ):
                        self.jump_probability[ site_label_1 ][ site_label_2 ][ coordination_1 ][ coordination_2 ] = self.relative_probability( site_label_1, site_label_2, coordination_1, coordination_2 )

    def generate_coordination_number_lookup_table( self ):
        self.jump_probability = {}
        for site_label_1 in self.connected_site_pairs:
            self.jump_probability[ site_label_1 ] = {}
            for site_label_2 in self.connected_site_pairs[ site_label_1 ]:
                self.jump_probability[ site_label_1 ][ site_label_2 ] = {}
                for neighbour_label_1 in self.site_specific_coordination_per_site[ site_label_1 ]:
                    print( neighbour_label_1 )
                    for coordination_1 in cn_combinations( self.site_specific_coordination_per_site[ site_label_1 ] ):
                        self.jump_probability[ site_label_1 ][ site_label_2 ][ coordination_1 ] = {}
                        for coordination_2 in cn_combinations( self.site_specific_coordination_per_site[ site_label_2 ] ):
                            self.jump_probability[ site_label_1 ][ site_label_2 ][ coordination_1 ][ coordination_2 ] = self.relative_probability( site_label_1, site_label_2, coordination_1, coordination_2 )
                            print( coordination_1, coordination_2 )
        print( self.jump_probability )




