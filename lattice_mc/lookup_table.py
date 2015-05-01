import math
import sys

k_boltzmann = 8.6173324e-5
temperature = 298.0
kT = k_boltzmann * temperature
rate_prefactor = 1e13 # standard vibrational frequency

class LookupTable:

    def __init__( self, lattice ):
        self.site_energies = lattice.site_energies
        self.nn_energy = lattice.nn_energy
        connected_site_pairs = lattice.connected_site_pairs()
        max_coordination_per_site = lattice.site_coordination_numbers()
        self.jump_probability = {}
        for site_label_1 in connected_site_pairs:
            self.jump_probability[ site_label_1 ] = {}
            for site_label_2 in connected_site_pairs[ site_label_1 ]:
                self.jump_probability[ site_label_1 ][ site_label_2 ] = {}
                for coordination_1 in range( max_coordination_per_site[ site_label_1 ] ):
                    self.jump_probability[ site_label_1 ][ site_label_2 ][ coordination_1 ] = {}
                    for coordination_2 in range( 1,max_coordination_per_site[ site_label_2 ] + 1 ):
                        self.jump_probability[ site_label_1 ][ site_label_2 ][ coordination_1 ][ coordination_2 ] = self.relative_probability( site_label_1, site_label_2, coordination_1, coordination_2 )

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
