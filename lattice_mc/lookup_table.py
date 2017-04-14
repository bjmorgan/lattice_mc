import math
import sys
import itertools as it
from lattice_mc.global_vars import kT, rate_prefactor

def metropolis( delta_E ):
    """
    Boltzmann probability factor for an event with an energy change `delta_E`, following the Metropolis algorithm.
        
    Args:
        delta_E (Float): The change in energy.

    Returns:
        (Float): Metropolis relative probability for this event.
        """
    if delta_E <= 0.0:
        return 1.0
    else:
       return math.exp( -delta_E / ( kT ) )

class LookupTable: #TODO if nearest-neighbour and coordination number dependent look-up tables have different data structures, they should each subclass this general class: the different implementations for setting these up and accessing the jump probabilities can then be self-contained
    """
    LookupTable class
    """

    def __init__( self, lattice, hamiltonian ):
        """
        Initialise a LookupTable object instance.

        Args:
            lattice (lattice_mc.Lattice): The lattice object, used to define the allowed jumps.
            hamiltonian (Str): The model Hamiltonian used to define the jump energies.
                Allowed values = `nearest-neigbour`
        
        Returns:
            None
        """
        expected_hamiltonian_values = [ 'nearest-neighbour' ]
        if hamiltonian not in expected_hamiltonian_values:
            raise ValueError( hamiltonian )
        self.site_energies = lattice.site_energies
        self.nn_energy = lattice.nn_energy
        self.cn_energy = lattice.cn_energies
        self.connected_site_pairs = lattice.connected_site_pairs()
        self.max_coordination_per_site = lattice.max_site_coordination_numbers()
        self.site_specific_coordination_per_site = lattice.site_specific_coordination_numbers()
        if hamiltonian == 'nearest-neighbour':
            self.generate_nearest_neighbour_lookup_table()
        
    def relative_probability( self, l1, l2, c1, c2 ):
        """
        The relative probability for a jump between two sites with specific site types and coordination numbers.

        Args:
            l1 (Str): Site label for the initial site.
            l2 (Str): Site label for the final site.
            c1 (Int): Coordination number for the initial site.
            c2 (Int): Coordination number for the final site.
   
        Returns:
            (Float): The relative probability of this jump occurring.
        """
        if self.site_energies:
            site_delta_E = self.site_energies[ l2 ] - self.site_energies[ l1 ]
        else:
            site_delta_E = 0.0
        if self.nn_energy:
            delta_nn = c2 - c1 - 1 # -1 because the hopping ion is not counted in the final site occupation number
            site_delta_E += delta_nn * self.nn_energy
        return metropolis( site_delta_E )

    def generate_nearest_neighbour_lookup_table( self ):
        """
        Construct a look-up table of relative jump probabilities for a nearest-neighbour interaction Hamiltonian.

        Args:
            None.

        Returns:
            None.
        """
        self.jump_probability = {}
        for site_label_1 in self.connected_site_pairs:
            self.jump_probability[ site_label_1 ] = {}
            for site_label_2 in self.connected_site_pairs[ site_label_1 ]:
                self.jump_probability[ site_label_1 ][ site_label_2 ] = {}
                for coordination_1 in range( self.max_coordination_per_site[ site_label_1 ] ):
                    self.jump_probability[ site_label_1 ][ site_label_2 ][ coordination_1 ] = {}
                    for coordination_2 in range( 1, self.max_coordination_per_site[ site_label_2 ] + 1 ):
                        self.jump_probability[ site_label_1 ][ site_label_2 ][ coordination_1 ][ coordination_2 ] = self.relative_probability( site_label_1, site_label_2, coordination_1, coordination_2 )

