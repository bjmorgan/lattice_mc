import math
from lattice_mc.global_vars import kT, rate_prefactor

"""
A Jump describes a possible move by a particle from one site to another site.
"""

# TODO (maybe) does it make more sense to subclass Jump (and other classes) depending on the type of Hamiltonian?

class Jump:

    def __init__( self, initial_site, final_site, nearest_neighbour_energy=False, coordination_number_energy=False, jump_lookup_table=None ):
        """
        Initialise a Jump instance.

        Args:
            initial_site (Site): Lattice site occupied before the jump.
            final_site (Site): Lattice site occupied after the jump.
            nearest_neighbour_energy (Bool, optional): Does the jump probability depend on a nearest-neighbour energy term? Defaults to False.
            coordination_number_energy (Bool, optional): Does the jump probability depend on a coordination-number dependent energy term? Defaults to False.
            jump_lookup_table (:obj:`LookupTable`, optional): If the jump relative probabilities have been precalculated and stored in a lookup-table, this table should be passed in here. If not, jump probabilities are calculated on the fly. Defaults to None.

        Returns:
            None
        """
        self.initial_site = initial_site
        self.final_site = final_site
        self.nearest_neighbour_energy = nearest_neighbour_energy
        self.coordination_number_energy = coordination_number_energy
        if jump_lookup_table:
            self._relative_probability = self.relative_probability_from_lookup_table( jump_lookup_table )
        else:
            self._relative_probability = self.boltzmann_factor()
        
    def rate( self ):
        """
        Average rate for this jump. Calculated as ( v_0 * P_jump ).

        Args:
            None

        Returns:
            (Float): The average rate for this jump.
        """
        return rate_prefactor * self.relative_probability
    
    def boltzmann_factor( self ):
        """
        Boltzmann probability factor for accepting this jump, following the Metropolis algorithm.
        
        Args:
            None

        Returns:
            (Float): Metropolis relative probability of accepting this jump.
        """
        if self.delta_E() <= 0.0:
            return 1.0
        else:
            return math.exp( -self.delta_E() / kT )

    def delta_E( self ):
        """
        The change in system energy if this jump were accepted.

        Args:
            None

        Returns:
            (Float): delta E
        """
        site_delta_E = self.final_site.energy - self.initial_site.energy
        if self.nearest_neighbour_energy:
            site_delta_E += self.nearest_neighbour_delta_E()
        if self.coordination_number_energy:
            site_delta_E += self.coordination_number_delta_E()
        return site_delta_E

    def nearest_neighbour_delta_E( self ):
        """
        Nearest-neighbour interaction contribution to the change in system energy if this jump were accepted.

        Args:
            None

        Returns:
            (Float): delta E (nearest-neighbour)
        """
        delta_nn = self.final_site.nn_occupation() - self.initial_site.nn_occupation() - 1 # -1 because the hopping ion is not counted in the final site occupation number
        return ( delta_nn * self.nearest_neighbour_energy )
    
    def coordination_number_delta_E( self ):
        """
        Coordination-number dependent energy conrtibution to the change in system energy if this jump were accepted.

        Args:
            None

        Returns:
            (Float): delta E (coordination-number)
        """
        initial_site_neighbours = [ s for s in self.initial_site.p_neighbours if s.is_occupied ] # excludes final site, since this is always unoccupied
        final_site_neighbours = [ s for s in self.final_site.p_neighbours if s.is_occupied and s is not self.initial_site ] # excludes initial site
        initial_cn_occupation_energy = ( self.initial_site.cn_occupation_energy() + 
            sum( [ site.cn_occupation_energy() for site in initial_site_neighbours ] ) +
            sum( [ site.cn_occupation_energy() for site in final_site_neighbours ] ) )
        final_cn_occupation_energy = ( self.final_site.cn_occupation_energy( delta_occupation = { self.initial_site.label : -1 } ) +
            sum( [ site.cn_occupation_energy( delta_occupation = { self.initial_site.label : -1 } ) for site in initial_site_neighbours ] ) +
            sum( [ site.cn_occupation_energy( delta_occupation = { self.final_site.label : +1 } ) for site in final_site_neighbours ] ) )
        return ( final_cn_occupation_energy - initial_cn_occupation_energy )

    def dr( self, cell_lengths ):
        """
        Particle displacement vector for this jump

        Args:
            cell_lengths (np.array(x,y,z)): Cell lengths for the orthogonal simulation cell.

        Returns
            (np.array(x,y,z)): dr
        """
        half_cell_lengths = cell_lengths / 2.0
        this_dr = self.final_site.r - self.initial_site.r
        for i in range( 3 ):
            if this_dr[ i ] > half_cell_lengths[ i ]:
                this_dr[ i ] -= cell_lengths[ i ]
            if this_dr[ i ] < -half_cell_lengths[ i ]:
                this_dr[ i ] += cell_lengths[ i ]
        return this_dr

    def relative_probability_from_lookup_table( self, jump_lookup_table ):
        """
        Relative probability of accepting this jump from a lookup-table.

        Args:
            jump_lookup_table (LookupTable): the lookup table to be used for this jump.

        Returns:
            (Float): relative probability of accepting this jump.
        """
        l1 = self.initial_site.label
        l2 = self.final_site.label
        c1 = self.initial_site.nn_occupation()
        c2 = self.final_site.nn_occupation()
        return jump_lookup_table.jump_probability[ l1 ][ l2 ][ c1 ][ c2 ]

    @property
    def relative_probability( self ):
        return self._relative_probability

    @relative_probability.setter
    def relative_probability( self, value ):
        self._relative_probability = value

