from collections import Counter

class Site:

    index = 0

    def __init__( self, number, coordinates, neighbours, energy, label, cn_energies = None ):
        self.number = number
        self.index = Site.index
        Site.index += 1
        self.r = coordinates
        self.neighbours = neighbours
        self.p_neighbours = None # pointer to neighbouring sites. initialised in Lattice.__init__
        self.energy = energy
        self.occupation = 0
        self.atom = None
        self.is_occupied = False
        self.label = label
        self.time_occupied = 0.0
        self.cn_occupation_energies = cn_energies

    def nn_occupation( self ):
        """The number of occupied nearest-neighbour sites

        Arguments: None

        Returns (int): The number of occupied nearest-neighbour sites
        """
        return sum( [ site.is_occupied for site in self.p_neighbours ] )

    def site_specific_nn_occupation( self ):
        to_return = { l : 0 for l in set( ( site.label for site in self.p_neighbours ) ) }
        for site in self.p_neighbours:
            if site.is_occupied:
             to_return[ site.label ] += 1
        return to_return

    def site_specific_neighbours( self ):
        return dict( Counter( ( site.label for site in self.p_neighbours ) ) )
        
    def set_cn_occupation_energies( self, cn_energies ):
        self.cn_occupation_energies = cn_energies
  
    def cn_occupation_energy( self, delta_occupation = None ):
        nn_occupations = self.site_specific_nn_occupation()
        if delta_occupation:
            for site in delta_occupation:
                assert( site in nn_occupations )
                nn_occupations[ site ] += delta_occupation[ site ]
        return sum( [ self.cn_occupation_energies[ s ][ n ] for s, n in nn_occupations.items() ] )
