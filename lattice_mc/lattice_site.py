class Site:

    index = 0

    def __init__( self, number, coordinates, neighbours, energy, label ):
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

    def nn_occupation( self ):
        return sum( [ site.is_occupied for site in self.p_neighbours ] )
