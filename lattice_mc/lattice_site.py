class Site:

    index = 0

    def __init__( self, number, coordinates, neighbours, energy ):
        self.number = number
        self.index = Site.index
        Site.index += 1
        self.r = coordinates
        self.neighbours = neighbours
        self.energy = energy
        self.occupation = 0
        self.atom = None
        self.is_occupied = False
