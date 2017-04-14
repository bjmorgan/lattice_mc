from collections import Counter

class Site:
    """
    Site class
    """

    index = 0

    def __init__( self, number, coordinates, neighbours, energy, label, cn_energies=None ):
        """
        Initialise a lattce Site object.

        Args:
            number (Int): An identifying number for this site.
            coordinates (np.array(x,y,z)): The coordinates of this site.
            neighbours (List(Int)): A list of the id numbers of the neighbouring sites.
            energy (Float): On-site occupation energy.
            label (Str): Label for classifying this as a specific site type.
            cn_energies (:obj:Dict(Int:Float), optional): Dictionary of coordination-number dependent energies, e.g. { 0 : 0.0, 1 : 0.5, 2 : 2.0 }. Defaults to None.

        Returns:
            None

        Notes:
            There should be a 1:1 mapping between sites and site numbers.
        """
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
        """
        The number of occupied nearest-neighbour sites.

        Args: 
            None

        Returns:
            (Int): The number of occupied nearest-neighbour sites.
        """
        return sum( [ site.is_occupied for site in self.p_neighbours ] )

    def site_specific_nn_occupation( self ):
        """
        Returns the number of occupied nearest neighbour sites, classified by site type.

        Args:
            None

        Returns:
            (Dict(Str:Int)): Dictionary of nearest-neighbour occupied site numbers, classified by site label, e.g. { 'A' : 2, 'B' : 1 }.
        """
        to_return = { l : 0 for l in set( ( site.label for site in self.p_neighbours ) ) }
        for site in self.p_neighbours:
            if site.is_occupied:
             to_return[ site.label ] += 1
        return to_return

    def site_specific_neighbours( self ):
        """
        Returns the number of neighbouring sites, classified by site type.

        Args:
            None

        Returns:
            (Dict(Str:Int)): Dictionary of neighboring sites, classified by site label, e.g. { 'A' : 1, 'B' : 1 }.
        """
        return dict( Counter( ( site.label for site in self.p_neighbours ) ) )
        
    def set_cn_occupation_energies( self, cn_energies ):
        """
        Set the coordination-number dependent energies for this site.

        Args:
            cn_energies (Dict(Int:Float)): Dictionary of coordination number dependent site energies, e.g. { 0 : 0.0, 1 : 0.5 }.
 
        Returns:
            None
        """
        self.cn_occupation_energies = cn_energies
  
    def cn_occupation_energy( self, delta_occupation=None ):
        """
        The coordination-number dependent energy for this site.

        Args:
            delta_occupation (:obj:Dict(Str:Int), optional): A dictionary of a change in (site-type specific) coordination number, e.g. { 'A' : 1, 'B' : -1 }.
                If this is not None, the coordination-number dependent energy is calculated including these changes in neighbour-site occupations. Defaults to None

        Returns:
            (Float): The coordination-number dependent energy for this site.
        """
        nn_occupations = self.site_specific_nn_occupation()
        if delta_occupation:
            for site in delta_occupation:
                assert( site in nn_occupations )
                nn_occupations[ site ] += delta_occupation[ site ]
        return sum( [ self.cn_occupation_energies[ s ][ n ] for s, n in nn_occupations.items() ] )
