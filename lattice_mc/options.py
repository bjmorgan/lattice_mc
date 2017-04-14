class Options:
    """
    Object for storing options for setting up and running a simulation.
    """

    def __init__( self ):
        """
        Initialise an Options object.

        Args:
            None

        Returns:
            None
        """
        self.site_energies = None
        self.nn_energy_scaling = None
        self.number_of_equilibration_jumps = 0

    def set_number_of_atoms( self, number_of_atoms ):
        """
        Set the number of atoms present in the simulation.

        Args:
            number_of_atoms (int):  The number of atoms.

        Returns:
            None

        """
        self.number_of_atoms = number_of_atoms

    def set_nn_energy_scaling( self, energy_scaling ):
        """
        Set a scaling factor for the nearest-neighbour interaction energy.

        Args:
            energy_scaling (Float): Nearest-neighbour energy scaling.

        Returns:
            None
        """
        self.nn_energy_scaling = energy_scaling

    def set_cn_energies( self, cn_energies ):
        """
        Set the coordination-number dependent energies.

        Args:
            cn_energies (Dict(Str:Dict(Int:Float))): Dict of Dicts specifying the coordination-number dependent energies for each site type.
                e.g. { 'A' : { 0 : 0.0, 1 : 1.0 }, 'B' : { 0 : 0.0, 1 : 2.0 } }
    
        Returns:
            None
        """
        self.cn_energies = cn_energies

    def set_cn_energy_scaling( self, cn_energy_scaling ):
        """
        Set a scaling factor for the coordination-number dependent energies.

        Args:
            cn_energy_scaling (Float): Coordination-number dependent energy scaling.

        Returns:
            None
        """
        self.cn_energy_scaling = cn_energy_scaling

    def set_site_energies( self, site_energies ):
        """
        Set the on-site energies for each site type.

        Args:
            site_energies (Dict(Str:Float)): On-site energies for each site type label.
                e.g. { 'A' : 1.0, 'B', -1.0 }
  
        Returns:
            None
        """
        self.site_energies = site_energies

    def set_number_of_jumps( self, number_of_jumps ):
        """
        Set the on-site energies for each site type.

        Args:
            site_energies (Dict(Str:Float)): On-site energies for each site type label.
                e.g. { 'A' : 1.0, 'B', -1.0 }
  
        Returns:
            None
        """
        self.number_of_jumps = number_of_jumps

    def set_number_of_equilibration_jumps( self, number_of_equilibration_jumps ):
        """
        Set the number of equilibration jumps for a simulation.

        Args:
            number_of_equilibration_jumps (Int): The number of equiibration jumps.
  
        Returns:
            None
        """
        self.number_of_equilibration_jumps = number_of_equilibration_jumps

    def read_lattice_from_file( self, filename ):
        """
        Set the filename for the sites file used to define the simulation lattice.

        Args:
            filename (Str): The filename for the sites file.
  
        Returns:
            None
        """
        self.lattice_site_file = filename

    def set_lattice_cell_lengths( self, cell_lengths ):
        """
        Set the lattice cell lengths for a simulation.

        Args:
            cell_lengths (np.array(x,y,z)): Lattice cell lengths.
  
        Returns:
            None
        """
        self.lattice_cell_lengths = cell_lengths

