from lattice_mc import init_lattice, species, lookup_table

class Simulation:
    """
    Simulation class
    """

    def __init__( self ):
        """
        Initialise a Simulation object.

        Args:
            None

        Returns:
            None
        
        Notes:
            Simulation parameters need to be set using their corresponding setter methods.
        """
        self.lattice = None
        self.number_of_atoms = None
        self.number_of_jumps = None
        self.for_time = None
        self.number_of_equilibration_jumps = 0
        self.atoms = None
        self.has_run = False

    def reset( self ):
        """
        Reset all counters for this simulation.

        Args:
            None

        Returns:
            None
        """
        self.lattice.reset()
        for atom in self.atoms.atoms:
            atom.reset()
 
    def set_number_of_atoms( self, n, selected_sites=None ):
        """
        Set the number of atoms for the simulation, and populate the simulation lattice.

        Args:
            n (Int): Number of atoms for this simulation.
            selected_sites (:obj:(List|Set|String), optional): Selects a subset of site types to be populated with atoms. Defaults to None.

        Returns:
            None
        """
        self.number_of_atoms = n
        self.atoms = species.Species( self.lattice.populate_sites( self.number_of_atoms, selected_sites=selected_sites ) )

    def set_number_of_jumps( self, n ):
        """
        Set the number of jumps for this simulation.

        Args:
            n (Int): number of jumps

        Returns:
            None
        """
        self.number_of_jumps = n 

    def set_number_of_equilibration_jumps( self, n ):
        """
        Set the number of equilibration jumps for this simulation.

        Args:
            n (Int): number of equilibration jumps

        Returns:
            None
        """
        self.number_of_equilibration_jumps = n 

    def define_lattice_from_file( self, filename, cell_lengths ):
        """
        Set up the simulation lattice from a file containing site data.
            Uses `init_lattice.lattice_from_sites_file`, which defines the site file spec.
  
        Args:
            filename (Str): sites file filename.
            cell_lengths (List(x,y,z)): cell lengths for the simulation cell.

        Returns:
            None
        """
        self.lattice = init_lattice.lattice_from_sites_file( filename, cell_lengths = cell_lengths )

    def set_nn_energy( self,  nn_energy ):
        """
        Set the nearest-neighbour energy for this simulation.

        Args:
            nn_energy (Float): nearest-neigbour energy.

        Returns:
            None
        """
        if nn_energy:
            self.lattice.set_nn_energy( nn_energy )

    def set_cn_energies( self, cn_energies ):
        """
        Set the coordination-number dependent energies for this simulation.

        Args:
            cn_energies	(Dict(Str:Dict(Int:Float))): Dict of Dicts specficying the coordination-number dependent energies.
                e.g. { 'A' : { 0 : 0.0, 1: 0.5 }, 'B' : { 0 : -0.5, 1 : -2.0 } }

        Returns:
            None	
        """
        if cn_energies:
            self.lattice.set_cn_energies( cn_energies )

    def set_site_energies( self, site_energies ):
        """
        Set the on-site energies for this simulation.

        Args:
            site_energies (Dict(Str:Float)): Dict of energies for each site type.
                e.g. { 'A' : 0.0,. 'B' : 1.0 }
        Returns:
            None
        """
        if site_energies:
            self.lattice.set_site_energies( site_energies )

    def is_initialised( self ):
        """
        Check whether the simulation has been initialised.

        Args:
            None

        Returns:
            None
        """
        if not self.lattice:
            raise AttributeError('Running a simulation needs the lattice to be initialised')
        if not self.atoms:
            raise AttributeError('Running a simulation needs the atoms to be initialised')
        if not self.number_of_jumps and not self.for_time:
            raise AttributeError('Running a simulation needs number_of_jumps or for_time to be set')
 
    def run( self, for_time=None ):
        """
        Run the simulation.

        Args:
            for_time (:obj:Float, optional): If `for_time` is set, then run the simulation until a set amount of time has passed. Otherwise, run the simulation for a set number of jumps. Defaults to None.

        Returns:
            None
        """
        self.for_time = for_time
        try:
            self.is_initialised()
        except AttributeError:
            raise
        if self.number_of_equilibration_jumps > 0:
            for step in range( self.number_of_equilibration_jumps ):
                self.lattice.jump()
            self.reset()
        if self.for_time:
            self.number_of_jumps = 0
            while self.lattice.time < self.for_time:
                self.lattice.jump()
                self.number_of_jumps += 1
        else: 
            for step in range( self.number_of_jumps ):
                self.lattice.jump()
        self.has_run = True

    @property
    def old_tracer_correlation( self ):
        """
        Deprecated tracer correlation factor for this simulation.
       
        Args:
            None

        Returns:
            (Float): The tracer correlation factor, f.

        Notes:
            This function assumes that the jump distance between sites has
            been normalised to a=1. If the jump distance is not equal to 1
            then the value returned by this function should be divided by a^2.
            Even better, use `self.tracer_correlation`.
        """
        if self.has_run:
            return self.atoms.sum_dr_squared() / float( self.number_of_jumps )
        else:
            return None
    
    @property
    def tracer_correlation( self ):
        """
        Tracer correlation factor, f.

        Args:
            None

        Returns:
            (Float): The tracer correlation factor, f.
        """
        if self.has_run:
            return self.atoms.tracer_correlation()
        else:
            return None

    @property
    def tracer_diffusion_coefficient( self ):
        """
        Tracer diffusion coefficient, D*.

        Args:
            None

        Returns:
            (Float): The tracer diffusion coefficient, D*.
        """
        if self.has_run:
            return self.atoms.sum_dr_squared() / ( 6.0 * float( self.number_of_atoms ) * self.lattice.time )
        else:
            return None

    @property
    def old_collective_correlation( self ):
        """
        Returns the collective correlation factor, f_I

        Args:
            None

        Returns:
            (Float): The collective correlation factor, f_I.

        Notes:
            This function assumes that the jump distance between sites has
            been normalised to a=1. If the jumps distance is not equal to 1
            then the value returned by this function should be divided by a^2.
            Even better, use self.collective_correlation
        """
        if self.has_run:
            return self.atoms.collective_dr_squared() / float( self.number_of_jumps )
        else:
            return None

    @property
    def collective_correlation( self ):
        """
        Returns the collective correlation factor, f_I

        Args: 
            None

        Returns: 
            (Float): The collective correlation factor, f_I.
        """
        if self.has_run:
            return self.atoms.collective_correlation()
        else:
            return None

    @property
    def collective_diffusion_coefficient( self ):
        """
        Returns the collective or "jump" diffusion coefficient, D_J.

        Args:
            None

        Returns:
            (Float): The collective diffusion coefficient, D_J.
        """
        if self.has_run:
            return self.atoms.collective_dr_squared() / ( 6.0 * self.lattice.time )
        else:
            return None

    @property
    def collective_diffusion_coefficient_per_atom( self ):
        """
        The collective diffusion coefficient per atom. D_J / n_atoms.

        Args:
            None

        Returns:
            (Float): The collective diffusion coefficient per atom, D_J / n_atoms.
        """
        if self.has_run:
            return self.collective_diffusion_coefficient / float( self.number_of_atoms )
        else:
            return None

    @property
    def average_site_occupations( self ):
        """
        Average site occupation numbers.

        Args:
            None

        Returns:
            (Dict(Str:Float)): Average site occupation numbers for each site label,
                e.g. { 'A' : 12.4, 'B' : 231.2 }
        """
        return self.lattice.site_occupation_statistics()

    def setup_lookup_table( self, hamiltonian='nearest-neighbour' ):
        """
        Create a jump-probability look-up table corresponding to the appropriate Hamiltonian.

        Args:
            hamiltonian (Str, optional): String specifying the simulation Hamiltonian.
                valid values are 'nearest-neighbour' (default) and 'coordination_number'.

        Returns:
            None
        """
        expected_hamiltonian_values = [ 'nearest-neighbour', 'coordination_number' ]
        if hamiltonian not in expected_hamiltonian_values:
            raise ValueError
        self.lattice.jump_lookup_table = lookup_table.LookupTable( self.lattice, hamiltonian )

