class Options:

    def __init__( self ):
        self.site_energies = None
        self.nn_energy_scaling = None
        self.number_of_equilibration_jumps = 0

    def set_number_of_atoms( self, number_of_atoms ):
        self.number_of_atoms = number_of_atoms

    def set_nn_energy_scaling( self, energy_scaling ):
        self.nn_energy_scaling = energy_scaling

    def set_cn_energies( self, cn_energies ):
        self.cn_energies = cn_energies

    def set_cn_energy_scaling( self, cn_energy_scaling ):
        self.cn_energy_scaling = cn_energy_scaling

    def set_site_energies( self, site_energies ):
        self.site_energies = site_energies

    def set_number_of_jumps( self, number_of_jumps ):
        self.number_of_jumps = number_of_jumps

    def set_number_of_equilibration_jumps( self, number_of_equilibration_jumps ):
        self.number_of_equilibration_jumps = number_of_equilibration_jumps

    def read_lattice_from_file( self, filename ):
        self.lattice_site_file = filename

    def set_lattice_cell_lengths( self, cell_lengths ):
        self.lattice_cell_lengths = cell_lengths

