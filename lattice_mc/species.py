import numpy as np

class Species:

    def __init__( self, atoms ):
        self.atoms = atoms

    def sites_occupied( self ):
        return [ atom.site.number for atom in self.atoms ]

    def sum_dr_squared( self ):
        return sum( [ atom.dr_squared() for atom in self.atoms ] )    

    def collective_dr_squared( self ):
        return sum( np.square( sum( [ atom.dr for atom in self.atoms ] ) ) )

    def occupations( self, site_label ):
        return sum( atom.site.label == site_label for atom in self.atoms )

    def summed_dr2( self ):
        return sum( [ atom.summed_dr2 for atom in self.atoms ] )

    def tracer_correlation( self ):
        return self.sum_dr_squared() / self.summed_dr2()

    def collective_correlation( self ):
        return self.collective_dr_squared / self.summed_dr2()
