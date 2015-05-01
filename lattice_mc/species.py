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

