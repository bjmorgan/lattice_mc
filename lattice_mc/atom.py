import numpy as np

class Atom:

    @classmethod
    def sites_occupied( cls, atoms ):
        return [ atom.site.number for atom in atoms ]
        
    atom_number = 0

    def __init__( self, initial_site ):
        self.site = initial_site
        # check this site is not already occupied
        Atom.atom_number += 1
        self.number = Atom.atom_number
        if self.site.occupation == 0:
            self.site.occupation = self.number
            self.site.is_occupied = True
            self.site.atom = self
        else:
            print( "Error: this site is already occupied" ) 
        self.number_of_hops = 0
        self.dr = np.array( [ 0.0, 0.0, 0.0 ] )
        self.summed_dr2 = 0.0
        self.sites_visited = [ self.site.number ]

    def dr_squared( self ):
        return( np.dot( self.dr, self.dr ) )

    def correlation_factor( self ):
        #print( self.dr )
        rk_2 = np.sum( self.dr * self.dr )
        return( rk_2 / self.summed_dr2 )

