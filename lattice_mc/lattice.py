import numpy as np
import random

from lattice_mc import atom, jump

class Lattice:

    def __init__( self, sites, cell_lengths ):
        self.cell_lengths = cell_lengths
        self.sites = sites
        self.enforce_periodic_boundary_conditions()
        self.initialise_site_lookup_table()

    def enforce_periodic_boundary_conditions( self ):
        for s in self.sites:
            for i in range(3):
                if s.r[i] < 0.0:
                    s.r[i] += self.cell_lengths[i]
                if s.r[i] > self.cell_lengths[i]:
                    s.r[i] -= self.cell_lengths[i]

    def initialise_site_lookup_table( self ):
        self.site_lookup = {}
        for site in self.sites:
            self.site_lookup[ site.number ] = site

    def site_with_id( self, number ):
        return self.site_lookup[ number ]

    def vacant_sites( self ):
        return [ site for site in self.sites if not site.is_occupied ]

    def vacant_site_numbers( self ):
        return [ site.number for site in self.sites if not site.is_occupied ]

    def occupied_sites( self ):
        return [ site for site in self.sites if site.is_occupied ]

    def occupied_site_numbers( self ):
        return [ site.number for site in self.sites if site.is_occupied ]
        
    def number_of_sites( self ):
        return len( self.sites )

    def potential_jumps( self ):
        '''return all possible nearest-neighbour jumps (from occupied to neighbouring unoccupied sites)'''
        jumps = []
        for vacant_site in self.vacant_sites():
            occupied_neighbours = [ site for site in [ self.site_with_id( n ) for n in vacant_site.neighbours ] if site.is_occupied ]
            for occupied_site in occupied_neighbours:
                jumps.append( jump.Jump( occupied_site, vacant_site ) )
        return jumps

    def update( self, jump ):
        atom = jump.initial_site.atom
        dr = jump.dr( self.cell_lengths )
        #print( "atom {} jumped from site {} to site {}".format( atom.number, jump.initial_site.number, jump.final_site.number ) )
        jump.final_site.occupation = atom.number
        jump.final_site.atom = atom
        jump.final_site.is_occupied = True
        jump.initial_site.occupation = 0
        jump.initial_site.atom = None
        jump.initial_site.is_occupied = False
        atom.site = jump.final_site
        atom.number_of_hops += 1
        atom.dr += dr
        atom.summed_dr2 += np.dot( dr, dr )

    def populate_sites( self, number_of_atoms ):
        atoms = [ atom.Atom( initial_site = site ) for site in random.sample( self.sites, number_of_atoms ) ]
        return atoms
