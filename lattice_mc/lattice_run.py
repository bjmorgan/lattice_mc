#! /usr/bin/env python3

import numpy as np
import math
import random
import bisect
from lattice_mc import atom, lattice_site, transitions, jump, init_lattice, lattice

def main():
    #lattice = square_lattice( 10, 10, 1.0 )
    lattice = init_lattice.cubic_lattice( 10, 10, 10, 1.0 )
    #lattice = lattice_from_sites_file( 'llzo_lattice_site_list.dat', cell_lengths = [ 49.0672361, 49.0672361, 49.0672361 ] )
    print( "lattice initialised" )
    #atoms = lattice.populate_sites( lattice.number_of_sites() - 1 )
    atoms = lattice.populate_sites( 10 )

    number_of_jumps = 500000
    for step in range( number_of_jumps ):
        all_transitions = transitions.Transitions( lattice.potential_jumps() )
        random_jump = all_transitions.random()
        jumping_atom = random_jump.initial_site.atom
        lattice.update( random_jump )

    # calculate correlation factor:
    sum_dr_squared = sum( [ atom.dr_squared() for atom in atoms ] )
    print( sum_dr_squared / float( number_of_jumps ) )

if __name__ == '__main__':
    main()
