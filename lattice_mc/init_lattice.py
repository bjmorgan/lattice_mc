import numpy as np
from lattice_mc import lattice, lattice_site
from math import sqrt
import re

"""
Functions for constructing lattices.
"""

def square_lattice( a, b, spacing ):
    """
    Generate a square lattice.

    Args:
        a (Int):         Number of lattice repeat units along x.
        b (Int):         Number of lattice repeat units along y.
        spacing (Float): Distance between lattice sites.

    Returns:
        (Lattice): The new lattice

    Notes:
        The returned lattice is 3D periodic, but all sites and edges lie in the xy plane.
    """
    grid = np.array( list( range( 1, a * b + 1 ) ) ).reshape( a, b, order='F' )
    it = np.nditer( grid, flags=['multi_index'] )
    sites = []
    while not it.finished:
        x, y = it.multi_index
        r = np.array( [ x * spacing, y * spacing, 0.0 ] )
        neighbours = [ np.roll( grid, +1, axis=0 )[x,y],
                       np.roll( grid, -1, axis=0 )[x,y],
                       np.roll( grid, +1, axis=1 )[x,y],
                       np.roll( grid, -1, axis=1 )[x,y] ]
        sites.append( lattice_site.Site( int( it[0] ), r, neighbours, 0.0, 'L' ) )
        it.iternext()
    return lattice.Lattice( sites, cell_lengths = np.array( [ a, b, 0.0 ] ) * spacing )

def honeycomb_lattice( a, b, spacing, alternating_sites=False ):
    """
    Generate a honeycomb lattice.

    Args:
        a (Int):         Number of lattice repeat units along x.
        b (Int):         Number of lattice repeat units along y.
        spacing (Float): Distance between lattice sites.
        alternating_sites (Bool, optional): Label alternating sites with 'A' and 'B'. Defaults to False.

    Returns:
        (Lattice): The new lattice

    Notes:
        The returned lattice is 3D periodic, but all sites and edges lie in the xy plane.
    """
    if alternating_sites:
        site_labels = [ 'A', 'B', 'A', 'B' ]
    else:
        site_labels = [ 'L', 'L', 'L', 'L' ]
    unit_cell_lengths = np.array( [ sqrt(3), 3.0, 0.0 ] ) * spacing
    cell_lengths = unit_cell_lengths * np.array( [ a, b, 1.0 ] )
    grid = np.array( list( range( 1, int( a * b * 4 + 1 ) ) ) ).reshape( a, b, 4, order='C' )
    sites = []
    for i in range( a ):
        for j in range( b ):
            # site 1
            r = np.array( [ i * sqrt(3) * spacing, j * 3 * spacing, 0.0 ] )
            neighbours = [ grid[ i, j, 1 ],
                           np.roll( grid, +1, axis=0 )[ i, j, 1 ],
                           np.roll( grid, +1, axis=1 )[ i, j, 3 ] ]
            sites.append( lattice_site.Site( grid[ i, j, 0 ], r, neighbours, 0.0, site_labels[0] ) )
            # site 2
            r = np.array( [ i * sqrt(3) * spacing + sqrt(3)/2 * spacing, ( j * 3 + 0.5 ) * spacing, 0.0 ] )
            neighbours = [ grid[ i, j, 0 ], 
                           grid[ i, j, 2 ], 
                           np.roll( grid, -1, axis=0 )[ i, j, 0 ] ]
            sites.append( lattice_site.Site( grid[ i, j, 1 ], r, neighbours, 0.0, site_labels[1] ) )
            # site 3
            r = np.array( [ i * sqrt(3) * spacing + sqrt(3)/2 * spacing, ( j * 3 + 1.5 ) * spacing, 0.0 ] )
            neighbours = [ grid[ i, j, 1 ],
                           grid[ i, j, 3 ],
                           np.roll( grid, -1, axis=0 )[ i, j, 3 ] ]
            sites.append( lattice_site.Site( grid[ i, j, 2 ], r, neighbours, 0.0, site_labels[2] ) )
            # site 4
            r = np.array( [ i * sqrt(3) * spacing, ( j * 3 + 2 ) * spacing, 0.0 ] )
            neighbours = [ grid[ i, j, 2 ], 
                           np.roll( grid, +1, axis=0 )[ i, j, 2 ],
                           np.roll( grid, -1, axis=1 )[ i, j, 0 ] ]
            sites.append( lattice_site.Site( grid[ i, j, 3 ], r, neighbours, 0.0, site_labels[3] ) )
    return lattice.Lattice( sites, cell_lengths=cell_lengths )

def cubic_lattice( a, b, c, spacing ):
    """
    Generate a cubic lattice.

    Args:
        a (Int):         Number of lattice repeat units along x.
        b (Int):         Number of lattice repeat units along y.
        c (Int):         Number of lattice repeat units along z.
        spacing (Float): Distance between lattice sites.

    Returns:
        (Lattice): The new lattice
    """
    grid = np.array( list( range( 1, a * b * c + 1 ) ) ).reshape( a, b, c, order='F' )
    it = np.nditer( grid, flags=[ 'multi_index' ] )
    sites = []
    while not it.finished:
        x, y, z = it.multi_index
        r = np.array( [ x, y, z ] ) * spacing
        neighbours = [ np.roll( grid, +1, axis=0 )[x,y,z],
                       np.roll( grid, -1, axis=0 )[x,y,z],
                       np.roll( grid, +1, axis=1 )[x,y,z],
                       np.roll( grid, -1, axis=1 )[x,y,z],
                       np.roll( grid, +1, axis=2 )[x,y,z],
                       np.roll( grid, -1, axis=2 )[x,y,z] ]
        sites.append( lattice_site.Site( int( it[0] ), r, neighbours, 0.0, 'L' ) )
        it.iternext()
    return lattice.Lattice( sites, cell_lengths = np.array( [ a, b, c ] ) * spacing )

#def lattice_from_sites_file( site_file, cell_lengths ):
#    sites = []
#    with open( site_file ) as f:
#        number_of_sites = int( f.readline() )
#        for i in range( number_of_sites ):
#            f.readline()
#            number = int( f.readline().split()[1] )
#            r = np.array( [ float( s ) for s in f.readline().split()[1:4] ] )
#            f.readline() # ignore vertices
#            neighbours = [ int( s ) for s in f.readline().split()[1:] ]
#            label = f.readline().split()[1]
#            sites.append( lattice_site.Site( number, r, neighbours, 0.0, label ) )
#        return lattice.Lattice( sites, cell_lengths = np.array( cell_lengths ) )

def lattice_from_sites_file( site_file, cell_lengths ):
    """
    Generate a lattice from a sites file.

    Args:
        site_file (Str): Filename for the file containing the site information.
        cell_lengths (List(Float,Float,Float)): A list containing the [ x, y, z ] cell lengths.

    Returns:
        (Lattice): The new lattice

    Notes:
        | The site information file format is:
        |     <number_of_sites> (Int).
        | Followed by blocks of data separated by double linebreaks; one block per site.
        |     site: <site number> (Int).
        |     centre: <x> <y> <z> (Float,Float,Float).
        |     neighbours: <list of site numbers of neighbouring sites> (List[Int]).
        |     label: <site group labal> (Str).
        |     energy: <site occupation energy> (Float).
        | The energy is optional, and will be set to 0.0 if not included.
        | Line order within each block is not meaningful.
        | British and American English spellings for centre|center and neighbour|neighbor are accepted.
        | An example file can be found in the examples directory. 
    """
    sites = []
    site_re = re.compile( 'site:\s+([-+]?\d+)' )
    r_re = re.compile( 'cent(?:er|re):\s+([-\d\.e]+)\s+([-\d\.e]+)\s+([-\d\.e]+)' )
    r_neighbours = re.compile( 'neighbou{0,1}rs:((\s+[-+]?\d+)+)' )
    r_label = re.compile( 'label:\s+(\S+)' )
    r_energy = re.compile( 'energy:\s([-+\d\.]+)' )
    with open( site_file ) as f:
        filein = f.read().split("\n\n")
    number_of_sites = int( filein[0] )
    for i in range( number_of_sites ):
        block = filein[ i+1 ]
        number = int( site_re.findall( block )[0] )
        r = np.array( [ float(s) for s in r_re.findall( block )[0] ] )
        neighbours = [ int( s ) for s in r_neighbours.findall( block )[0][0].split() ]
        label = r_label.findall( block )[0]
        energy = r_energy.findall( block )
        if energy:
            energy = float( energy[0] )
        else:
            energy = 0.0
        sites.append( lattice_site.Site( number, r, neighbours, energy, label ) )
    return lattice.Lattice( sites, cell_lengths = np.array( cell_lengths ) )
