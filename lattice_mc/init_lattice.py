import numpy as np
from lattice_mc import lattice, lattice_site

def square_lattice( a, b, spacing ):
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
        sites.append( lattice_site.Site( int( it[0] ), r, neighbours, 0.0 ) )
        it.iternext()
    return lattice.Lattice( sites, cell_lengths = np.array( [ a, b, 0.0 ] ) * spacing )

def cubic_lattice( a, b, c, spacing ):
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
        sites.append( lattice_site.Site( int( it[0] ), r, neighbours, 0.0 ) )
        it.iternext()
    return lattice.Lattice( sites, cell_lengths = np.array( [ a, b, c ] ) * spacing )

def lattice_from_sites_file( site_file, cell_lengths ):
    sites = []
    with open( site_file ) as f:
        number_of_sites = int( f.readline() )
        for i in range( number_of_sites ):
            f.readline()
            number = int( f.readline().split()[1] )
            r = np.array( [ float( s ) for s in f.readline().split()[1:4] ] )
            f.readline() # ignore vertices
            neighbours = [ int( s ) for s in f.readline().split()[1:] ]
            sites.append( lattice_site.Site( number, r, neighbours, 0.0 ) )
        return lattice.Lattice( sites, cell_lengths = np.array( cell_lengths ) )
