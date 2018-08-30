import numpy as np
import random
import itertools
import sys

from lattice_mc import atom, jump, transitions, cluster
from lattice_mc.error import BlockedLatticeError
from collections import Counter

class Lattice:
    """
    Lattice class
    """

    def __init__( self, sites, cell_lengths ):
        """
        Initialise a Lattice instance.

        Args:
            sites (List(Site)): List of sites contained in the lattice.
            cell_lengths (np.array(x,y,z)): Vector of cell lengths for the simulation cell.

        Returns:
            None
        """
        self.cell_lengths = cell_lengths
        self.sites = sites
        self.number_of_sites = len( self.sites )
        self.site_labels = set( [ site.label for site in self.sites ] )
        self.site_populations = Counter( [ site.label for site in self.sites ] )
        self.enforce_periodic_boundary_conditions()
        self.initialise_site_lookup_table()
        self.nn_energy = False
        self.cn_energies = False
        self.site_energies = False
        self.jump_lookup_table = False
        for site in self.sites:
            site.p_neighbours = [ self.site_with_id( i ) for i in site.neighbours ]
        self.reset()

    def enforce_periodic_boundary_conditions( self ):
        """
        Ensure that all lattice sites are within the central periodic image of the simulation cell.
        Sites that are outside the central simulation cell are mapped back into this cell.
        
        Args:
            None

        Returns:
            None
        """
        for s in self.sites:
            for i in range(3):
                if s.r[i] < 0.0:
                    s.r[i] += self.cell_lengths[i]
                if s.r[i] > self.cell_lengths[i]:
                    s.r[i] -= self.cell_lengths[i]

    def reset( self ):
        """
        Reset all time-dependent counters for this lattice and its constituent sites
   
        Args:
            None

        Returns:
            None
        """
        self.time = 0.0
        for site in self.sites:
            site.time_occupied = 0.0
          
    def initialise_site_lookup_table( self ):
        """
        Create a lookup table allowing sites in this lattice to be queried using `self.site_lookup[n]` where `n` is the identifying site numbe.

        Args:
            None

        Returns:
            None
        """
        self.site_lookup = {}
        for site in self.sites:
            self.site_lookup[ site.number ] = site

    def site_with_id( self, number ):
        """
        Select the site with a specific id number.

        Args:
            number (Int): The identifying number for a specific site.

        Returns:
            (Site): The site with id number equal to `number`
        """
        return self.site_lookup[ number ]

    def vacant_sites( self ):
        """
        The set of sites not occupied by atoms.

        Args:
            None

        Returns:
            List(Site): List of sites that are vacant.
        """
        return ( site for site in self.sites if not site.is_occupied )

    def occupied_sites( self ):
        """
        The set of sites occupied by atoms.

        Args:
            None

        Returns:
            List(Site): List of sites that are occupied.
        """
        return ( site for site in self.sites if site.is_occupied )

    def vacant_site_numbers( self ):
        """
        List of site id numbers for all sites that are vacant.

        Args:
            None
        
        Returns:
            List(Int): List of site id numbers for vacant sites.
        """
        return [ site.number for site in self.sites if not site.is_occupied ]

    def occupied_site_numbers( self ):
        """
        List of site id numbers for all sites that are occupied.

        Args:
            None
        
        Returns:
            List(Int): List of site id numbers for occupied sites.
        """
        return [ site.number for site in self.sites if site.is_occupied ]
        
    def potential_jumps( self ):
        """
        All nearest-neighbour jumps not blocked by volume exclusion 
        (i.e. from occupied to neighbouring unoccupied sites).

        Args:
            None

        Returns:
            (List(Jump)): List of possible jumps.
        """
        jumps = []
        if self.number_of_occupied_sites <= self.number_of_sites / 2:
            for occupied_site in self.occupied_sites():
                unoccupied_neighbours = [ site for site in [ self.site_with_id( n ) for n in occupied_site.neighbours ] if not site.is_occupied ]
                for vacant_site in unoccupied_neighbours:
                    jumps.append( jump.Jump( occupied_site, vacant_site, self.nn_energy, self.cn_energies, self.jump_lookup_table ) )
        else:
            for vacant_site in self.vacant_sites():
                occupied_neighbours = [ site for site in [ self.site_with_id( n ) for n in vacant_site.neighbours ] if site.is_occupied ]
                for occupied_site in occupied_neighbours:
                    jumps.append( jump.Jump( occupied_site, vacant_site, self.nn_energy, self.cn_energies, self.jump_lookup_table ) )
        return jumps

    def update( self, jump ):
        """
        Update the lattice state by accepting a specific jump

        Args:
            jump (Jump): The jump that has been accepted.

        Returns:
            None.
        """
        atom = jump.initial_site.atom
        dr = jump.dr( self.cell_lengths )
        #print( "atom {} jumped from site {} to site {}".format( atom.number, jump.initial_site.number, jump.final_site.number ) )
        jump.final_site.occupation = atom.number
        jump.final_site.atom = atom
        jump.final_site.is_occupied = True
        jump.initial_site.occupation = 0
        jump.initial_site.atom = None
        jump.initial_site.is_occupied = False
        # TODO: updating atom counters could be contained in an atom.move_to( site ) method
        atom.site = jump.final_site
        atom.number_of_hops += 1
        atom.dr += dr
        atom.summed_dr2 += np.dot( dr, dr )

    def populate_sites( self, number_of_atoms, selected_sites=None ):
        """
        Populate the lattice sites with a specific number of atoms.

        Args:
            number_of_atoms (Int): The number of atoms to populate the lattice sites with.
            selected_sites (:obj:List, optional): List of site labels if only some sites are to be occupied. Defaults to None.

        Returns:
            None
        """
        if number_of_atoms > self.number_of_sites:
            raise ValueError
        if selected_sites:
            atoms = [ atom.Atom( initial_site = site ) for site in random.sample( [ s for s in self.sites if s.label in selected_sites ], number_of_atoms ) ]
        else:
            atoms = [ atom.Atom( initial_site = site ) for site in random.sample( self.sites, number_of_atoms ) ]
        self.number_of_occupied_sites = number_of_atoms
        return atoms

    def jump( self ):
        """
        Select a jump at random from all potential jumps, then update the lattice state.

        Args:
            None

        Returns:
            None
        """
        potential_jumps = self.potential_jumps()
        if not potential_jumps:
            raise BlockedLatticeError('No moves are possible in this lattice')
        all_transitions = transitions.Transitions( self.potential_jumps() )
        random_jump = all_transitions.random()
        delta_t = all_transitions.time_to_jump()
        self.time += delta_t
        self.update_site_occupation_times( delta_t )
        self.update( random_jump )
        return( all_transitions.time_to_jump() )

    def update_site_occupation_times( self, delta_t ):
        """
        Increase the time occupied for all occupied sites by delta t

        Args:
            delta_t (Float): Timestep.

        Returns:
            None
        """
        for site in self.occupied_sites():
            site.time_occupied += delta_t

    def site_occupation_statistics( self ):
        """
        Average site occupation for each site type

        Args:
            None

        Returns:
            (Dict(Str:Float)): Dictionary of occupation statistics, e.g.::

                { 'A' : 2.5, 'B' : 25.3 } 
        """
        if self.time == 0.0:
            return None
        occupation_stats = { label : 0.0 for label in self.site_labels }
        for site in self.sites:
            occupation_stats[ site.label ] += site.time_occupied
        for label in self.site_labels:
            occupation_stats[ label ] /= self.time
        return occupation_stats
     
    def set_site_energies( self, energies ):
        """
        Set the energies for every site in the lattice according to the site labels.

        Args:
            energies (Dict(Str:Float): Dictionary of energies for each site label, e.g.::

                { 'A' : 1.0, 'B', 0.0 }

        Returns:
            None
        """
        self.site_energies = energies
        for site_label in energies:
            for site in self.sites:
                if site.label == site_label:
                    site.energy = energies[ site_label ]

    def set_nn_energy( self, delta_E ):
        """
        Set the lattice nearest-neighbour energy.

        Args:
            delta_E (Float): The nearest-neighbour energy E_nn.

        Returns:
            None
        """
        self.nn_energy = delta_E

    def set_cn_energies( self, cn_energies ):
        """
        Set the coordination number dependent energies for this lattice.

        Args:
            cn_energies (Dict(Str:Dict(Int:Float))): Dictionary of dictionaries specifying the coordination number dependent energies for each site type. e.g.::

                { 'A' : { 0 : 0.0, 1 : 1.0, 2 : 2.0 }, 'B' : { 0 : 0.0, 1 : 2.0 } }

        Returns:
            None
        """
        for site in self.sites:
            site.set_cn_occupation_energies( cn_energies[ site.label ] )
        self.cn_energies = cn_energies

    def site_coordination_numbers( self ):
        """
        Returns a dictionary of the coordination numbers for each site label. e.g.::
        
            { 'A' : { 4 }, 'B' : { 2, 4 } }
 
        Args:
            none

        Returns:
            coordination_numbers (Dict(Str:Set(Int))): dictionary of coordination
                                                       numbers for each site label.
        """
        coordination_numbers = {}
        for l in self.site_labels:
            coordination_numbers[ l ] = set( [ len( site.neighbours ) for site in self.sites if site.label is l ] ) 
        return coordination_numbers 

    def max_site_coordination_numbers( self ):
        """
        Returns a dictionary of the maximum coordination number for each site label.
        e.g.::
   
            { 'A' : 4, 'B' : 4 }

        Args:
            none

        Returns:
            max_coordination_numbers (Dict(Str:Int)): dictionary of maxmimum coordination
                                                      number for each site label.
        """
        return { l : max( c ) for l, c in self.site_coordination_numbers().items() }
       
    def site_specific_coordination_numbers( self ):
        """
        Returns a dictionary of coordination numbers for each site type.

        Args:
            None

        Returns:
            (Dict(Str:List(Int))) : Dictionary of coordination numbers for each site type, e.g.::

                { 'A' : [ 2, 4 ], 'B' : [ 2 ] }
        """
        specific_coordination_numbers = {}
        for site in self.sites:
            specific_coordination_numbers[ site.label ] = site.site_specific_neighbours()
        return specific_coordination_numbers

    def connected_site_pairs( self ):
        """
        Returns a dictionary of all connections between pair of sites (by site label).
        e.g. for a linear lattice A-B-C will return::
        
            { 'A' : [ 'B' ], 'B' : [ 'A', 'C' ], 'C' : [ 'B' ] }

        Args:
            None

        Returns:
            site_connections (Dict{Str List[Str]}): A dictionary of neighbouring site types in the lattice.
        """
        site_connections = {}
        for initial_site in self.sites:
            if not initial_site.label in site_connections:
                site_connections[ initial_site.label ] = []
            for final_site in initial_site.p_neighbours:
                if final_site.label not in site_connections[ initial_site.label ]:
                    site_connections[ initial_site.label ].append( final_site.label )
        return site_connections

    def transmute_sites( self, old_site_label, new_site_label, n_sites_to_change ):
        """
        Selects a random subset of sites with a specific label and gives them a different label.

        Args:
            old_site_label (String or List(String)): Site label(s) of the sites to be modified..
            new_site_label (String):                 Site label to be applied to the modified sites.
            n_sites_to_change (Int):                 Number of sites to modify.

        Returns:
            None
        """
        selected_sites = self.select_sites( old_site_label )
        for site in random.sample( selected_sites, n_sites_to_change ):
            site.label = new_site_label
        self.site_labels = set( [ site.label for site in self.sites ] )

    def connected_sites( self, site_labels=None ):
        """
        Searches the lattice to find sets of sites that are contiguously neighbouring.
        Mutually exclusive sets of contiguous sites are returned as Cluster objects.

        Args:
            site_labels (:obj:(List(Str)|Set(Str)|Str), optional): Labels for sites to be considered in the search.
                This can be a list::

                    [ 'A', 'B' ]

                a set::

                    ( 'A', 'B' )

                or a string::

                    'A'.

        Returns:
            (List(Cluster)): List of Cluster objects for groups of contiguous sites.
        """
        if site_labels:
           selected_sites = self.select_sites( site_labels )
        else:
           selected_sites = self.sites
        initial_clusters = [ cluster.Cluster( [ site ] ) for site in selected_sites ]
        if site_labels:
            blocking_sites = self.site_labels - set( site_labels )
            for c in initial_clusters:
                c.remove_sites_from_neighbours( blocking_sites )
        final_clusters = []
        while initial_clusters: # loop until initial_clusters is empty
            this_cluster = initial_clusters.pop(0)
            while this_cluster.neighbours:
                neighbouring_clusters = [ c for c in initial_clusters if this_cluster.is_neighbouring( c ) ] 
                for nc in neighbouring_clusters:
                    initial_clusters.remove( nc )
                    this_cluster = this_cluster.merge( nc ) 
            final_clusters.append( this_cluster )
        return final_clusters

    def select_sites( self, site_labels ):
        """
        Selects sites in the lattice with specified labels.

        Args:
            site_labels (List(Str)|Set(Str)|Str): Labels of sites to select.
                This can be a List [ 'A', 'B' ], a Set ( 'A', 'B' ), or a String 'A'.

        Returns:
            (List(Site)): List of sites with labels given by `site_labels`.
        """
        if type( site_labels ) in ( list, set ):
            selected_sites = [ s for s in self.sites if s.label in site_labels ]
        elif type( site_labels ) is str:
            selected_sites = [ s for s in self.sites if s.label is site_labels ]
        else:
            raise ValueError( str( site_labels ) )
        return selected_sites

    def detached_sites( self, site_labels=None ):
        """
        Returns all sites in the lattice (optionally from the set of sites with specific labels)
        that are not part of a percolating network.
        This is determined from clusters of connected sites that do not wrap round to
        themselves through a periodic boundary.

        Args:
            site_labels (String or List(String)): Lables of sites to be considered.

        Returns:
            (List(Site)): List of sites not in a periodic percolating network.
        """
        clusters = self.connected_sites( site_labels=site_labels )
        island_clusters = [ c for c in clusters if not any( c.is_periodically_contiguous() ) ]
        return list( itertools.chain.from_iterable( ( c.sites for c in island_clusters ) ) )

    def is_blocked( self ):
        """
        Check whether there are any possible jumps.

        Args:
            None

        Returns:
            (Bool): True if there are no possible jumps. Otherwise returns False.
        """
        if not self.potential_jumps():
            return True
        else:
            return False
