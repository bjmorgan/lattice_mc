class Cluster:
    """
    Clusters are sets of sites.
    """
    def __init__( self, sites ):
        """
        Initialise an Cluster instance.

        Args:
            sites (List(Site): The list of sites that make up the cluster.

        Returns:
            None
        """
        self.sites = set( sites )
        self.neighbours = set()
        for s in self.sites:
            self.neighbours.update( s.p_neighbours ) 
        self.neighbours = self.neighbours.difference( self.sites )

    def merge( self, other_cluster ):
        """
        Combine two clusters into a single cluster.

        Args:
            other_cluster (Cluster): The second cluster to combine.

        Returns:
            (Cluster):   The combination of both clusters.
        """
        new_cluster = Cluster( self.sites | other_cluster.sites )
        new_cluster.neighbours = ( self.neighbours | other_cluster.neighbours ).difference( new_cluster.sites )
        return new_cluster

    def is_neighbouring( self, other_cluster ):
        """
        logical check whether the neighbour list for cluster A includes any sites in cluster B
        
        Args:
            other_cluster (Cluster): the other cluster we are testing for neighbour connectivity

        Returns:
            (Bool): True if the other cluster neighbours this cluster.
        """
        return bool( self.neighbours & other_cluster.sites )

    def size( self ):
        """
        Number of sites in this cluster.

        Args:
            None

        Returns:
            (Int): The number of sites in this cluster.
        """
        return len( self.sites )

    def sites_at_edges( self ):
        """
        Finds the six sites with the maximum and minimum coordinates along x, y, and z.

        Args:
            None

        Returns:
             (List(List)): In the order [ +x, -x, +y, -y, +z, -z ]
        """
        min_x = min( [ s.r[0] for s in self.sites ] )
        max_x = max( [ s.r[0] for s in self.sites ] )
        min_y = min( [ s.r[1] for s in self.sites ] )
        max_y = max( [ s.r[1] for s in self.sites ] )
        min_z = min( [ s.r[2] for s in self.sites ] )
        max_z = max( [ s.r[2] for s in self.sites ] )
        x_max = [ s for s in self.sites if s.r[0] == min_x ]
        x_min = [ s for s in self.sites if s.r[0] == max_x ]
        y_max = [ s for s in self.sites if s.r[1] == min_y ]
        y_min = [ s for s in self.sites if s.r[1] == max_y ]
        z_max = [ s for s in self.sites if s.r[2] == min_z ]
        z_min = [ s for s in self.sites if s.r[2] == max_z ]
        return ( x_max, x_min, y_max, y_min, z_max, z_min )
   
    def is_periodically_contiguous( self ):
        """
        logical check whether a cluster connects with itself across the
        simulation periodic boundary conditions.

        Args:
            none

        Returns
            ( Bool, Bool, Bool ): Contiguity along the x, y, and z coordinate axes
        """
        edges = self.sites_at_edges()
        is_contiguous = [ False, False, False ]
        along_x = any( [ s2 in s1.p_neighbours for s1 in edges[0] for s2 in edges[1] ] )
        along_y = any( [ s2 in s1.p_neighbours for s1 in edges[2] for s2 in edges[3] ] )
        along_z = any( [ s2 in s1.p_neighbours for s1 in edges[4] for s2 in edges[5] ] )
        return ( along_x, along_y, along_z )

    def remove_sites_from_neighbours( self, remove_labels ):
        """
        Removes sites from the set of neighbouring sites if these have labels in remove_labels.

        Args:
            Remove_labels (List) or (Str): List of Site labels to be removed from the cluster neighbour set.

        Returns:
            None
        """
        if type( remove_labels ) is str:
            remove_labels = [ remove_labels ]
        self.neighbours = set( n for n in self.neighbours if n.label not in remove_labels ) 

