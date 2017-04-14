import numpy as np

class Species:
    """
    Species class.

    Contains methods that operate on sets of Atom objects
    """

    def __init__( self, atoms ):
        """
        Initialise a Species object.
  
        Args:
            atoms (List(Atom)): A list of Atom objects.

        Returns:
            None
        """
        self.atoms = atoms

    def sites_occupied( self ):
        """
        Returns a list of sites occupied by these atoms.
      
        Args:
            None

        Returns:
            (List): List of sites occupied by these atoms.
        """ 
        return [ atom.site.number for atom in self.atoms ]

    def sum_dr_squared( self ):
        """
        Sum of squared total displacements for these atoms.

        Args:
            None

        Returns:
            (Float): The sum of squared total displacements for these atoms.
        """
        return sum( [ atom.dr_squared() for atom in self.atoms ] )    

    def collective_dr_squared( self ):
        """
        Squared sum of total displacements for these atoms.

        Args:
            None

        Returns:
            (Float): The square of the summed total displacements for these atoms.
        """
        return sum( np.square( sum( [ atom.dr for atom in self.atoms ] ) ) )

    def occupations( self, site_label ):
        """
        Number of these atoms occupying a specific site type.

        Args:
            site_label (Str): Label for the site type being considered.

        Returns:
            (Int): Number of atoms occupying sites of type `site_label`.
        """
        return sum( atom.site.label == site_label for atom in self.atoms )

    def summed_dr2( self ):
        """
        Sum of squared individual displacements for these atoms.

        Args:
            None

        Returns:
            (Float): The sum of squared individual displacements for these atoms.
        """
        return sum( [ atom.summed_dr2 for atom in self.atoms ] )

    def tracer_correlation( self ):
        """
        Tracer correlation factor, f, for these atoms.

        Args:
            None

        Returns:
            (Float): The tracer correlation factor, f, for these atoms.
        """
        return self.sum_dr_squared() / self.summed_dr2()

    def collective_correlation( self ):
        """
        Collective correlation factor, f_I, for these atoms.

        Args:
            None

        Returns:
            (Float): The collective correlation factor, f_I, for these atoms.
        """
        return self.collective_dr_squared() / self.summed_dr2()
