from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from lattice_mc.atom import Atom


class Species:
    """
    Species class.

    Contains methods that operate on sets of Atom objects
    """

    def __init__(self, atoms: list[Atom]) -> None:
        """
        Initialise a Species object.

        Args:
            atoms (List(Atom)): A list of Atom objects.

        Returns:
            None
        """
        self.atoms: list[Atom] = atoms

    def sites_occupied(self) -> list[int]:
        """
        Returns a list of sites occupied by these atoms.

        Args:
            None

        Returns:
            (List): List of sites occupied by these atoms.
        """
        return [atom.site.number for atom in self.atoms]

    def sum_dr_squared(self) -> float:
        """
        Sum of squared total displacements for these atoms.

        Args:
            None

        Returns:
            (Float): The sum of squared total displacements for these atoms.
        """
        return float(sum([atom.dr_squared() for atom in self.atoms]))

    def collective_dr_squared(self) -> float:
        """
        Squared sum of total displacements for these atoms.

        Args:
            None

        Returns:
            (Float): The square of the summed total displacements for these atoms.
        """
        if not self.atoms:
            raise ValueError("Cannot compute collective displacement for empty Species.")
        return float(sum(np.square(sum([atom.dr for atom in self.atoms]))))

    def occupations(self, site_label: str) -> int:
        """
        Number of these atoms occupying a specific site type.

        Args:
            site_label (Str): Label for the site type being considered.

        Returns:
            (Int): Number of atoms occupying sites of type `site_label`.
        """
        return sum(atom.site.label == site_label for atom in self.atoms)

    def summed_dr2(self) -> float:
        """
        Sum of squared individual displacements for these atoms.

        Args:
            None

        Returns:
            (Float): The sum of squared individual displacements for these atoms.
        """
        return float(sum([atom.summed_dr2 for atom in self.atoms]))

    def tracer_correlation(self) -> float:
        """
        Tracer correlation factor, f, for these atoms.

        Args:
            None

        Returns:
            (Float): The tracer correlation factor, f, for these atoms.
        """
        return self.sum_dr_squared() / self.summed_dr2()

    def collective_correlation(self) -> float:
        """
        Collective correlation factor, f_I, for these atoms.

        Args:
            None

        Returns:
            (Float): The collective correlation factor, f_I, for these atoms.
        """
        return self.collective_dr_squared() / self.summed_dr2()
