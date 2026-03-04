from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:
    from lattice_mc.lattice_site import Site


class Atom:
    """
    Atoms are distinguishable particles, each occupying a specific lattice site.
    """

    atom_number: int = 0  # counter, so that every Atom instance is distinguishable

    def __init__(self, initial_site: Site) -> None:
        """
        Initialise an Atom instance.

        Args:
            initial_site (Site): Lattice site initially occupied by this Atom.

        Returns:
            None
        """
        Atom.atom_number += 1
        self.number: int = Atom.atom_number
        self._site = initial_site
        # check this site is not already occupied
        if self._site.occupation == 0:
            self._site.occupation = self.number
            self._site.is_occupied = True
            self._site.atom = self
        else:
            raise ValueError("This site is already occupied by atom {}".format(initial_site.occupation))
        self.reset()

    def reset(self) -> None:
        """
        Reinitialise the stored displacements, number of hops, and list of sites visited for this `Atom`.

        Args:
            None

        Returns:
            None
        """
        self.number_of_hops: int = 0
        self.dr: npt.NDArray[np.float64] = np.array([0.0, 0.0, 0.0])
        self.summed_dr2: float = 0.0
        self.sites_visited: list[int] = [self._site.number]

    def dr_squared(self) -> float:
        """
        :math:`|dr|^2`, where :math:`dr` is the total displacement vector for this `Atom`.

        Args:
            None

        Returns:
            dr_squared (float): :math:`|dr|^2`.
        """
        return float(np.dot(self.dr, self.dr))

    @property
    def site(self) -> Site:
        """
        Get or set `self.site` for this `Atom`.
        """
        return self._site

    @site.setter
    def site(self, value: Site) -> None:
        self._site = value
