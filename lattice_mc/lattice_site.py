from __future__ import annotations

from collections import Counter
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:
    from lattice_mc.atom import Atom


class Site:
    """
    Site class
    """

    index: int = 0

    def __init__(
        self,
        number: int,
        coordinates: npt.NDArray[np.float64],
        neighbours: list[int],
        energy: float,
        label: str,
        cn_energies: dict[str, dict[int, float]] | None = None,
    ) -> None:
        """
        Initialise a lattce Site object.

        Args:
            number (Int): An identifying number for this site.
            coordinates (np.array(x,y,z)): The coordinates of this site.
            neighbours (List(Int)): A list of the id numbers of the neighbouring sites.
            energy (Float): On-site occupation energy.
            label (Str): Label for classifying this as a specific site type.
            cn_energies (:obj:Dict(Int:Float), optional): Dictionary of coordination-number dependent energies, e.g. { 0 : 0.0, 1 : 0.5, 2 : 2.0 }. Defaults to None.

        Returns:
            None

        Notes:
            There should be a 1:1 mapping between sites and site numbers.
        """
        self.number: int = number
        self.index: int = Site.index
        Site.index += 1
        self.r: npt.NDArray[np.float64] = coordinates
        self.neighbours: list[int] = neighbours
        self.p_neighbours: list[Site] | None = None  # pointer to neighbouring sites. initialised in Lattice.__init__
        self.energy: float = energy
        self.occupation: int = 0
        self.atom: Atom | None = None
        self.is_occupied: bool = False
        self.label: str = label
        self.time_occupied: float = 0.0
        self.cn_occupation_energies: dict[str, dict[int, float]] | None = cn_energies

    def nn_occupation(self) -> int:
        """
        The number of occupied nearest-neighbour sites.

        Args:
            None

        Returns:
            (Int): The number of occupied nearest-neighbour sites.
        """
        assert self.p_neighbours is not None
        return sum([site.is_occupied for site in self.p_neighbours])

    def site_specific_nn_occupation(self) -> dict[str, int]:
        """
        Returns the number of occupied nearest neighbour sites, classified by site type.

        Args:
            None

        Returns:
            (Dict(Str:Int)): Dictionary of nearest-neighbour occupied site numbers, classified by site label, e.g. { 'A' : 2, 'B' : 1 }.
        """
        assert self.p_neighbours is not None
        to_return = {label: 0 for label in set((site.label for site in self.p_neighbours))}
        for site in self.p_neighbours:
            if site.is_occupied:
                to_return[site.label] += 1
        return to_return

    def site_specific_neighbours(self) -> dict[str, int]:
        """
        Returns the number of neighbouring sites, classified by site type.

        Args:
            None

        Returns:
            (Dict(Str:Int)): Dictionary of neighboring sites, classified by site label, e.g. { 'A' : 1, 'B' : 1 }.
        """
        assert self.p_neighbours is not None
        return dict(Counter((site.label for site in self.p_neighbours)))

    def set_cn_occupation_energies(self, cn_energies: dict[str, dict[int, float]]) -> None:
        """
        Set the coordination-number dependent energies for this site.

        Args:
            cn_energies (Dict(Int:Float)): Dictionary of coordination number dependent site energies, e.g. { 0 : 0.0, 1 : 0.5 }.

        Returns:
            None
        """
        self.cn_occupation_energies = cn_energies

    def cn_occupation_energy(self, delta_occupation: dict[str, int] | None = None) -> float:
        """
        The coordination-number dependent energy for this site.

        Args:
            delta_occupation (:obj:Dict(Str:Int), optional): A dictionary of a change in (site-type specific) coordination number, e.g. { 'A' : 1, 'B' : -1 }.
                If this is not None, the coordination-number dependent energy is calculated including these changes in neighbour-site occupations. Defaults to None

        Returns:
            (Float): The coordination-number dependent energy for this site.
        """
        nn_occupations = self.site_specific_nn_occupation()
        if delta_occupation:
            for site in delta_occupation:
                assert site in nn_occupations
                nn_occupations[site] += delta_occupation[site]
        assert self.cn_occupation_energies is not None
        return float(sum([self.cn_occupation_energies[s][n] for s, n in nn_occupations.items()]))
