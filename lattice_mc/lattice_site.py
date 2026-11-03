"""Individual lattice sites."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt


class Site:
    """
    A single lattice site: its identifying number, position, neighbours, occupation energy and label.
    """

    def __init__(
        self, number: int, coordinates: npt.NDArray[np.float64], neighbours: list[int], energy: float, label: str
    ) -> None:
        """
        Initialise a lattice Site object.

        Args:
            number: An identifying number for this site.
            coordinates: The coordinates of this site.
            neighbours: A list of the id numbers of the neighbouring sites.
            energy: On-site occupation energy, in eV.
            label: Label for classifying this as a specific site type.

        Notes:
            There should be a 1:1 mapping between sites and site numbers.
        """
        self.number: int = number
        self.r: npt.NDArray[np.float64] = coordinates
        self.neighbours: list[int] = neighbours
        self.energy: float = energy
        self.label: str = label
