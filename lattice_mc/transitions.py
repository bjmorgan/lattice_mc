from __future__ import annotations

import math
import random
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:
    from lattice_mc.jump import Jump
    from lattice_mc.simulation import SimulationParameters


class Transitions:
    """
    Transitions class

    Contains methods that operate on sets of Jumps.
    """

    def __init__(self, jumps: list[Jump], *, params: SimulationParameters) -> None:
        """
        Initialise a Transitions object.

        Args:
            jumps (List(Jump)): List of jumps to be contained in this Transitions object.
            params (SimulationParameters): Simulation parameters (temperature, rate prefactor).

        Returns:
            None
        """
        self.jumps: list[Jump] = jumps
        self.params: SimulationParameters = params
        self.p: npt.NDArray[np.float64] = np.array([jump.relative_probability for jump in self.jumps])

    def cumulative_probabilities(self) -> npt.NDArray[np.float64]:
        """
        Cumulative sum of the relative probabilities for all possible jumps.

        Args:
            None

        Returns:
            (np.array): Cumulative sum of relative jump probabilities.
        """
        partition_function = np.sum(self.p)
        return np.cumsum(self.p) / partition_function

    def random(self) -> Jump:
        """
        Select a jump at random with appropriate relative probabilities.

        Args:
            None

        Returns:
            (Jump): The randomly selected Jump.
        """
        j = np.searchsorted(self.cumulative_probabilities(), random.random())
        return self.jumps[j]

    def time_to_jump(self) -> float:
        """
        The timestep until the next jump.

        Args:
            None

        Returns:
            (Float): The timestep until the next jump.
        """
        k_tot = self.params.rate_prefactor * np.sum(self.p)
        return -(1.0 / k_tot) * math.log(random.random())
