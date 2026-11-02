"""Jumps of atoms between neighbouring lattice sites."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:
    from lattice_mc.lattice_site import Site
    from lattice_mc.lookup_table import LookupTable
    from lattice_mc.simulation import SimulationParameters

"""
A Jump describes a possible move by a particle from one site to another site.
"""

# TODO (maybe) does it make more sense to subclass Jump (and other classes) depending on the type of Hamiltonian?


class Jump:
    """
    A possible jump of an atom from an occupied site to a neighbouring vacant site.
    """

    def __init__(
        self,
        initial_site: Site,
        final_site: Site,
        nearest_neighbour_energy: float | None = None,
        coordination_number_energy: dict[str, dict[str, dict[int, float]]] | None = None,
        jump_lookup_table: LookupTable | None = None,
        *,
        params: SimulationParameters,
    ) -> None:
        """
        Initialise a Jump instance.

        Args:
            initial_site: Lattice site occupied before the jump.
            final_site: Lattice site occupied after the jump.
            nearest_neighbour_energy: Interaction energy, in eV, for each pair of occupied
                nearest-neighbour sites. Defaults to None.
            coordination_number_energy: Coordination-number dependent energies, in the form taken by
                `Lattice.set_cn_energies`. Defaults to None.
            jump_lookup_table: If the jump relative probabilities have been precalculated and stored in a lookup-table, this table should be passed in here. If not, jump probabilities are calculated on the fly. Defaults to None.
            params: Simulation parameters (temperature, rate prefactor).
        """
        self.initial_site: Site = initial_site
        self.final_site: Site = final_site
        self.nearest_neighbour_energy: float | None = nearest_neighbour_energy
        self.coordination_number_energy: dict[str, dict[str, dict[int, float]]] | None = coordination_number_energy
        self.params: SimulationParameters = params
        if jump_lookup_table:
            self._relative_probability: float = self.relative_probability_from_lookup_table(jump_lookup_table)
        else:
            self._relative_probability = self.boltzmann_factor()

    def rate(self) -> float:
        """
        Average rate for this jump. Calculated as v_0 * P_jump.

        Returns:
            The average rate for this jump.
        """
        return self.params.rate_prefactor * self.relative_probability

    def boltzmann_factor(self) -> float:
        """
        Boltzmann probability factor for accepting this jump, following the Metropolis algorithm.

        Returns:
            Metropolis relative probability of accepting this jump.
        """
        if self.delta_E() <= 0.0:
            return 1.0
        else:
            return math.exp(-self.delta_E() / self.params.kT)

    def delta_E(self) -> float:
        """
        The change in system energy if this jump were accepted.

        Returns:
            The change in energy, in eV.
        """
        site_delta_E = self.final_site.energy - self.initial_site.energy
        if self.nearest_neighbour_energy:
            site_delta_E += self.nearest_neighbour_delta_E()
        if self.coordination_number_energy:
            site_delta_E += self.coordination_number_delta_E()
        return site_delta_E

    def nearest_neighbour_delta_E(self) -> float:
        """
        Nearest-neighbour interaction contribution to the change in system energy if this jump were accepted.

        Returns:
            The nearest-neighbour contribution to the change in energy, in eV.
        """
        delta_nn = (
            self.final_site.nn_occupation() - self.initial_site.nn_occupation() - 1
        )  # -1 because the hopping ion is not counted in the final site occupation number
        assert self.nearest_neighbour_energy is not None
        return delta_nn * self.nearest_neighbour_energy

    def coordination_number_delta_E(self) -> float:
        """
        Coordination-number dependent energy contribution to the change in system energy if this jump were accepted.

        Returns:
            The coordination-number contribution to the change in energy, in eV.
        """
        assert self.initial_site.p_neighbours is not None
        assert self.final_site.p_neighbours is not None
        initial_site_neighbours = [
            s for s in self.initial_site.p_neighbours if s.is_occupied
        ]  # excludes final site, since this is always unoccupied
        final_site_neighbours = [
            s for s in self.final_site.p_neighbours if s.is_occupied and s is not self.initial_site
        ]  # excludes initial site
        initial_cn_occupation_energy = (
            self.initial_site.cn_occupation_energy()
            + sum([site.cn_occupation_energy() for site in initial_site_neighbours])
            + sum([site.cn_occupation_energy() for site in final_site_neighbours])
        )
        final_cn_occupation_energy = (
            self.final_site.cn_occupation_energy(delta_occupation={self.initial_site.label: -1})
            + sum([site.cn_occupation_energy(delta_occupation={self.initial_site.label: -1})
                   for site in initial_site_neighbours])
            + sum([site.cn_occupation_energy(delta_occupation={self.final_site.label: +1})
                   for site in final_site_neighbours])
        )
        return final_cn_occupation_energy - initial_cn_occupation_energy

    def dr(self, cell_lengths: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """
        Particle displacement vector for this jump.

        Args:
            cell_lengths: The [x, y, z] lengths of the orthorhombic simulation cell.

        Returns:
            The displacement vector of the jumping atom, using the minimum-image convention.
        """
        half_cell_lengths = cell_lengths / 2.0
        this_dr = self.final_site.r - self.initial_site.r
        for i in range(3):
            if this_dr[i] > half_cell_lengths[i]:
                this_dr[i] -= cell_lengths[i]
            if this_dr[i] < -half_cell_lengths[i]:
                this_dr[i] += cell_lengths[i]
        return this_dr

    def relative_probability_from_lookup_table(self, jump_lookup_table: LookupTable) -> float:
        """
        Relative probability of accepting this jump from a lookup-table.

        Args:
            jump_lookup_table: The lookup table to use for this jump.

        Returns:
            The relative probability of accepting this jump.
        """
        l1 = self.initial_site.label
        l2 = self.final_site.label
        c1 = self.initial_site.nn_occupation()
        c2 = self.final_site.nn_occupation()
        return jump_lookup_table.jump_probability[l1][l2][c1][c2]

    @property
    def relative_probability(self) -> float:
        """
        Relative probability of this jump.
        """
        return self._relative_probability

    @relative_probability.setter
    def relative_probability(self, value: float) -> None:
        self._relative_probability = value
