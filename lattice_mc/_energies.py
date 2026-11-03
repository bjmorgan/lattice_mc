"""Energy changes for proposed jumps, computed from site indices and an occupation array."""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Mapping
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:
    from lattice_mc.lattice import Lattice

CNEnergies = Mapping[str, Mapping[str, Mapping[int, float]]]
"""Coordination-number energies: site label -> neighbour label -> occupied count -> energy in eV."""


def metropolis(delta_E: float, kT: float) -> float:
    """
    Boltzmann probability factor for an event with an energy change `delta_E`, following the Metropolis algorithm.

    Args:
        delta_E: The change in energy, in eV.
        kT: Thermal energy kT in eV.

    Returns:
        Metropolis relative probability for this event.
    """
    if delta_E <= 0.0:
        return 1.0
    return math.exp(-delta_E / kT)


class EnergyModel:
    """
    The energy change of a jump, from site energies and nearest-neighbour and coordination-number interactions.

    The model reads the site energies, labels and neighbour lists of the lattice when it is
    created, so later changes to the lattice are not reflected.
    """

    def __init__(
        self,
        lattice: Lattice,
        *,
        nn_energy: float | None = None,
        cn_energies: CNEnergies | None = None,
    ) -> None:
        """
        Initialise an EnergyModel.

        Args:
            lattice: The lattice, which defines the site energies, labels and neighbours.
            nn_energy: Interaction energy, in eV, for each pair of occupied nearest-neighbour
                sites. Defaults to None, for no nearest-neighbour interaction.
            cn_energies: Coordination-number dependent energies of occupied sites. Defaults to
                None, for no coordination-number dependence. The energy of an occupied site
                is the sum, over the labels of its neighbours, of
                ``cn_energies[site_label][neighbour_label][number of occupied neighbours with that label]``.
                Any mapping given is validated and copied, so an empty mapping on a lattice
                with sites raises `ValueError`; pass None to turn the term off.

        Raises:
            ValueError: If `cn_energies` has no entry for a site label in the lattice, for a
                label of a neighbour of a site with that label, or for every count of occupied
                neighbours of that label from 0 up to the number of such neighbours.
        """
        self._site_energies: npt.NDArray[np.float64] = lattice.site_energies()
        self._labels: tuple[str, ...] = tuple(site.label for site in lattice.sites)
        self._neighbours: tuple[tuple[int, ...], ...] = lattice.neighbour_indices
        self._nn_energy: float | None = nn_energy
        self._cn_energies: dict[str, dict[str, dict[int, float]]] | None = None
        if cn_energies is not None:
            copied = {
                label: {neighbour_label: dict(by_count) for neighbour_label, by_count in by_neighbour.items()}
                for label, by_neighbour in cn_energies.items()
            }
            self._check_cn_energies(copied)
            self._cn_energies = copied

    @property
    def nn_energy(self) -> float | None:
        """Interaction energy, in eV, for each pair of occupied nearest-neighbour sites, or None for no interaction."""
        return self._nn_energy

    @property
    def has_cn_energies(self) -> bool:
        """Whether the model includes coordination-number dependent energies."""
        return self._cn_energies is not None

    def _check_cn_energies(self, cn_energies: CNEnergies) -> None:
        """
        Check that `cn_energies` covers every site label, neighbour label and occupied count.

        Args:
            cn_energies: The coordination-number energies to check.

        Raises:
            ValueError: If an entry is missing.
        """
        max_counts: dict[str, dict[str, int]] = {label: {} for label in self._labels}
        for index, label in enumerate(self._labels):
            neighbour_counts = Counter(self._labels[j] for j in self._neighbours[index])
            for neighbour_label, count in neighbour_counts.items():
                current = max_counts[label].get(neighbour_label, 0)
                max_counts[label][neighbour_label] = max(current, count)
        for label, by_neighbour in max_counts.items():
            if label not in cn_energies:
                raise ValueError(f"cn_energies has no entry for site label {label!r}")
            for neighbour_label, max_count in by_neighbour.items():
                if neighbour_label not in cn_energies[label]:
                    raise ValueError(
                        f"cn_energies for site label {label!r} has no entry for neighbour label {neighbour_label!r}"
                    )
                for count in range(max_count + 1):
                    if count not in cn_energies[label][neighbour_label]:
                        raise ValueError(
                            f"cn_energies for site label {label!r} and neighbour label {neighbour_label!r} "
                            f"has no entry for count {count}"
                        )

    def delta_E(self, occupied: npt.NDArray[np.bool_], initial: int, final: int) -> float:
        """
        The change in system energy if an atom jumped between two sites.

        The occupation array is not modified. If a site appears more than once in a neighbour
        list (only possible in very small built-in lattices), the nearest-neighbour term counts
        that neighbour once for each appearance.

        Args:
            occupied: Whether each site is occupied, before the jump.
            initial: Index of the occupied site the atom jumps from.
            final: Index of the vacant site the atom jumps to. It should be a neighbour of
                `initial`.

        Returns:
            The change in energy, in eV.
        """
        delta_E = float(self._site_energies[final] - self._site_energies[initial])
        if self._nn_energy is not None:
            n_initial = sum(1 for j in self._neighbours[initial] if occupied[j])
            n_final = sum(1 for j in self._neighbours[final] if occupied[j])
            # -1 because the jumping atom is not counted in the final site occupation number
            delta_E += self._nn_energy * (n_final - n_initial - 1)
        if self._cn_energies is not None:
            delta_E += self._coordination_number_delta_E(self._cn_energies, occupied, initial, final)
        return delta_E

    def _coordination_number_delta_E(
        self,
        cn_energies: CNEnergies,
        occupied: npt.NDArray[np.bool_],
        initial: int,
        final: int,
    ) -> float:
        """
        Coordination-number contribution to the change in energy for a jump.

        Only the two ends of the jump and their neighbours change energy. Each is visited once,
        so a site next to both ends gets both occupation changes at once.

        Args:
            cn_energies: The coordination-number energies.
            occupied: Whether each site is occupied, before the jump.
            initial: Index of the site the atom jumps from.
            final: Index of the site the atom jumps to.

        Returns:
            The coordination-number contribution to the change in energy, in eV.
        """
        affected = {initial, final, *self._neighbours[initial], *self._neighbours[final]}
        delta_E = 0.0
        for site in affected:
            if occupied[site]:
                delta_E -= self._cn_energy(site, cn_energies, occupied, None, None)
            if site == final or (occupied[site] and site != initial):
                delta_E += self._cn_energy(site, cn_energies, occupied, initial, final)
        return delta_E

    def _cn_energy(
        self,
        site: int,
        cn_energies: CNEnergies,
        occupied: npt.NDArray[np.bool_],
        vacated: int | None,
        filled: int | None,
    ) -> float:
        """
        Coordination-number energy of an occupied site.

        Args:
            site: Index of the site.
            cn_energies: The coordination-number energies.
            occupied: Whether each site is occupied.
            vacated: Index of a site to treat as vacant, or None for none.
            filled: Index of a site to treat as occupied, or None for none.

        Returns:
            The energy, in eV.
        """
        counts: dict[str, int] = {}
        for neighbour in self._neighbours[site]:
            label = self._labels[neighbour]
            is_occupied = neighbour == filled or (neighbour != vacated and bool(occupied[neighbour]))
            counts[label] = counts.get(label, 0) + is_occupied
        energies = cn_energies[self._labels[site]]
        return float(sum(energies[label][count] for label, count in counts.items()))
