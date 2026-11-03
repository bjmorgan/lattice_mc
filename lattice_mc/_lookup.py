"""Precomputed relative jump probabilities for the nearest-neighbour model."""

from __future__ import annotations

from itertools import product
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from lattice_mc._energies import EnergyModel, metropolis

if TYPE_CHECKING:
    from lattice_mc.lattice import Lattice


class LookupTable:
    """
    Relative jump probabilities, looked up by site labels and numbers of occupied neighbours.

    The table is built when it is created, so later changes to the lattice or the energy model
    are not reflected. It holds the same probabilities as ``metropolis(energy_model.delta_E(...), kT)``,
    for models with site energies and an optional nearest-neighbour interaction.
    """

    def __init__(self, lattice: Lattice, energy_model: EnergyModel, kT: float) -> None:
        """
        Initialise a LookupTable.

        Args:
            lattice: The lattice. The table reads its site energies, labels and neighbours from
                `energy_model`, so that it always agrees with the model.
            energy_model: The energy model. It must not have coordination-number energies.
            kT: Thermal energy kT in eV.

        Raises:
            ValueError: If `kT` is not positive, if `energy_model` has coordination-number
                energies, which depend on more than the numbers of occupied neighbours of the
                two sites, or if sites with the same label have different energies.
        """
        if kT <= 0.0:
            raise ValueError(f"kT must be positive, not {kT}.")
        if energy_model.has_cn_energies:
            raise ValueError("A LookupTable cannot be used with coordination-number energies.")
        nn_energy = energy_model.nn_energy
        self._labels: tuple[str, ...] = energy_model.labels
        site_energies: dict[str, float] = {}
        for label, energy in zip(self._labels, energy_model.site_energies, strict=True):
            if site_energies.setdefault(label, float(energy)) != energy:
                raise ValueError(
                    f"Sites with label {label} have different energies, {site_energies[label]} and {energy}."
                )
        self._neighbours: tuple[npt.NDArray[np.intp], ...] = tuple(
            np.array(neighbours, dtype=np.intp) for neighbours in energy_model.neighbour_indices
        )
        max_coordination: dict[str, int] = {}
        for label, neighbours in zip(self._labels, self._neighbours, strict=True):
            max_coordination[label] = max(max_coordination.get(label, 0), len(neighbours))
        self._probabilities: dict[tuple[str, str], npt.NDArray[np.float64]] = {}
        for initial_label, final_label in product(max_coordination, repeat=2):
            c_initial, c_final = np.indices((max_coordination[initial_label] + 1, max_coordination[final_label] + 1))
            delta_E = np.full(c_initial.shape, site_energies[final_label] - site_energies[initial_label])
            if nn_energy is not None:
                # -1 because the jumping atom is not counted in the final site occupation number
                delta_E += nn_energy * (c_final - c_initial - 1)
            self._probabilities[initial_label, final_label] = np.vectorize(metropolis)(delta_E, kT)

    def relative_probability(self, occupied: npt.NDArray[np.bool_], initial: int, final: int) -> float:
        """
        The relative probability of an atom jumping between two sites.

        Args:
            occupied: Whether each site is occupied, before the jump.
            initial: Index of the occupied site the atom jumps from.
            final: Index of the vacant site the atom jumps to.

        Returns:
            The relative probability of this jump occurring.
        """
        n_initial = int(occupied[self._neighbours[initial]].sum())
        n_final = int(occupied[self._neighbours[final]].sum())
        return float(self._probabilities[self._labels[initial], self._labels[final]][n_initial, n_final])
