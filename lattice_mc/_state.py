"""Mutable state of a running simulation, held as numpy arrays."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt


@dataclass(eq=False)
class SimulationState:
    """
    The mutable state of a running simulation.

    Sites are indexed 0..n_sites-1 and atoms 0..n_atoms-1.

    Attributes:
        site_atom: Index of the atom on each site, or -1 if the site is vacant.
        atom_site: Index of the site occupied by each atom.
        displacement: Cumulative displacement of each atom, with shape (n_atoms, 3).
        summed_jump_length_squared: Sum of the squared lengths of all jumps made.
        time: Total elapsed simulation time.
        number_of_jumps: Total number of jumps made.
        site_occupation_time: Time for which each site has been occupied.
    """

    site_atom: npt.NDArray[np.int_]
    atom_site: npt.NDArray[np.int_]
    displacement: npt.NDArray[np.float64]
    summed_jump_length_squared: float
    time: float
    number_of_jumps: int
    site_occupation_time: npt.NDArray[np.float64]

    @classmethod
    def from_occupied(cls, n_sites: int, occupied: Sequence[int]) -> SimulationState:
        """
        Create a state with one atom on each of the given sites.

        Atom indices are assigned in the order of `occupied`.

        Args:
            n_sites: The number of sites in the lattice.
            occupied: Indices of the occupied sites.

        Returns:
            A new state with all counters zeroed.

        Raises:
            ValueError: If `occupied` contains duplicate or out-of-range site indices.
        """
        atom_site = np.array(occupied, dtype=np.int_)
        if np.any((atom_site < 0) | (atom_site >= n_sites)):
            raise ValueError(f"Occupied site indices out of range for {n_sites} sites: {list(occupied)}")
        if len(np.unique(atom_site)) != len(atom_site):
            raise ValueError(f"Occupied site indices contain duplicate entries: {list(occupied)}")
        site_atom = np.full(n_sites, -1, dtype=np.int_)
        site_atom[atom_site] = np.arange(len(atom_site))
        return cls(
            site_atom=site_atom,
            atom_site=atom_site,
            displacement=np.zeros((len(atom_site), 3)),
            summed_jump_length_squared=0.0,
            time=0.0,
            number_of_jumps=0,
            site_occupation_time=np.zeros(n_sites),
        )

    @property
    def occupied(self) -> npt.NDArray[np.bool_]:
        """Boolean mask of occupied sites."""
        return self.site_atom >= 0

    def apply_jump(
        self,
        initial: int,
        final: int,
        dr: npt.NDArray[np.float64],
        delta_t: float,
    ) -> None:
        """
        Advance the clock and move one atom.

        Every site occupied before the move accrues `delta_t` of occupation time.

        Args:
            initial: Index of the site the atom jumps from.
            final: Index of the site the atom jumps to.
            dr: Displacement vector of the jump.
            delta_t: Time taken by the jump.

        Raises:
            ValueError: If the `initial` site is vacant or the `final` site is occupied.
                The state is left unchanged.
        """
        if self.site_atom[initial] < 0:
            raise ValueError(f"Cannot jump from site {initial}: site {initial} is vacant")
        if self.site_atom[final] >= 0:
            raise ValueError(f"Cannot jump to site {final}: site {final} is occupied")
        self.time += delta_t
        self.site_occupation_time[self.occupied] += delta_t
        atom = self.site_atom[initial]
        self.site_atom[initial] = -1
        self.site_atom[final] = atom
        self.atom_site[atom] = final
        self.displacement[atom] += dr
        self.summed_jump_length_squared += float(dr @ dr)
        self.number_of_jumps += 1

    def reset(self) -> None:
        """Zero all accumulated counters, leaving atom positions unchanged."""
        self.displacement[:] = 0.0
        self.summed_jump_length_squared = 0.0
        self.time = 0.0
        self.number_of_jumps = 0
        self.site_occupation_time[:] = 0.0
