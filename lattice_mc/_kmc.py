"""One rejection-free kinetic Monte Carlo step on index-based state."""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from lattice_mc.error import BlockedLatticeError

if TYPE_CHECKING:
    from lattice_mc._state import SimulationState
    from lattice_mc.lattice import Lattice
    from lattice_mc.simulation import SimulationParameters


def possible_jumps(lattice: Lattice, occupied: npt.NDArray[np.bool_]) -> list[tuple[int, int]]:
    """
    All jumps not blocked by volume exclusion.

    These are jumps from occupied sites to neighbouring vacant sites. The order of the
    jumps depends only on the occupation. Neighbour lists are assumed to be symmetric.

    Args:
        lattice: The lattice, which defines the neighbours of each site.
        occupied: Whether each site is occupied.

    Returns:
        The possible jumps, as (initial site index, final site index) pairs.
    """
    if np.count_nonzero(occupied) <= len(occupied) / 2:
        return [
            (initial, final)
            for initial in np.flatnonzero(occupied).tolist()
            for final in lattice.neighbour_indices[initial]
            if not occupied[final]
        ]
    # Mostly occupied, so it is quicker to start from the vacant sites.
    return [
        (initial, final)
        for final in np.flatnonzero(~occupied).tolist()
        for initial in lattice.neighbour_indices[final]
        if occupied[initial]
    ]


def jump_vector(lattice: Lattice, initial: int, final: int) -> npt.NDArray[np.float64]:
    """
    The displacement of an atom jumping between two sites.

    Args:
        lattice: The lattice, which defines the site positions and the cell lengths.
        initial: Index of the site the atom jumps from.
        final: Index of the site the atom jumps to.

    Returns:
        The displacement vector, using the minimum-image convention. Components along axes
        where the cell length is zero are zero.
    """
    cell_lengths = lattice.cell_lengths
    half_cell_lengths = cell_lengths / 2.0
    dr = lattice.sites[final].r - lattice.sites[initial].r
    dr = np.where(dr > half_cell_lengths, dr - cell_lengths, dr)
    dr = np.where(dr < -half_cell_lengths, dr + cell_lengths, dr)
    return np.where(cell_lengths == 0.0, 0.0, dr)


def step(
    lattice: Lattice,
    state: SimulationState,
    probability: Callable[[npt.NDArray[np.bool_], int, int], float],
    params: SimulationParameters,
    rng: np.random.Generator,
) -> None:
    """
    Make one kinetic Monte Carlo jump and advance the clock.

    A jump is chosen with probability proportional to its relative probability. The time
    step is -ln(u) / (rate_prefactor * sum of relative probabilities), with u drawn
    uniformly from (0, 1]. The state is left unchanged if an exception is raised.

    Args:
        lattice: The lattice, which defines the neighbours, positions and cell lengths.
        state: The state to update.
        probability: Returns the relative probability of a jump, given the occupation of
            every site, the initial site index and the final site index.
        params: The simulation parameters.
        rng: The random number generator, which is the only source of randomness.

    Raises:
        BlockedLatticeError: If no jump is possible, or no possible jump has a non-zero
            probability.
        ValueError: If any relative probability is negative or not finite.
    """
    occupied = state.occupied
    jumps = possible_jumps(lattice, occupied)
    if not jumps:
        raise BlockedLatticeError("No jumps are possible.")
    probabilities = np.array([probability(occupied, initial, final) for initial, final in jumps])
    if not np.isfinite(probabilities).all():
        raise ValueError("Jump probabilities must be finite.")
    if (probabilities < 0.0).any():
        raise ValueError("Jump probabilities must not be negative.")
    cumulative = np.cumsum(probabilities)
    total = float(cumulative[-1])
    if total <= 0.0:
        raise BlockedLatticeError("No jump has a non-zero probability.")
    chosen = int(np.searchsorted(cumulative, rng.random() * total, side="right"))
    delta_t = -math.log(1.0 - rng.random()) / (params.rate_prefactor * total)
    initial, final = jumps[chosen]
    state.apply_jump(initial, final, jump_vector(lattice, initial, final), delta_t)
