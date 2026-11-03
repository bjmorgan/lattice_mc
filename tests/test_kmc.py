"""Tests for lattice_mc._kmc."""

import numpy as np
import numpy.typing as npt
import pytest

from lattice_mc._kmc import jump_vector, possible_jumps, step
from lattice_mc._state import SimulationState
from lattice_mc.error import BlockedLatticeError
from lattice_mc.lattice import Lattice
from lattice_mc.lattice_site import Site
from lattice_mc.simulation import SimulationParameters

PARAMS = SimulationParameters(temperature=300.0, rate_prefactor=2.0)


def lattice_from_positions(
    positions: list[list[float]], neighbours: list[list[int]], cell_lengths: list[float]
) -> Lattice:
    """Build a lattice from site positions and zero-based neighbour indices."""
    sites = [
        Site(number + 1, np.array(position), [n + 1 for n in neighbours[number]], 0.0, "A")
        for number, position in enumerate(positions)
    ]
    return Lattice(sites, np.array(cell_lengths))


@pytest.fixture
def ring() -> Lattice:
    """A periodic one-dimensional ring of four sites at x = 0, 1, 2, 3."""
    return lattice_from_positions(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
        [[1, 3], [0, 2], [1, 3], [2, 0]],
        [4.0, 0.0, 0.0],
    )


def occupation(n_sites: int, occupied: list[int]) -> npt.NDArray[np.bool_]:
    mask = np.zeros(n_sites, dtype=np.bool_)
    mask[occupied] = True
    return mask


def test_possible_jumps_pairs_occupied_sites_with_vacant_neighbours(ring: Lattice) -> None:
    assert sorted(possible_jumps(ring, occupation(4, [0, 1]))) == [(0, 3), (1, 2)]


def test_possible_jumps_is_empty_when_every_neighbour_is_occupied(ring: Lattice) -> None:
    assert possible_jumps(ring, occupation(4, [0, 1, 2, 3])) == []


def test_possible_jumps_is_empty_when_no_site_is_occupied(ring: Lattice) -> None:
    assert possible_jumps(ring, occupation(4, [])) == []


def test_possible_jumps_on_a_mostly_occupied_lattice() -> None:
    lattice = lattice_from_positions(
        [[float(x), 0.0, 0.0] for x in range(6)],
        [[1, 5], [0, 2], [1, 3], [2, 4], [3, 5], [4, 0]],
        [6.0, 0.0, 0.0],
    )
    assert sorted(possible_jumps(lattice, occupation(6, [0, 1, 2, 3, 4]))) == [(0, 5), (4, 5)]
    assert sorted(possible_jumps(lattice, occupation(6, [0, 1, 3, 4, 5]))) == [(1, 2), (3, 2)]


def test_possible_jumps_gives_the_same_jumps_whichever_way_they_are_found() -> None:
    lattice = lattice_from_positions(
        [[float(x), 0.0, 0.0] for x in range(6)],
        [[1, 5], [0, 2], [1, 3], [2, 4], [3, 5], [4, 0]],
        [6.0, 0.0, 0.0],
    )
    for occupied in ([0, 1, 3], [0, 1, 2, 4], [1, 2, 3, 4, 5]):
        mask = occupation(6, occupied)
        expected = sorted((i, f) for i in occupied for f in lattice.neighbour_indices[i] if not mask[f])
        assert sorted(possible_jumps(lattice, mask)) == expected


def test_possible_jumps_order_depends_only_on_the_occupation(ring: Lattice) -> None:
    mask = occupation(4, [0, 1, 2])
    assert possible_jumps(ring, mask) == possible_jumps(ring, mask.copy())


def test_jump_vector_between_neighbours(ring: Lattice) -> None:
    assert jump_vector(ring, 0, 1).tolist() == [1.0, 0.0, 0.0]


def test_jump_vector_uses_the_minimum_image_across_the_boundary(ring: Lattice) -> None:
    assert jump_vector(ring, 3, 0).tolist() == [1.0, 0.0, 0.0]
    assert jump_vector(ring, 0, 3).tolist() == [-1.0, 0.0, 0.0]


def test_jump_vector_keeps_a_displacement_of_exactly_half_the_cell() -> None:
    lattice = lattice_from_positions([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], [[1], [0]], [4.0, 0.0, 0.0])
    assert jump_vector(lattice, 0, 1).tolist() == [2.0, 0.0, 0.0]
    assert jump_vector(lattice, 1, 0).tolist() == [-2.0, 0.0, 0.0]


def test_jump_vector_on_a_two_dimensional_lattice() -> None:
    lattice = lattice_from_positions(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        [[1, 2], [0], [0]],
        [3.0, 3.0, 0.0],
    )
    assert jump_vector(lattice, 0, 1).tolist() == [-1.0, 0.0, 0.0]
    assert jump_vector(lattice, 0, 2).tolist() == [0.0, 1.0, 0.0]


def test_jump_vector_is_zero_along_an_axis_with_zero_cell_length() -> None:
    lattice = lattice_from_positions(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.5]],
        [[1], [0]],
        [4.0, 0.0, 0.0],
    )
    assert jump_vector(lattice, 0, 1).tolist() == [1.0, 0.0, 0.0]


def test_step_moves_an_atom_and_advances_time(ring: Lattice) -> None:
    state = SimulationState.from_occupied(4, [0])
    step(ring, state, lambda occupied, initial, final: 1.0, PARAMS, np.random.default_rng(0))
    assert state.number_of_jumps == 1
    assert state.time > 0.0
    assert state.atom_site[0] in (1, 3)
    assert abs(state.displacement[0, 0]) == 1.0


def test_step_with_the_same_seed_is_reproducible(ring: Lattice) -> None:
    def run() -> SimulationState:
        state = SimulationState.from_occupied(4, [0, 1])
        rng = np.random.default_rng(42)
        for _ in range(20):
            step(ring, state, lambda occupied, initial, final: 1.0, PARAMS, rng)
        return state

    first, second = run(), run()
    assert first.time == second.time
    assert first.atom_site.tolist() == second.atom_site.tolist()
    assert first.displacement.tolist() == second.displacement.tolist()


def test_step_on_a_blocked_lattice_raises(ring: Lattice) -> None:
    state = SimulationState.from_occupied(4, [0, 1, 2, 3])
    with pytest.raises(BlockedLatticeError, match="No jumps are possible."):
        step(ring, state, lambda occupied, initial, final: 1.0, PARAMS, np.random.default_rng(0))


def test_step_raises_when_every_jump_has_zero_probability(ring: Lattice) -> None:
    state = SimulationState.from_occupied(4, [0])
    with pytest.raises(BlockedLatticeError, match="non-zero probability"):
        step(ring, state, lambda occupied, initial, final: 0.0, PARAMS, np.random.default_rng(0))
    assert state.number_of_jumps == 0
    assert state.time == 0.0


def test_step_never_chooses_a_jump_with_zero_probability(ring: Lattice) -> None:
    rng = np.random.default_rng(1)
    for _ in range(200):
        state = SimulationState.from_occupied(4, [0])
        step(ring, state, lambda occupied, initial, final: float(final == 3), PARAMS, rng)
        assert state.atom_site[0] == 3


def test_step_mean_time_step_is_the_inverse_of_the_total_rate(ring: Lattice) -> None:
    # A single atom on the ring always has one jump in each direction, with probabilities
    # 0.3 and 0.9, so the total rate is constant and differs from the number of jumps.
    def probability(occupied: npt.NDArray[np.bool_], initial: int, final: int) -> float:
        return 0.3 if final == (initial + 1) % 4 else 0.9

    state = SimulationState.from_occupied(4, [0])
    rng = np.random.default_rng(2024)
    n_steps = 20000
    times = np.empty(n_steps + 1)
    times[0] = state.time
    for i in range(n_steps):
        step(ring, state, probability, PARAMS, rng)
        times[i + 1] = state.time
    expected = 1.0 / (PARAMS.rate_prefactor * (0.3 + 0.9))
    # The standard error of the mean is expected / sqrt(n_steps), about 0.7 %.
    assert np.diff(times).mean() == pytest.approx(expected, rel=0.05)


def test_step_chooses_jumps_in_proportion_to_their_probabilities(ring: Lattice) -> None:
    rng = np.random.default_rng(7)
    n_steps = 5000
    n_to_site_1 = 0
    for _ in range(n_steps):
        state = SimulationState.from_occupied(4, [0])
        step(ring, state, lambda occupied, initial, final: 1.0 if final == 1 else 0.25, PARAMS, rng)
        n_to_site_1 += int(state.atom_site[0] == 1)
    # Expected fraction 0.8, with standard error sqrt(0.8 * 0.2 / n_steps), about 0.006.
    assert n_to_site_1 / n_steps == pytest.approx(0.8, abs=0.03)


def test_step_chooses_among_three_jumps_with_a_zero_weight_in_the_middle() -> None:
    star = lattice_from_positions(
        [[float(x), 0.0, 0.0] for x in range(4)],
        [[1, 2, 3], [0], [0], [0]],
        [10.0, 0.0, 0.0],
    )
    weights = {1: 1.0, 2: 0.0, 3: 3.0}
    rng = np.random.default_rng(11)
    n_steps = 4000
    counts = {1: 0, 2: 0, 3: 0}
    for _ in range(n_steps):
        state = SimulationState.from_occupied(4, [0])
        step(star, state, lambda occupied, initial, final: weights[final], PARAMS, rng)
        counts[int(state.atom_site[0])] += 1
    # Standard error at most sqrt(0.75 * 0.25 / n_steps), about 0.007.
    assert counts[2] == 0
    assert counts[1] / n_steps == pytest.approx(0.25, abs=0.03)
    assert counts[3] / n_steps == pytest.approx(0.75, abs=0.03)


@pytest.mark.parametrize(
    ("bad_value", "message"),
    [(-0.5, "negative"), (float("nan"), "finite"), (float("inf"), "finite"), (float("-inf"), "finite")],
)
def test_step_rejects_negative_or_non_finite_probabilities(ring: Lattice, bad_value: float, message: str) -> None:
    state = SimulationState.from_occupied(4, [0])
    with pytest.raises(ValueError, match=message):
        step(
            ring,
            state,
            lambda occupied, initial, final: 1.0 if final == 1 else bad_value,
            PARAMS,
            np.random.default_rng(0),
        )
    assert state.number_of_jumps == 0
    assert state.time == 0.0
    assert state.atom_site.tolist() == [0]
