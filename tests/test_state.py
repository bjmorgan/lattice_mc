"""Tests for lattice_mc._state."""

import numpy as np
import pytest

from lattice_mc._state import SimulationState


@pytest.fixture
def state() -> SimulationState:
    """A four-site state with atoms 0 and 1 on sites 2 and 0."""
    return SimulationState.from_occupied(4, [2, 0])


def test_from_occupied_assigns_atom_indices_in_order(state: SimulationState) -> None:
    assert state.atom_site.tolist() == [2, 0]
    assert state.site_atom.tolist() == [1, -1, 0, -1]


def test_from_occupied_zeroes_counters(state: SimulationState) -> None:
    assert state.displacement.shape == (2, 3)
    assert not state.displacement.any()
    assert state.summed_jump_length_squared == 0.0
    assert state.time == 0.0
    assert state.number_of_jumps == 0
    assert state.site_occupation_time.shape == (4,)
    assert not state.site_occupation_time.any()


def test_from_occupied_with_no_atoms() -> None:
    state = SimulationState.from_occupied(3, [])
    assert state.atom_site.shape == (0,)
    assert state.displacement.shape == (0, 3)
    assert state.site_atom.tolist() == [-1, -1, -1]


def test_occupied_property(state: SimulationState) -> None:
    assert state.occupied.tolist() == [True, False, True, False]


def test_from_occupied_raises_on_duplicate_sites() -> None:
    with pytest.raises(ValueError, match="duplicate"):
        SimulationState.from_occupied(4, [1, 1])


@pytest.mark.parametrize("bad_site", [-1, 4])
def test_from_occupied_raises_on_out_of_range_site(bad_site: int) -> None:
    with pytest.raises(ValueError, match="out of range"):
        SimulationState.from_occupied(4, [0, bad_site])


def test_apply_jump_moves_the_atom(state: SimulationState) -> None:
    state.apply_jump(2, 3, np.array([1.0, 0.0, 0.0]), 0.5)
    assert state.atom_site.tolist() == [3, 0]
    assert state.site_atom.tolist() == [1, -1, -1, 0]
    assert state.occupied.tolist() == [True, False, False, True]


def test_apply_jump_updates_counters_for_two_jumps(state: SimulationState) -> None:
    state.apply_jump(2, 3, np.array([1.0, 0.0, 0.0]), 0.5)
    state.apply_jump(3, 1, np.array([0.0, 2.0, 1.0]), 0.25)
    np.testing.assert_array_equal(state.displacement, [[1.0, 2.0, 1.0], [0.0, 0.0, 0.0]])
    assert state.summed_jump_length_squared == 1.0 + 5.0
    assert state.time == 0.75
    assert state.number_of_jumps == 2


def test_apply_jump_accrues_occupation_time_to_sites_occupied_before_the_move(
    state: SimulationState,
) -> None:
    state.apply_jump(2, 3, np.array([1.0, 0.0, 0.0]), 0.5)
    np.testing.assert_array_equal(state.site_occupation_time, [0.5, 0.0, 0.5, 0.0])
    state.apply_jump(3, 1, np.array([0.0, 2.0, 1.0]), 0.25)
    np.testing.assert_array_equal(state.site_occupation_time, [0.75, 0.0, 0.5, 0.25])


def test_apply_jump_does_not_alias_the_displacement_vector(state: SimulationState) -> None:
    dr = np.array([1.0, 0.0, 0.0])
    state.apply_jump(2, 3, dr, 0.5)
    dr[0] = 99.0
    assert state.displacement[0, 0] == 1.0


def test_reset_zeroes_counters_and_keeps_positions(state: SimulationState) -> None:
    state.apply_jump(2, 3, np.array([1.0, 0.0, 0.0]), 0.5)
    state.reset()
    assert not state.displacement.any()
    assert state.summed_jump_length_squared == 0.0
    assert state.time == 0.0
    assert state.number_of_jumps == 0
    assert not state.site_occupation_time.any()
    assert state.atom_site.tolist() == [3, 0]
    assert state.site_atom.tolist() == [1, -1, -1, 0]
