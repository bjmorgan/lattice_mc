"""Tests for lattice_mc.results."""

import copy
import dataclasses
import pickle
from collections.abc import Callable

import numpy as np
import pytest

from lattice_mc import init_lattice
from lattice_mc._state import SimulationState
from lattice_mc.lattice import Lattice
from lattice_mc.lattice_site import Site
from lattice_mc.results import Results


def labelled_ring(labels: tuple[str, ...] = ("A", "B", "A", "B")) -> Lattice:
    """A periodic one-dimensional ring of four sites with the given labels."""
    neighbours = [[2, 4], [1, 3], [2, 4], [3, 1]]
    sites = [
        Site(number + 1, np.array([float(number), 0.0, 0.0]), neighbours[number], 0.0, labels[number])
        for number in range(4)
    ]
    return Lattice(sites, np.array([4.0, 0.0, 0.0]))


def moved_state(n_sites: int) -> SimulationState:
    """
    A state of two atoms on the first two sites, as if after four jumps.

    The sum of squared displacements is 6, the squared summed displacement is 8,
    the summed squared jump length is 3 and the time is 2.
    """
    state = SimulationState.from_occupied(n_sites, [0, 1])
    state.displacement[:] = [[1.0, 0.0, 0.0], [1.0, 2.0, 0.0]]
    state.summed_jump_length_squared = 3.0
    state.time = 2.0
    state.number_of_jumps = 4
    state.site_occupation_time[:4] = [2.0, 1.0, 0.5, 0.0]
    return state


def results_from(state: SimulationState, lattice: Lattice) -> Results:
    """Results for a state in which jumps have been made."""
    result = Results.from_state(state, lattice)
    if result is None:
        pytest.fail("Results.from_state returned None for a state with jumps")
    return result


@pytest.fixture
def results() -> Results:
    return results_from(moved_state(4), labelled_ring())


def test_from_state_returns_none_before_any_jump() -> None:
    assert Results.from_state(SimulationState.from_occupied(4, [0, 1]), labelled_ring()) is None


def test_tracer_diffusion_coefficient(results: Results) -> None:
    # 6 / (2 * 1 * 2 atoms * 2 s)
    assert results.tracer_diffusion_coefficient == pytest.approx(0.75)


def test_collective_diffusion_coefficient_is_per_atom(results: Results) -> None:
    # 8 / (2 * 1 * 2 atoms * 2 s)
    assert results.collective_diffusion_coefficient == pytest.approx(1.0)


def test_reduced_ionic_conductivity(results: Results) -> None:
    # D_J * 2 atoms / 4 sites
    assert results.reduced_ionic_conductivity == pytest.approx(0.5)


def test_tracer_correlation(results: Results) -> None:
    assert results.tracer_correlation == pytest.approx(6.0 / 3.0)


def test_collective_correlation(results: Results) -> None:
    assert results.collective_correlation == pytest.approx(8.0 / 3.0)


def test_average_site_occupations_sum_over_sites_with_each_label(results: Results) -> None:
    assert results.average_site_occupations == pytest.approx({"A": 2.5 / 2.0, "B": 1.0 / 2.0})


def test_time_and_number_of_jumps_are_copied_from_the_state(results: Results) -> None:
    assert results.time == 2.0
    assert results.number_of_jumps == 4


@pytest.mark.parametrize(
    ("lattice", "tracer", "collective", "conductivity"),
    [
        (init_lattice.square_lattice(3, 3, 1.0), 0.375, 0.5, 1.0 / 9.0),
        (init_lattice.cubic_lattice(3, 3, 3, 1.0), 0.25, 1.0 / 3.0, 2.0 / 81.0),
    ],
    ids=["square", "cubic"],
)
def test_diffusion_coefficients_depend_on_the_dimensionality(
    lattice: Lattice, tracer: float, collective: float, conductivity: float
) -> None:
    result = results_from(moved_state(lattice.number_of_sites), lattice)
    assert result.tracer_diffusion_coefficient == pytest.approx(tracer)
    assert result.collective_diffusion_coefficient == pytest.approx(collective)
    assert result.reduced_ionic_conductivity == pytest.approx(conductivity)


def test_average_site_occupations_with_unequal_numbers_of_sites_per_label() -> None:
    lattice = labelled_ring(("A", "A", "A", "B"))
    result = results_from(moved_state(4), lattice)
    # A: (2.0 + 1.0 + 0.5) / 2 s, B: 0.0 / 2 s
    assert result.average_site_occupations == pytest.approx({"A": 1.75, "B": 0.0})


def test_average_site_occupations_are_read_only(results: Results) -> None:
    with pytest.raises(TypeError):
        results.average_site_occupations["A"] = 0.0  # type: ignore[index]


def test_results_survive_pickling(results: Results) -> None:
    assert pickle.loads(pickle.dumps(results)) == results


def test_results_survive_deepcopy(results: Results) -> None:
    assert copy.deepcopy(results) == results


@pytest.mark.parametrize("copier", [lambda r: pickle.loads(pickle.dumps(r)), copy.deepcopy], ids=["pickle", "deepcopy"])
def test_copied_average_site_occupations_are_read_only(results: Results, copier: Callable[[Results], Results]) -> None:
    copied = copier(results)
    with pytest.raises(TypeError):
        copied.average_site_occupations["A"] = 0.0


def test_equality_compares_the_average_site_occupations(results: Results) -> None:
    same = results_from(moved_state(4), labelled_ring())
    different = results_from(moved_state(4), labelled_ring(("A", "A", "A", "B")))
    assert same == results
    assert different != results


def test_results_are_frozen(results: Results) -> None:
    with pytest.raises(dataclasses.FrozenInstanceError):
        results.time = 3.0  # type: ignore[misc]
