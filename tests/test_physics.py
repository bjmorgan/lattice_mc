"""Seeded checks of simulation results against known physical values."""

import pytest

from lattice_mc import Results, Simulation, SimulationParameters, cubic_lattice, honeycomb_lattice, square_lattice

PARAMS = SimulationParameters(temperature=298.0, rate_prefactor=1.0)


def segment_results(simulation: Simulation, n_segments: int, n_jumps: int) -> list[Results]:
    """
    The results of consecutive segments of one run, with a reset after each segment.

    The segments are statistically independent samples, which is quicker than one
    simulation for each seed.
    """
    results = []
    for _ in range(n_segments):
        simulation.run(n_jumps=n_jumps)
        segment = simulation.results
        if segment is None:
            pytest.fail("Simulation.results is None after jumps were made")
        results.append(segment)
        simulation.reset()
    return results


@pytest.mark.parametrize(
    "build_lattice", [lambda: square_lattice(10, 10, 1.0), lambda: cubic_lattice(6, 6, 6, 1.0)], ids=["square", "cubic"]
)
def test_single_atom_tracer_diffusion_coefficient_is_rate_times_spacing_squared(build_lattice):
    # D* = nu a^2 for one atom, with rate_prefactor nu = 1 and spacing a = 1. Weighting each
    # segment by its time gives the summed squared displacement over the total time, without
    # the bias of averaging 1 / t over short segments.
    simulation = Simulation(build_lattice(), PARAMS, n_atoms=1, seed=0)
    results = segment_results(simulation, n_segments=1500, n_jumps=10)
    d = sum(r.tracer_diffusion_coefficient * r.time for r in results) / sum(r.time for r in results)
    assert d == pytest.approx(1.0, rel=0.1)


def test_honeycomb_single_vacancy_tracer_correlation_is_one_third():
    # Friauf (1962): f = 1/3 for vacancy-mediated diffusion on the honeycomb lattice. This
    # finite lattice and run length give about 0.35.
    simulation = Simulation(honeycomb_lattice(4, 3, 1.0), PARAMS, n_atoms=47, seed=0)
    results = segment_results(simulation, n_segments=40, n_jumps=2000)
    f = sum(r.tracer_correlation for r in results) / len(results)
    assert f == pytest.approx(1 / 3, abs=0.05)


def test_non_interacting_atoms_have_unit_collective_correlation():
    # The tracer correlation of these atoms is about 0.76, so this also checks that the
    # displacements are summed before squaring.
    simulation = Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=12, seed=0)
    results = segment_results(simulation, n_segments=1000, n_jumps=10)
    f = sum(r.collective_correlation for r in results) / len(results)
    assert f == pytest.approx(1.0, rel=0.1)


def test_lookup_table_and_on_the_fly_probabilities_give_the_same_results():
    lattice = honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
    lattice.set_site_energies(by_label={"A": 0.05, "B": -0.03})
    simulations = [
        Simulation(lattice, PARAMS, n_atoms=8, nn_energy=0.04, lookup_table=lookup_table, seed=0)
        for lookup_table in (False, True)
    ]
    for simulation in simulations:
        simulation.run(n_jumps=500)
    assert simulations[0].results == simulations[1].results
