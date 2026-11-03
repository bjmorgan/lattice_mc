"""Tests for lattice_mc.simulation.Simulation."""

import pathlib

import pytest

from lattice_mc import (
    BlockedLatticeError,
    Results,
    Simulation,
    SimulationParameters,
    cubic_lattice,
    honeycomb_lattice,
    lattice_from_sites_file,
    square_lattice,
)

PARAMS = SimulationParameters(temperature=298.0, rate_prefactor=1.0)
LLZO = pathlib.Path(__file__).resolve().parents[1] / "examples" / "llzo_lattice_site_list.dat"


def results_of(simulation: Simulation) -> Results:
    """The results of a simulation that has made jumps."""
    results = simulation.results
    if results is None:
        pytest.fail("Simulation.results is None after jumps were made")
    return results


class TestConstruction:
    def test_needs_n_atoms_or_occupied_sites(self):
        with pytest.raises(ValueError):
            Simulation(square_lattice(4, 4, 1.0), PARAMS)

    def test_rejects_both_n_atoms_and_occupied_sites(self):
        with pytest.raises(ValueError):
            Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=1, occupied_sites=[1])

    def test_rejects_zero_atoms(self):
        with pytest.raises(ValueError):
            Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=0)

    def test_rejects_more_atoms_than_sites(self):
        with pytest.raises(ValueError):
            Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=17)

    def test_places_n_atoms(self):
        simulation = Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=5, seed=0)
        assert len(simulation.occupied_sites) == 5

    def test_selected_labels_restrict_placement(self):
        lattice = honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
        simulation = Simulation(lattice, PARAMS, n_atoms=8, selected_labels=["A"], seed=0)
        a_sites = {site.number for site in lattice.sites if site.label == "A"}
        assert set(simulation.occupied_sites) == a_sites

    def test_rejects_more_atoms_than_selected_sites(self):
        lattice = honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
        with pytest.raises(ValueError):
            Simulation(lattice, PARAMS, n_atoms=9, selected_labels=["A"])

    def test_rejects_unknown_selected_labels(self):
        lattice = honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
        with pytest.raises(ValueError):
            Simulation(lattice, PARAMS, n_atoms=1, selected_labels=["A", "C"])

    def test_rejects_empty_selected_labels(self):
        lattice = honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
        with pytest.raises(ValueError, match="selected_labels"):
            Simulation(lattice, PARAMS, n_atoms=1, selected_labels=[])

    def test_rejects_selected_labels_with_occupied_sites(self):
        lattice = honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
        with pytest.raises(ValueError):
            Simulation(lattice, PARAMS, occupied_sites=[1], selected_labels=["A"])

    def test_occupied_sites_places_atoms_on_those_sites(self):
        simulation = Simulation(square_lattice(4, 4, 1.0), PARAMS, occupied_sites=[7, 2])
        assert sorted(simulation.occupied_sites) == [2, 7]

    def test_rejects_unknown_occupied_sites(self):
        with pytest.raises(ValueError):
            Simulation(square_lattice(4, 4, 1.0), PARAMS, occupied_sites=[1, 99])

    def test_rejects_empty_occupied_sites(self):
        with pytest.raises(ValueError):
            Simulation(square_lattice(4, 4, 1.0), PARAMS, occupied_sites=[])

    def test_rejects_incomplete_cn_energies(self):
        cn_energies = {"L": {"L": {count: 0.0 for count in range(4)}}}
        with pytest.raises(ValueError):
            Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=1, cn_energies=cn_energies)

    def test_rejects_a_non_finite_nn_energy(self):
        with pytest.raises(ValueError, match="nn_energy must be finite"):
            Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=1, nn_energy=float("nan"))

    def test_rejects_non_finite_cn_energies(self):
        cn_energies = {"L": {"L": {count: float("inf") for count in range(5)}}}
        with pytest.raises(ValueError, match="must be finite"):
            Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=1, cn_energies=cn_energies)

    def test_lookup_table_rejects_cn_energies(self):
        cn_energies = {"L": {"L": {count: 0.0 for count in range(5)}}}
        with pytest.raises(ValueError):
            Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=1, cn_energies=cn_energies, lookup_table=True)

    def test_lookup_table_rejects_different_energies_for_one_label(self):
        lattice = square_lattice(4, 4, 1.0)
        lattice.set_site_energies(by_site={1: 0.1})
        with pytest.raises(ValueError):
            Simulation(lattice, PARAMS, n_atoms=1, lookup_table=True)


class TestRun:
    @pytest.fixture
    def simulation(self) -> Simulation:
        return Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=4, seed=0)

    def test_zero_jumps_is_a_no_op(self, simulation):
        simulation.run(n_jumps=0)
        assert simulation.results is None

    def test_zero_time_is_a_no_op(self, simulation):
        simulation.run(for_time=0.0)
        assert simulation.results is None

    def test_run_rejects_negative_n_jumps(self, simulation):
        with pytest.raises(ValueError, match="n_jumps"):
            simulation.run(n_jumps=-1)

    @pytest.mark.parametrize("for_time", [-1.0, float("nan"), float("inf")])
    def test_run_rejects_negative_or_non_finite_for_time(self, simulation, for_time):
        with pytest.raises(ValueError, match="for_time"):
            simulation.run(for_time=for_time)

    def test_results_are_none_before_any_run(self, simulation):
        assert simulation.results is None

    def test_run_n_jumps_makes_that_many_jumps(self, simulation):
        simulation.run(n_jumps=10)
        assert results_of(simulation).number_of_jumps == 10

    def test_run_needs_an_argument(self, simulation):
        with pytest.raises(ValueError):
            simulation.run()

    def test_run_rejects_both_arguments(self, simulation):
        with pytest.raises(ValueError):
            simulation.run(n_jumps=10, for_time=1.0)

    def test_run_for_time_reaches_that_time(self, simulation):
        simulation.run(for_time=5.0)
        assert results_of(simulation).time >= 5.0

    def test_second_run_for_the_same_time_makes_no_jumps(self, simulation):
        simulation.run(for_time=5.0)
        number_of_jumps = results_of(simulation).number_of_jumps
        simulation.run(for_time=5.0)
        assert results_of(simulation).number_of_jumps == number_of_jumps

    def test_reset_clears_results(self, simulation):
        simulation.run(n_jumps=10)
        simulation.reset()
        assert simulation.results is None

    def test_reset_keeps_atom_positions(self, simulation):
        simulation.run(n_jumps=10)
        occupied_sites = simulation.occupied_sites
        simulation.reset()
        assert simulation.occupied_sites == occupied_sites

    def test_same_seed_gives_the_same_results(self):
        simulations = [Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=4, seed=3) for _ in range(2)]
        for simulation in simulations:
            simulation.run(n_jumps=100)
        assert results_of(simulations[0]) == results_of(simulations[1])

    def test_split_run_equals_one_run(self):
        split = Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=4, seed=3)
        split.run(n_jumps=300)
        split.run(n_jumps=700)
        whole = Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=4, seed=3)
        whole.run(n_jumps=1000)
        assert results_of(split) == results_of(whole)


class TestBlocked:
    def test_partly_occupied_lattice_is_not_blocked(self):
        assert Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=4).is_blocked() is False

    def test_fully_occupied_lattice_is_blocked_before_any_run(self):
        assert Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=16).is_blocked() is True

    def test_run_on_a_blocked_lattice_raises_and_leaves_no_results(self):
        simulation = Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=16)
        with pytest.raises(BlockedLatticeError):
            simulation.run(n_jumps=1)
        assert simulation.results is None


class TestLatticeCopy:
    def test_later_changes_to_the_callers_lattice_do_not_change_results(self):
        reference = Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=4, seed=5)
        reference.run(n_jumps=100)
        lattice = square_lattice(4, 4, 1.0)
        simulation = Simulation(lattice, PARAMS, n_atoms=4, seed=5)
        lattice.set_site_energies(by_site={number: 0.5 * (number % 2) for number in range(1, 17)})
        for site in lattice.sites:
            site.r = site.r * 0.5
        simulation.run(n_jumps=100)
        assert results_of(simulation) == results_of(reference)

    def test_lattice_is_a_copy(self):
        lattice = square_lattice(4, 4, 1.0)
        simulation = Simulation(lattice, PARAMS, n_atoms=4)
        assert simulation.lattice is not lattice
        assert simulation.lattice.number_of_sites == lattice.number_of_sites

    def test_lattice_is_a_new_copy_on_each_access(self):
        simulation = Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=4)
        assert simulation.lattice is not simulation.lattice

    def test_changes_to_the_returned_lattice_do_not_change_results(self):
        reference = Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=4, seed=5)
        reference.run(n_jumps=100)
        simulation = Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=4, seed=5)
        lattice = simulation.lattice
        lattice.transmute_sites("L", "M", 8, seed=0)
        lattice.set_site_energies(by_label={"M": 0.5})
        for site in lattice.sites:
            site.r = site.r * 0.5
        simulation.run(n_jumps=100)
        assert results_of(simulation) == results_of(reference)


class TestLattices:
    def test_sites_file_with_negative_unsorted_numbers_runs(self):
        lattice = lattice_from_sites_file(str(LLZO), cell_lengths=[49.0672361] * 3)
        simulation = Simulation(lattice, PARAMS, n_atoms=50, seed=0)
        simulation.run(n_jumps=20)
        site_numbers = {site.number for site in lattice.sites}
        assert len(set(simulation.occupied_sites)) == 50
        assert set(simulation.occupied_sites) <= site_numbers

    def test_two_dimensional_lattice_gives_results(self):
        simulation = Simulation(square_lattice(4, 4, 1.0), PARAMS, n_atoms=4, seed=0)
        simulation.run(n_jumps=10)
        assert simulation.lattice.dimensionality == 2
        assert isinstance(simulation.results, Results)

    def test_cn_energies_run_on_the_fly(self):
        cn_energies = {"L": {"L": {count: -0.01 * count for count in range(7)}}}
        simulation = Simulation(cubic_lattice(3, 3, 3, 1.0), PARAMS, n_atoms=5, cn_energies=cn_energies, seed=0)
        simulation.run(n_jumps=10)
        assert results_of(simulation).number_of_jumps == 10
