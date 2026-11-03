"""Set-up, running, and results of lattice-gas Monte Carlo simulations."""

from __future__ import annotations

import copy
import math
from collections.abc import Callable, Iterable
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from lattice_mc._energies import CNEnergies, EnergyModel, metropolis
from lattice_mc._kmc import possible_jumps, step
from lattice_mc._lookup import LookupTable
from lattice_mc._state import SimulationState
from lattice_mc.constants import k_boltzmann
from lattice_mc.lattice import Lattice
from lattice_mc.results import Results


@dataclass(frozen=True)
class SimulationParameters:
    """
    Immutable container for the physical parameters of a simulation.

    Raises:
        ValueError: If `temperature` or `rate_prefactor` is not finite and positive.
    """

    temperature: float
    """Temperature, in K. Must be finite and positive."""

    rate_prefactor: float
    """Jump rate prefactor, nu_0, in s^-1. Must be finite and positive."""

    def __post_init__(self) -> None:
        """
        Check that the temperature and rate prefactor are finite and positive.
        """
        if not (math.isfinite(self.temperature) and self.temperature > 0):
            raise ValueError(f"temperature must be finite and positive; got {self.temperature!r}.")
        if not (math.isfinite(self.rate_prefactor) and self.rate_prefactor > 0):
            raise ValueError(f"rate_prefactor must be finite and positive; got {self.rate_prefactor!r}.")

    @property
    def kT(self) -> float:
        """Thermal energy kT in eV."""
        return k_boltzmann * self.temperature


class Simulation:
    """
    A lattice-gas kinetic Monte Carlo simulation.

    The simulation takes its own copy of the lattice when it is created, and reads the site
    energies from that copy, so later changes to the caller's lattice have no effect.
    """

    def __init__(
        self,
        lattice: Lattice,
        params: SimulationParameters,
        *,
        n_atoms: int | None = None,
        occupied_sites: Iterable[int] | None = None,
        selected_labels: str | Iterable[str] | None = None,
        nn_energy: float | None = None,
        cn_energies: CNEnergies | None = None,
        lookup_table: bool = False,
        seed: int | np.random.Generator | None = None,
    ) -> None:
        """
        Initialise a Simulation.

        Args:
            lattice: The lattice to simulate. The simulation uses its own copy.
            params: Physical parameters for the simulation (temperature, rate prefactor).
            n_atoms: Number of atoms to place on randomly chosen sites. Give either this or
                `occupied_sites`.
            occupied_sites: Numbers of the sites to place one atom on each. Give either this
                or `n_atoms`.
            selected_labels: Site labels, or a single site label, to restrict the random
                placement under `n_atoms` to. Defaults to None, for sites with any label.
            nn_energy: Interaction energy, in eV, for each pair of occupied nearest-neighbour
                sites. A positive energy is a repulsion. Defaults to None, for no interaction.
            cn_energies: Energies, in eV, for each occupied site, nested by site label, then by
                the label of its neighbouring sites, then by the number of those neighbours that
                are occupied. The energy of an occupied site is the sum of the energies for each
                neighbouring label; vacant sites contribute no energy. The mapping must include
                every site label, every label of a neighbouring site, and every count from 0 to
                the number of those neighbours. For example, on a honeycomb lattice built with
                ``alternating_sites=True``, every 'A' site has three 'B' neighbours::

                    cn = {0: 0.0, 1: 0.1, 2: 0.2, 3: 0.3}
                    {'A': {'B': cn}, 'B': {'A': cn}}

                Defaults to None, for no coordination-number dependence.
            lookup_table: If True, precompute the jump probabilities, keyed by site label and
                numbers of occupied neighbours, instead of calculating each as it is needed.
                The table supports site energies and `nn_energy` only. Defaults to False.
            seed: Seed, or random number generator, for all randomness in the simulation:
                the initial placement of atoms and the choice and timing of jumps. A generator
                is used directly, not copied, so it is advanced in place as the simulation runs.
                The same integer seed gives the same results. Defaults to None, for an
                unseeded generator.

        Raises:
            ValueError: If not exactly one of `n_atoms` and `occupied_sites` is given; if there
                are no atoms, or more atoms than (selected) sites; if `selected_labels` is given
                with `occupied_sites`, is empty, or names a label not in the lattice; if `occupied_sites`
                has a duplicate or unknown site number; if `nn_energy` or an energy in
                `cn_energies` is not finite; if `cn_energies` is incomplete; or if
                `lookup_table` is True with `cn_energies`, or with sites of one label that have
                different energies.
        """
        self._lattice: Lattice = copy.deepcopy(lattice)
        self._params: SimulationParameters = params
        self._kT: float = params.kT
        self._rng: np.random.Generator = np.random.default_rng(seed)
        self._energy_model: EnergyModel = EnergyModel(self._lattice, nn_energy=nn_energy, cn_energies=cn_energies)
        self._probability: Callable[[npt.NDArray[np.bool_], int, int], float] = (
            LookupTable(self._energy_model, self._kT).relative_probability
            if lookup_table
            else self._metropolis_probability
        )
        if n_atoms is not None and occupied_sites is None:
            occupied = self._random_sites(n_atoms, selected_labels)
        elif occupied_sites is not None and n_atoms is None:
            if selected_labels is not None:
                raise ValueError("selected_labels can only be given with n_atoms.")
            occupied = self._site_indices(occupied_sites)
        else:
            raise ValueError("Give exactly one of n_atoms and occupied_sites.")
        self._state: SimulationState = SimulationState.from_occupied(self._lattice.number_of_sites, occupied)

    def _random_sites(self, n_atoms: int, selected_labels: str | Iterable[str] | None) -> list[int]:
        """
        Choose distinct sites at random for the atoms.

        Args:
            n_atoms: The number of atoms.
            selected_labels: Site labels, or a single site label, to choose from, or None for
                any label.

        Returns:
            The indices of the chosen sites.

        Raises:
            ValueError: If `n_atoms` is not positive or is more than the number of sites to
                choose from, or if `selected_labels` is empty or has a label not in the lattice.
        """
        if selected_labels is None:
            labels = self._lattice.site_labels
        else:
            labels = {selected_labels} if isinstance(selected_labels, str) else set(selected_labels)
        if not labels:
            raise ValueError("selected_labels must contain at least one site label.")
        unknown_labels = labels - self._lattice.site_labels
        if unknown_labels:
            raise ValueError(f"Unknown site labels: {sorted(unknown_labels)}.")
        candidates = [index for index, site in enumerate(self._lattice.sites) if site.label in labels]
        if not 0 < n_atoms <= len(candidates):
            raise ValueError(f"n_atoms must be between 1 and {len(candidates)}; got {n_atoms}.")
        return [int(index) for index in self._rng.choice(candidates, size=n_atoms, replace=False)]

    def _site_indices(self, site_numbers: Iterable[int]) -> list[int]:
        """
        The indices of the sites with the given numbers.

        Args:
            site_numbers: The site numbers.

        Returns:
            The index of each site.

        Raises:
            ValueError: If there are no site numbers, or a site number is not in the lattice.
        """
        numbers = list(site_numbers)
        if not numbers:
            raise ValueError("occupied_sites must contain at least one site number.")
        try:
            return [self._lattice.site_index(number) for number in numbers]
        except KeyError as error:
            raise ValueError(f"Unknown site number: {error.args[0]}.") from None

    def _metropolis_probability(self, occupied: npt.NDArray[np.bool_], initial: int, final: int) -> float:
        """
        The relative probability of a jump, calculated from the energy model.

        Args:
            occupied: Whether each site is occupied, before the jump.
            initial: Index of the site the atom jumps from.
            final: Index of the site the atom jumps to.

        Returns:
            The Metropolis relative probability of the jump.
        """
        return metropolis(self._energy_model.delta_E(occupied, initial, final), self._kT)

    @property
    def lattice(self) -> Lattice:
        """
        A copy of the simulation's lattice.

        Each access returns a new copy, so changes to it have no effect on the simulation.
        """
        return copy.deepcopy(self._lattice)

    @property
    def occupied_sites(self) -> tuple[int, ...]:
        """The numbers of the occupied sites, in the order of `Lattice.sites`."""
        return tuple(int(self._lattice.sites[index].number) for index in np.flatnonzero(self._state.occupied))

    @property
    def results(self) -> Results | None:
        """The results accumulated since construction or the last reset, or None if no jump has been made."""
        return Results.from_state(self._state, self._lattice)

    def is_blocked(self) -> bool:
        """
        Check whether no jump is possible.

        Returns:
            True if no occupied site has a vacant neighbour, otherwise False.
        """
        return not possible_jumps(self._lattice, self._state.occupied)

    def run(self, *, n_jumps: int | None = None, for_time: float | None = None) -> None:
        """
        Advance the simulation.

        Results accumulate across runs, so with the same integer seed ``run(n_jumps=300)`` then
        ``run(n_jumps=700)`` gives the same results as ``run(n_jumps=1000)``.

        Args:
            n_jumps: Number of jumps to make. Zero makes no jumps.
            for_time: Run until the time elapsed since construction or the last reset
                reaches this value. If it has already been reached, including when it is
                zero, no jumps are made.

        Raises:
            ValueError: If not exactly one of `n_jumps` and `for_time` is given, if `n_jumps`
                is negative, or if `for_time` is negative or not finite.
            BlockedLatticeError: If no jump is possible. Jumps made before this are kept.
        """
        if n_jumps is not None and for_time is None:
            if n_jumps < 0:
                raise ValueError(f"n_jumps must not be negative; got {n_jumps}.")
            for _ in range(n_jumps):
                self._step()
        elif for_time is not None and n_jumps is None:
            if not (math.isfinite(for_time) and for_time >= 0.0):
                raise ValueError(f"for_time must be finite and not negative; got {for_time}.")
            while self._state.time < for_time:
                self._step()
        else:
            raise ValueError("Give exactly one of n_jumps and for_time.")

    def _step(self) -> None:
        """Make one jump."""
        step(self._lattice, self._state, self._probability, self._params, self._rng)

    def reset(self) -> None:
        """Zero the elapsed time, displacements and other accumulated counters, leaving the atoms where they are."""
        self._state.reset()
