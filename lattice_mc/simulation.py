"""Set-up, running, and results of lattice-gas Monte Carlo simulations."""

from __future__ import annotations

from dataclasses import dataclass

from lattice_mc import init_lattice, lookup_table, species
from lattice_mc.constants import k_boltzmann
from lattice_mc.lattice import Lattice


@dataclass(frozen=True)
class SimulationParameters:
    """
    Immutable container for the physical parameters of a simulation.

    Raises:
        ValueError: If `temperature` or `rate_prefactor` is not positive.
    """

    temperature: float
    """Temperature, in K. Must be positive."""

    rate_prefactor: float
    """Jump rate prefactor, nu_0, in s^-1. Must be positive."""

    def __post_init__(self) -> None:
        """
        Check that the temperature and rate prefactor are positive.
        """
        if self.temperature <= 0:
            raise ValueError(
                f"temperature must be positive; got {self.temperature!r}."
            )
        if self.rate_prefactor <= 0:
            raise ValueError(
                f"rate_prefactor must be positive; got {self.rate_prefactor!r}."
            )

    @property
    def kT(self) -> float:
        """Thermal energy kT in eV."""
        return k_boltzmann * self.temperature


class Simulation:
    """
    A lattice-gas Monte Carlo simulation.
    """

    def __init__(self, params: SimulationParameters) -> None:
        """
        Initialise a Simulation object.

        Args:
            params: Physical parameters for the simulation
                (temperature, rate prefactor).
        """
        self.params: SimulationParameters = params
        self.lattice: Lattice | None = None
        self.number_of_atoms: int | None = None
        self.number_of_jumps: int | None = None
        self.for_time: float | None = None
        self.number_of_equilibration_jumps: int = 0
        self.atoms: species.Species | None = None
        self.has_run: bool = False

    def reset(self) -> None:
        """
        Reset all counters for this simulation.

        Notes:
            Reading the diffusion coefficients or correlation factors after a reset, before
            running again, raises a ZeroDivisionError. See
            https://github.com/bjmorgan/lattice_mc/issues/32.
        """
        assert self.lattice is not None
        assert self.atoms is not None
        self.lattice.reset()
        for atom in self.atoms.atoms:
            atom.reset()

    def set_number_of_atoms(self, n: int, selected_sites: list[str] | None = None) -> None:
        """
        Set the number of atoms for the simulation, and populate the simulation lattice.

        Args:
            n: Number of atoms for this simulation.
            selected_sites: Site labels to populate with atoms. If None, atoms are placed on
                randomly chosen sites of any label. Defaults to None.
        """
        assert self.lattice is not None
        self.number_of_atoms = n
        self.atoms = species.Species(self.lattice.populate_sites(self.number_of_atoms, selected_sites=selected_sites))

    def set_number_of_jumps(self, n: int) -> None:
        """
        Set the number of jumps for this simulation.

        Args:
            n: Number of jumps.
        """
        self.number_of_jumps = n

    def set_number_of_equilibration_jumps(self, n: int) -> None:
        """
        Set the number of equilibration jumps for this simulation.

        Args:
            n: Number of equilibration jumps, run (and then discarded) before the jumps that are
                analysed.
        """
        self.number_of_equilibration_jumps = n

    def define_lattice_from_file(self, filename: str, cell_lengths: list[float]) -> None:
        """
        Set up the simulation lattice from a file containing site data.

        Uses `init_lattice.lattice_from_sites_file`, which defines the site file spec.

        Args:
            filename: Path to the sites file.
            cell_lengths: The [x, y, z] lengths of the orthorhombic simulation cell.
                Set the z cell length to 0 for a two-dimensional lattice.
        """
        self.lattice = init_lattice.lattice_from_sites_file(filename, cell_lengths=cell_lengths)

    def set_nn_energy(self, nn_energy: float | None) -> None:
        """
        Set the nearest-neighbour energy for this simulation.

        Args:
            nn_energy: Interaction energy, in eV, for each pair of occupied nearest-neighbour
                sites. A positive energy is a repulsion. If None, the lattice is left unchanged.
        """
        assert self.lattice is not None
        if nn_energy is not None:
            self.lattice.set_nn_energy(nn_energy)

    def set_cn_energies(self, cn_energies: dict[str, dict[str, dict[int, float]]] | None) -> None:
        """
        Set the coordination-number dependent energies for this simulation.

        Args:
            cn_energies: Energies, in eV, for each occupied site, nested by site label, then by
                the label of its neighbouring sites, then by the number of those neighbours that
                are occupied. The energy of an occupied site is the sum of the energies for each
                neighbouring label; vacant sites contribute no energy. If None, the lattice is left
                unchanged. The dict must be complete:
                every site label, every label of a neighbouring site, and every count from 0 to the
                number of those neighbours, or the simulation raises a KeyError. For example, on a
                honeycomb lattice built with ``alternating_sites=True``, every 'A' site has three
                'B' neighbours::

                    cn = {0: 0.0, 1: 0.1, 2: 0.2, 3: 0.3}
                    {'A': {'B': cn}, 'B': {'A': cn}}
        """
        assert self.lattice is not None
        if cn_energies is not None:
            self.lattice.set_cn_energies(cn_energies)

    def set_site_energies(self, site_energies: dict[str, float] | None) -> None:
        """
        Set the on-site energies for this simulation.

        Args:
            site_energies: Occupation energy, in eV, for each site label, for example
                ``{'A': 0.0, 'B': 1.0}``. If None, the lattice is left unchanged.
        """
        assert self.lattice is not None
        if site_energies is not None:
            self.lattice.set_site_energies(by_label=site_energies)

    def is_initialised(self) -> None:
        """
        Check whether the simulation has been initialised.

        Raises:
            AttributeError: If the lattice or the atoms have not been set, or if neither a number
                of jumps nor a run time has been set.
        """
        if not self.lattice:
            raise AttributeError("Running a simulation needs the lattice to be initialised")
        if not self.atoms:
            raise AttributeError("Running a simulation needs the atoms to be initialised")
        if not self.number_of_jumps and not self.for_time:
            raise AttributeError("Running a simulation needs number_of_jumps or for_time to be set")

    def run(self, for_time: float | None = None) -> None:
        """
        Run the simulation.

        Args:
            for_time: If given, run until the simulation time reaches this value, instead of for
                the number of jumps set by `set_number_of_jumps`. Any equilibration jumps are run
                first, and the time is counted from the end of equilibration. After the run,
                `number_of_jumps` is the number of jumps made. Defaults to None.

        Raises:
            AttributeError: If the simulation has not been initialised.
            BlockedLatticeError: If no jumps are possible on the lattice.
        """
        self.for_time = for_time
        self.is_initialised()
        assert self.lattice is not None
        assert self.atoms is not None
        self.lattice.params = self.params
        if self.number_of_equilibration_jumps > 0:
            for step in range(self.number_of_equilibration_jumps):
                self.lattice.jump()
            self.reset()
        if self.for_time:
            self.number_of_jumps = 0
            while self.lattice.time < self.for_time:
                self.lattice.jump()
                self.number_of_jumps += 1
        else:
            assert self.number_of_jumps is not None
            for step in range(self.number_of_jumps):
                self.lattice.jump()
        self.has_run = True

    @property
    def tracer_correlation(self) -> float | None:
        """
        Tracer correlation factor, f.

        Returns:
            The tracer correlation factor, f, or None if the simulation has not been run.
        """
        if self.has_run:
            assert self.atoms is not None
            return self.atoms.tracer_correlation()
        else:
            return None

    @property
    def tracer_diffusion_coefficient(self) -> float | None:
        """
        Tracer diffusion coefficient, D*.

        Calculated from the Einstein relation, <r^2> = 2 d D* t, where <r^2> is
        the mean squared displacement over all atoms and d is the dimensionality
        of the lattice.

        Returns:
            The tracer diffusion coefficient, D*, or None if the simulation has not been run.
        """
        if self.has_run:
            assert self.atoms is not None
            assert self.lattice is not None
            assert self.number_of_atoms is not None
            return self.atoms.sum_dr_squared() / (
                2.0 * self.lattice.dimensionality * float(self.number_of_atoms) * self.lattice.time
            )
        else:
            return None

    @property
    def collective_correlation(self) -> float | None:
        """
        Collective correlation factor, f_I.

        Returns:
            The collective correlation factor, f_I, or None if the simulation has not been run.
        """
        if self.has_run:
            assert self.atoms is not None
            return self.atoms.collective_correlation()
        else:
            return None

    @property
    def collective_diffusion_coefficient(self) -> float | None:
        """
        Collective, or "jump", diffusion coefficient, D_J.

        Calculated from the Einstein relation, <R^2> = 2 d D_J t, where
        R = sum_i dr_i is the summed displacement of all atoms (not the centre of
        mass displacement) and d is the dimensionality of the lattice. No 1/N
        normalisation is applied; see `collective_diffusion_coefficient_per_atom`.

        Returns:
            The collective diffusion coefficient, D_J, or None if the simulation has not been run.
        """
        if self.has_run:
            assert self.atoms is not None
            assert self.lattice is not None
            return self.atoms.collective_dr_squared() / (2.0 * self.lattice.dimensionality * self.lattice.time)
        else:
            return None

    @property
    def collective_diffusion_coefficient_per_atom(self) -> float | None:
        """
        Collective diffusion coefficient per atom, D_J / n_atoms.

        Returns:
            The collective diffusion coefficient per atom, or None if the simulation has not been run.
        """
        if self.has_run:
            assert self.number_of_atoms is not None
            assert self.collective_diffusion_coefficient is not None
            return self.collective_diffusion_coefficient / float(self.number_of_atoms)
        else:
            return None

    @property
    def average_site_occupations(self) -> dict[str, float] | None:
        """
        Average site occupation numbers.

        Returns:
            Average site occupation numbers for each site label, for example
            ``{'A': 12.4, 'B': 231.2}``, or None if no simulation time has elapsed.
        """
        assert self.lattice is not None
        return self.lattice.site_occupation_statistics()

    def setup_lookup_table(self, hamiltonian: str = "nearest-neighbour") -> None:
        """
        Create a look-up table of jump probabilities, to speed up the simulation.

        Without a table, the probability of each jump is calculated as it is needed. The table
        includes the site energies, whether read from a sites file or set with
        `set_site_energies`, and the nearest-neighbour energy, as they are when the table is
        created, so set those first. It supports nearest-neighbour interactions only, so do not
        use a table with coordination-number energies.

        Args:
            hamiltonian: The model Hamiltonian. The only supported value is 'nearest-neighbour'.
                Defaults to 'nearest-neighbour'.

        Raises:
            ValueError: If `hamiltonian` is not supported, or if sites with the same label have
                different energies.
        """
        expected_hamiltonian_values = ["nearest-neighbour"]
        if hamiltonian not in expected_hamiltonian_values:
            raise ValueError(
                f"Unsupported hamiltonian {hamiltonian!r}. "
                f"Expected one of {expected_hamiltonian_values!r}."
            )
        assert self.lattice is not None
        self.lattice.params = self.params
        self.lattice.jump_lookup_table = lookup_table.LookupTable(self.lattice, hamiltonian)
