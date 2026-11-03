"""The simulation lattice: its sites, the atoms occupying them, and the possible jumps."""

from __future__ import annotations

import itertools
import random
from collections import Counter
from collections.abc import Iterator
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from lattice_mc import atom, cluster, jump, transitions
from lattice_mc.error import BlockedLatticeError

if TYPE_CHECKING:
    from lattice_mc.atom import Atom
    from lattice_mc.jump import Jump
    from lattice_mc.lattice_site import Site
    from lattice_mc.lookup_table import LookupTable
    from lattice_mc.simulation import SimulationParameters


class Lattice:
    """
    A periodic lattice of sites, and the atoms that occupy them.
    """

    def __init__(self, sites: list[Site], cell_lengths: npt.NDArray[np.float64]) -> None:
        """
        Initialise a Lattice instance.

        Args:
            sites: The sites contained in the lattice.
            cell_lengths: The [x, y, z] lengths of the orthorhombic simulation cell.
                A cell length of 0 means the lattice has no extent along that axis.
        """
        self.cell_lengths: npt.NDArray[np.float64] = cell_lengths
        self.sites: list[Site] = sites
        self.number_of_sites: int = len(self.sites)
        self.site_labels: set[str] = set([site.label for site in self.sites])
        self.site_populations: Counter[str] = Counter([site.label for site in self.sites])
        self.enforce_periodic_boundary_conditions()
        self.initialise_site_lookup_table()
        self.params: SimulationParameters | None = None
        self.nn_energy: float | None = None
        self.cn_energies: dict[str, dict[str, dict[int, float]]] | None = None
        self.jump_lookup_table: LookupTable | None = None
        self.number_of_occupied_sites: int = 0
        for site in self.sites:
            site.p_neighbours = [self.site_with_id(i) for i in site.neighbours]
        self.reset()

    @property
    def dimensionality(self) -> int:
        """
        Number of dimensions spanned by the lattice.

        A cell length of zero means the lattice has no extent along that axis,
        so a two-dimensional lattice is defined with a z cell length of 0.

        Returns:
            The number of non-zero cell lengths.
        """
        return int(np.count_nonzero(self.cell_lengths))

    def enforce_periodic_boundary_conditions(self) -> None:
        """
        Ensure that all lattice sites are within the central periodic image of the simulation cell.

        Sites that are outside the central simulation cell are mapped back into this cell.
        """
        for s in self.sites:
            for i in range(3):
                if s.r[i] < 0.0:
                    s.r[i] += self.cell_lengths[i]
                if s.r[i] > self.cell_lengths[i]:
                    s.r[i] -= self.cell_lengths[i]

    def reset(self) -> None:
        """
        Reset all time-dependent counters for this lattice and its constituent sites.
        """
        self.time: float = 0.0
        for site in self.sites:
            site.time_occupied = 0.0

    def initialise_site_lookup_table(self) -> None:
        """
        Create a lookup table of sites by number.

        After this, `self.site_lookup[n]` is the site whose identifying number is `n`.
        """
        self.site_lookup: dict[int, Site] = {}
        for site in self.sites:
            self.site_lookup[site.number] = site

    def site_with_id(self, number: int) -> Site:
        """
        Select the site with a specific id number.

        Args:
            number: The identifying number for a specific site.

        Returns:
            The site with id number equal to `number`.
        """
        return self.site_lookup[number]

    def vacant_sites(self) -> Iterator[Site]:
        """
        The sites not occupied by atoms.

        Returns:
            An iterator over the vacant sites.
        """
        return (site for site in self.sites if not site.is_occupied)

    def occupied_sites(self) -> Iterator[Site]:
        """
        The sites occupied by atoms.

        Returns:
            An iterator over the occupied sites.
        """
        return (site for site in self.sites if site.is_occupied)

    def vacant_site_numbers(self) -> list[int]:
        """
        List of site id numbers for all sites that are vacant.

        Returns:
            The id numbers of the vacant sites.
        """
        return [site.number for site in self.sites if not site.is_occupied]

    def occupied_site_numbers(self) -> list[int]:
        """
        List of site id numbers for all sites that are occupied.

        Returns:
            The id numbers of the occupied sites.
        """
        return [site.number for site in self.sites if site.is_occupied]

    def potential_jumps(self) -> list[Jump]:
        """
        All nearest-neighbour jumps not blocked by volume exclusion.

        These are jumps from occupied sites to neighbouring unoccupied sites.

        Returns:
            The possible jumps.

        Raises:
            RuntimeError: If `self.params` has not been set.
        """
        if self.params is None:
            raise RuntimeError("Lattice.params must be set before computing jumps")
        jumps = []
        if self.number_of_occupied_sites <= self.number_of_sites / 2:
            for occupied_site in self.occupied_sites():
                unoccupied_neighbours = [
                    site for site in [self.site_with_id(n) for n in occupied_site.neighbours] if not site.is_occupied
                ]
                for vacant_site in unoccupied_neighbours:
                    jumps.append(
                        jump.Jump(
                            occupied_site,
                            vacant_site,
                            self.nn_energy,
                            self.cn_energies,
                            self.jump_lookup_table,
                            params=self.params,
                        )
                    )
        else:
            for vacant_site in self.vacant_sites():
                occupied_neighbours = [
                    site for site in [self.site_with_id(n) for n in vacant_site.neighbours] if site.is_occupied
                ]
                for occupied_site in occupied_neighbours:
                    jumps.append(
                        jump.Jump(
                            occupied_site,
                            vacant_site,
                            self.nn_energy,
                            self.cn_energies,
                            self.jump_lookup_table,
                            params=self.params,
                        )
                    )
        return jumps

    def update(self, accepted_jump: Jump) -> None:
        """
        Update the lattice state by accepting a specific jump.

        Args:
            accepted_jump: The jump that has been accepted.
        """
        jumping_atom = accepted_jump.initial_site.atom
        assert jumping_atom is not None
        dr = accepted_jump.dr(self.cell_lengths)
        accepted_jump.final_site.occupation = jumping_atom.number
        accepted_jump.final_site.atom = jumping_atom
        accepted_jump.final_site.is_occupied = True
        accepted_jump.initial_site.occupation = 0
        accepted_jump.initial_site.atom = None
        accepted_jump.initial_site.is_occupied = False
        # TODO: updating atom counters could be contained in an atom.move_to( site ) method
        jumping_atom.site = accepted_jump.final_site
        jumping_atom.number_of_hops += 1
        jumping_atom.dr += dr
        jumping_atom.summed_dr2 += np.dot(dr, dr)

    def populate_sites(self, number_of_atoms: int, selected_sites: list[str] | None = None) -> list[Atom]:
        """
        Populate the lattice sites with a specific number of atoms.

        Args:
            number_of_atoms: The number of atoms to populate the lattice sites with.
            selected_sites: Site labels to populate with atoms. If None, atoms are placed on
                randomly chosen sites of any label. Defaults to None.

        Returns:
            The new atoms.

        Raises:
            ValueError: If there are more atoms than sites, or than sites with the labels in
                `selected_sites`.
        """
        if number_of_atoms > self.number_of_sites:
            raise ValueError
        if selected_sites:
            atoms = [
                atom.Atom(initial_site=site)
                for site in random.sample([s for s in self.sites if s.label in selected_sites], number_of_atoms)
            ]
        else:
            atoms = [atom.Atom(initial_site=site) for site in random.sample(self.sites, number_of_atoms)]
        self.number_of_occupied_sites = number_of_atoms
        return atoms

    def jump(self) -> None:
        """
        Select a jump at random from all potential jumps, then update the lattice state.

        Raises:
            BlockedLatticeError: If no jumps are possible.
            RuntimeError: If `self.params` has not been set.
        """
        potential_jumps = self.potential_jumps()
        if not potential_jumps:
            raise BlockedLatticeError("No moves are possible in this lattice")
        assert self.params is not None
        all_transitions = transitions.Transitions(potential_jumps, params=self.params)
        random_jump = all_transitions.random()
        delta_t = all_transitions.time_to_jump()
        self.time += delta_t
        self.update_site_occupation_times(delta_t)
        self.update(random_jump)

    def update_site_occupation_times(self, delta_t: float) -> None:
        """
        Increase the time occupied for all occupied sites by delta t.

        Args:
            delta_t: The time step.
        """
        for site in self.occupied_sites():
            site.time_occupied += delta_t

    def site_occupation_statistics(self) -> dict[str, float] | None:
        """
        Average site occupation for each site label.

        Returns:
            The time-averaged number of atoms on the sites with each label, for example
            ``{'A': 2.5, 'B': 25.3}``, or None if no time has elapsed.
        """
        if self.time == 0.0:
            return None
        occupation_stats = {label: 0.0 for label in self.site_labels}
        for site in self.sites:
            occupation_stats[site.label] += site.time_occupied
        for label in self.site_labels:
            occupation_stats[label] /= self.time
        return occupation_stats

    def set_site_energies(self, energies: dict[str, float]) -> None:
        """
        Set the energies for every site in the lattice according to the site labels.

        Args:
            energies: Occupation energy, in eV, for each site label, for example
                ``{'A': 1.0, 'B': 0.0}``.
        """
        for site_label in energies:
            for site in self.sites:
                if site.label == site_label:
                    site.energy = energies[site_label]

    def site_energies_by_label(self) -> dict[str, float]:
        """
        The occupation energy of the sites with each label.

        Sites with the same label can have different energies, but the lookup table
        is keyed by label, so it needs every site with a given label to have the same
        energy.

        Returns:
            The occupation energy, in eV, for each site label.

        Raises:
            ValueError: If sites with the same label have different energies.
        """
        energies: dict[str, float] = {}
        for site in self.sites:
            if energies.setdefault(site.label, site.energy) != site.energy:
                raise ValueError(
                    f"Sites with label {site.label} have different energies, "
                    f"{energies[site.label]} and {site.energy}."
                )
        return energies

    def set_nn_energy(self, delta_E: float) -> None:
        """
        Set the lattice nearest-neighbour energy.

        Args:
            delta_E: Interaction energy, in eV, for each pair of occupied nearest-neighbour sites.
                A positive energy is a repulsion.
        """
        self.nn_energy = delta_E

    def set_cn_energies(self, cn_energies: dict[str, dict[str, dict[int, float]]]) -> None:
        """
        Set the coordination-number dependent energies for this lattice.

        Args:
            cn_energies: Energies, in eV, for each occupied site, nested by site label, then by
                the label of its neighbouring sites, then by the number of those neighbours that
                are occupied. The energy of an occupied site is the sum of the energies for each
                neighbouring label; vacant sites contribute no energy. The dict must be complete:
                every site label, every label of a neighbouring site, and every count from 0 to the
                number of those neighbours, or the simulation raises a KeyError. For example, on a
                honeycomb lattice built with ``alternating_sites=True``, every 'A' site has three
                'B' neighbours::

                    cn = {0: 0.0, 1: 0.1, 2: 0.2, 3: 0.3}
                    {'A': {'B': cn}, 'B': {'A': cn}}
        """
        for site in self.sites:
            site.set_cn_occupation_energies(cn_energies[site.label])
        self.cn_energies = cn_energies

    def site_coordination_numbers(self) -> dict[str, set[int]]:
        """
        Returns a dictionary of the coordination numbers for each site label.

        For example::

            { 'A' : { 4 }, 'B' : { 2, 4 } }

        Returns:
            The set of coordination numbers found for the sites with each label.
        """
        coordination_numbers = {}
        for label in self.site_labels:
            coordination_numbers[label] = set([len(site.neighbours) for site in self.sites if site.label == label])
        return coordination_numbers

    def max_site_coordination_numbers(self) -> dict[str, int]:
        """
        Returns a dictionary of the maximum coordination number for each site label.

        For example::

            { 'A' : 4, 'B' : 4 }

        Returns:
            The maximum coordination number for the sites with each label.
        """
        return {label: max(c) for label, c in self.site_coordination_numbers().items()}

    def site_specific_coordination_numbers(self) -> dict[str, dict[str, int]]:
        """
        Returns the number of neighbours with each label, for each site label.

        The counts for each site label are taken from one site with that label, so this
        assumes that all sites with the same label have the same neighbours.

        Returns:
            For each site label, the number of neighbouring sites with each label, for example
            ``{'A': {'B': 3}, 'B': {'A': 3}}``.
        """
        specific_coordination_numbers = {}
        for site in self.sites:
            specific_coordination_numbers[site.label] = site.site_specific_neighbours()
        return specific_coordination_numbers

    def connected_site_pairs(self) -> dict[str, list[str]]:
        """
        Returns a dictionary of all connections between pairs of sites (by site label).

        For example, for a linear lattice A-B-C this returns::

            { 'A' : [ 'B' ], 'B' : [ 'A', 'C' ], 'C' : [ 'B' ] }

        Returns:
            For each site label, the labels of its neighbouring sites.
        """
        site_connections: dict[str, list[str]] = {}
        for initial_site in self.sites:
            if initial_site.label not in site_connections:
                site_connections[initial_site.label] = []
            assert initial_site.p_neighbours is not None
            for final_site in initial_site.p_neighbours:
                if final_site.label not in site_connections[initial_site.label]:
                    site_connections[initial_site.label].append(final_site.label)
        return site_connections

    def transmute_sites(self, old_site_label: str | list[str], new_site_label: str, n_sites_to_change: int) -> None:
        """
        Selects a random subset of sites with a specific label and gives them a different label.

        Args:
            old_site_label: Site label, or labels, of the sites that can be modified.
            new_site_label: Site label to apply to the modified sites.
            n_sites_to_change: Number of sites to modify.
        """
        selected_sites = self.select_sites(old_site_label)
        for site in random.sample(selected_sites, n_sites_to_change):
            site.label = new_site_label
        self.site_labels = set([site.label for site in self.sites])

    def connected_sites(self, site_labels: list[str] | set[str] | str | None = None) -> list[cluster.Cluster]:
        """
        Searches the lattice to find sets of sites that are contiguously neighbouring.

        Mutually exclusive sets of contiguous sites are returned as Cluster objects.

        Args:
            site_labels: Labels of the sites to consider in the search: a list or set of labels,
                such as ``['A', 'B']``, or a single label, such as ``'A'``. If None, all sites
                are considered. Defaults to None.

        Returns:
            The clusters of contiguous sites.
        """
        if site_labels:
            selected_sites = self.select_sites(site_labels)
        else:
            selected_sites = self.sites
        initial_clusters = [cluster.Cluster([site]) for site in selected_sites]
        if site_labels:
            blocking_sites = self.site_labels - set(site_labels)
            for c in initial_clusters:
                c.remove_sites_from_neighbours(blocking_sites)
        final_clusters = []
        while initial_clusters:  # loop until initial_clusters is empty
            this_cluster = initial_clusters.pop(0)
            while this_cluster.neighbours:
                neighbouring_clusters = [c for c in initial_clusters if this_cluster.is_neighbouring(c)]
                for nc in neighbouring_clusters:
                    initial_clusters.remove(nc)
                    this_cluster = this_cluster.merge(nc)
            final_clusters.append(this_cluster)
        return final_clusters

    def select_sites(self, site_labels: list[str] | set[str] | str) -> list[Site]:
        """
        Selects sites in the lattice with specified labels.

        Args:
            site_labels: Labels of the sites to select: a list or set of labels, such as
                ``['A', 'B']``, or a single label, such as ``'A'``.

        Returns:
            The sites with labels given by `site_labels`.

        Raises:
            ValueError: If `site_labels` is not a list, set or string.
        """
        if isinstance(site_labels, (list, set)):
            selected_sites = [s for s in self.sites if s.label in site_labels]
        elif isinstance(site_labels, str):
            selected_sites = [s for s in self.sites if s.label == site_labels]
        else:
            raise ValueError(str(site_labels))
        return selected_sites

    def detached_sites(self, site_labels: list[str] | set[str] | str | None = None) -> list[Site]:
        """
        Returns all sites in the lattice that are not part of a percolating network.

        These are the sites in clusters of connected sites that do not wrap round to themselves
        through a periodic boundary. The search can optionally be restricted to sites with
        specific labels.

        Args:
            site_labels: Labels of the sites to consider. If None, all sites are considered.
                Defaults to None.

        Returns:
            List of sites not in a periodic percolating network.

        Notes:
            Every bond must span less than half the cell length along each axis
            (see `Cluster.is_periodically_contiguous`). In smaller cells, clusters
            that wind round the cell can be reported as detached.
        """
        clusters = self.connected_sites(site_labels=site_labels)
        island_clusters = [c for c in clusters if not any(c.is_periodically_contiguous(self.cell_lengths))]
        return list(itertools.chain.from_iterable((c.sites for c in island_clusters)))

    def is_blocked(self) -> bool:
        """
        Check whether there are any possible jumps.

        Returns:
            True if there are no possible jumps, otherwise False.

        Raises:
            RuntimeError: If `self.params` has not been set. `Simulation.run` and
                `Simulation.setup_lookup_table` set it.

        Notes:
            Checking for a blocked lattice should not need the simulation parameters. See
            https://github.com/bjmorgan/lattice_mc/issues/29.
        """
        if not self.potential_jumps():
            return True
        else:
            return False
