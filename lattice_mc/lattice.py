"""The simulation lattice: its sites, their positions, neighbours, labels and energies."""

from __future__ import annotations

import itertools
from collections import Counter
from collections.abc import Mapping
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from lattice_mc import cluster

if TYPE_CHECKING:
    from lattice_mc.lattice_site import Site


class Lattice:
    """
    A periodic lattice of sites, with their positions, neighbours, labels and occupation energies.
    """

    def __init__(self, sites: list[Site], cell_lengths: npt.NDArray[np.float64]) -> None:
        """
        Initialise a Lattice instance.

        Args:
            sites: The sites contained in the lattice.
            cell_lengths: The [x, y, z] lengths of the orthorhombic simulation cell.
                A cell length of 0 means the lattice has no extent along that axis.

        Raises:
            ValueError: If a site lists itself as a neighbour, or lists the same neighbour
                more than once.
        """
        for site in sites:
            if site.number in site.neighbours:
                raise ValueError(f"Site {site.number} lists itself as a neighbour.")
            repeated = sorted(number for number, count in Counter(site.neighbours).items() if count > 1)
            if repeated:
                raise ValueError(f"Site {site.number} lists site {repeated[0]} more than once as a neighbour.")
        self.cell_lengths: npt.NDArray[np.float64] = cell_lengths
        self.sites: list[Site] = sites
        self.number_of_sites: int = len(self.sites)
        self.enforce_periodic_boundary_conditions()
        self._site_indices: dict[int, int] = {site.number: index for index, site in enumerate(self.sites)}
        self.neighbour_indices: tuple[tuple[int, ...], ...] = tuple(
            tuple(self._site_indices[number] for number in site.neighbours) for site in self.sites
        )

    @property
    def site_labels(self) -> set[str]:
        """The labels of the sites, computed from the sites when read."""
        return {site.label for site in self.sites}

    @property
    def site_populations(self) -> Counter[str]:
        """The number of sites with each label, computed from the sites when read."""
        return Counter(site.label for site in self.sites)

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

    def site_with_id(self, number: int) -> Site:
        """
        Select the site with a specific id number.

        Args:
            number: The identifying number for a specific site.

        Returns:
            The site with id number equal to `number`.

        Raises:
            KeyError: If no site has the id number `number`.
        """
        return self.sites[self._site_indices[number]]

    def site_index(self, number: int) -> int:
        """
        The position of a site in `Lattice.sites`.

        Site numbers are arbitrary identifiers, so this is not related to `number` in general.

        Args:
            number: The identifying number for a specific site.

        Returns:
            The index of the site with id number equal to `number`.

        Raises:
            KeyError: If no site has the id number `number`.
        """
        return self._site_indices[number]

    def neighbours_of(self, site: Site) -> tuple[Site, ...]:
        """
        The neighbouring sites of a site.

        Args:
            site: A site in this lattice.

        Returns:
            The neighbours of `site`, in the order of its `neighbours` list.

        Raises:
            KeyError: If no site in this lattice has the number of `site`.
        """
        return tuple(self.sites[index] for index in self.neighbour_indices[self._site_indices[site.number]])

    def set_site_energies(
        self,
        *,
        by_label: Mapping[str, float] | None = None,
        by_site: Mapping[int, float] | None = None,
    ) -> None:
        """
        Set the occupation energies of the sites in the lattice.

        Energies given by label are applied first, then energies given by site number, so
        a site energy overrides the energy for that site's label. The lattice is unchanged
        if any label or site number is unknown.

        Args:
            by_label: Occupation energy, in eV, for every site with each label, for example
                ``{'A': 1.0, 'B': 0.0}``. Defaults to None.
            by_site: Occupation energy, in eV, for individual sites, keyed by site number.
                Defaults to None.

        Raises:
            ValueError: If neither `by_label` nor `by_site` is given, or if a label
                or site number is not in the lattice.
        """
        if by_label is None and by_site is None:
            raise ValueError("One of by_label or by_site must be given.")
        by_label = by_label or {}
        by_site = by_site or {}
        unknown_labels = set(by_label) - self.site_labels
        if unknown_labels:
            raise ValueError(f"Unknown site labels: {sorted(unknown_labels)}.")
        unknown_numbers = set(by_site) - self._site_indices.keys()
        if unknown_numbers:
            raise ValueError(f"Unknown site numbers: {sorted(unknown_numbers)}.")
        for site in self.sites:
            if site.label in by_label:
                site.energy = by_label[site.label]
        for number, energy in by_site.items():
            self.sites[self._site_indices[number]].energy = energy

    def site_energies(self) -> npt.NDArray[np.float64]:
        """
        The occupation energy of every site.

        Returns:
            The occupation energy, in eV, for each site, in the order of `Lattice.sites`.
        """
        return np.array([site.energy for site in self.sites], dtype=np.float64)

    def site_energies_by_label(self) -> dict[str, float]:
        """
        The occupation energy of the sites with each label.

        Returns:
            The occupation energy, in eV, for each site label.

        Raises:
            ValueError: If sites with the same label have different energies.
        """
        energies: dict[str, float] = {}
        for site in self.sites:
            if energies.setdefault(site.label, site.energy) != site.energy:
                raise ValueError(
                    f"Sites with label {site.label} have different energies, {energies[site.label]} and {site.energy}."
                )
        return energies

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
            specific_coordination_numbers[site.label] = dict(Counter(n.label for n in self.neighbours_of(site)))
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
            for final_site in self.neighbours_of(initial_site):
                if final_site.label not in site_connections[initial_site.label]:
                    site_connections[initial_site.label].append(final_site.label)
        return site_connections

    def transmute_sites(
        self,
        old_site_label: str | list[str],
        new_site_label: str,
        n_sites_to_change: int,
        *,
        seed: int | np.random.Generator | None = None,
    ) -> None:
        """
        Selects a random subset of sites with a specific label and gives them a different label.

        Args:
            old_site_label: Site label, or labels, of the sites that can be modified.
            new_site_label: Site label to apply to the modified sites.
            n_sites_to_change: Number of sites to modify.
            seed: Seed, or random number generator, for choosing the sites. A generator is
                used directly, and is advanced in place. The same integer seed chooses the
                same sites. Defaults to None, for an unseeded generator.

        Raises:
            ValueError: If `n_sites_to_change` is negative or more than the number of sites
                that can be modified.
        """
        selected_sites = self.select_sites(old_site_label)
        rng = np.random.default_rng(seed)
        for index in rng.choice(len(selected_sites), size=n_sites_to_change, replace=False):
            selected_sites[index].label = new_site_label

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
        # A single label is wrapped, because set("Li") would be {"L", "i"}.
        labels = {site_labels} if isinstance(site_labels, str) and site_labels else set(site_labels or ())
        if labels:
            selected_sites = self.select_sites(labels)
        else:
            selected_sites = self.sites
        initial_clusters = [cluster.Cluster([site], self) for site in selected_sites]
        if labels:
            blocking_sites = self.site_labels - labels
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
        island_clusters = [c for c in clusters if not any(c.is_periodically_contiguous())]
        return list(itertools.chain.from_iterable((c.sites for c in island_clusters)))
