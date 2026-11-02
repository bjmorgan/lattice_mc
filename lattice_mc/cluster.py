"""Clusters of contiguously connected lattice sites."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:
    from lattice_mc.lattice_site import Site


class Cluster:
    """
    Clusters are sets of sites.
    """

    def __init__(self, sites: Iterable[Site]) -> None:
        """
        Initialise a Cluster instance.

        Args:
            sites: The sites that make up the cluster.
        """
        self.sites: set[Site] = set(sites)
        self.neighbours: set[Site] = set()
        for s in self.sites:
            assert s.p_neighbours is not None
            self.neighbours.update(s.p_neighbours)
        self.neighbours = self.neighbours.difference(self.sites)

    def merge(self, other_cluster: Cluster) -> Cluster:
        """
        Combine two clusters into a single cluster.

        Args:
            other_cluster: The second cluster to combine.

        Returns:
            The combination of both clusters.
        """
        new_cluster = Cluster(self.sites | other_cluster.sites)
        new_cluster.neighbours = (self.neighbours | other_cluster.neighbours).difference(new_cluster.sites)
        return new_cluster

    def is_neighbouring(self, other_cluster: Cluster) -> bool:
        """
        Logical check whether the neighbour list for cluster A includes any sites in cluster B.

        Args:
            other_cluster: The other cluster to test for neighbour connectivity.

        Returns:
            True if the other cluster neighbours this cluster.
        """
        return bool(self.neighbours & other_cluster.sites)

    def size(self) -> int:
        """
        Number of sites in this cluster.

        Returns:
            The number of sites in this cluster.
        """
        return len(self.sites)

    def is_periodically_contiguous(self, cell_lengths: npt.NDArray[np.float64]) -> tuple[bool, bool, bool]:
        """
        Check whether this cluster winds round the periodic boundaries of the simulation cell.

        A cluster that winds round a boundary connects to its own periodic images. Sites
        are assigned integer periodic-image vectors by traversing the cluster from an
        arbitrary root site, using the minimum-image convention to find the wrap count of
        each bond. The cluster is periodically contiguous along an axis if and only if a
        bond reaches an already-visited site with a different image vector on that axis
        (a non-zero winding number).

        Args:
            cell_lengths: Cell lengths for the orthorhombic simulation cell.
                An axis with a zero cell length (e.g. z for a two-dimensional
                lattice) is not periodic, and is never reported as contiguous.

        Returns:
            Contiguity along the x, y, and z axes.

        Notes:
            Assumes the minimum-image convention: every bond must span less
            than half the cell length along each axis. In practice this means
            winding can only be detected along an axis where the cluster is at
            least three sites across the cell.
        """
        root = next(iter(self.sites))
        images: dict[Site, npt.NDArray[np.int_]] = {root: np.zeros(3, dtype=int)}
        stack = [root]
        winds = np.zeros(3, dtype=bool)
        while stack:
            site = stack.pop()
            assert site.p_neighbours is not None
            for neighbour in site.p_neighbours:
                if neighbour not in self.sites:
                    continue
                bond_wrap = np.round(
                    np.divide(neighbour.r - site.r, cell_lengths, out=np.zeros(3), where=cell_lengths != 0)
                ).astype(int)
                neighbour_image = images[site] - bond_wrap
                if neighbour in images:
                    winds |= images[neighbour] != neighbour_image
                else:
                    images[neighbour] = neighbour_image
                    stack.append(neighbour)
        return (bool(winds[0]), bool(winds[1]), bool(winds[2]))

    def remove_sites_from_neighbours(self, remove_labels: list[str] | set[str] | str) -> None:
        """
        Removes sites from the set of neighbouring sites if these have labels in remove_labels.

        Args:
            remove_labels: Site labels, or a single site label, to remove from the cluster neighbour set.
        """
        if type(remove_labels) is str:
            remove_labels = [remove_labels]
        self.neighbours = set(n for n in self.neighbours if n.label not in remove_labels)
