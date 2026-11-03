"""Functions that build lattices, either programmatically or from a sites file."""

from __future__ import annotations

import re
from dataclasses import dataclass
from math import sqrt

import numpy as np
import numpy.typing as npt

from lattice_mc import lattice, lattice_site


def square_lattice(a: int, b: int, spacing: float) -> lattice.Lattice:
    """
    Generate a square lattice.

    Args:
        a: Number of lattice repeat units along x.
        b: Number of lattice repeat units along y.
        spacing: Distance between lattice sites.

    Returns:
        The new lattice.

    Notes:
        The returned lattice is two-dimensional: all sites lie in the xy plane and the z cell length is 0.
    """
    grid = np.array(list(range(1, a * b + 1))).reshape(a, b, order="F")
    it = np.nditer(grid, flags=["multi_index"])
    sites: list[lattice_site.Site] = []
    while not it.finished:
        x, y = it.multi_index
        r = np.array([x * spacing, y * spacing, 0.0])
        neighbours = [
            np.roll(grid, +1, axis=0)[x, y],
            np.roll(grid, -1, axis=0)[x, y],
            np.roll(grid, +1, axis=1)[x, y],
            np.roll(grid, -1, axis=1)[x, y],
        ]
        sites.append(lattice_site.Site(int(it[0]), r, neighbours, 0.0, "L"))
        it.iternext()
    return lattice.Lattice(sites, cell_lengths=np.array([a, b, 0.0]) * spacing)


def honeycomb_lattice(a: int, b: int, spacing: float, alternating_sites: bool = False) -> lattice.Lattice:
    """
    Generate a honeycomb lattice.

    Args:
        a: Number of lattice repeat units along x.
        b: Number of lattice repeat units along y.
        spacing: Distance between lattice sites.
        alternating_sites: If True, label alternating sites 'A' and 'B', so that every 'A' site
            has only 'B' neighbours. Otherwise every site is labelled 'L'. Defaults to False.

    Returns:
        The new lattice.

    Notes:
        The returned lattice is two-dimensional: all sites lie in the xy plane and the z cell length is 0.
    """
    if alternating_sites:
        site_labels = ["A", "B", "A", "B"]
    else:
        site_labels = ["L", "L", "L", "L"]
    unit_cell_lengths = np.array([sqrt(3), 3.0, 0.0]) * spacing
    cell_lengths = unit_cell_lengths * np.array([a, b, 1.0])
    grid = np.array(list(range(1, int(a * b * 4 + 1)))).reshape(a, b, 4, order="C")
    sites: list[lattice_site.Site] = []
    for i in range(a):
        for j in range(b):
            # site 1
            r = np.array([i * sqrt(3) * spacing, j * 3 * spacing, 0.0])
            neighbours = [grid[i, j, 1], np.roll(grid, +1, axis=0)[i, j, 1], np.roll(grid, +1, axis=1)[i, j, 3]]
            sites.append(lattice_site.Site(grid[i, j, 0], r, neighbours, 0.0, site_labels[0]))
            # site 2
            r = np.array([i * sqrt(3) * spacing + sqrt(3) / 2 * spacing, (j * 3 + 0.5) * spacing, 0.0])
            neighbours = [grid[i, j, 0], grid[i, j, 2], np.roll(grid, -1, axis=0)[i, j, 0]]
            sites.append(lattice_site.Site(grid[i, j, 1], r, neighbours, 0.0, site_labels[1]))
            # site 3
            r = np.array([i * sqrt(3) * spacing + sqrt(3) / 2 * spacing, (j * 3 + 1.5) * spacing, 0.0])
            neighbours = [grid[i, j, 1], grid[i, j, 3], np.roll(grid, -1, axis=0)[i, j, 3]]
            sites.append(lattice_site.Site(grid[i, j, 2], r, neighbours, 0.0, site_labels[2]))
            # site 4
            r = np.array([i * sqrt(3) * spacing, (j * 3 + 2) * spacing, 0.0])
            neighbours = [grid[i, j, 2], np.roll(grid, +1, axis=0)[i, j, 2], np.roll(grid, -1, axis=1)[i, j, 0]]
            sites.append(lattice_site.Site(grid[i, j, 3], r, neighbours, 0.0, site_labels[3]))
    return lattice.Lattice(sites, cell_lengths=cell_lengths)


def cubic_lattice(a: int, b: int, c: int, spacing: float) -> lattice.Lattice:
    """
    Generate a cubic lattice.

    Args:
        a: Number of lattice repeat units along x.
        b: Number of lattice repeat units along y.
        c: Number of lattice repeat units along z.
        spacing: Distance between lattice sites.

    Returns:
        The new lattice.
    """
    grid = np.array(list(range(1, a * b * c + 1))).reshape(a, b, c, order="F")
    it = np.nditer(grid, flags=["multi_index"])
    sites: list[lattice_site.Site] = []
    while not it.finished:
        x, y, z = it.multi_index
        r = np.array([x, y, z]) * spacing
        neighbours = [
            np.roll(grid, +1, axis=0)[x, y, z],
            np.roll(grid, -1, axis=0)[x, y, z],
            np.roll(grid, +1, axis=1)[x, y, z],
            np.roll(grid, -1, axis=1)[x, y, z],
            np.roll(grid, +1, axis=2)[x, y, z],
            np.roll(grid, -1, axis=2)[x, y, z],
        ]
        sites.append(lattice_site.Site(int(it[0]), r, neighbours, 0.0, "L"))
        it.iternext()
    return lattice.Lattice(sites, cell_lengths=np.array([a, b, c]) * spacing)


def lattice_from_sites_file(site_file: str, cell_lengths: list[float]) -> lattice.Lattice:
    """
    Generate a lattice from a sites file.

    Args:
        site_file: Path to the file containing the site information.
        cell_lengths: The [x, y, z] lengths of the orthorhombic simulation cell.
            Set the z cell length to 0 for a two-dimensional lattice.

    Returns:
        The new lattice.

    Raises:
        ValueError: If the file is malformed or inconsistent: the site count does
            not match the number of site blocks, a key is unknown, repeated or
            missing, a site number is repeated, a neighbour does not exist, two
            sites do not list each other as neighbours, or sites with the same
            label have different energies.

    Notes:
        The file starts with the number of sites, followed by one block per site.
        Blocks are separated by blank lines. Each block gives::

            site: <site number>
            centre: <x> <y> <z>
            neighbours: <site numbers of neighbouring sites>
            label: <site label>
            energy: <site occupation energy>

        Site numbers are arbitrary integers. The energy line is optional, and
        defaults to 0.0. Lines within a block can be in any order. Both
        centre/center and neighbours/neighbors are accepted, and legacy vertices
        lines are ignored. A full example is
        ``examples/llzo_lattice_site_list.dat`` in the lattice_mc repository.

        For `Lattice.detached_sites`, every bond must span less than half the cell
        length along each axis, so a cell of a single conventional unit cell is
        usually too small.
    """
    with open(site_file) as f:
        records = _read_site_records(f.read())
    sites = [lattice_site.Site(r.number, r.centre, r.neighbours, r.energy, r.label) for r in records]
    return lattice.Lattice(sites, cell_lengths=np.array(cell_lengths))


_KEY_ALIASES = {"center": "centre", "neighbors": "neighbours"}
_REQUIRED_KEYS = ("site", "centre", "neighbours", "label")
_OPTIONAL_KEYS = ("energy", "vertices")


@dataclass(frozen=True)
class _SiteRecord:
    """The contents of one site block in a sites file."""

    number: int
    centre: npt.NDArray[np.float64]
    neighbours: list[int]
    label: str
    energy: float


def _read_site_records(text: str) -> list[_SiteRecord]:
    """
    Parse and check the contents of a sites file.

    Args:
        text: The contents of the sites file.

    Returns:
        One record for each site block, in file order.

    Raises:
        ValueError: If the file is malformed or inconsistent.
    """
    header, *blocks = re.split(r"\n[ \t]*\n", text.strip())
    try:
        number_of_sites = int(header)
    except ValueError:
        raise ValueError(f"The first line of a sites file must be the number of sites, not {header!r}.") from None
    if len(blocks) != number_of_sites:
        raise ValueError(f"The file header gives {number_of_sites} sites, but the file has {len(blocks)} site blocks.")
    records = [_parse_site_block(block, index) for index, block in enumerate(blocks, start=1)]
    _check_site_records(records)
    return records


def _parse_site_block(block: str, index: int) -> _SiteRecord:
    """
    Parse one site block.

    Args:
        block: The lines of the block.
        index: The position of the block in the file, used in error messages.

    Returns:
        The parsed site.

    Raises:
        ValueError: If a line is malformed, a key is unknown or repeated, or a required key is missing.
    """
    fields: dict[str, str] = {}
    for line in block.splitlines():
        if not line.strip():
            continue
        key, colon, value = line.partition(":")
        key = _KEY_ALIASES.get(key.strip(), key.strip())
        if not colon:
            raise ValueError(f"Line {line.strip()!r} in site block {index} has no 'key: value' form.")
        if key not in _REQUIRED_KEYS + _OPTIONAL_KEYS:
            raise ValueError(f"Unknown key {key!r} in site block {index}.")
        if key in fields:
            raise ValueError(f"Key {key!r} appears more than once in site block {index}.")
        fields[key] = value.strip()
    name = f"site {fields['site']}" if "site" in fields else f"site block {index}"
    for key in _REQUIRED_KEYS:
        if key not in fields:
            raise ValueError(f"The block for {name} has no {key} line.")
    try:
        number = int(fields["site"])
        centre = np.array([float(x) for x in fields["centre"].split()])
        neighbours = [int(x) for x in fields["neighbours"].split()]
        energy = float(fields.get("energy", "0.0"))
    except ValueError as error:
        raise ValueError(f"Could not read the block for {name}: {error}") from None
    if centre.shape != (3,):
        raise ValueError(f"The centre of {name} must have three coordinates.")
    if len(fields["label"].split()) != 1:
        raise ValueError(f"The label of {name} must be a single word.")
    return _SiteRecord(number, centre, neighbours, fields["label"], energy)


def _check_site_records(records: list[_SiteRecord]) -> None:
    """
    Check that a set of site records is consistent.

    Args:
        records: The parsed site records.

    Raises:
        ValueError: If a site number is repeated, a neighbour does not exist, a
            neighbour pair is one-way, or sites with the same label have
            different energies.
    """
    neighbours: dict[int, set[int]] = {}
    for record in records:
        if record.number in neighbours:
            raise ValueError(f"The file defines site {record.number} more than once.")
        neighbours[record.number] = set(record.neighbours)
    for number, site_neighbours in neighbours.items():
        for other in sorted(site_neighbours):
            if other not in neighbours:
                raise ValueError(f"Site {number} lists site {other} as a neighbour, but there is no site {other}.")
            if number not in neighbours[other]:
                raise ValueError(
                    f"Site {number} lists site {other} as a neighbour, but site {other} does not list site {number}."
                )
    energies: dict[str, set[float]] = {}
    for record in records:
        energies.setdefault(record.label, set()).add(record.energy)
    for label, label_energies in energies.items():
        if len(label_energies) > 1:
            raise ValueError(
                f"Sites with label {label} have different energies, {sorted(label_energies)}. "
                "Sites with different energies need different labels."
            )
