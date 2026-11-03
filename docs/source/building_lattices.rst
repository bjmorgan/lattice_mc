Building a lattice
==================

A simulation runs on a :class:`~lattice_mc.lattice.Lattice`: a set of sites, each
with a position, a label, and a list of neighbouring sites that atoms can jump to.
Lattices are periodic, inside an orthorhombic simulation cell. Simple lattices can be
built programmatically, and lattices of any geometry can be read from a sites file.

Simple lattices
---------------

:func:`~lattice_mc.init_lattice.square_lattice`,
:func:`~lattice_mc.init_lattice.honeycomb_lattice` and
:func:`~lattice_mc.init_lattice.cubic_lattice` build square, honeycomb and cubic
lattices from a number of repeat units along each axis and a site spacing. They are
defined in :mod:`lattice_mc.init_lattice`, and can be imported from the top-level
package.

A square lattice of 4 x 4 sites:

.. testcode::

    import lattice_mc as lmc

    lattice = lmc.square_lattice(4, 4, 1.0)
    print(lattice.number_of_sites, lattice.dimensionality, sorted(lattice.site_labels))

.. testoutput::

    16 2 ['L']

A cubic lattice of 3 x 3 x 3 sites:

.. testcode::

    lattice = lmc.cubic_lattice(3, 3, 3, 1.0)
    print(lattice.number_of_sites, lattice.dimensionality)

.. testoutput::

    27 3

Use at least 3 repeat units along each axis of square and cubic lattices, and at
least 2 along each axis of honeycomb lattices. With fewer, a site's neighbours on
either side along an axis are the same site, so it is listed as a neighbour twice.

Each repeat unit of a honeycomb lattice contains 4 sites. Every site in these
lattices is labelled ``'L'``. A honeycomb lattice can instead
have alternating ``'A'`` and ``'B'`` sites, so that every ``'A'`` site has only
``'B'`` neighbours. This is useful for giving the two sublattices different energies
(see `Site energies`_).

.. testcode::

    lattice = lmc.honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
    print(lattice.number_of_sites, lattice.dimensionality, sorted(lattice.site_labels))

.. testoutput::

    16 2 ['A', 'B']

Two-dimensional lattices
------------------------

Square and honeycomb lattices lie in the xy plane, and have a z cell length of 0.
:attr:`~lattice_mc.lattice.Lattice.dimensionality` is the number of non-zero cell
lengths, and is used as :math:`d` in the Einstein relations for the diffusion
coefficients (see :doc:`introduction`). A cubic lattice with one repeat unit along
z still has a non-zero z cell length, so it counts as three-dimensional.

Lattices from a sites file
--------------------------

For other geometries, such as lattices based on crystallographic data, the sites and
their connectivity are read from a file. The file starts with the number of sites,
followed by one block per site, separated by blank lines:

.. literalinclude:: data/four_site_ring.dat
    :caption: four_site_ring.dat

Each block gives:

``site``
    An integer that identifies the site. Site numbers are arbitrary, and can be
    negative.
``centre``
    The x, y and z coordinates of the site. ``center`` is also accepted.
``neighbours``
    The numbers of the sites that an atom on this site can jump to. ``neighbors`` is
    also accepted.
``label``
    A label for the site type.
``energy``
    The occupation energy of the site, in eV. This line is optional, and defaults to
    0.0. Sites with the same label can have different energies, but a lookup table
    then cannot be used (see :doc:`interaction_models`).

Lines within a block can be in any order. Site numbers must be unique. Each listed
neighbour must exist, be listed only once, and not be the site itself, and neighbours
must be listed both ways: if site 1 lists site 2, site 2 must list site 1. A file that
breaks any of these rules raises a ``ValueError``.

The cell lengths are given separately, and must match the periodicity built into the
neighbour lists. Set the z cell length to 0 for a two-dimensional lattice. In this
example, ``DATA_DIR`` is the directory that holds the file; replace the path with the
location of your own sites file:

.. testcode::

    lattice = lmc.lattice_from_sites_file(
        str(DATA_DIR / "four_site_ring.dat"), cell_lengths=[4.0, 4.0, 0.0]
    )
    print(lattice.number_of_sites, lattice.dimensionality, sorted(lattice.site_labels))

.. testoutput::

    4 2 ['A', 'B']

A full real-world example, the lithium sublattice of the garnet LLZO studied in
[Morgan2017]_, is `llzo_lattice_site_list.dat
<https://github.com/bjmorgan/lattice_mc/blob/main/examples/llzo_lattice_site_list.dat>`_
in the ``examples`` directory of the repository.

Site energies
-------------

:meth:`~lattice_mc.lattice.Lattice.set_site_energies` sets the occupation energies,
in eV, of the sites in a lattice. ``by_label`` gives an energy to every site with each
label, and ``by_site`` gives energies to individual sites, keyed by site number.
Energies given by label are applied first, so an energy given by site number
overrides the energy for that site's label. Both arguments are keyword-only. Giving
neither raises a ``ValueError``. An unknown label or site number raises a
``ValueError``, and leaves the lattice unchanged.

.. testcode::

    lattice = lmc.honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
    lattice.set_site_energies(by_label={"A": 0.0, "B": 0.1}, by_site={1: 0.2})
    for number in (1, 2, 3):
        site = lattice.site_with_id(number)
        print(site.number, site.label, site.energy)

.. testoutput::

    1 A 0.2
    2 B 0.1
    3 A 0.0

``set_site_energies`` overwrites any energies read from a sites file for the labels
and sites it is given. A :class:`~lattice_mc.simulation.Simulation` reads the site
energies when it is created, so set them before creating the simulation.
