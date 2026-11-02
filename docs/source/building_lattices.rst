Building a lattice
==================

A simulation runs on a :class:`~lattice_mc.lattice.Lattice`: a set of sites, each
with a position, a label, and a list of neighbouring sites that atoms can jump to.
Lattices are periodic, inside an orthorhombic simulation cell. Simple lattices can be
built programmatically, and lattices of any geometry can be read from a sites file.

Simple lattices
---------------

The functions in :mod:`lattice_mc.init_lattice` build square, honeycomb and cubic
lattices from a number of repeat units along each axis and a site spacing.

A square lattice of 4 x 4 sites:

.. testcode::

    from lattice_mc import init_lattice

    lattice = init_lattice.square_lattice(4, 4, 1.0)
    print(lattice.number_of_sites, lattice.dimensionality, sorted(lattice.site_labels))

.. testoutput::

    16 2 ['L']

A cubic lattice of 2 x 2 x 2 sites:

.. testcode::

    lattice = init_lattice.cubic_lattice(2, 2, 2, 1.0)
    print(lattice.number_of_sites, lattice.dimensionality)

.. testoutput::

    8 3

Every site in these lattices is labelled ``'L'``. A honeycomb lattice can instead
have alternating ``'A'`` and ``'B'`` sites, so that every ``'A'`` site has only
``'B'`` neighbours. This is useful for giving the two sublattices different energies
(see :doc:`interaction_models`).

.. testcode::

    lattice = init_lattice.honeycomb_lattice(2, 2, 1.0, alternating_sites=True)
    print(lattice.number_of_sites, lattice.dimensionality, sorted(lattice.site_labels))

.. testoutput::

    16 2 ['A', 'B']

Two-dimensional lattices
------------------------

Square and honeycomb lattices lie in the xy plane, and have a z cell length of 0.
:attr:`~lattice_mc.lattice.Lattice.dimensionality` is the number of non-zero cell
lengths, and is used as :math:`d` in the Einstein relations for the diffusion
coefficients (see :doc:`introduction`).

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
    0.0.

Lines within a block can be in any order. The cell lengths are given separately. Set
the z cell length to 0 for a two-dimensional lattice:

.. testcode::

    lattice = init_lattice.lattice_from_sites_file(
        str(DATA_DIR / "four_site_ring.dat"), cell_lengths=[4.0, 4.0, 0.0]
    )
    print(lattice.number_of_sites, lattice.dimensionality, sorted(lattice.site_labels))

.. testoutput::

    4 2 ['A', 'B']

:meth:`Simulation.define_lattice_from_file
<lattice_mc.simulation.Simulation.define_lattice_from_file>` does the same, and
assigns the lattice to a simulation.

A full real-world example, the lithium sublattice of the garnet LLZO studied in
[Morgan2017]_, is `llzo_lattice_site_list.dat
<https://github.com/bjmorgan/lattice_mc/blob/master/examples/llzo_lattice_site_list.dat>`_
in the ``examples`` directory of the repository.
