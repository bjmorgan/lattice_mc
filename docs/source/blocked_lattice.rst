Blocked Lattices
================

It is possible to construct a simulation lattice where there are no possible jumps; such a lattice is classified as "blocked".
To test for this, a ``Lattice`` object can be queried using the ``is_blocked()`` method:

.. code:: python

    >>> simulation.lattice.is_blocked()
    True

If a simulation with a blocked lattice is run, a ``BlockedLatticeError`` exception is raised:

.. code::

    Traceback (most recent call last):
      …
      self.lattice.jump()
    File "/Users/bjm42/source/lattice_mc/lattice_mc/lattice.py", line 228, in jump
      raise BlockedLatticeError('No moves are possible in this lattice')
    lattice_mc.error.BlockedLatticeError: No moves are possible in this lattice
