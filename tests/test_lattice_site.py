import numpy as np

from lattice_mc.lattice_site import Site


def test_site_is_initialised():
    coordinates = np.array([1.0, 2.0, 3.0])
    site = Site(7, coordinates, [2, 3], 0.4, "A")
    assert site.number == 7
    assert site.r is coordinates
    assert site.neighbours == [2, 3]
    assert site.energy == 0.4
    assert site.label == "A"
