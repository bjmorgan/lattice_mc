from unittest.mock import Mock, patch

import numpy as np
import pytest

from lattice_mc.atom import Atom
from lattice_mc.lattice_site import Site


class TestAtom:
    """Test for Atom class"""

    @pytest.fixture(autouse=True)
    def set_up(self):
        self.mock_site = Mock(spec=Site)
        self.mock_site.occupation = 0
        self.mock_site.is_occupied = False

    @patch("lattice_mc.atom.Atom.reset")
    def test_atom_is_initialised(self, mock_reset):
        atom = Atom(self.mock_site)
        assert atom.number == Atom.atom_number
        assert atom._site == self.mock_site
        assert self.mock_site.is_occupied is True
        assert self.mock_site.occupation == atom.number
        assert self.mock_site.atom is atom
        assert mock_reset.called

    @patch("lattice_mc.atom.Atom.reset")
    def test_atom_initialised_with_occupied_site_raises_ValueError(self, mock_reset):
        self.mock_site.occupation = 3
        self.mock_site.is_occupied = True
        with pytest.raises(ValueError):
            Atom(self.mock_site)

    def test_reset(self):
        with patch("lattice_mc.atom.Atom.reset"):
            self.mock_site.number = 5
            atom = Atom(self.mock_site)  # does not call Atom.reset() in __init__()
        atom.reset()
        assert atom.number_of_hops == 0
        np.testing.assert_array_equal(atom.dr, np.array([0.0, 0.0, 0.0]))
        assert atom.summed_dr2 == 0.0
        assert atom.sites_visited == [self.mock_site.number]

    @patch("lattice_mc.atom.Atom.reset")
    def test_dr_squared(self, mock_reset):
        atom = Atom(self.mock_site)
        atom.dr = np.array([1.0, 2.0, 3.0])
        assert atom.dr_squared() == 14.0

    @patch("lattice_mc.atom.Atom.reset")
    def test_site(self, mock_reset):
        atom = Atom(self.mock_site)
        assert atom.site == self.mock_site

    @patch("lattice_mc.atom.Atom.reset")
    def test_site_setter(self, mock_reset):
        atom = Atom(self.mock_site)
        new_mock_site = Mock(spec=Site)
        atom.site = new_mock_site
        assert atom._site == new_mock_site
