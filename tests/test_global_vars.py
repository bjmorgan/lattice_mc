import unittest
from lattice_mc.global_vars import kT


class GlobalVarsTestCase(unittest.TestCase):
    """Tests for global_vars constants."""

    def test_kT_at_298K_is_approximately_correct(self):
        """Validate kT at 298 K is approximately 0.02569 eV."""
        self.assertAlmostEqual(kT, 0.02569, places=4)


if __name__ == "__main__":
    unittest.main()
