import unittest

from lattice_mc.constants import k_boltzmann
from lattice_mc.simulation import SimulationParameters


class ConstantsTestCase(unittest.TestCase):
    """Tests for physical constants."""

    def test_k_boltzmann_is_approximately_correct(self):
        """Validate k_boltzmann is approximately 8.617e-5 eV/K."""
        self.assertAlmostEqual(k_boltzmann, 8.617e-5, places=8)


class SimulationParametersTestCase(unittest.TestCase):
    """Tests for SimulationParameters."""

    def test_kT_at_298K(self):
        """Validate kT at 298 K is approximately 0.02569 eV."""
        params = SimulationParameters(temperature=298.0, rate_prefactor=1e13)
        self.assertAlmostEqual(params.kT, 0.02569, places=4)

    def test_frozen(self):
        """SimulationParameters should be immutable."""
        params = SimulationParameters(temperature=298.0, rate_prefactor=1e13)
        with self.assertRaises(AttributeError):
            params.temperature = 500.0


if __name__ == "__main__":
    unittest.main()
