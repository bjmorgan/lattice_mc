import pytest

from lattice_mc.constants import k_boltzmann
from lattice_mc.simulation import SimulationParameters


class TestConstants:
    """Tests for physical constants."""

    def test_k_boltzmann_is_approximately_correct(self):
        """Validate k_boltzmann is approximately 8.617e-5 eV/K."""
        assert k_boltzmann == pytest.approx(8.617e-5, abs=5e-9)


class TestSimulationParameters:
    """Tests for SimulationParameters."""

    def test_kT_at_298K(self):
        """Validate kT at 298 K is approximately 0.02569 eV."""
        params = SimulationParameters(temperature=298.0, rate_prefactor=1e13)
        assert params.kT == pytest.approx(0.02569, abs=5e-5)

    def test_frozen(self):
        """SimulationParameters should be immutable."""
        params = SimulationParameters(temperature=298.0, rate_prefactor=1e13)
        with pytest.raises(AttributeError):
            params.temperature = 500.0

    def test_negative_temperature_raises_ValueError(self):
        with pytest.raises(ValueError):
            SimulationParameters(temperature=-100.0, rate_prefactor=1e13)

    def test_zero_temperature_raises_ValueError(self):
        with pytest.raises(ValueError):
            SimulationParameters(temperature=0.0, rate_prefactor=1e13)

    def test_negative_rate_prefactor_raises_ValueError(self):
        with pytest.raises(ValueError):
            SimulationParameters(temperature=298.0, rate_prefactor=-1e13)

    def test_zero_rate_prefactor_raises_ValueError(self):
        with pytest.raises(ValueError):
            SimulationParameters(temperature=298.0, rate_prefactor=0.0)
