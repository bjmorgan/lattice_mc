from unittest.mock import Mock, patch

import numpy as np
import pytest

from lattice_mc.jump import Jump
from lattice_mc.simulation import SimulationParameters
from lattice_mc.transitions import Transitions

PARAMS = SimulationParameters(temperature=298.0, rate_prefactor=1e13)


class TestTransitions:
    """Test for Transitions class"""

    @pytest.fixture(autouse=True)
    def set_up(self):
        self.jumps = [Mock(spec=Jump, relative_probability=0.25) for i in range(4)]
        self.transitions = Transitions(self.jumps, params=PARAMS)

    def test_transitions_is_initialised(self):
        assert self.transitions.jumps == self.jumps
        np.testing.assert_array_equal(self.transitions.p, np.array([0.25, 0.25, 0.25, 0.25]))

    def test_cumulative_probabilities(self):
        self.transitions.p = np.array([0.1, 0.2, 0.3, 0.4])
        np.testing.assert_allclose(self.transitions.cumulative_probabilities(), np.array([0.1, 0.3, 0.6, 1.0]))

    @patch("random.random")
    def test_random(self, mock_random):
        mock_random.return_value = 0.05
        assert self.transitions.random() is self.jumps[0]
        mock_random.return_value = 0.3
        assert self.transitions.random() is self.jumps[1]

    @patch("random.random")
    def test_time_to_jump(self, mock_random):
        self.transitions.p = np.array([0.1, 0.2, 0.3, 0.4])
        mock_random.return_value = 0.15
        assert self.transitions.time_to_jump() == pytest.approx(1.89711998489e-13, rel=1e-6, abs=0)
