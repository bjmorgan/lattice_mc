import warnings
from importlib.metadata import PackageNotFoundError, version

from lattice_mc.simulation import Simulation as Simulation
from lattice_mc.simulation import SimulationParameters as SimulationParameters

try:
    __version__ = version("lattice_mc")
except PackageNotFoundError:
    __version__ = "0.0.0.dev0"
    warnings.warn("lattice_mc package is not installed. __version__ is set to '0.0.0.dev0'.", stacklevel=2)
