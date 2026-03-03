from importlib.metadata import version, PackageNotFoundError

from lattice_mc.simulation import Simulation
from lattice_mc.options import Options

try:
    __version__ = version("lattice_mc")
except PackageNotFoundError:
    __version__ = "unknown"
