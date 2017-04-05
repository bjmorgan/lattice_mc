from scipy.constants import physical_constants

k_boltzmann = physical_constants[ 'Boltzmann constant in eV/K' ][0]
temperature = 298.0 # arbitrary temperature if we always use energies in scaled units of kT.
rate_prefactor = 1e13 # standard vibrational frequency

kT = k_boltzmann * temperature

