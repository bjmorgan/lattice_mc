"""The quantities reported by a simulation."""

from __future__ import annotations

import dataclasses
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from lattice_mc._state import SimulationState
    from lattice_mc.lattice import Lattice


def _rebuild_results(fields: dict[str, Any]) -> Results:
    """
    Rebuild a `Results` from its fields, as plain values, when unpickled or copied.

    Args:
        fields: The field values, with the site occupations as a plain dict.

    Returns:
        The rebuilt results.
    """
    return Results(**{**fields, "average_site_occupations": MappingProxyType(fields["average_site_occupations"])})


@dataclass(frozen=True)
class Results:
    """
    Everything a simulation reports, accumulated since the last reset.

    Lengths are in the length unit of the site coordinates and cell lengths, and times
    are in seconds.

    Attributes:
        tracer_diffusion_coefficient: Tracer diffusion coefficient, in (length unit)^2 s^-1.
            Defined as the sum of the squared displacements of the atoms divided by
            `2 d N t`, for `d` dimensions, `N` atoms and elapsed time `t`.
        collective_diffusion_coefficient: Collective (jump) diffusion coefficient per atom,
            in (length unit)^2 s^-1. Defined as the squared length of the summed displacement
            of all atoms divided by `2 d N t`.
        reduced_ionic_conductivity: Reduced ionic conductivity, in (length unit)^2 s^-1.
            Defined as the collective diffusion coefficient multiplied by the number of atoms
            and divided by the number of sites.
        tracer_correlation: Tracer correlation factor, dimensionless. Defined as the sum of
            the squared displacements of the atoms divided by the sum of the squared jump lengths.
        collective_correlation: Collective correlation factor, dimensionless.
            Defined as the squared length of the summed displacement of all atoms divided by
            the sum of the squared jump lengths.
        average_site_occupations: Time-averaged number of atoms on the sites with each label,
            keyed by site label, in sorted label order. The mapping is read-only.
        time: Simulation time elapsed since the last reset, in s.
        number_of_jumps: Number of jumps made since the last reset.
    """

    tracer_diffusion_coefficient: float
    collective_diffusion_coefficient: float
    reduced_ionic_conductivity: float
    tracer_correlation: float
    collective_correlation: float
    average_site_occupations: Mapping[str, float]
    time: float
    number_of_jumps: int

    def __reduce__(self) -> tuple[Callable[[dict[str, Any]], Results], tuple[dict[str, Any]]]:
        """
        Support pickling and copying, which cannot handle the read-only mapping directly.

        Returns:
            The function and arguments that rebuild these results.
        """
        fields = {field.name: getattr(self, field.name) for field in dataclasses.fields(self)}
        fields["average_site_occupations"] = dict(self.average_site_occupations)
        return _rebuild_results, (fields,)

    @classmethod
    def from_state(cls, state: SimulationState, lattice: Lattice) -> Results | None:
        """
        Compute the results from the state of a simulation.

        A state with at least one jump has at least one atom, and a positive elapsed time,
        so no division by zero arises and the case of zero atoms needs no handling.

        Args:
            state: The state of the simulation.
            lattice: The lattice the simulation ran on.

        Returns:
            The results, or None if no jump has been made since the last reset.
        """
        if state.number_of_jumps == 0:
            return None
        n_atoms = len(state.atom_site)
        denominator = 2 * lattice.dimensionality * n_atoms * state.time
        sum_squared_displacement = float(np.sum(state.displacement**2))
        squared_summed_displacement = float(np.sum(state.displacement.sum(axis=0) ** 2))
        collective_diffusion_coefficient = squared_summed_displacement / denominator

        site_labels = np.array([site.label for site in lattice.sites])
        labels, label_index = np.unique(site_labels, return_inverse=True)
        occupation_time = np.bincount(label_index, weights=state.site_occupation_time, minlength=len(labels))
        return cls(
            tracer_diffusion_coefficient=sum_squared_displacement / denominator,
            collective_diffusion_coefficient=collective_diffusion_coefficient,
            reduced_ionic_conductivity=collective_diffusion_coefficient * n_atoms / lattice.number_of_sites,
            tracer_correlation=sum_squared_displacement / state.summed_jump_length_squared,
            collective_correlation=squared_summed_displacement / state.summed_jump_length_squared,
            average_site_occupations=MappingProxyType(
                dict(zip(labels.tolist(), (occupation_time / state.time).tolist(), strict=True))
            ),
            time=state.time,
            number_of_jumps=state.number_of_jumps,
        )
