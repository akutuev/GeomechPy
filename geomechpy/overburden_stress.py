import math
from dataclasses import dataclass


@dataclass(frozen=True)
class OverburdenStress:
    """Calculation of the Overburden Stress"""

    overburden_stress: float
    overburden_stress_mw: float


class OverburdenStressCalculation:
    """Computation of the Overburden Stress using various methods based on gradient, density

    Reference:
       Zhang, Jon Jincai. Applied petroleum geomechanics. Vol. 1. Cambridge: Gulf Professional Publishing, 2019. Chapter 6.1

    """

    @staticmethod
    def overburden_stress_onshore(tvd, litho_gradient='1.05', Air_gap='0') -> OverburdenStress:
        """
        Calculates overburden stress (vertical stress) from tvd and lithogradient.

        Args:
            tvd (array-like): Array of depths. Unit: Depth Unit [ft]
            litho_gradient float): Overburden stress depth gradient. Unit: Depth Gradient Unit [psi/ft]

        Returns:
            array-like: Array of overburden stress values. Unit: Depth Unit [psi]
        """

        air_gradient = 0.0004
        Air_gap = float(Air_gap)
        air_pressure = air_gradient * Air_gap
        if tvd >= 0 and tvd < Air_gap:
            overburden_stress = air_gradient * tvd
        else:
            litho_gradient = float(litho_gradient)
            overburden_stress = air_pressure + litho_gradient * (tvd - Air_gap)

        overburden_stress_mw = overburden_stress / (19.25 * tvd)

        return OverburdenStress(overburden_stress, overburden_stress_mw)

    @staticmethod
    def overburden_stress_offshore(tvd, litho_gradient='1.05', Air_gap='0', Water_depth='0', Sea_water_pressure_gradient='0.47') -> OverburdenStress:
        """
        Calculates overburden stress (vertical stress) from tvd and lithogradient.

        Args:
            tvd (array-like): Array of depths. Unit: Depth Unit [ft]
            litho_gradient float): Overburden stress depth gradient. Unit: Depth Gradient Unit [psi/ft]

        Returns:
            array-like: Array of overburden stress values. Unit: Depth Unit [psi]
        """

        air_gradient = 0.0004
        Air_gap = float(Air_gap)
        air_pressure = air_gradient * Air_gap

        Sea_water_pressure_gradient = float(Sea_water_pressure_gradient)
        Water_depth = float(Water_depth)
        water_pressure = Sea_water_pressure_gradient * Water_depth

        if tvd >= 0 and tvd < Air_gap:
            overburden_stress = air_gradient * tvd
        elif tvd >= Air_gap and tvd <= (Air_gap + Water_depth):
            overburden_stress = air_pressure + Sea_water_pressure_gradient * (tvd - Air_gap)
        else:
            overburden_stress = air_pressure + water_pressure + float(litho_gradient) * (tvd - Water_depth - Air_gap)

        overburden_stress_mw = overburden_stress / (19.25 * tvd)

        return OverburdenStress(overburden_stress, overburden_stress_mw)
