import math


class OverburdenStressCalculation:
    """Computation of the Overburden Stress using various methods based on gradient, density.

    Reference:
       Zhang, Jon Jincai. Applied petroleum geomechanics. Vol. 1. Cambridge: Gulf Professional Publishing, 2019. Chapter 6.1
    """

    @staticmethod
    def calculate_overburden_stress_onshore(tvd: float, lithostatic_gradient: float = 1.05, air_gap: float = 0.0) -> float:
        """Calculates overburden stress (vertical stress) from tvd and lithogradient.

        Args:
            tvd (float): TODO Depth value. Unit: Depth Unit [ft]
            lithostatic_gradient (float): Overburden stress depth gradient. Unit: Depth Gradient Unit [psi/ft]. Defaults to 1.05
            air_gap (float): TODO. Defaults to 0.0

        Returns:
            Overburden stress onshore value.
            Unit: Depth Unit [psi]
        """
        air_gradient = 0.0004
        air_pressure = air_gradient * air_gap

        if tvd >= 0 and tvd < air_gap:
            overburden_stress = air_gradient * tvd
        else:
            overburden_stress = air_pressure + lithostatic_gradient * (tvd - air_gap)

        return overburden_stress

    @staticmethod
    def calculate_overburden_stress_offshore(tvd: float, lithostatic_gradient: float = 1.05, air_gap: float = 0.0, water_depth: float = 0.0, sea_water_pressure_gradient: float = 0.47) -> float:
        """Calculates overburden stress (vertical stress) from tvd and lithogradient. {TODO: Difference in description from previous?}

        Args:
            tvd (float): TODO Depth value. Unit: Depth Unit [ft]
            lithostatic_gradient (float): Overburden stress depth gradient. Unit: Depth Gradient Unit [psi/ft]. Defaults to 1.05
            air_gap (float): TODO. Defaults to 0.0
            water_depth (float):  TODO. Defaults to 0.0
            sea_water_pressure_gradient (float): TODO. Defaults to 0.47

        Returns:
            Overburden stress offshore value.
            Unit: Depth Unit [psi]
        """
        air_gradient = 0.0004
        air_pressure = air_gradient * air_gap

        water_pressure = sea_water_pressure_gradient * water_depth

        if tvd >= 0 and tvd < air_gap:
            overburden_stress = air_gradient * tvd
        elif tvd >= air_gap and tvd <= (air_gap + water_depth):
            overburden_stress = air_pressure + sea_water_pressure_gradient * (tvd - air_gap)
        else:
            overburden_stress = air_pressure + water_pressure + lithostatic_gradient * (tvd - water_depth - air_gap)

        return overburden_stress
