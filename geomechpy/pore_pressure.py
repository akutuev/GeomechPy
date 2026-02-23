from dataclasses import dataclass


class PorePressureCalculation:
    """Computation of the Pore Pressure using gradient based methods.

    Reference:
       Zhang, Jon Jincai. Applied petroleum geomechanics. Vol. 1. Cambridge: Gulf Professional Publishing, 2019. Chapter 6.1
    """

    @staticmethod
    def calculate_pore_pressure_onshore(tvd: float, formation_pore_pressure_gradient: float = 0.47, air_gap: float = 0.0) -> float:
        """Calculates pore pressure from tvd and pore pressure gradient.

        Args:
            tvd (float): TODO
            formation_pore_pressure_gradient (float): TODO. Defaults to 0.47
            air_gap: TODO Defaults to 0.0

        Returns:
            Overburden stress onshore value.
            Unit: Depth Unit [psi]
        """
        air_gradient = 0.0004
        air_pressure = air_gradient * air_gap
        if tvd >= 0 and tvd < air_gap:
            pore_pressure = air_gradient * tvd
        else:
            pore_pressure = air_pressure + formation_pore_pressure_gradient * (tvd - air_gap)

        return pore_pressure

    @staticmethod
    def calculate_pore_pressure_offshore(tvd: float, formation_pore_pressure_gradient: float = 0.47, air_gap: float = 0.0, water_depth: float = 0.0, sea_water_pressure_gradient: float = 0.47) -> float:
        """Calculates overburden stress (vertical stress) from tvd and lithogradient.

        Args:
            tvd (float): TODO
            formation_pore_pressure_gradient (float): TODO. Defaults to 0.47
            air_gap: TODO. Defaults to 0.0
            water_depth: TODO Defaults to 0.0
            sea_water_pressure_gradient: TODO Defaults to 0.47

        Returns:
            Overburden stress offshore value.
            Unit: Depth Unit [psi]
        """
        air_gradient = 0.0004
        air_pressure = air_gradient * air_gap

        water_pressure = sea_water_pressure_gradient * water_depth

        if tvd >= 0 and tvd < air_gap:
            pore_pressure = air_gradient * tvd
        elif tvd >= air_gap and tvd <= (air_gap + water_depth):
            pore_pressure = air_pressure + sea_water_pressure_gradient * (tvd - air_gap)
        else:
            pore_pressure = air_pressure + water_pressure + formation_pore_pressure_gradient * (tvd - water_depth - air_gap)

        return pore_pressure
