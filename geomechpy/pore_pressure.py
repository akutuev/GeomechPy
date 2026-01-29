PorePressureimport math
from dataclasses import dataclass


@dataclass(frozen=True)
class PorePressure:
    """Calculation of the Overburden Stress

    """

    PorePressure: float
    PorePressur_MW: float


class PorePressureCalculation:
    """Computation of the Pore Pressure using gradient based methods 

    Reference:
       Zhang, Jon Jincai. Applied petroleum geomechanics. Vol. 1. Cambridge: Gulf Professional Publishing, 2019. Chapter 6.1

    """
    @staticmethod
    def pore_pressure_onshore(tvd, Formation_pore_pressure_gradient = '0.47', Air_gap = '0') -> PorePressure:
        """
        Calculates pore pressure from tvd and pore pressure gradient.

        Args:
            tvd (array-like): Array of depths. Unit: Depth Unit [ft]
            litho_gradient float): Overburden stress depth gradient. Unit: Depth Gradient Unit [psi/ft]

        Returns:
            array-like: Array of overburden stress values. Unit: Depth Unit [psi]
        """
        air_gradient = 0.0004
        Air_gap = float(Air_gap)
        air_pressure = air_gradient*Air_gap
        if tvd >= 0 and tvd < Air_gap:
            PorePressure = air_gradient*tvd
        else:
            pp_gradient = float(pp_gradient)
            PorePressure = air_pressure + Formation_pore_pressure_gradient*(tvd - Air_gap)
            
        PorePressure_MW = PorePressure/(19.25*tvd)


        return PorePressure(PorePressure, PorePressure_MW)

    @staticmethod
    def pore_pressure_offshore(tvd, Formation_pore_pressure_gradient = '0.47', Air_gap = '0', Water_depth = '0', Sea_water_pressure_gradient = '0.47') -> PorePressure:
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
        air_pressure = air_gradient*Air_gap

        Sea_water_pressure_gradient= float(Sea_water_pressure_gradient)
        Water_depth = float(Water_depth)
        water_pressure = Sea_water_pressure_gradient*Water_depth

        Formation_pore_pressure_gradient = float(Formation_pore_pressure_gradient)

        if tvd >= 0 and tvd < Air_gap:
            air_gradient = 0.0004
            PorePressure = air_gradient*tvd
        elif tvd >= Air_gap and tvd <= (Air_gap + Water_depth):      
            PorePressure = air_pressure + Sea_water_pressure_gradient * (tvd - Air_gap)
        else:   
            PorePressure =  air_pressure + water_pressure + Formation_pore_pressure_gradient*(tvd - Water_depth - Air_gap)
        
        PorePressure_MW = PorePressure/(19.25*tvd)
            
        return PorePressure(PorePressure, PorePressure_MW)