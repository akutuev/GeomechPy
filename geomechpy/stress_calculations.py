import math
from dataclasses import dataclass


@dataclass(frozen=True)
class HorizontalStresses:
    """Calculation of stress components and properties.

    Attributes:
        Shmin(float): Minimum horizontal stress magnitude. Unit: Pressure
        SHmax(float): Maximum horizontal stress magnitude. Unit: Pressure

    """

    Shmin: float
    SHmax: float


class HorizontalStressesCalculation:
    """Calculation of stress components and properties.

    Reference:
       Plumb hlkh
    """

    @staticmethod
    def calculate_poroelastic_horizontal_stresses(overburden_stress: float, pore_pressure: float, poisson_ratio: float, youngs_modulus: float, biot_coefficient: float = 1, EX: float = 0.0001, EY: float = 0.009):
        """
        Calculates horizontal stress using Poroelastic horizontal stress equation

        Args:
            overburden_stress (array-like): Array of overburden stress values.
            pore_pressure (array-like): Array of pore pressure values.
            biot_coefficient (float): Biot's coefficient.
            poisson_ratio (float): Poisson's ratio.

        Returns:
            array-like: Array of horizontal stress values.
        """
        biot_coefficient = int(biot_coefficient)
        EX = float(EX) / 1e-3
        EY = float(EY) / 1e-3
        A = poisson_ratio / (1 - poisson_ratio)
        B = youngs_modulus / (1 - poisson_ratio * poisson_ratio)
        C = (poisson_ratio * youngs_modulus) / (1 - poisson_ratio * poisson_ratio)

        shmin_phs = A * overburden_stress + (1 - A) * biot_coefficient * pore_pressure + B * EX + C * EY
        shmax_phs = A * overburden_stress + (1 - A) * biot_coefficient * pore_pressure + B * EY + C * EX

        return shmin_phs, shmax_phs

    @staticmethod
    def calculate_shmax_multiplier(shmin: float, shmax_multiplier: float = 1.1) -> float:
        """
        Calculates maximum horizontal stress from minimum horizontal strss using multiploer

        Args:
            shmin (float): TODO
            shmax_multiplier (float): TODO Defaults to 1.1

        Returns:
            shmax (float): TODO
        """
        shmax_multiplier = float(shmax_multiplier)
        shmax = shmin * shmax_multiplier

        return shmax

    @staticmethod
    def calculate_stress_regime_q_factor(sigv: float, shmin: float, shmax: float) -> float:
        """
        Calculates maximum horizontal stress from minimum horizontal strss using multiploer

        Args:
            sigv (float): TODO
            shmin (float): TODO
            shmax (float): TODO

        Returns:
            shmax (float): TODO
        """
        if sigv > shmax and shmax >= shmin:
            q_factor = (shmax - shmin) / (sigv - shmin)
        elif shmin < sigv and sigv <= shmax:
            q_factor = 2 - (sigv - shmin) / (shmax - shmin)
        elif sigv <= shmin and shmin < shmax:
            q_factor = 2 + (shmin - sigv) / (shmax - sigv)
        else:
            q_factor = 4

        return q_factor

    @staticmethod
    def calculate_horizontal_stress_ratio(shmax: float, shmin: float) -> float:
        """
        Calculates maximum horizontal stress from minimum horizontal strss using multiploer

        Args:
            shmax (float): TODO
            shmin (float): TODO

        Returns:
            shmax_shmin_ratio (float): TODO
        """
        shmax_shmin_ratio = shmax / shmin
        return shmax_shmin_ratio
