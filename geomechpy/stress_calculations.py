import math
from dataclasses import dataclass


@dataclass(frozen=True)
class HorizontalStresses:
    """Calculation of stress components and properties.

    Attributes:
        Shmin(float): Minimum horizontal stress magnitude. Unit: Pressure
        SHmax(float): Maximum horizontal stress magnitude. Unit: Pressure

    """

    shmin: float
    shmax: float
    q_factor: float
    shmax_shmin_ratio: float


class HorizontalStressesCalculation:
    """Calculation of stress components and properties.

    Reference:
       Plumb hlkh
    """

    @staticmethod
    def calculate_poroelastic_horizontal_stresses(overburden_stress: float, pore_pressure: float, poisson_ratio: float, youngs_modulus: float, biot_coefficient: float = 1.0, EX: float = 0.0001, EY: float = 0.009) -> HorizontalStresses:
        """Calculates horizontal stress using Poroelastic horizontal stress equation.

        Args:
            overburden_stress (float): Array of overburden stress values.
            pore_pressure (float): Array of pore pressure values.
            poisson_ratio (float): Poisson's ratio.
            youngs_modulus (float): TODO
            biot_coefficient (float): Biot's coefficient. Defaults to 1.0
            EX (float): Defaults to 0.0001
            EY (float): Defaults to 0.009

        Returns:
            HorizontalStresses: Dataclass containing computed horizontal stresses properties. See `HorizontalStresses` for details
            Output unit: consistent with input unit
        """
        EX = float(EX) / 1e-3
        EY = float(EY) / 1e-3
        A = poisson_ratio / (1 - poisson_ratio)
        B = youngs_modulus / (1 - poisson_ratio * poisson_ratio)
        C = (poisson_ratio * youngs_modulus) / (1 - poisson_ratio * poisson_ratio)

        shmin = A * overburden_stress + (1 - A) * biot_coefficient * pore_pressure + B * EX + C * EY
        shmax = A * overburden_stress + (1 - A) * biot_coefficient * pore_pressure + B * EY + C * EX

        q_factor = HorizontalStressesCalculation.calculate_stress_regime_q_factor(0.0, shmax, shmin)
        shmax_shmin_ratio = HorizontalStressesCalculation.calculate_horizontal_stress_ratio(shmax, shmin)

        return HorizontalStresses(shmin, shmax, q_factor, shmax_shmin_ratio)

    @staticmethod
    def calculate_shmax_multiplier(shmin: float, shmax_multiplier: float = 1.1) -> float:
        """Calculates maximum horizontal stress from minimum horizontal strss using multiplier.

        Args:
            shmin (float): TODO
            shmax_multiplier (float): TODO Defaults to 1.1

        Returns:
            shmax (float): TODO
        """
        shmax = shmin * shmax_multiplier

        return shmax

    @staticmethod
    def calculate_stress_regime_q_factor(sigv: float, shmax: float, shmin: float) -> float:
        """Calculates maximum horizontal stress from minimum horizontal strss using multiplier.

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
        """Calculates maximum horizontal stress from minimum horizontal strss using multiplier.

        Args:
            shmax (float): TODO
            shmin (float): TODO

        Returns:
            shmax_shmin_ratio (float): TODO
        """
        shmax_shmin_ratio = shmax / shmin
        return shmax_shmin_ratio
