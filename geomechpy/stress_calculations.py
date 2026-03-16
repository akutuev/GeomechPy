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
       Zhang, Jon Jincai. Applied petroleum geomechanics. Vol. 1. Cambridge: Gulf Professional Publishing, 2019 Chapter 6.
    """

    @staticmethod
    def calculate_poroelastic_horizontal_stresses(overburden_stress: float, pore_pressure: float, poisson_ratio: float, youngs_modulus: float, biot_coefficient: float = 1.0, EX: float = 0.0001, EY: float = 0.009) -> HorizontalStresses:
        """Calculates horizontal stress using Poroelastic horizontal stress equation.

        Reference: Thiercelin, Marc Jean, and Richard A. Plumb. "Core-based prediction of lithologic stress contrasts in East Texas formations." SPE Formation Evaluation 9.04 (1994): 251-258.

        Args:
            overburden_stress (float): Array of overburden stress values.
            pore_pressure (float): Array of pore pressure values pore_pressure in pressure unit  Unit: Pressure Unit [psi].
            poisson_ratio (float): Static Poisson's ratio. Unit: unitless.
            youngs_modulus (float): Static Young's modulus Unit:Elastic Modulus Unit [Mpsi].
            biot_coefficient (float): Biot's coefficient. Defaults to 1.0
            EX (float): Tectonic strain term Unit: unitless Defaults to 0.0001.
            EY (float): Tectonic strain term Unit: unitless Defaults to 0.009.

        Returns:
            shmin (float): Minimum horizontal stress magnitude Unit [psi].
            shmax (float): MAximum horizontal stress magnitude Unit [psi].
            q_factor (float): Stress Regime Indicator  Unit [unitless].
            shmax_shmin_ratio (float): Stress ratio Unit [unitless].

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
            shmin (float): Minimum horizontal stress magnitude Unit: [psi]
            shmax_multiplier (float): A unitless multiplier representing the stress anisotropy Defaults to 1.1

        Returns:
            shmax (float): Maximum horizontal stress magnitude Unit [psi].
        """
        shmax = shmin * shmax_multiplier

        return shmax

    @staticmethod
    def calculate_stress_regime_q_factor(sigv: float, shmax: float, shmin: float) -> float:
        """Calculates q factor represnting the str.

        Reference: Prats, M., Effect of Burial History on the Subsurface Horizontal Stresses of Formations Having Different Material Properties, SPE 9017-PA, 1981.

        Args:
            sigv (float): Vertical stress magnitude Unit [psi].
            shmin (float): Minimum horizontal stress magnitude Unit [psi].
            shmax (float): Maximum horizontal stress magnitude Unit [psi].

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
        """Calculates the ratio between maximum and minimum horizontal stress magnitudes.

        Args:
            shmin (float): Minimum horizontal stress magnitude Unit [psi].
            shmax (float): Maximum horizontal stress magnitude Unit [psi].

        Returns:
            shmax_shmin_ratio (float): stress ratio between maximum and minimum horizontal stress magnitudes Unit [unitless].
            Value needs to be equal or bigger than 1
        """
        shmax_shmin_ratio = shmax / shmin
        return shmax_shmin_ratio
