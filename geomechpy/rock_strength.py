import math
from dataclasses import dataclass


@dataclass(frozen=True)
class RockStrengthPropertiesConverter:
    """Compute rock strength properties from sonic slownesses, dynamic or static elastic properties
       Rock strength properties are:
         - Unconfined compressive Strength (UCS)
         - Tensile Strength (TSTR)
         - Friction Angle (FANG)

    Reference:
       Zhang, Yuliang, et al. "Extracting static elastic moduli of rock through elastic wave velocities." Acta Geophysica 72.2 (2024): 915-931.
    """

    @staticmethod
    def convert_yme_sta_to_ucs_plumb(yme_sta: float) -> float:
        """
        Convert static Young's modulus to UCS using Plumb Generic correlation

        Equation type: Linear law (y = a*x)

        Applicable for: generic.

        Reference: Plumb, R. A., 1994. Influence of composition and texture on the failure properties of clastic rocks. Eurock 94, Rock Mechanics in Petroleum Engineering conference, Delft, Netherlands, pp. 13-20.

        Args:
           yme_sta (float): Static Young's modulus magnitude Unit: Mpsi

        Returns:
           ucs_sta_plumb_generic: Unconfined compressive strength (UCS) from Plumb generic Young's modulus correlation Unit: psi
        """

        multiplier = 0.210306770614015
        ucs_plumb_generic = multiplier * yme_sta

        return float(ucs_plumb_generic)

    @staticmethod
    def convert_ucs_to_tstr(ucs: float, multiplier: float = 0.15) -> float:
        """
        Convert UCS to tensile strength using a constant multiploer

        Equation type: Linear law (y = a*x)

        Applicable for: Generic

        Reference: N/A
        Args:
           ucs (float): Unconfined Compressive Strength Unit: psi
           multiplier (float): constant multiplier Default set to 0.15.
        Returns:
           tstr (float): Tensile Stength Unit: psi
        """

        tstr = multiplier * ucs

        return float(tstr)

    @staticmethod
    def convert_friction_angle_lal(dtco: float) -> float:
        """
        Convert compressional slowness to friction angle using Lal correlation

        Equation type: trigonometric correlation

        Applicable for: Shale

        Reference: Lal, M., 1999. Shale stability: drilling fluid interaction and shale strength. SPE Latin American and Caribbean Petroleum Engineering Conference held in Caracas, Venezuela.

        Args:
           dtco (float): Compressional Slowness Unit: us/ft

        Returns:
          fang_lal (float): Friction angle from Lal correlation Unit: dega
        """

        fang_lal = (180 / 3.141592) * math.asin((304800 - 1000 * dtco) / (304800 + 1000 * dtco))

        return float(fang_lal)
