import math
from dataclasses import dataclass


@dataclass(frozen=True)
class ElasticProperties:
    """Elastic properties model per depth. Can be populated manually or via converter

    Attributes:
        Bulk_modulus (float): The bulk modulus (K or B or k) of a substance is a measure of the resistance of a substance to bulk compression.
        Youngs_modulus (float): Young modulus is a mechanical property of solid materials that measures the tensile or compressive stiffness when the force is applied lengthwise.
        Lame_parameter (float): Lame parameters are two material-dependent quantities denoted by λ and μ that arise in strain-stress relationships.
        Shear_modulus (float): Shear modulus is a measure of the elastic shear stiffness of a material and is defined as the ratio of shear stress to the shear strain:
        Poissons_ratio (float): Poisson's ratio is a measure of the Poisson effect, the deformation of a material in directions perpendicular to the specific direction of loading.
        P_wave_modulus (float): P-wave modulus is one of the elastic moduli available to describe isotropic homogeneous materials.
    """

    bulk_modulus: float
    youngs_modulus: float
    lames_first_parameter: float
    shear_modulus: float
    poissons_ratio: float
    p_wave_modulus: float


class ElasticPropertiesConverter:
    """Convert any pair of two elastic properties into all other types of elastic property notations

    Reference:
        https://en.wikipedia.org/wiki/Elastic_modulus

    """

    @staticmethod
    def from_bulk_and_youngs_modulus(bulk_modulus: float, youngs_modulus: float) -> ElasticProperties:
        """
        Convert Bulk and Youngs modulus to  other elastic property types

         Args:
             bulk_modulus (float): Bulk modulus magnitude TODO: [Unitless]
             youngs_modulus (float): Young modulus magnitude TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """

        lames_first_parameter = 3 * bulk_modulus * (3 * bulk_modulus - youngs_modulus) / (9 * bulk_modulus - youngs_modulus)
        shear_modulus = 3 * bulk_modulus * youngs_modulus / (9 * bulk_modulus - youngs_modulus)
        poissons_ratio = (3 * bulk_modulus - youngs_modulus) / (6 * bulk_modulus)
        p_wave_modulus = 3 * bulk_modulus * (3 * bulk_modulus + youngs_modulus) / (9 * bulk_modulus - youngs_modulus)

        return ElasticProperties(bulk_modulus, youngs_modulus, lames_first_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_bulk_and_lames_modulus(bulk_modulus: float, lames_first_parameter: float) -> ElasticProperties:
        """
        Convert xxx and xxx to  other elastic property types
        Pressure unit between inputs needs to be consistent

         Args:
            xxx (float): Bulk modulus magnitude Unit: Pressure Unit
            xxx (float): Young modulus magnitude Unit: Pressure Unit

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """

        youngs_modulus = 9 * bulk_modulus * (bulk_modulus - lames_first_parameter) / (3 * bulk_modulus - lames_first_parameter)
        shear_modulus = 3 * (bulk_modulus - lames_first_parameter) / 2
        poissons_ratio = lames_first_parameter / (3 * bulk_modulus - lames_first_parameter)
        p_wave_modulus = 3 * bulk_modulus - 2 * lames_first_parameter

        return ElasticProperties(bulk_modulus, youngs_modulus, lames_first_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_bulk_and_shear_modulus(bulk_modulus: float, shear_modulus: float) -> ElasticProperties:
        """
        Convert xxx and xxx to  other elastic property types

         Args:
             xxx (float): Bulk modulus magnitude TODO: [Unitless]
            xxx (float): Young modulus magnitude TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """
        youngs_modulus = 9 * bulk_modulus * shear_modulus / (3 * bulk_modulus + shear_modulus)
        lames_first_parameter = bulk_modulus - 2 * shear_modulus / 3
        poissons_ratio = (3 * bulk_modulus - 2 * shear_modulus) / (2 * (3 * bulk_modulus + shear_modulus))
        p_wave_modulus = bulk_modulus + 4 * shear_modulus / 3

        return ElasticProperties(bulk_modulus, youngs_modulus, lames_first_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_bulk_and_poissons_ratio_modulus(bulk_modulus: float, poissons_ratio: float) -> ElasticProperties:
        """
        Convert xxx and xxx to  other elastic property types

         Args:
             xxx (float): Bulk modulus magnitude TODO: [Unitless]
            xxx (float): Young modulus magnitude TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """
        youngs_modulus = 3 * bulk_modulus * (1 - 2 * poissons_ratio)
        lames_first_parameter = 3 * bulk_modulus * poissons_ratio / (1 + poissons_ratio)
        shear_modulus = 3 * bulk_modulus * (1 - 2 * poissons_ratio) / (2 * (1 + poissons_ratio))
        p_wave_modulus = 3 * bulk_modulus * (1 - poissons_ratio) / (1 + poissons_ratio)

        return ElasticProperties(bulk_modulus, youngs_modulus, lames_first_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_bulk_and_p_wave_modulus(bulk_modulus: float, p_wave_modulus: float) -> ElasticProperties:
        """
        Convert xxx and xxx to  other elastic property types

         Args:
             xxx (float): Bulk modulus magnitude TODO: [Unitless]
            xxx (float): Young modulus magnitude TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """
        youngs_modulus = 9 * bulk_modulus * (p_wave_modulus - bulk_modulus) / (3 * bulk_modulus + p_wave_modulus)
        lames_first_parameter = (3 * bulk_modulus - p_wave_modulus) / 2
        shear_modulus = (3 * (p_wave_modulus - bulk_modulus)) / 4
        poissons_ratio = (3 * bulk_modulus - p_wave_modulus) / (3 * bulk_modulus + p_wave_modulus)

        return ElasticProperties(bulk_modulus, youngs_modulus, lames_first_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_youngs_and_lames_modulus(youngs_modulus: float, lames_first_parameter: float) -> ElasticProperties:
        """
        Convert xxx and xxx to  other elastic property types

         Args:
             xxx (float): Bulk modulus magnitude TODO: [Unitless]
            xxx (float): Young modulus magnitude TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """
        R = math.sqrt(youngs_modulus**2 + 9 * lames_first_parameter**2 + 2 * youngs_modulus * lames_first_parameter)

        bulk_modulus = (youngs_modulus + 3 * lames_first_parameter + R) / 6
        shear_modulus = (youngs_modulus - 3 * lames_first_parameter + R) / 4
        poissons_ratio = 2 * lames_first_parameter / (youngs_modulus + lames_first_parameter + R)
        p_wave_modulus = (youngs_modulus - lames_first_parameter + R) / 2

        return ElasticProperties(bulk_modulus, youngs_modulus, lames_first_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_youngs_and_shear_modulus(youngs_modulus: float, shear_modulus: float) -> ElasticProperties:
        """
        Convert xxx and xxx to  other elastic property types

         Args:
             xxx (float): Bulk modulus magnitude TODO: [Unitless]
            xxx (float): Young modulus magnitude TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """
        bulk_modulus = youngs_modulus * shear_modulus / (3 * (3 * shear_modulus - youngs_modulus))
        lames_first_parameter = (shear_modulus * (youngs_modulus - 2 * shear_modulus)) / (3 * shear_modulus - youngs_modulus)
        poissons_ratio = 0.5 * (youngs_modulus / shear_modulus) - 1
        p_wave_modulus = (shear_modulus * (4 * shear_modulus - youngs_modulus)) / (3 * shear_modulus - youngs_modulus)

        return ElasticProperties(bulk_modulus, youngs_modulus, lames_first_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_youngs_and_poissons_ratio_modulus(youngs_modulus: float, poissons_ratio: float) -> ElasticProperties:
        """
        Convert xxx and xxx to  other elastic property types

         Args:
             xxx (float): Bulk modulus magnitude TODO: [Unitless]
            xxx (float): Young modulus magnitude TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """
        bulk_modulus = youngs_modulus / (3 * (1 - 2 * poissons_ratio))
        lames_first_parameter = youngs_modulus * poissons_ratio / ((1 + poissons_ratio) * (1 - 2 * poissons_ratio))
        shear_modulus = youngs_modulus / (2 * (1 + poissons_ratio))
        p_wave_modulus = (youngs_modulus * (1 - poissons_ratio)) / ((1 + poissons_ratio) * (1 - 2 * poissons_ratio))

        return ElasticProperties(bulk_modulus, youngs_modulus, lames_first_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_youngs_and_p_wave_modulus(youngs_modulus: float, p_wave_modulus: float) -> ElasticProperties:
        """
        Convert xxx and xxx to  other elastic property types

         Args:
             xxx (float): Bulk modulus magnitude TODO: [Unitless]
            xxx (float): Young modulus magnitude TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """
        S = math.sqrt(youngs_modulus**2 + 9 * p_wave_modulus**2 - 10 * youngs_modulus * p_wave_modulus)

        bulk_modulus = (3 * p_wave_modulus - youngs_modulus + S) / 6
        lames_first_parameter = (p_wave_modulus - youngs_modulus + S) / 4
        shear_modulus = (3 * p_wave_modulus + youngs_modulus - S) / 8
        poissons_ratio = (youngs_modulus - p_wave_modulus + S) / (4 * p_wave_modulus)

        return ElasticProperties(bulk_modulus, youngs_modulus, lames_first_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_lames_and_shear_modulus(lames_first_parameter: float, shear_modulus: float) -> ElasticProperties:
        """
        Convert xxx and xxx to  other elastic property types

         Args:
             xxx (float): Bulk modulus magnitude TODO: [Unitless]
            xxx (float): Young modulus magnitude TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """
        bulk_modulus = lames_first_parameter + (2 * shear_modulus / 3)
        youngs_modulus = shear_modulus * (3 * lames_first_parameter + 2 * shear_modulus) / (lames_first_parameter + shear_modulus)
        poissons_ratio = lames_first_parameter / (2 * (lames_first_parameter + shear_modulus))
        p_wave_modulus = lames_first_parameter + 2 * shear_modulus

        return ElasticProperties(bulk_modulus, youngs_modulus, lames_first_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_lames_and_poissons_ratio_modulus(lames_first_parameter: float, poissons_ratio: float) -> ElasticProperties:
        """
        Convert xxx and xxx to  other elastic property types

         Args:
             xxx (float): Bulk modulus magnitude TODO: [Unitless]
            xxx (float): Young modulus magnitude TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """
        bulk_modulus = lames_first_parameter * (1 + poissons_ratio) / (3 * poissons_ratio)
        youngs_modulus = (lames_first_parameter * (1 + poissons_ratio) * (1 - 2 * poissons_ratio)) / poissons_ratio
        shear_modulus = lames_first_parameter * (1 - 2 * poissons_ratio) / (2 * poissons_ratio)
        p_wave_modulus = (lames_first_parameter * (1 - poissons_ratio)) / poissons_ratio

        return ElasticProperties(bulk_modulus, youngs_modulus, lames_first_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_lames_and_p_wave_modulus(lames_first_parameter: float, p_wave_modulus: float) -> ElasticProperties:
        """
        Convert xxx and xxx to  other elastic property types

         Args:
             xxx (float): Bulk modulus magnitude TODO: [Unitless]
            xxx (float): Young modulus magnitude TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """
        bulk_modulus = (p_wave_modulus + 2 * lames_first_parameter) / 3
        youngs_modulus = ((p_wave_modulus - lames_first_parameter) * (p_wave_modulus + 2 * lames_first_parameter)) / (p_wave_modulus + lames_first_parameter)
        shear_modulus = (p_wave_modulus - lames_first_parameter) / 2
        poissons_ratio = lames_first_parameter / (p_wave_modulus + lames_first_parameter)

        return ElasticProperties(bulk_modulus, youngs_modulus, lames_first_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_shear_and_poissons_ratio_modulus(shear_modulus: float, poissons_ratio: float) -> ElasticProperties:
        """
        Convert xxx and xxx to  other elastic property types

         Args:
             xxx (float): Bulk modulus magnitude TODO: [Unitless]
            xxx (float): Young modulus magnitude TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """
        bulk_modulus = (2 * shear_modulus * (1 + poissons_ratio)) / (3 * (1 - 2 * poissons_ratio))
        youngs_modulus = 2 * shear_modulus * (1 + poissons_ratio)
        lames_first_parameter = 2 * shear_modulus * poissons_ratio / (1 - 2 * poissons_ratio)
        p_wave_modulus = (2 * shear_modulus * (1 - poissons_ratio)) / (1 - 2 * poissons_ratio)

        return ElasticProperties(bulk_modulus, youngs_modulus, lames_first_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_shear_and_p_wave_modulus(shear_modulus: float, p_wave_modulus: float) -> ElasticProperties:
        """
        Convert xxx and xxx to  other elastic property types

         Args:
             xxx (float): Bulk modulus magnitude TODO: [Unitless]
            xxx (float): Young modulus magnitude TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """
        bulk_modulus = p_wave_modulus - (4 / 3) * shear_modulus
        youngs_modulus = (shear_modulus * (3 * p_wave_modulus - 4 * shear_modulus)) / (p_wave_modulus - shear_modulus)
        lames_first_parameter = p_wave_modulus - 2 * shear_modulus
        poissons_ratio = (p_wave_modulus - 2 * shear_modulus) / (2 * p_wave_modulus - 2 * shear_modulus)

        return ElasticProperties(bulk_modulus, youngs_modulus, lames_first_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_poissons_ratio_and_p_wave_modulus(poissons_ratio: float, p_wave_modulus: float) -> ElasticProperties:
        """
        Convert xxx and xxx to  other elastic property types

         Args:
             xxx (float): Bulk modulus magnitude TODO: [Unitless]
            xxx (float): Young modulus magnitude TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """
        bulk_modulus = (p_wave_modulus * (1 + poissons_ratio)) / (3 * (1 - poissons_ratio))
        youngs_modulus = (p_wave_modulus * (1 + poissons_ratio) * (1 - 2 * poissons_ratio)) / (1 - poissons_ratio)
        lames_first_parameter = p_wave_modulus * poissons_ratio / (1 - poissons_ratio)
        shear_modulus = (p_wave_modulus * (1 - 2 * poissons_ratio)) / (2 * (1 - poissons_ratio))

        return ElasticProperties(bulk_modulus, youngs_modulus, lames_first_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_velocity(p_wave_velocity: float, s_wave_velocity: float, density: float) -> ElasticProperties:
        """
        Convert Convert P and S Wave Velocity and Bulk Density to all elastic property types

         Args:
            p_wave_velocity (float): Compressional Velocity TODO: [Unitless]
            s_wave_velocity (float): Shear Velocity TODO: [Unitless]
            density (float): Bulk Density TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_shear_and_p_wave_modulus(shear_modulus, p_wave_modulus)
             ElasticProperties(bulk_modulus=40.56, youngs_modulus=32.55, lames_first_parameter=32.60, shear_modulus=11.91, poissons_ratio=0.366, p_wave_modulus==56.43)
        """
        shear_modulus = density * s_wave_velocity**2
        p_wave_modulus = density * p_wave_velocity**2

        return ElasticPropertiesConverter.from_shear_and_p_wave_modulus(shear_modulus, p_wave_modulus)

    @staticmethod
    def from_slowness(p_wave_slowness: float, s_wave_slowness: float, density: float) -> ElasticProperties:
        """
        Convert Convert P and S Wave Velocity and Bulk Density to all elastic property types

         Args:
            p_wave_slowness (float): Compressional Velocity TODO: [Unitless]
            s_wave_slowness (float): Shear Velocity TODO: [Unitless]
            density (float): Bulk Density TODO: [Unitless]

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> from_slowness(150, 200, 2500)
        """
        shear_modulus = density * (304800 / s_wave_slowness) ** 2
        p_wave_modulus = density * (304800 / p_wave_slowness) ** 2

        return ElasticPropertiesConverter.from_shear_and_p_wave_modulus(shear_modulus, p_wave_modulus)
