import math
from dataclasses import dataclass


@dataclass(frozen=True)
class ElasticProperties:
    """Elastic properties model per depth. Can be populated manually or via converter

    Attributes:
        Bulk_modulus (float): The bulk modulus (K or B or k) of a substance is a measure of the resistance of a substance to bulk compression. Unit: Pressure
        Youngs_modulus (float): Young modulus is a mechanical property of solid materials that measures the tensile or compressive stiffness when the force is applied lengthwise. Unit: Pressure
        Lame_parameter (float): Lame parameters is a material-dependent quantity denoted by λ which arises in strain-stress relationships. Unit: Pressure
        Shear_modulus (float): Shear modulus is a measure of the elastic shear stiffness of a material and is defined as the ratio of shear stress to the shear strain. Unit: Pressure
        Poissons_ratio (float): Poisson's ratio is a measure of the Poisson effect, the deformation of a material in directions perpendicular to the specific direction of loading. Unit: unitless
        P_wave_modulus (float): P-wave modulus is one of the elastic moduli available to describe isotropic homogeneous materials. Unit: Pressure
    """

    bulk_modulus: float
    youngs_modulus: float
    lame_parameter: float
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
        InputUnit: Pressure Unit of any type (GPa, psi, Mpsi etc)
        Pressure unit of both inputs needs to be consistent

         Args:
             bulk_modulus (float): Bulk modulus magnitude Unit: Pressure Unit
             youngs_modulus (float): Young modulus magnitude Unit: Pressure Unit

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example for input in:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(19.7184,33.2748)
            ElasticProperties(bulk_modulus=19.7184, youngs_modulus=33.2748, lame_parameter=10.617599999999998, shear_modulus=13.6512, poissons_ratio=0.21874999999999997, p_wave_modulus=37.919999999999995)
        """

        lame_parameter = 3 * bulk_modulus * (3 * bulk_modulus - youngs_modulus) / (9 * bulk_modulus - youngs_modulus)
        shear_modulus = 3 * bulk_modulus * youngs_modulus / (9 * bulk_modulus - youngs_modulus)
        poissons_ratio = (3 * bulk_modulus - youngs_modulus) / (6 * bulk_modulus)
        p_wave_modulus = 3 * bulk_modulus * (3 * bulk_modulus + youngs_modulus) / (9 * bulk_modulus - youngs_modulus)

        return ElasticProperties(bulk_modulus, youngs_modulus, lame_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_bulk_and_lame_modulus(bulk_modulus: float, lame_parameter: float) -> ElasticProperties:
        """
        Convert Bulk modulus and Lame parameter to  other elastic property types
        InputUnit: Pressure Unit of any type (GPa, psi, Mpsi etc)
        Pressure unit of both inputs needs to be consistent

         Args:
            bulk_modulus (float): Bulk modulus magnitude Unit: Pressure Unit
            lame_parameter (float): Lame parameter magnitude Unit: Pressure Unit

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_lame_modulus(19.7184, 10.61851)
             ElasticProperties(bulk_modulus=19.7184, youngs_modulus=33.2748, lame_parameter=10.617599999999998, shear_modulus=13.6512, poissons_ratio=0.21874999999999997, p_wave_modulus=37.919999999999995)
        """

        youngs_modulus = 9 * bulk_modulus * (bulk_modulus - lame_parameter) / (3 * bulk_modulus - lame_parameter)
        shear_modulus = 3 * (bulk_modulus - lame_parameter) / 2
        poissons_ratio = lame_parameter / (3 * bulk_modulus - lame_parameter)
        p_wave_modulus = 3 * bulk_modulus - 2 * lame_parameter

        return ElasticProperties(bulk_modulus, youngs_modulus, lame_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_bulk_and_shear_modulus(bulk_modulus: float, shear_modulus: float) -> ElasticProperties:
        """
        Convert Bulk modulus and Shear modulus to  other elastic property types
        InputUnit: Pressure Unit of any type (GPa, psi, Mpsi etc)
        Pressure unit of both inputs needs to be consistent


         Args:
             bulk_modulus (float): Bulk modulus magnitude Unit: Pressure Unit
             shear_modulus (float): Shear modulus magnitude Unit: Pressure Unit

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_shear_modulus(19.7184, 13.6512)
             ElasticProperties(bulk_modulus=19.7184, youngs_modulus=33.2748, lame_parameter=10.617599999999998, shear_modulus=13.6512, poissons_ratio=0.21874999999999997, p_wave_modulus=37.919999999999995)
        """
        youngs_modulus = 9 * bulk_modulus * shear_modulus / (3 * bulk_modulus + shear_modulus)
        lame_parameter = bulk_modulus - 2 * shear_modulus / 3
        poissons_ratio = (3 * bulk_modulus - 2 * shear_modulus) / (2 * (3 * bulk_modulus + shear_modulus))
        p_wave_modulus = bulk_modulus + 4 * shear_modulus / 3

        return ElasticProperties(bulk_modulus, youngs_modulus, lame_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_bulk_and_poissons_ratio_modulus(bulk_modulus: float, poissons_ratio: float) -> ElasticProperties:
        """
        Convert xxx and xxx to  other elastic property types
        InputUnit: Pressure Unit of any type (GPa, psi, Mpsi etc)
        Pressure unit of both inputs needs to be consistent
         Args:
             bulk_modulus (float): Bulk modulus magnitude Unit: Pressure Unit
             poissons_ratio (float): Young modulus magnitude Unit: Unitless

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_poissons_ratio_modulus(40.56, 32.55)
             ElasticProperties(bulk_modulus=19.7184, youngs_modulus=33.2748, lame_parameter=10.617599999999998, shear_modulus=13.6512, poissons_ratio=0.21874999999999997, p_wave_modulus=37.919999999999995)
        """
        youngs_modulus = 3 * bulk_modulus * (1 - 2 * poissons_ratio)
        lame_parameter = 3 * bulk_modulus * poissons_ratio / (1 + poissons_ratio)
        shear_modulus = 3 * bulk_modulus * (1 - 2 * poissons_ratio) / (2 * (1 + poissons_ratio))
        p_wave_modulus = 3 * bulk_modulus * (1 - poissons_ratio) / (1 + poissons_ratio)

        return ElasticProperties(bulk_modulus, youngs_modulus, lame_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_bulk_and_p_wave_modulus(bulk_modulus: float, p_wave_modulus: float) -> ElasticProperties:
        """
        Convert Bulk and P-wave modulus to  other elastic property types
        InputUnit: Pressure Unit of any type (GPa, psi, Mpsi etc)
        Pressure unit of both inputs needs to be consistent
         Args:
             bulk_modulus (float): Bulk modulus magnitude Unit: Pressure Unit
             p_wave_modulus (float): P-wave modulus magnitude Unit: Pressure Unit

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_p_wave_modulus(19.7184, 37.92)
             ElasticProperties(bulk_modulus=19.7184, youngs_modulus=33.2748, lame_parameter=10.617599999999998, shear_modulus=13.6512, poissons_ratio=0.21874999999999997, p_wave_modulus=37.919999999999995)
        """
        youngs_modulus = 9 * bulk_modulus * (p_wave_modulus - bulk_modulus) / (3 * bulk_modulus + p_wave_modulus)
        lame_parameter = (3 * bulk_modulus - p_wave_modulus) / 2
        shear_modulus = (3 * (p_wave_modulus - bulk_modulus)) / 4
        poissons_ratio = (3 * bulk_modulus - p_wave_modulus) / (3 * bulk_modulus + p_wave_modulus)

        return ElasticProperties(bulk_modulus, youngs_modulus, lame_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_youngs_and_lame_modulus(youngs_modulus: float, lame_parameter: float) -> ElasticProperties:
        """
        Convert Young's modulus and Lame parameter to  other elastic property types
        InputUnit: Pressure Unit of any type (GPa, psi, Mpsi etc)
        Pressure unit of both inputs needs to be consistent
         Args:
             youngs_modulus (float): Young's modulus magnitude Unit: Pressure Unit
             lame_parameter (float): Lame parameter magnitude Unit: Pressure Unit

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_youngs_and_lame_modulus(33.2748, 10.61851)
             ElasticProperties(bulk_modulus=19.7184, youngs_modulus=33.2748, lame_parameter=10.617599999999998, shear_modulus=13.6512, poissons_ratio=0.21874999999999997, p_wave_modulus=37.919999999999995)
        """
        R = math.sqrt(youngs_modulus**2 + 9 * lame_parameter**2 + 2 * youngs_modulus * lame_parameter)

        bulk_modulus = (youngs_modulus + 3 * lame_parameter + R) / 6
        shear_modulus = (youngs_modulus - 3 * lame_parameter + R) / 4
        poissons_ratio = 2 * lame_parameter / (youngs_modulus + lame_parameter + R)
        p_wave_modulus = (youngs_modulus - lame_parameter + R) / 2

        return ElasticProperties(bulk_modulus, youngs_modulus, lame_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_youngs_and_shear_modulus(youngs_modulus: float, shear_modulus: float) -> ElasticProperties:
        """
        Convert Young's and Shear modulus to  other elastic property types
        InputUnit: Pressure Unit of any type (GPa, psi, Mpsi etc)
        Pressure unit of both inputs needs to be consistent
         Args:
             youngs_modulus (float): Young's modulus magnitude Unit: Pressure Unit
             shear_modulus (float):  Shear modulus magnitude Unit: Pressure Unit

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_youngs_and_shear_modulus(33.2748, 13.6512)
             ElasticProperties(bulk_modulus=19.7184, youngs_modulus=33.2748, lame_parameter=10.617599999999998, shear_modulus=13.6512, poissons_ratio=0.21874999999999997, p_wave_modulus=37.919999999999995)
        """
        bulk_modulus = youngs_modulus * shear_modulus / (3 * (3 * shear_modulus - youngs_modulus))
        lame_parameter = (shear_modulus * (youngs_modulus - 2 * shear_modulus)) / (3 * shear_modulus - youngs_modulus)
        poissons_ratio = 0.5 * (youngs_modulus / shear_modulus) - 1
        p_wave_modulus = (shear_modulus * (4 * shear_modulus - youngs_modulus)) / (3 * shear_modulus - youngs_modulus)

        return ElasticProperties(bulk_modulus, youngs_modulus, lame_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_youngs_and_poissons_ratio_modulus(youngs_modulus: float, poissons_ratio: float) -> ElasticProperties:
        """
        Convert Young's modulus and Poisson's ratio to  other elastic property types
        InputUnit: Young's modulus - pressure unit  | Poisson's ratio - unitless

         Args:
             youngs_modulus (float): Bulk modulus magnitude Unit: Pressure Unit
             poissons_ratio (float): Poisson's ratio magnitude Unit: unitless

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_youngs_and_poissons_ratio_modulus(33.2748, 0.21875)
             ElasticProperties(bulk_modulus=19.7184, youngs_modulus=33.2748, lame_parameter=10.617599999999998, shear_modulus=13.6512, poissons_ratio=0.21874999999999997, p_wave_modulus=37.919999999999995)
        """
        bulk_modulus = youngs_modulus / (3 * (1 - 2 * poissons_ratio))
        lame_parameter = youngs_modulus * poissons_ratio / ((1 + poissons_ratio) * (1 - 2 * poissons_ratio))
        shear_modulus = youngs_modulus / (2 * (1 + poissons_ratio))
        p_wave_modulus = (youngs_modulus * (1 - poissons_ratio)) / ((1 + poissons_ratio) * (1 - 2 * poissons_ratio))

        return ElasticProperties(bulk_modulus, youngs_modulus, lame_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_youngs_and_p_wave_modulus(youngs_modulus: float, p_wave_modulus: float) -> ElasticProperties:
        """
        Convert Young's modulus and P-wave modulus to  other elastic property types
        InputUnit: Pressure Unit of any type (GPa, psi, Mpsi etc)
        Pressure unit of both inputs needs to be consistent

         Args:
             youngs_modulus (float): Young's modulus magnitude Unit: Pressure Unit
             p_wave_modulus (float): P-wave modulus magnitude Unit: Pressure Unit

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_youngs_and_p_wave_modulus(33.2748, 37.92)
             ElasticProperties(bulk_modulus=19.7184, youngs_modulus=33.2748, lame_parameter=10.617599999999998, shear_modulus=13.6512, poissons_ratio=0.21874999999999997, p_wave_modulus=37.919999999999995)
        """
        S = math.sqrt(youngs_modulus**2 + 9 * p_wave_modulus**2 - 10 * youngs_modulus * p_wave_modulus)

        bulk_modulus = (3 * p_wave_modulus - youngs_modulus + S) / 6
        lame_parameter = (p_wave_modulus - youngs_modulus + S) / 4
        shear_modulus = (3 * p_wave_modulus + youngs_modulus - S) / 8
        poissons_ratio = (youngs_modulus - p_wave_modulus + S) / (4 * p_wave_modulus)

        return ElasticProperties(bulk_modulus, youngs_modulus, lame_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_lame_and_shear_modulus(lame_parameter: float, shear_modulus: float) -> ElasticProperties:
        """
        Convert Lame parameter and Shear modulus to  other elastic property types
        InputUnit: Pressure Unit of any type (GPa, psi, Mpsi etc)
        Pressure unit of both inputs needs to be consistent

         Args:
             lame_parameter (float): Lame parameter magnitude Unit: Pressure Unit
             shear_modulus (float): Shear modulus magnitude Unit: Pressure Unit

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_lame_and_shear_modulus(10.61851, 13.6512)
             ElasticProperties(bulk_modulus=19.7184, youngs_modulus=33.2748, lame_parameter=10.617599999999998, shear_modulus=13.6512, poissons_ratio=0.21874999999999997, p_wave_modulus=37.919999999999995)
        """
        bulk_modulus = lame_parameter + (2 * shear_modulus / 3)
        youngs_modulus = shear_modulus * (3 * lame_parameter + 2 * shear_modulus) / (lame_parameter + shear_modulus)
        poissons_ratio = lame_parameter / (2 * (lame_parameter + shear_modulus))
        p_wave_modulus = lame_parameter + 2 * shear_modulus

        return ElasticProperties(bulk_modulus, youngs_modulus, lame_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_lame_and_poissons_ratio_modulus(lame_parameter: float, poissons_ratio: float) -> ElasticProperties:
        """
        Convert Lame parameter and Poisson's ratio to  other elastic property types
        InputUnit: Lame parameter - pressure unit  | Poisson's ratio - unitless

         Args:
             lame_parameter (float): Lame parameter magnitude Unit: Pressure Unit
             poissons_ratio (float): Poisson's ratio magnitude Unit: unitless

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_lame_and_poissons_ratio_modulus(10.61851, 0.21875)
             ElasticProperties(bulk_modulus=19.7184, youngs_modulus=33.2748, lame_parameter=10.617599999999998, shear_modulus=13.6512, poissons_ratio=0.21874999999999997, p_wave_modulus=37.919999999999995)
        """
        bulk_modulus = lame_parameter * (1 + poissons_ratio) / (3 * poissons_ratio)
        youngs_modulus = (lame_parameter * (1 + poissons_ratio) * (1 - 2 * poissons_ratio)) / poissons_ratio
        shear_modulus = lame_parameter * (1 - 2 * poissons_ratio) / (2 * poissons_ratio)
        p_wave_modulus = (lame_parameter * (1 - poissons_ratio)) / poissons_ratio

        return ElasticProperties(bulk_modulus, youngs_modulus, lame_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_lame_and_p_wave_modulus(lame_parameter: float, p_wave_modulus: float) -> ElasticProperties:
        """
        Convert Lame parameter and P-wave modulus to  other elastic property types
        InputUnit: Pressure Unit of any type (GPa, psi, Mpsi etc)
        Pressure unit of both inputs needs to be consistent

         Args:
             lame_parameter (float): Lame parameter magnitude Unit: Pressure Unit
             p_wave_modulus (float): P-wave modulus magnitude Unit: Pressure Unit

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_lame_and_p_wave_modulus(10.61851, 37.92)
             ElasticProperties(bulk_modulus=19.7184, youngs_modulus=33.2748, lame_parameter=10.617599999999998, shear_modulus=13.6512, poissons_ratio=0.21874999999999997, p_wave_modulus=37.919999999999995)
        """
        bulk_modulus = (p_wave_modulus + 2 * lame_parameter) / 3
        youngs_modulus = ((p_wave_modulus - lame_parameter) * (p_wave_modulus + 2 * lame_parameter)) / (p_wave_modulus + lame_parameter)
        shear_modulus = (p_wave_modulus - lame_parameter) / 2
        poissons_ratio = lame_parameter / (p_wave_modulus + lame_parameter)

        return ElasticProperties(bulk_modulus, youngs_modulus, lame_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_shear_and_poissons_ratio_modulus(shear_modulus: float, poissons_ratio: float) -> ElasticProperties:
        """
        Convert Shear modulus and Poisson's ratio to  other elastic property types
        InputUnit: Young's modulus - pressure unit  | Poisson's ratio - unitless

         Args:
             shear_modulus (float): Shear modulus magnitude Unit: Pressure Unit
             poissons_ratio (float): Poisson's ratio magnitude Unit: unitless

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_shear_and_poissons_ratio_modulus(13.6512, 0.21875)
             ElasticProperties(bulk_modulus=19.7184, youngs_modulus=33.2748, lame_parameter=10.617599999999998, shear_modulus=13.6512, poissons_ratio=0.21874999999999997, p_wave_modulus=37.919999999999995)
        """
        bulk_modulus = (2 * shear_modulus * (1 + poissons_ratio)) / (3 * (1 - 2 * poissons_ratio))
        youngs_modulus = 2 * shear_modulus * (1 + poissons_ratio)
        lame_parameter = 2 * shear_modulus * poissons_ratio / (1 - 2 * poissons_ratio)
        p_wave_modulus = (2 * shear_modulus * (1 - poissons_ratio)) / (1 - 2 * poissons_ratio)

        return ElasticProperties(bulk_modulus, youngs_modulus, lame_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_shear_and_p_wave_modulus(shear_modulus: float, p_wave_modulus: float) -> ElasticProperties:
        """
        Convert Shear modulus and P-wave modulus to  other elastic property types
        InputUnit: Pressure Unit of any type (GPa, psi, Mpsi etc)
        Pressure unit of both inputs needs to be consistent

         Args:
             shear_modulus (float): Shear modulus magnitude Unit: Pressure Unit
             p_wave_modulus (float): P-wave modulus magnitude Unit: Pressure Unit

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_shear_and_p_wave_modulus(13.6512, 37.92)
             ElasticProperties(bulk_modulus=19.7184, youngs_modulus=33.2748, lame_parameter=10.617599999999998, shear_modulus=13.6512, poissons_ratio=0.21874999999999997, p_wave_modulus=37.919999999999995)
        """
        bulk_modulus = p_wave_modulus - (4 / 3) * shear_modulus
        youngs_modulus = (shear_modulus * (3 * p_wave_modulus - 4 * shear_modulus)) / (p_wave_modulus - shear_modulus)
        lame_parameter = p_wave_modulus - 2 * shear_modulus
        poissons_ratio = (p_wave_modulus - 2 * shear_modulus) / (2 * p_wave_modulus - 2 * shear_modulus)

        return ElasticProperties(bulk_modulus, youngs_modulus, lame_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_poissons_ratio_and_p_wave_modulus(poissons_ratio: float, p_wave_modulus: float) -> ElasticProperties:
        """
        Convert Poisson's ratio and P-wave modulus to  other elastic property types
        InputUnit: Poisson's ratio - unitless | P-wave modulus - pressure unit

         Args:
             poissons_ratio (float): Poisson's ratio magnitude Unit: unitless
             p_wave_modulus (float): P-wave modulus magnitude Unit: Pressure Unit

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> ElasticPropertiesConverter.from_bulk_and_youngs_modulus(0.21875, 37.92)
             ElasticProperties(bulk_modulus=19.7184, youngs_modulus=33.2748, lame_parameter=10.617599999999998, shear_modulus=13.6512, poissons_ratio=0.21874999999999997, p_wave_modulus=37.919999999999995)
        """
        bulk_modulus = (p_wave_modulus * (1 + poissons_ratio)) / (3 * (1 - poissons_ratio))
        youngs_modulus = (p_wave_modulus * (1 + poissons_ratio) * (1 - 2 * poissons_ratio)) / (1 - poissons_ratio)
        lame_parameter = p_wave_modulus * poissons_ratio / (1 - poissons_ratio)
        shear_modulus = (p_wave_modulus * (1 - 2 * poissons_ratio)) / (2 * (1 - poissons_ratio))

        return ElasticProperties(bulk_modulus, youngs_modulus, lame_parameter, shear_modulus, poissons_ratio, p_wave_modulus)

    @staticmethod
    def from_velocity(p_wave_velocity: float, s_wave_velocity: float, density: float) -> ElasticProperties:
        """
        Convert Convert P and S Wave Velocity and Bulk Density to all elastic property types

         Args:
            p_wave_velocity (float): Compressional Velocity Unit: Velocity Unit m/s
            s_wave_velocity (float): Shear Velocity Unit: Velocity Unit m/s
            density (float): Bulk Density Unit: Density unit kg/m3

         Returns:
             ElasticProperties: Dataclass containing computed elastic. See `ElasticProperties` for details
            Output Unit is in Pascal [Pa]
         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> from_velocity(4000, 2400, 2370)
             ElasticProperties(bulk_modulus=19718400000.0, youngs_modulus=33274800000.0, lame_parameter=10617600000, shear_modulus=13651200000, poissons_ratio=0.21875, p_wave_modulus=37920000000)
        """
        shear_modulus = density * s_wave_velocity**2
        p_wave_modulus = density * p_wave_velocity**2

        return ElasticPropertiesConverter.from_shear_and_p_wave_modulus(shear_modulus, p_wave_modulus)

    @staticmethod
    def from_slowness(p_wave_slowness: float, s_wave_slowness: float, density: float) -> ElasticProperties:
        """
        Convert Convert P and S Wave Slowness and Bulk Density to all elastic property types
        Output Unit is in Pascal [Pa]

         Args:
            p_wave_slowness (float):
                Compressional Slowness Unit: us/ft
            s_wave_slowness (float):
                Shear Slowness Unit: us/ft
            density (float):
                Bulk Density Unit: kg/m3

         Returns:
             ElasticProperties: Dataclass containing computed elastic properties. See `ElasticProperties` for details

         Raises:
             ValueError: In case if we need to check values

         Example:
             >>> from_slowness(76.2, 127, 2370)
             ElasticProperties(bulk_modulus=19718400000.0, youngs_modulus=33274800000.0, lame_parameter=10617600000.0, shear_modulus=13651200000.0, poissons_ratio=0.21875, p_wave_modulus=37920000000.0)
        """
        shear_modulus = density * (304800 / s_wave_slowness) ** 2
        p_wave_modulus = density * (304800 / p_wave_slowness) ** 2

        return ElasticPropertiesConverter.from_shear_and_p_wave_modulus(shear_modulus, p_wave_modulus)
