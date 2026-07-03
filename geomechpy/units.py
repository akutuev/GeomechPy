import math
from enum import Enum


class PressureUnit(Enum):
    """Pressure/stress unit family, used for elastic moduli, stresses and rock strength."""

    PA = "Pa"
    KPA = "kPa"
    MPA = "MPa"
    GPA = "GPa"
    PSI = "psi"
    KPSI = "kpsi"
    MPSI = "Mpsi"


class DensityUnit(Enum):
    """Density unit family, used for bulk density."""

    KG_M3 = "kg/m3"
    G_CC = "g/cc"


class LengthUnit(Enum):
    """Length unit family, used for depth."""

    FT = "ft"
    M = "m"


class AngleUnit(Enum):
    """Angle unit family, used for azimuths, deviations and friction angle."""

    DEG = "deg"
    RAD = "rad"


_PRESSURE_TO_PA = {
    PressureUnit.PA: 1.0,
    PressureUnit.KPA: 1e3,
    PressureUnit.MPA: 1e6,
    PressureUnit.GPA: 1e9,
    PressureUnit.PSI: 6894.757293168,
    PressureUnit.KPSI: 6894.757293168e3,
    PressureUnit.MPSI: 6894.757293168e6,
}

_DENSITY_TO_KG_M3 = {
    DensityUnit.KG_M3: 1.0,
    DensityUnit.G_CC: 1000.0,
}

_LENGTH_TO_M = {
    LengthUnit.FT: 0.3048,
    LengthUnit.M: 1.0,
}

_ANGLE_TO_RAD = {
    AngleUnit.DEG: math.pi / 180,
    AngleUnit.RAD: 1.0,
}

_FAMILY_TO_BASE_TABLE: dict[type, dict[Enum, float]] = {
    PressureUnit: _PRESSURE_TO_PA,
    DensityUnit: _DENSITY_TO_KG_M3,
    LengthUnit: _LENGTH_TO_M,
    AngleUnit: _ANGLE_TO_RAD,
}


class UnitConverter:
    """Convert values between units within the same unit family (e.g. psi -> GPa, ft -> m).

    Conversion is only permitted between units of the same family (same Enum type). Converting between units of
    different families (e.g. a `PressureUnit` and a `DensityUnit`) raises a `ValueError`, since such a conversion
    is physically meaningless."""

    @staticmethod
    def convert(value: float, from_unit: Enum, to_unit: Enum) -> float:
        """Convert a value from one unit to another unit of the same family.

        Args:
            value (float): Magnitude to convert. Unit: `from_unit`
            from_unit (Enum): Unit of `value`. Must be a member of `PressureUnit`, `DensityUnit`, `LengthUnit` or `AngleUnit`.
            to_unit (Enum): Unit to convert to. Must be a member of the same family (Enum type) as `from_unit`.

        Returns:
            float: Converted magnitude. Unit: `to_unit`

        Raises:
            ValueError: If `from_unit` and `to_unit` belong to different unit families."""
        if type(from_unit) is not type(to_unit):
            raise ValueError(
                f"Cannot convert between different unit families: {type(from_unit).__name__} and {type(to_unit).__name__}"
            )

        conversion_table = _FAMILY_TO_BASE_TABLE[type(from_unit)]
        base_value = value * conversion_table[from_unit]

        return base_value / conversion_table[to_unit]

    @staticmethod
    def convert_array(values: list[float], from_unit: Enum, to_unit: Enum) -> list[float]:
        """Convert an array of values from one unit to another unit of the same family.

        Args:
            values (list[float]): Magnitudes to convert. Unit: `from_unit`
            from_unit (Enum): Unit of `values`. Must be a member of `PressureUnit`, `DensityUnit`, `LengthUnit` or `AngleUnit`.
            to_unit (Enum): Unit to convert to. Must be a member of the same family (Enum type) as `from_unit`.

        Returns:
            list[float]: Converted magnitudes. Unit: `to_unit`

        Raises:
            ValueError: If `from_unit` and `to_unit` belong to different unit families."""
        return [UnitConverter.convert(value, from_unit, to_unit) for value in values]
