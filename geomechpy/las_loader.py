from dataclasses import dataclass
from pathlib import Path

import lasio


@dataclass(frozen=True)
class WellLogData:
    """Well log curves loaded from a LAS file.

    Attributes:
        well_name (str): Well name as recorded in the LAS well information section.
        depth (list[float]): Index curve values (measured or true vertical depth). Unit: as recorded in `depth_unit`
        depth_unit (str): Unit of the depth index curve as recorded in the LAS file (e.g. 'FT', 'M').
        curves (dict[str, list[float]]): Curve values keyed by mnemonic (e.g. 'DTCO', 'RHOB'). Missing/NULL samples are `NaN`.
        curve_units (dict[str, str]): Unit string per curve mnemonic as recorded in the LAS file."""

    well_name: str
    depth: list[float]
    depth_unit: str
    curves: dict[str, list[float]]
    curve_units: dict[str, str]


class LasFileLoader:
    """Load well log curves from LAS (Log ASCII Standard) files for use with the geomechpy converters.

    Reference:
        https://www.cwls.org/products/#products-las"""

    @staticmethod
    def load_las_file(filepath: str | Path) -> WellLogData:
        """Load a LAS file into a `WellLogData` instance.

        NULL samples (as declared by the LAS file's NULL header item) are represented as `NaN` rather than being dropped,
        so gaps in the log remain visible to downstream calculations.

        Args:
            filepath (str | Path): Path to the LAS file.

        Returns:
            WellLogData: Well name, depth index and curve values/units read from the file. See `WellLogData` for details."""
        las = lasio.read(str(filepath))

        well_name = str(las.well.WELL.value) if "WELL" in las.well else ""
        depth_mnemonic = las.curves[0].mnemonic

        curves: dict[str, list[float]] = {}
        curve_units: dict[str, str] = {}
        for curve in las.curves:
            if curve.mnemonic == depth_mnemonic:
                continue
            curves[curve.mnemonic] = curve.data.tolist()
            curve_units[curve.mnemonic] = curve.unit

        return WellLogData(
            well_name=well_name,
            depth=las.index.tolist(),
            depth_unit=las.index_unit or "",
            curves=curves,
            curve_units=curve_units,
        )

    @staticmethod
    def get_curve(well_log_data: WellLogData, mnemonic: str) -> list[float]:
        """Retrieve a curve's values from loaded well log data by mnemonic.

        Args:
            well_log_data (WellLogData): Well log data loaded via `load_las_file`.
            mnemonic (str): Curve mnemonic to retrieve (e.g. 'DTCO').

        Returns:
            list[float]: Curve values. Unit: as recorded in `well_log_data.curve_units[mnemonic]`

        Raises:
            KeyError: If the mnemonic is not present in the well log data. The error lists the available mnemonics."""
        if mnemonic not in well_log_data.curves:
            available = ", ".join(sorted(well_log_data.curves)) or "none"
            raise KeyError(f"Curve mnemonic '{mnemonic}' not found. Available mnemonics: {available}")

        return well_log_data.curves[mnemonic]

    @staticmethod
    def convert_density_g_cc_to_kg_m3(density_g_cc: list[float]) -> list[float]:
        """Convert bulk density from grams per cubic centimeter (a common LAS RHOB unit) to kilograms per cubic meter.

        Args:
            density_g_cc (list[float]): Bulk density values. Unit: g/cc

        Returns:
            list[float]: Bulk density values. Unit: kg/m3"""
        return [value * 1000 for value in density_g_cc]
