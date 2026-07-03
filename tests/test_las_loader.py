import math
from pathlib import Path

import pytest

from geomechpy.las_loader import LasFileLoader, WellLogData

LAS_CONTENT = """~VERSION INFORMATION
 VERS.                  2.0 : CWLS LOG ASCII STANDARD -VERSION 2.0
 WRAP.                  NO  : ONE LINE PER DEPTH STEP
~WELL INFORMATION
 STRT.FT              100.0 : START DEPTH
 STOP.FT               102.0 : STOP DEPTH
 STEP.FT                 1.0 : STEP
 NULL.                -999.25 : NULL VALUE
 WELL.               TEST WELL 1 : WELL NAME
~CURVE INFORMATION
 DEPT.FT                  : DEPTH
 DTCO.US/F                : COMPRESSIONAL SLOWNESS
 DTSM.US/F                : SHEAR SLOWNESS
 RHOB.G/C3                : BULK DENSITY
~ASCII
100.0  90.0  180.0  2.40
101.0  -999.25  181.0  2.41
102.0  89.5  179.5  2.42
"""


@pytest.fixture
def las_file(tmp_path: Path) -> Path:
    path = tmp_path / "test_well.las"
    path.write_text(LAS_CONTENT)
    return path


class TestLoadLasFile:
    def test_reads_well_name(self, las_file: Path) -> None:
        result = LasFileLoader.load_las_file(las_file)
        assert result.well_name == "TEST WELL 1"

    def test_reads_depth_index_and_unit(self, las_file: Path) -> None:
        result = LasFileLoader.load_las_file(las_file)
        assert result.depth == pytest.approx([100.0, 101.0, 102.0])
        assert result.depth_unit == "FT"

    def test_depth_curve_excluded_from_curves(self, las_file: Path) -> None:
        result = LasFileLoader.load_las_file(las_file)
        assert "DEPT" not in result.curves

    def test_reads_curve_values(self, las_file: Path) -> None:
        result = LasFileLoader.load_las_file(las_file)
        assert result.curves["DTSM"] == pytest.approx([180.0, 181.0, 179.5])
        assert result.curves["RHOB"] == pytest.approx([2.40, 2.41, 2.42])

    def test_reads_curve_units(self, las_file: Path) -> None:
        result = LasFileLoader.load_las_file(las_file)
        assert result.curve_units["DTCO"] == "US/F"
        assert result.curve_units["RHOB"] == "G/C3"

    def test_null_value_becomes_nan(self, las_file: Path) -> None:
        result = LasFileLoader.load_las_file(las_file)
        assert math.isnan(result.curves["DTCO"][1])


class TestGetCurve:
    def test_returns_curve_values(self) -> None:
        well_log_data = WellLogData(
            well_name="TEST",
            depth=[100.0, 101.0],
            depth_unit="FT",
            curves={"DTCO": [90.0, 91.0]},
            curve_units={"DTCO": "US/F"},
        )
        assert LasFileLoader.get_curve(well_log_data, "DTCO") == [90.0, 91.0]

    def test_missing_mnemonic_raises_key_error_listing_available(self) -> None:
        well_log_data = WellLogData(
            well_name="TEST",
            depth=[100.0],
            depth_unit="FT",
            curves={"DTCO": [90.0], "RHOB": [2.4]},
            curve_units={"DTCO": "US/F", "RHOB": "G/C3"},
        )
        with pytest.raises(KeyError, match="DTCO, RHOB"):
            LasFileLoader.get_curve(well_log_data, "DTSM")


class TestConvertDensityGCcToKgM3:
    def test_converts_single_value(self) -> None:
        assert LasFileLoader.convert_density_g_cc_to_kg_m3([2.4]) == pytest.approx([2400.0])

    def test_converts_multiple_values(self) -> None:
        result = LasFileLoader.convert_density_g_cc_to_kg_m3([2.4, 2.6, 2.65])
        assert result == pytest.approx([2400.0, 2600.0, 2650.0])
