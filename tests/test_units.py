import math

import pytest

from geomechpy.units import AngleUnit, DensityUnit, LengthUnit, PressureUnit, UnitConverter

TOLERANCE = 1e-6


class TestConvertPressure:
    def test_psi_to_pa(self) -> None:
        result = UnitConverter.convert(1.0, PressureUnit.PSI, PressureUnit.PA)
        assert result == pytest.approx(6894.757293168, rel=TOLERANCE)

    def test_mpsi_to_psi(self) -> None:
        result = UnitConverter.convert(2.0, PressureUnit.MPSI, PressureUnit.PSI)
        assert result == pytest.approx(2_000_000.0, rel=TOLERANCE)

    def test_gpa_to_psi_roundtrip(self) -> None:
        gpa_value = 10.0
        psi_value = UnitConverter.convert(gpa_value, PressureUnit.GPA, PressureUnit.PSI)
        roundtrip = UnitConverter.convert(psi_value, PressureUnit.PSI, PressureUnit.GPA)
        assert roundtrip == pytest.approx(gpa_value, rel=TOLERANCE)

    def test_same_unit_is_identity(self) -> None:
        result = UnitConverter.convert(42.0, PressureUnit.MPA, PressureUnit.MPA)
        assert result == pytest.approx(42.0, rel=TOLERANCE)


class TestConvertDensity:
    def test_g_cc_to_kg_m3(self) -> None:
        result = UnitConverter.convert(2.4, DensityUnit.G_CC, DensityUnit.KG_M3)
        assert result == pytest.approx(2400.0, rel=TOLERANCE)

    def test_kg_m3_to_g_cc(self) -> None:
        result = UnitConverter.convert(2400.0, DensityUnit.KG_M3, DensityUnit.G_CC)
        assert result == pytest.approx(2.4, rel=TOLERANCE)


class TestConvertLength:
    def test_ft_to_m(self) -> None:
        result = UnitConverter.convert(1.0, LengthUnit.FT, LengthUnit.M)
        assert result == pytest.approx(0.3048, rel=TOLERANCE)

    def test_m_to_ft(self) -> None:
        result = UnitConverter.convert(1.0, LengthUnit.M, LengthUnit.FT)
        assert result == pytest.approx(1 / 0.3048, rel=TOLERANCE)


class TestConvertAngle:
    def test_deg_to_rad(self) -> None:
        result = UnitConverter.convert(180.0, AngleUnit.DEG, AngleUnit.RAD)
        assert result == pytest.approx(math.pi, rel=TOLERANCE)

    def test_rad_to_deg(self) -> None:
        result = UnitConverter.convert(math.pi, AngleUnit.RAD, AngleUnit.DEG)
        assert result == pytest.approx(180.0, rel=TOLERANCE)


class TestCrossFamilyConversionRejected:
    def test_pressure_to_density_raises(self) -> None:
        with pytest.raises(ValueError, match="different unit families"):
            UnitConverter.convert(1.0, PressureUnit.PSI, DensityUnit.KG_M3)

    def test_length_to_angle_raises(self) -> None:
        with pytest.raises(ValueError, match="different unit families"):
            UnitConverter.convert(1.0, LengthUnit.FT, AngleUnit.DEG)


class TestConvertArray:
    def test_converts_each_value(self) -> None:
        result = UnitConverter.convert_array([1.0, 2.0, 3.0], PressureUnit.MPSI, PressureUnit.PSI)
        assert result == pytest.approx([1_000_000.0, 2_000_000.0, 3_000_000.0], rel=TOLERANCE)

    def test_cross_family_raises(self) -> None:
        with pytest.raises(ValueError, match="different unit families"):
            UnitConverter.convert_array([1.0], PressureUnit.PSI, DensityUnit.KG_M3)
