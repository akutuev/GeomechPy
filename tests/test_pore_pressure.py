import pytest
from geomechpy.pore_pressure import PorePressureCalculation

TOLERANCE = 1e-6


class TestPorePressureOnshore:
    def test_default_gradient_no_air_gap(self) -> None:
        result = PorePressureCalculation.calculate_pore_pressure_onshore(tvd=10000)
        assert result == pytest.approx(4700.0, rel=TOLERANCE)

    def test_below_air_gap_uses_air_gradient(self) -> None:
        result = PorePressureCalculation.calculate_pore_pressure_onshore(tvd=50, air_gap=100)
        assert result == pytest.approx(0.0004 * 50, rel=TOLERANCE)

    def test_above_air_gap_includes_air_pressure(self) -> None:
        result = PorePressureCalculation.calculate_pore_pressure_onshore(
            tvd=10000, formation_pore_pressure_gradient=0.47, air_gap=100
        )
        expected = 0.0004 * 100 + 0.47 * (10000 - 100)
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_zero_tvd(self) -> None:
        result = PorePressureCalculation.calculate_pore_pressure_onshore(tvd=0, air_gap=100)
        assert result == pytest.approx(0.0, abs=TOLERANCE)

    def test_overpressured_gradient(self) -> None:
        result = PorePressureCalculation.calculate_pore_pressure_onshore(
            tvd=8000, formation_pore_pressure_gradient=0.7
        )
        assert result == pytest.approx(5600.0, rel=TOLERANCE)


class TestPorePressureOffshore:
    def test_below_air_gap(self) -> None:
        result = PorePressureCalculation.calculate_pore_pressure_offshore(tvd=50, air_gap=100, water_depth=1000)
        assert result == pytest.approx(0.0004 * 50, rel=TOLERANCE)

    def test_within_water_column(self) -> None:
        result = PorePressureCalculation.calculate_pore_pressure_offshore(tvd=500, air_gap=100, water_depth=1000)
        expected = 0.0004 * 100 + 0.47 * (500 - 100)
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_below_seabed(self) -> None:
        result = PorePressureCalculation.calculate_pore_pressure_offshore(
            tvd=10000, air_gap=100, water_depth=1000
        )
        expected = 0.0004 * 100 + 0.47 * 1000 + 0.47 * (10000 - 1000 - 100)
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_no_water_no_air_matches_onshore(self) -> None:
        offshore = PorePressureCalculation.calculate_pore_pressure_offshore(tvd=8000)
        onshore = PorePressureCalculation.calculate_pore_pressure_onshore(tvd=8000)
        assert offshore == pytest.approx(onshore, rel=TOLERANCE)
