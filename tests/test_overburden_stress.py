import pytest
from geomechpy.overburden_stress import OverburdenStressCalculation

TOLERANCE = 1e-6


class TestOverburdenStressOnshore:
    def test_default_gradient_no_air_gap(self) -> None:
        result = OverburdenStressCalculation.calculate_overburden_stress_onshore(tvd=10000)
        assert result == pytest.approx(10500.0, rel=TOLERANCE)

    def test_below_air_gap_uses_air_gradient(self) -> None:
        result = OverburdenStressCalculation.calculate_overburden_stress_onshore(tvd=50, air_gap=100)
        assert result == pytest.approx(0.0004 * 50, rel=TOLERANCE)

    def test_above_air_gap_includes_air_pressure(self) -> None:
        result = OverburdenStressCalculation.calculate_overburden_stress_onshore(tvd=10000, lithostatic_gradient=1.05, air_gap=100)
        expected = 0.0004 * 100 + 1.05 * (10000 - 100)
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_zero_tvd(self) -> None:
        result = OverburdenStressCalculation.calculate_overburden_stress_onshore(tvd=0, air_gap=100)
        assert result == pytest.approx(0.0, abs=TOLERANCE)

    def test_custom_lithostatic_gradient(self) -> None:
        result = OverburdenStressCalculation.calculate_overburden_stress_onshore(tvd=5000, lithostatic_gradient=0.9)
        assert result == pytest.approx(4500.0, rel=TOLERANCE)


class TestOverburdenStressOffshore:
    def test_below_air_gap(self) -> None:
        result = OverburdenStressCalculation.calculate_overburden_stress_offshore(tvd=50, air_gap=100, water_depth=1000)
        assert result == pytest.approx(0.0004 * 50, rel=TOLERANCE)

    def test_within_water_column(self) -> None:
        result = OverburdenStressCalculation.calculate_overburden_stress_offshore(tvd=500, air_gap=100, water_depth=1000)
        expected = 0.0004 * 100 + 0.47 * (500 - 100)
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_below_seabed(self) -> None:
        result = OverburdenStressCalculation.calculate_overburden_stress_offshore(tvd=10000, air_gap=100, water_depth=1000)
        expected = 0.0004 * 100 + 0.47 * 1000 + 1.05 * (10000 - 1000 - 100)
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_no_water_no_air_matches_onshore(self) -> None:
        offshore = OverburdenStressCalculation.calculate_overburden_stress_offshore(tvd=8000)
        onshore = OverburdenStressCalculation.calculate_overburden_stress_onshore(tvd=8000)
        assert offshore == pytest.approx(onshore, rel=TOLERANCE)

    def test_custom_sea_water_gradient(self) -> None:
        result = OverburdenStressCalculation.calculate_overburden_stress_offshore(
            tvd=500, air_gap=0, water_depth=1000, sea_water_pressure_gradient=0.45
        )
        assert result == pytest.approx(0.45 * 500, rel=TOLERANCE)
