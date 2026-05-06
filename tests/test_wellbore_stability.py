import math

import pytest
from geomechpy.wellbore_stability import WellboreStabilityCalculation

TOLERANCE = 1e-6


class TestBreakdownPressure:
    def test_known_value(self) -> None:
        result = WellboreStabilityCalculation.calculate_breakdown_calculation_vertical_well_analytical(
            shmax=12000, shmin=10000, pprs=5000, tstr=500
        )
        # 3*shmin - shmax - pprs + tstr = 30000 - 12000 - 5000 + 500
        assert result == pytest.approx(13500.0, rel=TOLERANCE)

    def test_returns_float(self) -> None:
        result = WellboreStabilityCalculation.calculate_breakdown_calculation_vertical_well_analytical(
            shmax=12000, shmin=10000, pprs=5000, tstr=500
        )
        assert isinstance(result, float)

    def test_zero_tensile_strength(self) -> None:
        result = WellboreStabilityCalculation.calculate_breakdown_calculation_vertical_well_analytical(
            shmax=12000, shmin=10000, pprs=5000, tstr=0
        )
        assert result == pytest.approx(13000.0, rel=TOLERANCE)

    def test_breakdown_increases_with_tensile_strength(self) -> None:
        low_tstr_breakdown = WellboreStabilityCalculation.calculate_breakdown_calculation_vertical_well_analytical(
            shmax=12000, shmin=10000, pprs=5000, tstr=100
        )
        high_tstr_breakdown = WellboreStabilityCalculation.calculate_breakdown_calculation_vertical_well_analytical(
            shmax=12000, shmin=10000, pprs=5000, tstr=1000
        )
        assert high_tstr_breakdown > low_tstr_breakdown


class TestBreakoutPressure:
    def test_known_value(self) -> None:
        # shmax=12000, shmin=10000, pprs=5000, overburden=13000, ucs=8000, fang=30, pr_sta=0.25
        # q = tan(60)^2 = 3
        # CC = 8000 - 5000*(3-1) = -2000
        # A = 3*12000 - 10000 = 26000
        # B = 13000 + 2*0.25*(12000-10000) = 14000
        # Pw_z_t_r = (14000 - (-2000))/3 = 16000/3
        # Pw_t_z_r = (26000 - (-2000))/4 = 7000
        # Pw_t_r_z = 26000 - (-2000) - 3*14000 = -14000
        # max -> 7000
        result = WellboreStabilityCalculation.calculate_breakout_calculation_vertical_well_mohr_coulomb_analytical(
            shmax=12000, shmin=10000, pprs=5000, overburden_stress=13000, ucs=8000, fang=30, pr_sta=0.25
        )
        assert result == pytest.approx(7000.0, rel=TOLERANCE)

    def test_returns_float(self) -> None:
        result = WellboreStabilityCalculation.calculate_breakout_calculation_vertical_well_mohr_coulomb_analytical(
            shmax=12000, shmin=10000, pprs=5000, overburden_stress=13000, ucs=8000, fang=30, pr_sta=0.25
        )
        assert isinstance(result, float)

    def test_breakout_decreases_with_higher_ucs(self) -> None:
        low_ucs_breakout = WellboreStabilityCalculation.calculate_breakout_calculation_vertical_well_mohr_coulomb_analytical(
            shmax=12000, shmin=10000, pprs=5000, overburden_stress=13000, ucs=2000, fang=30, pr_sta=0.25
        )
        high_ucs_breakout = WellboreStabilityCalculation.calculate_breakout_calculation_vertical_well_mohr_coulomb_analytical(
            shmax=12000, shmin=10000, pprs=5000, overburden_stress=13000, ucs=15000, fang=30, pr_sta=0.25
        )
        assert low_ucs_breakout > high_ucs_breakout

    def test_picks_maximum_of_three_scenarios(self) -> None:
        # Manually compute the three scenarios and confirm the function returns the max
        shmax = 12000
        shmin = 10000
        pore_pressure = 5000
        overburden_stress = 13000
        ucs = 8000
        friction_angle = 30
        poisson_ratio_static = 0.25

        passive_earth_pressure_coefficient = math.tan(math.radians(45) + math.radians(friction_angle / 2)) ** 2
        cohesion_term = ucs - pore_pressure * (passive_earth_pressure_coefficient - 1)
        tangential_stress_term = 3 * shmax - shmin
        axial_stress_term = overburden_stress + 2 * poisson_ratio_static * (shmax - shmin)

        scenario_pressures = [
            (axial_stress_term - cohesion_term) / passive_earth_pressure_coefficient,
            (tangential_stress_term - cohesion_term) / (1 + passive_earth_pressure_coefficient),
            tangential_stress_term - cohesion_term - passive_earth_pressure_coefficient * axial_stress_term,
        ]
        expected_breakout = max(scenario_pressures)

        result = WellboreStabilityCalculation.calculate_breakout_calculation_vertical_well_mohr_coulomb_analytical(
            shmax=shmax,
            shmin=shmin,
            pprs=pore_pressure,
            overburden_stress=overburden_stress,
            ucs=ucs,
            fang=friction_angle,
            pr_sta=poisson_ratio_static,
        )
        assert result == pytest.approx(expected_breakout, rel=TOLERANCE)
