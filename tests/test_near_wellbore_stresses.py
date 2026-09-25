import math

import numpy as np
import pytest

from geomechpy.near_wellbore_stresses import (
    BoreholeGeneralStresses,
    BoreholeWallStresses,
    NearWellboreStressesCalculation,
    PrincipalStressesBoreholewall,
)

TOLERANCE = 1e-6


class TestKirschBoreholeWallStresses:
    def test_returns_borehole_wall_stresses_dataclass(self) -> None:
        result = NearWellboreStressesCalculation.calculate_kirsch_borehole_wall_stresses(
            shmin=10000,
            shmax=12000,
            svert=13000,
            pore_pressure=5000,
            shmax_azimuth=0,
            mud_pressure=5000,
            theta=np.linspace(0, 360, 37),
            poisson_ratio_static=0.25,
            borehole_deviation=0,
            borehole_azimuth=0,
        )
        assert isinstance(result, BoreholeWallStresses)

    def test_output_arrays_have_correct_length(self) -> None:
        theta = np.linspace(0, 360, 37)
        result = NearWellboreStressesCalculation.calculate_kirsch_borehole_wall_stresses(
            shmin=10000,
            shmax=12000,
            svert=13000,
            pore_pressure=5000,
            shmax_azimuth=0,
            mud_pressure=5000,
            theta=theta,
            poisson_ratio_static=0.25,
            borehole_deviation=0,
            borehole_azimuth=0,
        )
        assert len(result.sigma_boreholewall_rr) == len(theta)
        assert len(result.sigma_boreholewall_tt) == len(theta)
        assert len(result.sigma_boreholewall_zz) == len(theta)
        assert len(result.sigma_boreholewall_tz) == len(theta)
        assert len(result.sigma_boreholewall_rt) == len(theta)
        assert len(result.sigma_boreholewall_rz) == len(theta)

    def test_sigma_rr_equals_mud_minus_pore_pressure(self) -> None:
        # sigma_rr is the borehole wall boundary condition; it is constant around the circumference
        mud_pressure = 6000
        pore_pressure = 5000
        theta = np.linspace(0, 360, 37)
        result = NearWellboreStressesCalculation.calculate_kirsch_borehole_wall_stresses(
            shmin=10000,
            shmax=12000,
            svert=13000,
            pore_pressure=pore_pressure,
            shmax_azimuth=0,
            mud_pressure=mud_pressure,
            theta=theta,
            poisson_ratio_static=0.25,
            borehole_deviation=0,
            borehole_azimuth=0,
        )
        assert result.sigma_boreholewall_rr == pytest.approx(np.full(len(theta), mud_pressure - pore_pressure), rel=TOLERANCE)

    def test_sigma_rt_is_zero_on_borehole_wall(self) -> None:
        theta = np.linspace(0, 360, 37)
        result = NearWellboreStressesCalculation.calculate_kirsch_borehole_wall_stresses(
            shmin=10000,
            shmax=12000,
            svert=13000,
            pore_pressure=5000,
            shmax_azimuth=0,
            mud_pressure=5000,
            theta=theta,
            poisson_ratio_static=0.25,
            borehole_deviation=0,
            borehole_azimuth=0,
        )
        assert result.sigma_boreholewall_rt == pytest.approx(np.zeros(len(theta)), abs=TOLERANCE)

    def test_sigma_rz_is_zero_on_borehole_wall(self) -> None:
        theta = np.linspace(0, 360, 37)
        result = NearWellboreStressesCalculation.calculate_kirsch_borehole_wall_stresses(
            shmin=10000,
            shmax=12000,
            svert=13000,
            pore_pressure=5000,
            shmax_azimuth=0,
            mud_pressure=5000,
            theta=theta,
            poisson_ratio_static=0.25,
            borehole_deviation=0,
            borehole_azimuth=0,
        )
        assert result.sigma_boreholewall_rz == pytest.approx(np.zeros(len(theta)), abs=TOLERANCE)

    def test_known_value_sigma_tt_vertical_well_balanced_mud(self) -> None:
        # Vertical borehole, shmax aligned North, balanced mud (mud=pprs), no shear (sxy=sxz=syz=0)
        # At theta=0 (top of hole): sigma_tt = 3*shmin - shmax = 3*10000 - 12000 = 18000
        # At theta=90:              sigma_tt = 3*shmax - shmin = 3*12000 - 10000 = 26000
        result = NearWellboreStressesCalculation.calculate_kirsch_borehole_wall_stresses(
            shmin=10000,
            shmax=12000,
            svert=13000,
            pore_pressure=5000,
            shmax_azimuth=0,
            mud_pressure=5000,
            theta=np.array([0.0, 90.0]),
            poisson_ratio_static=0.25,
            borehole_deviation=0,
            borehole_azimuth=0,
        )
        assert result.sigma_boreholewall_tt[0] == pytest.approx(18000, rel=TOLERANCE)
        assert result.sigma_boreholewall_tt[1] == pytest.approx(26000, rel=TOLERANCE)

    def test_known_value_sigma_zz_vertical_well_balanced_mud(self) -> None:
        # sigma_zz at theta=0: svert - pr*(2*(shmax-shmin)) = 13000 - 0.25*4000 = 12000
        # sigma_zz at theta=90: svert + pr*(2*(shmax-shmin)) = 13000 + 0.25*4000 = 14000
        result = NearWellboreStressesCalculation.calculate_kirsch_borehole_wall_stresses(
            shmin=10000,
            shmax=12000,
            svert=13000,
            pore_pressure=5000,
            shmax_azimuth=0,
            mud_pressure=5000,
            theta=np.array([0.0, 90.0]),
            poisson_ratio_static=0.25,
            borehole_deviation=0,
            borehole_azimuth=0,
        )
        assert result.sigma_boreholewall_zz[0] == pytest.approx(12000, rel=TOLERANCE)
        assert result.sigma_boreholewall_zz[1] == pytest.approx(14000, rel=TOLERANCE)


class TestKirschBoreholeGeneralStresses:
    def test_returns_borehole_general_stresses_dataclass(self) -> None:
        bitsize = 8.5
        result = NearWellboreStressesCalculation.calculate_kirsch_borehole_stresses_general(
            shmin=10000,
            shmax=12000,
            svert=13000,
            pore_pressure=5000,
            shmax_azimuth=0,
            mud_pressure=5000,
            theta=np.linspace(0, 360, 37),
            radius=np.full(37, bitsize / 2),
            poisson_ratio_static=0.25,
            borehole_deviation=0,
            borehole_azimuth=0,
            bitsize=bitsize,
        )
        assert isinstance(result, BoreholeGeneralStresses)

    def test_output_arrays_have_correct_length(self) -> None:
        bitsize = 8.5
        theta = np.linspace(0, 360, 37)
        result = NearWellboreStressesCalculation.calculate_kirsch_borehole_stresses_general(
            shmin=10000,
            shmax=12000,
            svert=13000,
            pore_pressure=5000,
            shmax_azimuth=0,
            mud_pressure=5000,
            theta=theta,
            radius=np.full(len(theta), bitsize / 2),
            poisson_ratio_static=0.25,
            borehole_deviation=0,
            borehole_azimuth=0,
            bitsize=bitsize,
        )
        assert len(result.sigma_general_rr) == len(theta)
        assert len(result.sigma_general_tt) == len(theta)
        assert len(result.sigma_general_zz) == len(theta)
        assert len(result.sigma_general_tz) == len(theta)
        assert len(result.sigma_general_rt) == len(theta)
        assert len(result.sigma_general_rz) == len(theta)

    def test_general_solution_at_borehole_wall_matches_wall_solution(self) -> None:
        # At r = bh_radius the general Kirsch solution must reduce to the borehole-wall solution
        shmin, shmax, svert = 10000, 12000, 13000
        pore_pressure, mud_pressure = 5000, 6000
        shmax_azimuth, poisson_ratio_static = 45, 0.25
        borehole_deviation, borehole_azimuth = 30, 60
        bitsize = 8.5
        theta = np.linspace(0, 360, 37)

        wall = NearWellboreStressesCalculation.calculate_kirsch_borehole_wall_stresses(
            shmin=shmin,
            shmax=shmax,
            svert=svert,
            pore_pressure=pore_pressure,
            shmax_azimuth=shmax_azimuth,
            mud_pressure=mud_pressure,
            theta=theta,
            poisson_ratio_static=poisson_ratio_static,
            borehole_deviation=borehole_deviation,
            borehole_azimuth=borehole_azimuth,
        )
        general = NearWellboreStressesCalculation.calculate_kirsch_borehole_stresses_general(
            shmin=shmin,
            shmax=shmax,
            svert=svert,
            pore_pressure=pore_pressure,
            shmax_azimuth=shmax_azimuth,
            mud_pressure=mud_pressure,
            theta=theta,
            radius=np.full(len(theta), bitsize / 2),
            poisson_ratio_static=poisson_ratio_static,
            borehole_deviation=borehole_deviation,
            borehole_azimuth=borehole_azimuth,
            bitsize=bitsize,
        )
        assert general.sigma_general_rr == pytest.approx(wall.sigma_boreholewall_rr, rel=TOLERANCE)
        assert general.sigma_general_tt == pytest.approx(wall.sigma_boreholewall_tt, rel=TOLERANCE)
        assert general.sigma_general_zz == pytest.approx(wall.sigma_boreholewall_zz, rel=TOLERANCE)
        assert general.sigma_general_tz == pytest.approx(wall.sigma_boreholewall_tz, rel=TOLERANCE)
        assert general.sigma_general_rt == pytest.approx(wall.sigma_boreholewall_rt, abs=TOLERANCE)
        assert general.sigma_general_rz == pytest.approx(wall.sigma_boreholewall_rz, abs=TOLERANCE)

    def test_sigma_rr_at_borehole_wall_equals_mud_minus_pore_pressure(self) -> None:
        # At r = bh_radius the radial stress boundary condition must equal mud - pore pressure
        bitsize = 8.5
        mud_pressure = 6000
        pore_pressure = 5000
        theta = np.linspace(0, 360, 37)
        result = NearWellboreStressesCalculation.calculate_kirsch_borehole_stresses_general(
            shmin=10000,
            shmax=12000,
            svert=13000,
            pore_pressure=pore_pressure,
            shmax_azimuth=0,
            mud_pressure=mud_pressure,
            theta=theta,
            radius=np.full(len(theta), bitsize / 2),
            poisson_ratio_static=0.25,
            borehole_deviation=0,
            borehole_azimuth=0,
            bitsize=bitsize,
        )
        assert result.sigma_general_rr == pytest.approx(np.full(len(theta), mud_pressure - pore_pressure), rel=TOLERANCE)

    def test_sigma_rt_and_rz_are_zero_at_borehole_wall(self) -> None:
        # At r = bh_radius the shear components sigma_rt and sigma_rz must vanish
        bitsize = 8.5
        theta = np.linspace(0, 360, 37)
        result = NearWellboreStressesCalculation.calculate_kirsch_borehole_stresses_general(
            shmin=10000,
            shmax=12000,
            svert=13000,
            pore_pressure=5000,
            shmax_azimuth=0,
            mud_pressure=5000,
            theta=theta,
            radius=np.full(len(theta), bitsize / 2),
            poisson_ratio_static=0.25,
            borehole_deviation=0,
            borehole_azimuth=0,
            bitsize=bitsize,
        )
        assert result.sigma_general_rt == pytest.approx(np.zeros(len(theta)), abs=TOLERANCE)
        assert result.sigma_general_rz == pytest.approx(np.zeros(len(theta)), abs=TOLERANCE)

    def test_stresses_converge_to_far_field_at_large_radius(self) -> None:
        # At r >> bh_radius all terms with a/r vanish and the stress reduces to far-field values.
        # For a vertical well with no shear (sxy=sxz=syz=0), balanced mud, and theta=0:
        #   sigma_rr -> 0.5*(sx0+sy0) + 0.5*(sx0-sy0)*cos(0) = sx0 = shmax
        #   sigma_tt -> 0.5*(sx0+sy0) - 0.5*(sx0-sy0)*cos(0) = sy0 = shmin
        bitsize = 8.5
        far_radius = 1e6  # effectively infinite
        theta = np.array([0.0])
        result = NearWellboreStressesCalculation.calculate_kirsch_borehole_stresses_general(
            shmin=10000,
            shmax=12000,
            svert=13000,
            pore_pressure=5000,
            shmax_azimuth=0,
            mud_pressure=5000,
            theta=theta,
            radius=np.full(len(theta), far_radius),
            poisson_ratio_static=0.25,
            borehole_deviation=0,
            borehole_azimuth=0,
            bitsize=bitsize,
        )
        assert result.sigma_general_rr[0] == pytest.approx(12000, rel=1e-3)
        assert result.sigma_general_tt[0] == pytest.approx(10000, rel=1e-3)


class TestPrincipalStressesAnalytical:
    def test_returns_principal_stresses_dataclass(self) -> None:
        result = NearWellboreStressesCalculation.calculate_principal_stresses_analytical_boreholewall(
            sigma_boreholewall_rr=np.array([5000.0]),
            sigma_boreholewall_tt=np.array([10000.0]),
            sigma_boreholewall_zz=np.array([8000.0]),
            sigma_boreholewall_tz=np.array([0.0]),
        )
        assert isinstance(result, PrincipalStressesBoreholewall)

    def test_sigma_1_is_greater_than_or_equal_to_sigma_2(self) -> None:
        result = NearWellboreStressesCalculation.calculate_principal_stresses_analytical_boreholewall(
            sigma_boreholewall_rr=np.array([5000.0, 5000.0]),
            sigma_boreholewall_tt=np.array([10000.0, 8000.0]),
            sigma_boreholewall_zz=np.array([8000.0, 10000.0]),
            sigma_boreholewall_tz=np.array([1000.0, 1500.0]),
        )
        assert all(result.sigma_boreholewall_1 >= result.sigma_boreholewall_2)

    def test_sigma_2_is_greater_than_or_equal_to_sigma_3(self) -> None:
        result = NearWellboreStressesCalculation.calculate_principal_stresses_analytical_boreholewall(
            sigma_boreholewall_rr=np.array([5000.0, 5000.0]),
            sigma_boreholewall_tt=np.array([10000.0, 8000.0]),
            sigma_boreholewall_zz=np.array([8000.0, 10000.0]),
            sigma_boreholewall_tz=np.array([1000.0, 1500.0]),
        )
        assert all(result.sigma_boreholewall_2 >= result.sigma_boreholewall_3)

    def test_no_shear_stress_principal_stresses_equal_inputs(self) -> None:
        # When sigma_tz=0 the principal stresses reduce to the input components directly
        result = NearWellboreStressesCalculation.calculate_principal_stresses_analytical_boreholewall(
            sigma_boreholewall_rr=np.array([5000.0]),
            sigma_boreholewall_tt=np.array([10000.0]),
            sigma_boreholewall_zz=np.array([8000.0]),
            sigma_boreholewall_tz=np.array([0.0]),
        )
        assert result.sigma_boreholewall_1[0] == pytest.approx(10000.0, rel=TOLERANCE)
        assert result.sigma_boreholewall_2[0] == pytest.approx(8000.0, rel=TOLERANCE)
        assert result.sigma_boreholewall_3[0] == pytest.approx(5000.0, rel=TOLERANCE)

    def test_known_values(self) -> None:
        # sigma_tt=10000, sigma_zz=8000, sigma_tz=1000, sigma_rr=5000
        # sigma_a = 9000 + 0.5*sqrt((2000)^2 + 4*1000^2) = 9000 + 1000*sqrt(2) ≈ 10414.214
        # sigma_b = 9000 - 1000*sqrt(2) ≈ 7585.786
        # Sorted descending: sigma_1=10414.2, sigma_2=7585.8, sigma_3=5000
        # theta_tortuosity = 0.5 * degrees(arctan(2000/2000)) = 0.5 * 45 = 22.5
        result = NearWellboreStressesCalculation.calculate_principal_stresses_analytical_boreholewall(
            sigma_boreholewall_rr=np.array([5000.0]),
            sigma_boreholewall_tt=np.array([10000.0]),
            sigma_boreholewall_zz=np.array([8000.0]),
            sigma_boreholewall_tz=np.array([1000.0]),
        )
        assert result.sigma_boreholewall_1[0] == pytest.approx(9000 + 1000 * math.sqrt(2), rel=TOLERANCE)
        assert result.sigma_boreholewall_2[0] == pytest.approx(9000 - 1000 * math.sqrt(2), rel=TOLERANCE)
        assert result.sigma_boreholewall_3[0] == pytest.approx(5000.0, rel=TOLERANCE)
        assert result.theta_boreholewall_tortuosity[0] == pytest.approx(22.5, rel=TOLERANCE)
