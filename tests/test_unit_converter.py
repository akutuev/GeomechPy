import pytest
from geomechpy.unit_converter import UnitConverter

TOLERANCE = 1e-6

class TestConvertMudWeight:
    def test_ppg_to_sg_fresh_water(self) -> None:
        # 8.3454 ppg / 8.3454 (ppg/SG) * 1.0 = 1.0 SG
        result = UnitConverter.convert_mud_weight(8.3454, 'ppg', 'sg')
        assert result == pytest.approx(1.0, rel=TOLERANCE)

    def test_sg_to_kg_per_m3_fresh_water(self) -> None:
        # 1.0 SG / 1.0 * 1000.0 (kg/m3 per SG) = 1000.0 kg/m3
        result = UnitConverter.convert_mud_weight(1.0, 'sg', 'kg/m3')
        assert result == pytest.approx(1000.0, rel=TOLERANCE)

    def test_sg_to_ppg_drilling_mud(self) -> None:
        # 1.2 SG / 1.0 * 8.3454 (ppg/SG) = 10.01448 ppg
        result = UnitConverter.convert_mud_weight(1.2, 'sg', 'ppg')
        expected = 1.2 * 8.3454
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_kg_per_m3_to_ppg_seawater(self) -> None:
        # 1025 kg/m3 / 1000 (kg/m3 per SG) * 8.3454 (ppg/SG) = 8.5540 ppg
        result = UnitConverter.convert_mud_weight(1025.0, 'kg/m3', 'ppg')
        expected = (1025.0 / 1000.0) * 8.3454
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_g_per_cm3_to_lb_per_ft3_fresh_water(self) -> None:
        # 1.0 g/cm3 / 1.0 * 62.42796 (lb/ft3 per SG) = 62.42796 lb/ft3
        result = UnitConverter.convert_mud_weight(1.0, 'g/cm3', 'lb/ft3')
        assert result == pytest.approx(62.42796, rel=TOLERANCE)

    def test_psi_per_ft_to_sg_fresh_water_gradient(self) -> None:
        # 0.433515 psi/ft is the hydrostatic gradient of fresh water (1.0 SG)
        # 0.433515 / 0.433515 * 1.0 = 1.0 SG
        result = UnitConverter.convert_mud_weight(0.433515, 'psi/ft', 'sg')
        assert result == pytest.approx(1.0, rel=TOLERANCE)

    def test_pcf_to_kg_per_m3(self) -> None:
        # 65.0 pcf / 62.42796 (pcf/SG) * 1000 (kg/m3 per SG) = 1041.2001 kg/m3
        result = UnitConverter.convert_mud_weight(65.0, 'pcf', 'kg/m3')
        expected = (65.0 / 62.42796) * 1000.0
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_lb_gal_alias_is_identical_to_ppg(self) -> None:
        # 'lb/gal' and 'ppg' share factor 8.3454; converting between them is identity
        result = UnitConverter.convert_mud_weight(10.0, 'lb/gal', 'ppg')
        assert result == pytest.approx(10.0, rel=TOLERANCE)

    def test_pcf_alias_is_identical_to_lb_per_ft3(self) -> None:
        # 'pcf' and 'lb/ft3' share factor 62.42796; converting between them is identity
        result = UnitConverter.convert_mud_weight(65.0, 'pcf', 'lb/ft3')
        assert result == pytest.approx(65.0, rel=TOLERANCE)

    def test_same_unit_returns_input_unchanged(self) -> None:
        result = UnitConverter.convert_mud_weight(9.5, 'ppg', 'ppg')
        assert result == 9.5

    def test_round_trip_ppg_to_sg_to_ppg(self) -> None:
        # Forward and reverse conversion must recover the original value
        intermediate = UnitConverter.convert_mud_weight(10.5, 'ppg', 'sg')
        result = UnitConverter.convert_mud_weight(intermediate, 'sg', 'ppg')
        assert result == pytest.approx(10.5, rel=TOLERANCE)

    def test_invalid_output_unit_raises_value_error(self) -> None:
        with pytest.raises(ValueError):
            UnitConverter.convert_mud_weight(1.0, 'ppg', 'psi')

    def test_invalid_input_unit_raises_value_error(self) -> None:
        with pytest.raises(ValueError):
            UnitConverter.convert_mud_weight(1.0, 'invalid', 'sg')

    def test_returns_float_type(self) -> None:
        assert isinstance(UnitConverter.convert_mud_weight(1.0, 'sg', 'ppg'), float)


class TestConvertMudWeightArray:
    def test_each_element_matches_independent_scalar_result(self) -> None:
        # 8.3454 ppg = 1.0 SG; 16.6908 ppg = 2.0 SG  (both exact by definition)
        result = UnitConverter.convert_mud_weight_array([8.3454, 16.6908], 'ppg', 'sg')
        assert result[0] == pytest.approx(1.0, rel=TOLERANCE)
        assert result[1] == pytest.approx(2.0, rel=TOLERANCE)

    def test_returns_list_type(self) -> None:
        result = UnitConverter.convert_mud_weight_array([1.0, 1.2], 'sg', 'ppg')
        assert isinstance(result, list)

    def test_empty_input_returns_empty_list(self) -> None:
        assert UnitConverter.convert_mud_weight_array([], 'sg', 'ppg') == []


class TestConvertPressure:
    def test_mpa_to_psi(self) -> None:
        # 1 MPa * 1e6 Pa/MPa / 6894.757293168 Pa/psi = 145.03773773 psi
        result = UnitConverter.convert_pressure(1.0, 'mpa', 'psi')
        expected = 1_000_000.0 / 6894.757293168
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_kpa_to_bar(self) -> None:
        # 1000 kPa * 1000 Pa/kPa / 100000 Pa/bar = 10.0 bar (exact)
        result = UnitConverter.convert_pressure(1000.0, 'kpa', 'bar')
        assert result == pytest.approx(10.0, rel=TOLERANCE)

    def test_atm_to_pa(self) -> None:
        # 1 atm * 101325 Pa/atm = 101325 Pa (exact by ISO 31-4)
        result = UnitConverter.convert_pressure(1.0, 'atm', 'pa')
        assert result == pytest.approx(101325.0, rel=TOLERANCE)

    def test_psi_to_pa(self) -> None:
        # 1 psi * 6894.757293168 Pa/psi = 6894.757293168 Pa (NIST definition)
        result = UnitConverter.convert_pressure(1.0, 'psi', 'pa')
        assert result == pytest.approx(6894.757293168, rel=TOLERANCE)

    def test_760_mmhg_to_kpa_standard_atmosphere(self) -> None:
        # 760 mmHg * 133.322387415 Pa/mmHg / 1000 Pa/kPa = 101.325014 kPa
        # This is the definition of standard atmospheric pressure.
        result = UnitConverter.convert_pressure(760.0, 'mmhg', 'kpa')
        expected = (760.0 * 133.322387415) / 1000.0
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_bar_to_mbar(self) -> None:
        # 1 bar * 100000 Pa/bar / 100 Pa/mbar = 1000.0 mbar (exact)
        result = UnitConverter.convert_pressure(1.0, 'bar', 'mbar')
        assert result == pytest.approx(1000.0, rel=TOLERANCE)

    def test_gpa_to_mpa(self) -> None:
        # 1 GPa * 1e9 Pa/GPa / 1e6 Pa/MPa = 1000.0 MPa (exact)
        result = UnitConverter.convert_pressure(1.0, 'gpa', 'mpa')
        assert result == pytest.approx(1000.0, rel=TOLERANCE)

    def test_atm_to_psi(self) -> None:
        # 1 atm = 101325 Pa / 6894.757293168 Pa/psi = 14.69594878 psi
        result = UnitConverter.convert_pressure(1.0, 'atm', 'psi')
        expected = 101325.0 / 6894.757293168
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_atm_to_inhg_standard_atmosphere(self) -> None:
        # 1 atm = 101325 Pa / 3386.389 Pa/inHg = 29.92125240 inHg
        result = UnitConverter.convert_pressure(1.0, 'atm', 'inhg')
        expected = 101325.0 / 3386.389
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_same_unit_returns_input_unchanged(self) -> None:
        result = UnitConverter.convert_pressure(5000.0, 'psi', 'psi')
        assert result == 5000.0

    def test_round_trip_psi_to_mpa_to_psi(self) -> None:
        # Forward and reverse conversion must recover the original value
        intermediate = UnitConverter.convert_pressure(5000.0, 'psi', 'mpa')
        result = UnitConverter.convert_pressure(intermediate, 'mpa', 'psi')
        assert result == pytest.approx(5000.0, rel=TOLERANCE)

    def test_invalid_output_unit_raises_value_error(self) -> None:
        with pytest.raises(ValueError):
            UnitConverter.convert_pressure(1.0, 'mpa', 'ppg')

    def test_invalid_input_unit_raises_value_error(self) -> None:
        with pytest.raises(ValueError):
            UnitConverter.convert_pressure(1.0, 'invalid', 'psi')

    def test_returns_float_type(self) -> None:
        assert isinstance(UnitConverter.convert_pressure(1.0, 'mpa', 'psi'), float)


class TestConvertPressureArray:
    def test_each_element_matches_independent_scalar_result(self) -> None:
        # 1 MPa = 1e6/6894.757293168 psi; 2 MPa = twice that value
        result = UnitConverter.convert_pressure_array([1.0, 2.0], 'mpa', 'psi')
        expected_one_mpa_in_psi = 1_000_000.0 / 6894.757293168
        assert result[0] == pytest.approx(expected_one_mpa_in_psi, rel=TOLERANCE)
        assert result[1] == pytest.approx(2.0 * expected_one_mpa_in_psi, rel=TOLERANCE)

    def test_returns_list_type(self) -> None:
        result = UnitConverter.convert_pressure_array([1.0, 2.0], 'bar', 'kpa')
        assert isinstance(result, list)

    def test_empty_input_returns_empty_list(self) -> None:
        assert UnitConverter.convert_pressure_array([], 'psi', 'mpa') == []



class TestConvertLength:
    def test_ft_to_m(self) -> None:
        # 1 ft * 0.3048 m/ft = 0.3048 m (exact definition of the international foot)
        result = UnitConverter.convert_length(1.0, 'ft', 'm')
        assert result == pytest.approx(0.3048, rel=TOLERANCE)

    def test_ft_to_inches(self) -> None:
        # 1 ft * 0.3048 m/ft / 0.0254 m/in = 12.0 in (exact)
        result = UnitConverter.convert_length(1.0, 'ft', 'in')
        assert result == pytest.approx(12.0, rel=TOLERANCE)

    def test_one_mile_to_feet(self) -> None:
        # 1 mi * 1609.344 m/mi / 0.3048 m/ft = 5280.0 ft (exact)
        result = UnitConverter.convert_length(1.0, 'mi', 'ft')
        assert result == pytest.approx(5280.0, rel=TOLERANCE)

    def test_one_yard_to_feet(self) -> None:
        # 1 yd * 0.9144 m/yd / 0.3048 m/ft = 3.0 ft (exact)
        result = UnitConverter.convert_length(1.0, 'yd', 'ft')
        assert result == pytest.approx(3.0, rel=TOLERANCE)

    def test_100_cm_to_m(self) -> None:
        # 100 cm * 0.01 m/cm = 1.0 m (exact)
        result = UnitConverter.convert_length(100.0, 'cm', 'm')
        assert result == pytest.approx(1.0, rel=TOLERANCE)

    def test_1000_mm_to_m(self) -> None:
        # 1000 mm * 0.001 m/mm = 1.0 m (exact)
        result = UnitConverter.convert_length(1000.0, 'mm', 'm')
        assert result == pytest.approx(1.0, rel=TOLERANCE)

    def test_km_to_miles(self) -> None:
        # 1 km * 1000 m/km / 1609.344 m/mi = 0.62137119 mi
        result = UnitConverter.convert_length(1.0, 'km', 'mi')
        expected = 1000.0 / 1609.344
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_m_to_inches(self) -> None:
        # 1 m * 1.0 m/m / 0.0254 m/in = 39.37007874 in
        result = UnitConverter.convert_length(1.0, 'm', 'in')
        expected = 1.0 / 0.0254
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_foot_alias_is_identical_to_ft(self) -> None:
        # 'foot' and 'ft' share factor 0.3048; converting between them is identity
        result = UnitConverter.convert_length(1.0, 'foot', 'ft')
        assert result == pytest.approx(1.0, rel=TOLERANCE)

    def test_feet_alias_is_identical_to_ft(self) -> None:
        # 'feet' and 'ft' share factor 0.3048; converting between them is identity
        result = UnitConverter.convert_length(1.0, 'feet', 'ft')
        assert result == pytest.approx(1.0, rel=TOLERANCE)

    def test_same_unit_returns_input_unchanged(self) -> None:
        result = UnitConverter.convert_length(5000.0, 'ft', 'ft')
        assert result == 5000.0

    def test_round_trip_ft_to_m_to_ft(self) -> None:
        # Forward and reverse conversion must recover the original value
        intermediate = UnitConverter.convert_length(10000.0, 'ft', 'm')
        result = UnitConverter.convert_length(intermediate, 'm', 'ft')
        assert result == pytest.approx(10000.0, rel=TOLERANCE)

    def test_invalid_output_unit_raises_value_error(self) -> None:
        with pytest.raises(ValueError):
            UnitConverter.convert_length(1.0, 'ft', 'nautical_mile')

    def test_invalid_input_unit_raises_value_error(self) -> None:
        with pytest.raises(ValueError):
            UnitConverter.convert_length(1.0, 'fathom', 'm')

    def test_returns_float_type(self) -> None:
        assert isinstance(UnitConverter.convert_length(1.0, 'ft', 'm'), float)


class TestConvertLengthArray:
    def test_each_element_matches_independent_scalar_result(self) -> None:
        # 1 ft = 0.3048 m; 2 ft = 0.6096 m (exact by definition)
        result = UnitConverter.convert_length_array([1.0, 2.0], 'ft', 'm')
        assert result[0] == pytest.approx(0.3048, rel=TOLERANCE)
        assert result[1] == pytest.approx(0.6096, rel=TOLERANCE)

    def test_returns_list_type(self) -> None:
        result = UnitConverter.convert_length_array([100.0, 200.0], 'cm', 'm')
        assert isinstance(result, list)

    def test_empty_input_returns_empty_list(self) -> None:
        assert UnitConverter.convert_length_array([], 'ft', 'm') == []

class TestConvertAcousticVelocityToVelocity:
    def test_ms_to_fts(self) -> None:
        # 1000 m/s * 1.0 m/s per ms / 0.3048 m/s per fts = 3280.8399 ft/s
        result = UnitConverter.convert_acoustic(1000.0, 'ms', 'fts')
        expected = 1000.0 / 0.3048
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_mph_to_ms(self) -> None:
        # 60 mph * 0.44704 m/s per mph = 26.8224 m/s (exact SI conversion)
        result = UnitConverter.convert_acoustic(60.0, 'mph', 'ms')
        assert result == pytest.approx(26.8224, rel=TOLERANCE)

    def test_kfts_to_fts(self) -> None:
        # 1 kft/s * 304.8 m/s per kfts / 0.3048 m/s per fts = 1000.0 ft/s (exact)
        result = UnitConverter.convert_acoustic(1.0, 'kfts', 'fts')
        assert result == pytest.approx(1000.0, rel=TOLERANCE)

    def test_ftsec_alias_is_identical_to_fts(self) -> None:
        # 'ftsec' and 'fts' share factor 0.3048; converting between them is identity
        result = UnitConverter.convert_acoustic(5000.0, 'ftsec', 'fts')
        assert result == pytest.approx(5000.0, rel=TOLERANCE)


class TestConvertAcousticSlownessToSlowness:
    def test_usft_to_usm(self) -> None:
        # 305 µs/ft * 3.280839895 µs/m per µs/ft = 1000.6562 µs/m
        result = UnitConverter.convert_acoustic(305.0, 'usft', 'usm')
        expected = 305.0 * 3.280839895
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_usm_to_sm(self) -> None:
        # 100 µs/m * 1e-6 s/m per µs/m / 1.0 s/m per s/m = 1e-4 s/m
        result = UnitConverter.convert_acoustic(100.0, 'usm', 'sm')
        assert result == pytest.approx(1e-4, rel=TOLERANCE)

    def test_usfoot_alias_is_identical_to_usft(self) -> None:
        # 'usfoot' and 'usft' share factor 3.280839895e-6; converting is identity
        result = UnitConverter.convert_acoustic(70.0, 'usfoot', 'usft')
        assert result == pytest.approx(70.0, rel=TOLERANCE)

    def test_usmeter_alias_is_identical_to_usm(self) -> None:
        # 'usmeter' and 'usm' share factor 1e-6; converting is identity
        result = UnitConverter.convert_acoustic(120.0, 'usmeter', 'usm')
        assert result == pytest.approx(120.0, rel=TOLERANCE)


class TestConvertAcousticVelocityToSlowness:
    def test_ms_to_usft(self) -> None:
        # 5000 m/s → s_sm = 1/5000 = 2e-4 s/m → usft = 2e-4 / 3.280839895e-6 = 60.96 µs/ft
        result = UnitConverter.convert_acoustic(5000.0, 'ms', 'usft')
        slowness_seconds_per_meter = 1.0 / 5000.0
        expected = slowness_seconds_per_meter / 3.280839895e-6
        assert result == pytest.approx(expected, rel=TOLERANCE)

    def test_fts_to_usft_round_number(self) -> None:
        # 2000 ft/s → v_ms = 2000 * 0.3048 = 609.6 m/s
        # → s_sm = 1/609.6 s/m → usft = 1/609.6 / 3.280839895e-6 = 500.0 µs/ft (exact)
        result = UnitConverter.convert_acoustic(2000.0, 'fts', 'usft')
        assert result == pytest.approx(500.0, rel=TOLERANCE)


class TestConvertAcousticSlownessToVelocity:
    def test_usft_to_ms(self) -> None:
        # 200 µs/ft → s_sm = 200 * 3.280839895e-6 = 6.56168e-4 s/m
        # → v_ms = 1/6.56168e-4 = 1524.0 m/s  (= 5000 ft/s * 0.3048 m/s per ft/s, exact)
        result = UnitConverter.convert_acoustic(200.0, 'usft', 'ms')
        assert result == pytest.approx(1524.0, rel=TOLERANCE)

    def test_usft_to_fts_round_number(self) -> None:
        # 500 µs/ft → s_ft = 5e-4 s/ft → v = 1/5e-4 = 2000 ft/s (exact)
        result = UnitConverter.convert_acoustic(500.0, 'usft', 'fts')
        assert result == pytest.approx(2000.0, rel=TOLERANCE)


class TestConvertAcousticEdgeCasesAndErrors:
    def test_velocity_slowness_round_trip(self) -> None:
        # 70 µs/ft → ms → µs/ft: must recover the original value
        intermediate = UnitConverter.convert_acoustic(70.0, 'usft', 'ms')
        result = UnitConverter.convert_acoustic(intermediate, 'ms', 'usft')
        assert result == pytest.approx(70.0, rel=TOLERANCE)

    def test_same_unit_returns_input_unchanged(self) -> None:
        result = UnitConverter.convert_acoustic(80.0, 'usft', 'usft')
        assert result == 80.0

    def test_zero_value_raises_value_error(self) -> None:
        with pytest.raises(ValueError):
            UnitConverter.convert_acoustic(0.0, 'ms', 'fts')

    def test_negative_value_raises_value_error(self) -> None:
        with pytest.raises(ValueError):
            UnitConverter.convert_acoustic(-100.0, 'ms', 'fts')

    def test_invalid_unit_raises_value_error(self) -> None:
        with pytest.raises(ValueError):
            UnitConverter.convert_acoustic(1000.0, 'ms', 'hands_per_fortnight')

    def test_returns_float_type(self) -> None:
        assert isinstance(UnitConverter.convert_acoustic(1000.0, 'ms', 'fts'), float)


class TestConvertAcousticArray:
    def test_each_element_matches_independent_scalar_result(self) -> None:
        # 1000 m/s = 1000/0.3048 ft/s; 2000 m/s = 2000/0.3048 ft/s
        result = UnitConverter.convert_acoustic_array([1000.0, 2000.0], 'ms', 'fts')
        assert result[0] == pytest.approx(1000.0 / 0.3048, rel=TOLERANCE)
        assert result[1] == pytest.approx(2000.0 / 0.3048, rel=TOLERANCE)

    def test_returns_list_type(self) -> None:
        result = UnitConverter.convert_acoustic_array([70.0, 80.0], 'usft', 'ms')
        assert isinstance(result, list)

    def test_empty_input_returns_empty_list(self) -> None:
        assert UnitConverter.convert_acoustic_array([], 'ms', 'usft') == []
