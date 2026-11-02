from dataclasses import dataclass


@dataclass(frozen=True)
class UnitConverter:
    """Convert physical quantities between standard industry and SI units.

    Supports mud weight, pressure, length, and acoustic (velocity/slowness) unit
    conversions used in petroleum geomechanics workflows. All conversions route
    through a well-defined SI baseline value (SG, Pascal, metre, m/s or s/m) and
    back, ensuring a consistent and extendable conversion pattern.

    Reference: N/A
    """

    @staticmethod
    def convert_mud_weight(value: float, unit_in: str, unit_out: str) -> float:
        """Convert mud weight between standard drilling-engineering units via an SG baseline.

        Equation type: Ratio scaling via fresh-water SG baseline.

        Applicable for: All standard mud-weight unit pairs used in drilling engineering.

        Reference: N/A

        Supported units: 'sg', 'g/cm3', 'ppg', 'lb/gal', 'kg/m3', 'pcf', 'lb/ft3', 'psi/ft'

        Args:
            value (float): Mud weight magnitude expressed in the input unit. Unit: unit_in
            unit_in (str): Source unit key, e.g. 'ppg', 'sg', 'kg/m3'. Unit: dimensionless
            unit_out (str): Target unit key, e.g. 'ppg', 'sg', 'kg/m3'. Unit: dimensionless

        Returns:
            converted_mud_weight (float): Mud weight magnitude in the output unit. Unit: unit_out
        """
        if unit_in == unit_out:
            return float(value)

        # Conversion factors relative to 1.0 SG (fresh-water baseline).
        # value_in_unit / factor = equivalent SG value.
        sg_conversion_factors = {
            'sg': 1.0,
            'g/cm3': 1.0,
            'ppg': 8.3454,
            'lb/gal': 8.3454,
            'kg/m3': 1000.0,
            'pcf': 62.42796,
            'lb/ft3': 62.42796,
            'psi/ft': 0.433515,
        }

        if unit_in not in sg_conversion_factors or unit_out not in sg_conversion_factors:
            raise ValueError(f"Unsupported unit conversion from '{unit_in}' to '{unit_out}'")

        # Convert the input value to the SG baseline
        value_in_sg = value / sg_conversion_factors[unit_in]

        # Convert the SG baseline value to the requested output unit
        converted_mud_weight = value_in_sg * sg_conversion_factors[unit_out]

        return float(converted_mud_weight)

    @staticmethod
    def convert_mud_weight_array(
        values: list[float],
        unit_in: str,
        unit_out: str,
    ) -> list[float]:
        """Convert a list of mud weight values between standard drilling-engineering units.

        Args:
            values (list[float]): Mud weight magnitudes in the input unit. Unit: unit_in
            unit_in (str): Source unit key, e.g. 'ppg', 'sg', 'kg/m3'. Unit: dimensionless
            unit_out (str): Target unit key, e.g. 'ppg', 'sg', 'kg/m3'. Unit: dimensionless

        Returns:
            converted_mud_weights (list[float]): Mud weight magnitudes in the output unit. Unit: unit_out
        """
        return [
            UnitConverter.convert_mud_weight(value=value, unit_in=unit_in, unit_out=unit_out)
            for value in values
        ]

    @staticmethod
    def convert_pressure(value: float, unit_in: str, unit_out: str) -> float:
        """Convert pressure between standard SI, metric, and oilfield units via a Pascal baseline.

        Equation type: Ratio scaling via Pascal baseline.

        Applicable for: All standard pressure unit pairs used in petroleum engineering
        and geomechanics.

        Reference: N/A

        Supported units: 'pa', 'kpa', 'mpa', 'gpa', 'psi', 'mpsi', 'bar', 'mbar',
                         'atm', 'torr', 'mmhg', 'inhg'

        Args:
            value (float): Pressure magnitude expressed in the input unit. Unit: unit_in
            unit_in (str): Source unit key, e.g. 'psi', 'mpa', 'bar'. Unit: dimensionless
            unit_out (str): Target unit key, e.g. 'psi', 'mpa', 'bar'. Unit: dimensionless

        Returns:
            converted_pressure (float): Pressure magnitude in the output unit. Unit: unit_out
        """
        if unit_in == unit_out:
            return float(value)

        # Conversion factors relative to 1.0 Pascal (Pa baseline).
        # value_in_unit * factor = equivalent Pascal value.
        pa_conversion_factors = {
            'pa': 1.0,
            'kpa': 1000.0,
            'mpa': 1000000.0,
            'gpa': 1000000000.0,
            'psi': 6894.757293168,
            'mpsi': 6894757293.168,
            'bar': 100000.0,
            'mbar': 100.0,
            'atm': 101325.0,
            'torr': 133.322368421,
            'mmhg': 133.322387415,
            'inhg': 3386.389,
        }

        if unit_in not in pa_conversion_factors or unit_out not in pa_conversion_factors:
            raise ValueError(f"Unsupported unit conversion from '{unit_in}' to '{unit_out}'")

        # Convert the input value to the Pascal baseline
        value_in_pascal = value * pa_conversion_factors[unit_in]

        # Convert the Pascal baseline value to the requested output unit
        converted_pressure = value_in_pascal / pa_conversion_factors[unit_out]

        return float(converted_pressure)

    @staticmethod
    def convert_pressure_array(
        values: list[float],
        unit_in: str,
        unit_out: str,
    ) -> list[float]:
        """Convert a list of pressure values between standard SI, metric, and oilfield units.

        Args:
            values (list[float]): Pressure magnitudes in the input unit. Unit: unit_in
            unit_in (str): Source unit key, e.g. 'psi', 'mpa', 'bar'. Unit: dimensionless
            unit_out (str): Target unit key, e.g. 'psi', 'mpa', 'bar'. Unit: dimensionless

        Returns:
            converted_pressures (list[float]): Pressure magnitudes in the output unit. Unit: unit_out
        """
        return [
            UnitConverter.convert_pressure(value=value, unit_in=unit_in, unit_out=unit_out)
            for value in values
        ]

    @staticmethod
    def convert_length(value: float, unit_in: str, unit_out: str) -> float:
        """Convert length between metric, imperial, and oilfield units via a metre baseline.

        Equation type: Ratio scaling via metre baseline.

        Applicable for: All standard length unit pairs used in petroleum engineering
        and geomechanics.

        Reference: N/A

        Supported units: 'm', 'meter', 'meters', 'mm', 'cm', 'km', 'in', 'inch', 'inches',
                         'ft', 'foot', 'feet', 'yd', 'yard', 'mi', 'mile'

        Args:
            value (float): Length magnitude expressed in the input unit. Unit: unit_in
            unit_in (str): Source unit key, e.g. 'ft', 'm', 'km'. Unit: dimensionless
            unit_out (str): Target unit key, e.g. 'ft', 'm', 'km'. Unit: dimensionless

        Returns:
            converted_length (float): Length magnitude in the output unit. Unit: unit_out
        """
        if unit_in == unit_out:
            return float(value)

        # Conversion factors relative to 1.0 metre (m baseline).
        # value_in_unit * factor = equivalent metre value.
        m_conversion_factors = {
            'm': 1.0,
            'meter': 1.0,
            'meters': 1.0,
            'mm': 0.001,
            'cm': 0.01,
            'km': 1000.0,
            'in': 0.0254,
            'inch': 0.0254,
            'inches': 0.0254,
            'ft': 0.3048,
            'foot': 0.3048,
            'feet': 0.3048,
            'yd': 0.9144,
            'yard': 0.9144,
            'mi': 1609.344,
            'mile': 1609.344,
        }

        if unit_in not in m_conversion_factors or unit_out not in m_conversion_factors:
            raise ValueError(f"Unsupported unit conversion from '{unit_in}' to '{unit_out}'")

        # Convert the input value to the metre baseline
        value_in_metres = value * m_conversion_factors[unit_in]

        # Convert the metre baseline value to the requested output unit
        converted_length = value_in_metres / m_conversion_factors[unit_out]

        return float(converted_length)

    @staticmethod
    def convert_length_array(
        values: list[float],
        unit_in: str,
        unit_out: str,
    ) -> list[float]:
        """Convert a list of length values between metric, imperial, and oilfield units.

        Args:
            values (list[float]): Length magnitudes in the input unit. Unit: unit_in
            unit_in (str): Source unit key, e.g. 'ft', 'm', 'km'. Unit: dimensionless
            unit_out (str): Target unit key, e.g. 'ft', 'm', 'km'. Unit: dimensionless

        Returns:
            converted_lengths (list[float]): Length magnitudes in the output unit. Unit: unit_out
        """
        return [
            UnitConverter.convert_length(value=value, unit_in=unit_in, unit_out=unit_out)
            for value in values
        ]

    @staticmethod
    def convert_acoustic(value: float, unit_in: str, unit_out: str) -> float:
        """Convert acoustic velocity or slowness between standard oilfield and SI units.

        Handles same-domain conversions (velocity-to-velocity, slowness-to-slowness) and
        cross-domain conversions (velocity-to-slowness and slowness-to-velocity) using
        the reciprocal relationship: slowness [s/m] = 1 / velocity [m/s].

        Equation type: Ratio scaling with reciprocal cross-domain inversion.

        Applicable for: Compressional and shear sonic log data in petroleum geomechanics
        workflows.

        Reference: N/A

        Supported velocity units: 'ms' (m/s), 'fts' (ft/s), 'ftsec' (ft/s alias),
                                   'kmh' (km/h), 'mph' (miles/h), 'kfts' (kft/s)
        Supported slowness units: 'sm' (s/m), 'sft' (s/ft), 'usm' (µs/m),
                                   'usmeter' (µs/m alias), 'usft' (µs/ft),
                                   'usfoot' (µs/ft alias)

        Args:
            value (float): Acoustic value expressed in the input unit. Must be > 0. Unit: unit_in
            unit_in (str): Source unit key, e.g. 'ms', 'usft'. Unit: dimensionless
            unit_out (str): Target unit key, e.g. 'ms', 'usft'. Unit: dimensionless

        Returns:
            converted_acoustic (float): Acoustic value in the output unit. Unit: unit_out
        """
        if value <= 0:
            raise ValueError("Acoustic values must be strictly greater than zero.")

        if unit_in == unit_out:
            return float(value)

        # Conversion factors: metres per second per one unit of the given velocity key.
        velocity_unit_to_meters_per_second = {
            'ms': 1.0,
            'fts': 0.3048,
            'ftsec': 0.3048,
            'kmh': 1.0 / 3.6,
            'mph': 0.44704,
            'kfts': 304.8,
        }

        # Conversion factors: seconds per metre per one unit of the given slowness key.
        slowness_unit_to_seconds_per_meter = {
            'sm': 1.0,
            'sft': 3.280839895,
            'usm': 1e-6,
            'usmeter': 1e-6,
            'usft': 3.280839895e-6,
            'usfoot': 3.280839895e-6,
        }

        # Classify input and output units by domain (velocity or slowness)
        input_is_velocity = unit_in in velocity_unit_to_meters_per_second
        input_is_slowness = unit_in in slowness_unit_to_seconds_per_meter
        output_is_velocity = unit_out in velocity_unit_to_meters_per_second
        output_is_slowness = unit_out in slowness_unit_to_seconds_per_meter

        if not (input_is_velocity or input_is_slowness) or not (output_is_velocity or output_is_slowness):
            raise ValueError(f"Unsupported unit conversion from '{unit_in}' to '{unit_out}'.")

        # Velocity to Velocity: scale via m/s baseline
        if input_is_velocity and output_is_velocity:
            velocity_meters_per_second = value * velocity_unit_to_meters_per_second[unit_in]
            converted_acoustic = velocity_meters_per_second / velocity_unit_to_meters_per_second[unit_out]

        # Slowness to Slowness: scale via s/m baseline
        elif input_is_slowness and output_is_slowness:
            slowness_seconds_per_meter = value * slowness_unit_to_seconds_per_meter[unit_in]
            converted_acoustic = slowness_seconds_per_meter / slowness_unit_to_seconds_per_meter[unit_out]

        # Velocity to Slowness: convert to m/s baseline, then invert to s/m, then scale
        elif input_is_velocity and output_is_slowness:
            velocity_meters_per_second = value * velocity_unit_to_meters_per_second[unit_in]
            slowness_seconds_per_meter = 1.0 / velocity_meters_per_second
            converted_acoustic = slowness_seconds_per_meter / slowness_unit_to_seconds_per_meter[unit_out]

        # Slowness to Velocity: convert to s/m baseline, then invert to m/s, then scale
        else:
            slowness_seconds_per_meter = value * slowness_unit_to_seconds_per_meter[unit_in]
            velocity_meters_per_second = 1.0 / slowness_seconds_per_meter
            converted_acoustic = velocity_meters_per_second / velocity_unit_to_meters_per_second[unit_out]

        return float(converted_acoustic)

    @staticmethod
    def convert_acoustic_array(
        values: list[float],
        unit_in: str,
        unit_out: str,
    ) -> list[float]:
        """Convert a list of acoustic velocity or slowness values between standard units.

        Args:
            values (list[float]): Acoustic values in the input unit. All values must be > 0. Unit: unit_in
            unit_in (str): Source unit key, e.g. 'ms', 'usft'. Unit: dimensionless
            unit_out (str): Target unit key, e.g. 'ms', 'usft'. Unit: dimensionless

        Returns:
            converted_acoustics (list[float]): Acoustic values in the output unit. Unit: unit_out
        """
        return [
            UnitConverter.convert_acoustic(value=value, unit_in=unit_in, unit_out=unit_out)
            for value in values
        ]
