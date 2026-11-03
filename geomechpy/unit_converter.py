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
    # Core constants relative to a 1.0 SG (Fresh Water Baseline)
    # Density / Gradient values divided by these factors = SG
    SG_CONVERSION_FACTORS = {
        'sg': 1.0,
        'g/cm3': 1.0,
        'ppg': 8.345404,
        'lb/gal': 8.345404,
        'kg/m3': 1000.0,
        'pcf': 62.42796,
        'lb/ft3': 62.42796,
        'mpa/m': 0.00980665,
        'psi/ft': 0.433515,
        'bar/m': 0.0980665,
        'kpa/m': 9.80665
    }
    # Conversion factors relative to 1.0 Pascal 
    # Value multiplied by this factor = Pascal value
    PA_CONVERSION_FACTORS = {
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
            'inhg': 3386.389
        }   
    # Conversion factors relative to 1.0 Meter (m Baseline)
    # Value multiplied by this factor = Meter value
    METER_CONVERSION_FACTORS = {
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
            'mile': 1609.344
        }    

    # Velocity scaling map to m/s baseline
    METER_PER_SECOND_CONVERSION_FACTORS = {
        'ms': 1.0,
        'fts': 0.3048,
        'ftsec': 0.3048,
        'kmh': 1.0 / 3.6,
        'mph': 0.44704,
        'kfts': 304.8
    }

    # Slowness scaling map to s/m baseline
    SECONDS_PER_METER_CONVERSION_FACTORS = {
        'sm': 1.0,
        'sft': 3.280839895,
        'usm': 1e-6,
        'usmeter': 1e-6,
        'usft': 3.280839895e-6,
        'usfoot': 3.280839895e-6
    }
    GRAVITY = 0.00980665  # 1 SG fluid down 1 meter TVD exerts 0.00980665 MPa
    FIELD_UNITS = {'ppg', 'lb/gal', 'pcf', 'lb/ft3', 'psi/ft', 'psi'}

    @staticmethod
    def convert_equivalent_mud_weight(value: float, unit_in: str, unit_out: str) -> float:
        """Convert mud weight between standard drilling-engineering units via an SG baseline.

        Equation type: Ratio scaling via fresh-water SG baseline.

        Applicable for: All standard mud-weight unit pairs used in drilling engineering.

        Reference: N/A

        Supported units: 'sg', 'g/cm3', 'ppg', 'lb/gal', 'kg/m3', 'pcf', 'lb/ft3', 'psi/ft'

        Args:
            value (float): Mud weight magnitude expressed in the input unit. Unit: unit_in
            unit_in (str): Source unit key, e.g. 'ppg', 'sg', 'kg/m3'. Unit: mud weight or pressure gradient unit
            unit_out (str): Target unit key, e.g. 'ppg', 'sg', 'kg/m3'. mud weight or pressure gradient unit

        Returns:
            converted_equivalent_mud_weight (float): Mud weight magnitude in the output unit. Unit: unit_out
        """
        if unit_in == unit_out:
            return float(value)

        if unit_in not in UnitConverter.SG_CONVERSION_FACTORS or unit_out not in UnitConverter.SG_CONVERSION_FACTORS:
            raise ValueError(f"Unsupported unit conversion from '{unit_in}' to '{unit_out}'")

        # Convert the input value to the SG baseline
        value_in_sg = value / UnitConverter.SG_CONVERSION_FACTORS[unit_in]

        # Convert the SG baseline value to the requested output unit
        converted_equivalent_mud_weight = value_in_sg * UnitConverter.SG_CONVERSION_FACTORS[unit_out]

        return float(converted_equivalent_mud_weight)

    @staticmethod
    def convert_equivalent_mud_weight_array(
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
            converted_equivalent_mud_weight (list[float]): Mud weight magnitudes in the output unit. Unit: unit_out
        """
        return [
            UnitConverter.convert_equivalent_mud_weight(value=value, unit_in=unit_in, unit_out=unit_out)
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

        if unit_in not in UnitConverter.PA_CONVERSION_FACTORS or unit_out not in UnitConverter.PA_CONVERSION_FACTORS:
            raise ValueError(f"Unsupported unit conversion from '{unit_in}' to '{unit_out}'")

        # Convert the input value to the Pascal baseline
        value_in_pascal = value * UnitConverter.PA_CONVERSION_FACTORS[unit_in]

        # Convert the Pascal baseline value to the requested output unit
        converted_pressure = value_in_pascal / UnitConverter.PA_CONVERSION_FACTORS[unit_out]

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

        if unit_in not in UnitConverter.METER_CONVERSION_FACTORS or unit_out not in UnitConverter.METER_CONVERSION_FACTORS:
            raise ValueError(f"Unsupported unit conversion from '{unit_in}' to '{unit_out}'")

        # Convert the input value to the metre baseline
        value_in_metres = value * UnitConverter.METER_CONVERSION_FACTORS[unit_in]

        # Convert the metre baseline value to the requested output unit
        converted_length = value_in_metres / UnitConverter.METER_CONVERSION_FACTORS[unit_out]

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

        # Classify input and output units by domain (velocity or slowness)
        input_is_velocity = unit_in in UnitConverter.METER_PER_SECOND_CONVERSION_FACTORS
        input_is_slowness = unit_in in UnitConverter.SECONDS_PER_METER_CONVERSION_FACTORS
        output_is_velocity = unit_out in UnitConverter.METER_PER_SECOND_CONVERSION_FACTORS
        output_is_slowness = unit_out in UnitConverter.SECONDS_PER_METER_CONVERSION_FACTORS

        if not (input_is_velocity or input_is_slowness) or not (output_is_velocity or output_is_slowness):
            raise ValueError(f"Unsupported unit conversion from '{unit_in}' to '{unit_out}'.")

        # Velocity to Velocity: scale via m/s baseline
        if input_is_velocity and output_is_velocity:
            velocity_meters_per_second = value * UnitConverter.METER_PER_SECOND_CONVERSION_FACTORS[unit_in]
            converted_acoustic = velocity_meters_per_second / UnitConverter.METER_PER_SECOND_CONVERSION_FACTORS[unit_out]

        # Slowness to Slowness: scale via s/m baseline
        elif input_is_slowness and output_is_slowness:
            slowness_seconds_per_meter = value * UnitConverter.SECONDS_PER_METER_CONVERSION_FACTORS[unit_in]
            converted_acoustic = slowness_seconds_per_meter / UnitConverter.SECONDS_PER_METER_CONVERSION_FACTORS[unit_out]

        # Velocity to Slowness: convert to m/s baseline, then invert to s/m, then scale
        elif input_is_velocity and output_is_slowness:
            velocity_meters_per_second = value * UnitConverter.METER_PER_SECOND_CONVERSION_FACTORS[unit_in]
            slowness_seconds_per_meter = 1.0 / velocity_meters_per_second
            converted_acoustic = slowness_seconds_per_meter / UnitConverter.SECONDS_PER_METER_CONVERSION_FACTORS[unit_out]

        # Slowness to Velocity: convert to s/m baseline, then invert to m/s, then scale
        else:
            slowness_seconds_per_meter = value * UnitConverter.SECONDS_PER_METER_CONVERSION_FACTORS[unit_in]
            velocity_meters_per_second = 1.0 / slowness_seconds_per_meter
            converted_acoustic = velocity_meters_per_second / UnitConverter.METER_PER_SECOND_CONVERSION_FACTORS[unit_out]

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

    @staticmethod
    def convert_from_pressure_to_equivalent_mud_weight(pressure_value: float, unit_in: str, unit_out: str, tvd: float) -> float:
        """
        Convert an absolute fluid pressure (+ depth) into an equivalent mud weight or gradient.
        The 'tvd' unit is determined by the input pressure unit (psi = feet, others = meters).

        Args:
            pressure_value (float): Geomechanical property in pressure unit e.g. Pore Pressure, Stress Magnitudes. 
            tvd (float): True vertical depth. 
            unit_in (str): SI units: [mpa,bar,kpa] Field Unit: dimensionless

            Note: If unit_in is SI unit --> TVD is automatically identified as meters 
                  if unit_in is Field unit --> TVD is automatically identified as feet
            unit_out (str): equivalent mud weight units: [sg,g/cm3,ppg,lb/gal,kg/m3,pcf,lb/ft3,mpa/m,psi/ft,bar/m,kpa/m] Unit: dimensionless
            
        Returns:
           equivalent_mud_weight in the unit provided in unit_out
        """
        
        if tvd <= 0:
            raise ValueError("True Vertical Depth (tvd) must be greater than zero.")
        if unit_out not in UnitConverter.SG_CONVERSION_FACTORS:
            raise ValueError(f"Unsupported output density/gradient unit: '{unit_out}'")

        # 1. Resolve TVD to meters based on the input pressure system
        tvd_m = tvd * 0.3048 if unit_in in UnitConverter.FIELD_UNITS else float(tvd)

        # 2. Standardize input pressure value to baseline pa
        if unit_in == 'pa':
            pa_pressure = pressure_value
        elif unit_in == 'mpa':
            pa_pressure = pressure_value / UnitConverter.convert_pressure(unit_in)
        elif unit_in == 'psi':
            pa_pressure = pressure_value / UnitConverter.convert_pressure(unit_in)
        elif unit_in == 'bar':
            pa_pressure = pressure_value / UnitConverter.convert_pressure(unit_in)
        elif unit_in == 'kpa':
            pa_pressure = pressure_value / UnitConverter.convert_pressure(unit_in)
        else:
            raise ValueError(f"Unsupported input pressure unit: '{unit_in}'")

        # 3. Calculate baseline Specific Gravity (SG)
        sg_value = pa_pressure * UnitConverter.convert_pressure('mpa') / (UnitConverter.GRAVITY * tvd_m)

        # 4. Scale baseline SG out to the target density/gradient unit
        converted_equivalent_mud_weight =sg_value * UnitConverter.SG_CONVERSION_FACTORS[unit_out]

        return float(converted_equivalent_mud_weight) 

    @staticmethod
    def from_equivalent_mud_weight_to_pressure(mw_value: float, unit_in: str, unit_out: str, tvd: float) -> float:
        """
        Convert a fluid density/mud weight (+ depth) into an absolute hydrostatic pressure.
        The 'tvd' unit is determined by the output pressure unit (psi = feet, others = meters).
        """
        unit_in = unit_in.lower().strip()
        unit_out = unit_out.lower().strip()

        if tvd <= 0:
            raise ValueError("True Vertical Depth (tvd) must be greater than zero.")
        if unit_in not in UnitConverter.SG_CONVERSION_FACTORS:
            raise ValueError(f"Unsupported input density/gradient unit: '{unit_in}'")

        # 1. Resolve TVD to meters based on the output pressure system
        tvd_m = tvd * 0.3048 if unit_out in UnitConverter.FIELD_UNITS else float(tvd)

        # 2. Normalize input density/gradient to baseline Specific Gravity (SG)
        sg_value = mw_value / UnitConverter.SG_CONVERSION_FACTORS[unit_in]

        # 3. Compute absolute baseline hydrostatic pressure in Pa
        target_pa = sg_value * UnitConverter.GRAVITY * tvd_m * UnitConverter.convert_pressure('mpa')

        # 4. Convert baseline pa pressure to the desired output pressure scale
        if unit_out == 'pa':
            converted_pressure =  target_pa
        elif unit_out == 'psi':
            converted_pressure =   target_pa / UnitConverter.convert_pressure(unit_out)
        elif unit_out == 'bar':
            converted_pressure =   target_pa / UnitConverter.convert_pressure(unit_out)
        elif unit_out == 'kpa':
            converted_pressure =   target_pa / UnitConverter.convert_pressure(unit_out)
        elif unit_out == 'mpa':
            converted_pressure =   target_pa / UnitConverter.convert_pressure(unit_out)
        else:
            raise ValueError(f"Unsupported output pressure unit: '{unit_out}'")    

        return float(converted_pressure)