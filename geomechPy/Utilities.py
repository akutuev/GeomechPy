class Utilities:

	def convert_velocity(self,vel, unit_in, unit_out):
		# Conversion factors to meters per second
		conversion_to_ms = {
			"m/s": 1,
			"km/h": 1 / 3.6,
			"km/s": 1000,
			"ft/s": 0.3048,
			"cm/s": 0.01,
			"in/s": 0.0254
		}

		# Validate units
		if unit_in not in conversion_to_ms or unit_out not in conversion_to_ms:
			print('Input or Output units are not compatible')

		# Convert from source unit to m/s
		vel_ms = vel * conversion_to_ms[unit_in]

		# Convert from m/s to target unit
		vel_out = vel_ms / conversion_to_ms[unit_out]

		return vel_out
	
	def convert_pressure(self,P, unit_in, unit_out):
		# Conversion factors to pascal
		conversion_to_pascal = {
			'Pa': 1,
			'MPa': 1000000,
			'GPa': 1000000000,
			'bar': 100000,
			'psi': 6894.76,
			'Mpsi': 6894760000,
			'atm': 101325,
			'torr': 133.322
    	}


		# Validate units
		if unit_in not in conversion_to_pascal or unit_out not in conversion_to_pascal:
			print('Input or Output units are not compatible')

		# Convert from source unit to m/s
		P_pascal = P * conversion_to_pascal[unit_in]

		# Convert from m/s to target unit
		P_out = P_pascal / conversion_to_pascal[unit_out]

		return P_out
	




	def velocity2slowness(self,vp,vs,unit_in):
		"""
       Converts velocity to slowness 

        Args:
            vp: Compressional Velocity.
            vs: Shear velocity.
			Velocity Units: "m/s","km/h","km/s","ft/s","cm/s","in/s"

        Returns:
            dt: Compressional Slowness.
            ds: Shear Slowness.
			Slowness Unit: "us/ft"
        """
		if vp > 0 and vs > 0 and vp/vs > 1.414:
			dt = 304800 / self.convert_velocity(vp, unit_in, "m/s")
			ds = 304800 / self.convert_velocity(vs, unit_in, "m/s")
		elif vp <= 0 or vs <= 0:
			dt = None
			ds = None
		return dt,ds


	def slowness2velocity(self,dt,ds,unit_out):
		"""
       Converts slowness with Unit [us/ft] unit to velocity in [m/s]

        Args:
            dt: Compressional Slowness.
            ds: Shear Slowness.
			Slowness Unit: "us/ft"
		Conditions:
			dt and ds must bothe be bigger than 0
			dt must have smaller value

        Returns:
            vp: Compressional Velocity.
            vs: Shear velocity.
			Velocity Units: "m/s","km/h","km/s","ft/s","cm/s","in/s"
        """
		if dt > 0 and ds > 0 and ds/dt > 1.414:
			vp = self.convert_velocity(304800 / dt, "m/s", unit_out)
			vs = self.convert_velocity(304800 / ds, "m/s", unit_out)
		elif dt <= 0 or ds <= 0:
			vp = None
			vs = None
		return vp,vs
	

	def convert_mud_weight(value, unit_in, unit_out):
		"""
		Convert mud weight between units:
		- ppg: pounds per gallon
		- sg: specific gravity
		- kg/m3: kilograms per cubic meter
		- psi/ft: pressure gradient in psi per foot

		Parameters:
		- value: numeric value of mud weight
		- from_unit: source unit ('ppg', 'sg', 'kg/m3', 'psi/ft')
		- to_unit: target unit ('ppg', 'sg', 'kg/m3', 'psi/ft')

		Returns:
		- Converted value
		"""
		# Conversion factors to base unit (ppg)
		to_ppg = {
			'ppg': 1.0,
			'sg': 8.33,
			'kg/m3': 8.33 / 1000 * 1000,  # 1 sg = 1000 kg/m3
			'lb/ft3': 8.33 / 62.4,  # 1 sg = 62.4 lb/ft3
			'psi/ft': 1 / 0.052,
			'Pa/m': 1 / 10253.4,  # 1 ppg ≈ 10253.4 Pa/m
			'MPa/m':  0.0976  # 1 MPa/m≈0.0976 ppg
		}

		if unit_in not in to_ppg or unit_out not in to_ppg:
			raise ValueError("Supported units: 'ppg', 'sg', 'kg/m3', 'lb/ft3', 'psi/ft','Pa/m','MPa/m'")

		# Convert to ppg
		value_in_ppg = value * to_ppg[unit_in]

		# Convert from ppg to target unit
		converted_value = value_in_ppg / to_ppg[unit_out]

		return converted_value





