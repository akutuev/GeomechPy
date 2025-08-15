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
	







