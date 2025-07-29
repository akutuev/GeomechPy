class Stress:
    def overburden_stress(depth, ovb_gradient = '1.05'):
        """
        Calculates overburden stress (vertical stress) from density log.

        Args:
            depth (array-like): Array of depths.
            density_log (array-like): Array of density values corresponding to depths.

        Returns:
            array-like: Array of overburden stress values.
        """
        # Assuming density_log is in g/cm^3 and depth is in meters
        # Convert density to kg/m^3 and depth to meters
        ovb_gradient = float(ovb_gradient)
        Vertical_Stress = depth*ovb_gradient 

        return Vertical_Stress


    def pore_pressure_from_gradient(depth, pore_pressure_gradient = '0.47'):
        """
        Calculates pore pressure gradient.

        Args:
            depth (array-like): Array of depths.
            pore_pressure (array-like): Array of pore pressure values.

        Returns:
            array-like: Array of pore pressure gradient values.
        """
        pore_pressure_gradient = float(pore_pressure_gradient)
        pore_pressure = depth*pore_pressure_gradient    
        return pore_pressure 

    def poroelastic_horizontal_stresses(overburden_stress, pore_pressure, poisson_ratio,youngs_modulus, biot_coefficient ='1',EX ='0.0001', EY='0.009'):
        """
        Calculates horizontal stress using Poroelastic horizontal stress equation

        Args:
            overburden_stress (array-like): Array of overburden stress values.
            pore_pressure (array-like): Array of pore pressure values.
            biot_coefficient (float): Biot's coefficient.
            poisson_ratio (float): Poisson's ratio.

        Returns:
            array-like: Array of horizontal stress values.
        """
        biot_coefficient = int(biot_coefficient)
        EX = float(EX)/1e-3
        EY = float(EY)/1e-3
        A = poisson_ratio/(1-poisson_ratio)
        B = youngs_modulus/(1-poisson_ratio*poisson_ratio)
        C = (poisson_ratio*youngs_modulus)/(1-poisson_ratio*poisson_ratio)
                
        shmin_phs = A*overburden_stress + (1 - A)*biot_coefficient*pore_pressure + B*EX + C*EY
        shmax_phs = A*overburden_stress + (1 - A)*biot_coefficient*pore_pressure + B*EY + C*EX
        return shmin_phs, shmax_phs

    def shmax_multiplier(shmin, shmax_multiplier='1.1'):
        """
        Calculates maximum horizontal stress from minimum horizontal strss using multiploer

        Args:
            shmin

        Returns:
            shmax
        """
        shmax_multiplier = float(shmax_multiplier)
        shmax = shmin*shmax_multiplier
        return shmax

    def Stress_Regime_Q_Factor(sigv,shmin,shmax):
        """
        Calculates maximum horizontal stress from minimum horizontal strss using multiploer

        Args:
            shmin

        Returns:
            shmax
        """
        shmax_shmin_ratio = shmax / shmin
                    
        if sigv >shmax and shmax >= shmin:
        Q = (shmax - shmin)/(sigv - shmin)
        elif shmin < sigv and sigv <= shmax:
            Q = 2 - (sigv - shmin)/(shmax- shmin)
        elif sigv <= shmin and shmin < shmax:
            Q = 2 + (shmin - sigv)/(shmax - sigv)
        else:
            Q = 4
        return Q

    def Horizontal_Stress_Ratio(shmax,shmin):
        """
        Calculates maximum horizontal stress from minimum horizontal strss using multiploer

        Args:
            shmin

        Returns:
            shmax
        """
        shmax_shmin_ratio = shmax / shmin
        return shmax_shmin_ratio