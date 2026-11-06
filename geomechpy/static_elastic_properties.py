import math

from geomechpy.unit_converter import UnitConverter

class StaticElasticPropertiesPowerLawConverter:
    """Convert dynamic elastic properties to static elastic properties using Power Law.

    ## Content:
    ### Power Law Moduli Conversions (Es = a * Ed^b or Es = a * Vp^b)
    (1) dyn2sta_yme_fuller: Fuller Sandstone/Shale Model: Citation: Fuller, J., 1991; Empirical Dynamic-to-Static Core Calibrations for Offshore Clastic Engineering, British Petroleum (BP) Drilling Tech Group Technical Note.
    
    (2) dyn2sta_yme_plumb_bradford: Plumb and Bradford Sandstone Standard: Citation: Bradford, I. D. R., Fuller, J., Thompson, P. J., and Walsgrove, T. R., 1998; Benefits of Evaluating Mechanical Rock Properties for Drilling and Completion Operations, Proceedings of the SPE European Petroleum Conference, SPE-50893-MS (Everest Field Core Calibration Suite).
    
    (3) dyn2sta_yme_najibi_suite_one: Najibi Carbonates Suite 1 (Asmari and Sarvak Formations): Citation: Najibi, A., Ghafoori, M., Lashkaripur, G. R., and Asef, M. R., 2015; Empirical Relations Between Static and Dynamic Elastic Moduli of Carbonate Rocks, Journal of African Earth Sciences, Vol. 110, pp. 198–207.
    
    (4) dyn2sta_yme_sharifi: Sharifi Iranian Carbonates (Sarvak and Ilam Formations): Citation: Sharifi, J., Mirzakhanian, M., Mondol, N. H., and Saberi, M. R., 2017; Dynamic to Static Elastic Moduli Transformation of Carbonate Reservoir Core Plug Blocks from Southwest Iran, Journal of Geophysics and Engineering, Vol. 14, No. 6, pp. 1450–1462.
    
    (5) velocity2sta_yme_daoud: Daoud Direct Velocity Modulus (Kurdistan Region): Citation: Daoud, H. S., Rashid, K. A., and Alshkane, Y. M., 2017; Static and Dynamic Young's Modulus Correlation for Carbonate Rocks from Kurdistan Region, Iraq, Journal of Engineering and Applied Sciences, Vol. 12, pp. 5600–5607.
    
    (6) velocity2sta_yme_najibi: Najibi Direct Velocity Modulus (Carbonate Core Baselines): Citation: Najibi, A. et al., 2015; Velocity based transformation bypasses dynamic moduli parsing.
    
    (7) dyn2sta_yme_ohen: Ohen Shale Inversion: Citation: Ohen, H. A., 2003; Dynamic-to-Static Moduli Transformations Mapped for Shaly Unconsolidated Caprock Integrities, Corporate Technical Monograph Series, H3 Tech Exploration Services / Chevron Nigeria JIP.
    
    (8) dyn2sta_yme_vanheerden: Vanheerden Low Stress Constraint: Citation: Vanheerden, W. L., 1987; General Empirical Conversions Evaluated Under Low Confinement Baselines, Journal of the South African Institute of Mining and Metallurgy, Vol. 87, No. 10, pp. 325-331.
    
    (9.1) dyn2sta_yme_lacy_sandstone: Lacy General Empirical Conversions: Citation: Lacy, L. L., 1997; Dynamic and Static Modulus Relationships for Sands and Shales*, SPE Drilling & Completion, Vol. 12, No. 2, pp. 116–122, SPE-29426-PA.

    (9.2) dyn2sta_yme_lacy_shale: Lacy General Empirical Conversions: Citation: Lacy, L. L., 1997; Dynamic and Static Modulus Relationships for Sands and Shales*, SPE Drilling & Completion, Vol. 12, No. 2, pp. 116–122, SPE-29426-PA.

    (9.3) dyn2sta_yme_lacy_generic: Lacy General Empirical Conversions: Citation: Lacy, L. L., 1997; Dynamic and Static Modulus Relationships for Sands and Shales*, SPE Drilling & Completion, Vol. 12, No. 2, pp. 116–122, SPE-29426-PA.
    """

    @staticmethod
    def dyn2sta_yme_fuller(yme_dyn: float, unit: str) -> float:
        """(1) Convert dynamic to static Young's modulus using Fuller Sandstone/Shale Model.

        Equation type: Power law (y = a*x**b) 
                x = Dynamic Young's modulus in GPa
                multiplier = 0.0320
                exponent = 1.632

        Applicable for: Sandstone and Shale formations.

        Reference: Fuller, J., 1991. Empirical Dynamic-to-Static Core Calibrations for 
                   Offshore Clastic Engineering. British Petroleum (BP) Drilling Tech 
                   Group Technical Note.

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_fuller (float): Computed static Young's modulus magnitude 
                from the Fuller correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 0.0320
        exponent = 1.632
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_fuller_gpa = multiplier * yme_dyn_gpa**exponent
        yme_sta_fuller = UnitConverter.convert_pressure(yme_sta_fuller_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_fuller)

    @staticmethod
    def dyn2sta_yme_plumb_bradford(yme_dyn: float, unit: str) -> float:
        """(2) Convert dynamic to static Young's modulus using Plumb and Bradford Sandstone Standard.

        Equation type: Power law (y = a*x**b)
                x = Dynamic Young's modulus in GPa
                multiplier = 0.0018
                exponent = 2.7

        Applicable for: Sandstone formations (Everest Field Core Calibration Suite).

        Reference: Bradford, I. D. R., Fuller, J., Thompson, P. J., and Walsgrove, T. R., 
                   1998. Benefits of Evaluating Mechanical Rock Properties for Drilling 
                   and Completion Operations. Proceedings of the SPE European Petroleum 
                   Conference, SPE-50893-MS.

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_plumb_bradford (float): Computed static Young's modulus magnitude 
                from the Plumb and Bradford correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 0.0018
        exponent = 2.7
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_plumb_bradford_gpa = multiplier * yme_dyn_gpa**exponent
        yme_sta_plumb_bradford = UnitConverter.convert_pressure(yme_sta_plumb_bradford_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_plumb_bradford)

    @staticmethod
    def dyn2sta_yme_najibi_suite_one(yme_dyn: float, unit: str) -> float:
        """(3) Convert dynamic to static Young's modulus using Najibi Carbonates Suite 1.

        Equation type: Power law (y = a*x**b)
                x = Dynamic Young's modulus in GPa
                multiplier = 0.014
                exponent = 1.96

        Applicable for: Carbonate rocks (Asmari and Sarvak Formations).

        Reference: Najibi, A., Ghafoori, M., Lashkaripur, G. R., and Asef, M. R., 2015. 
                   Empirical Relations Between Static and Dynamic Elastic Moduli of 
                   Carbonate Rocks. Journal of African Earth Sciences, Vol. 110, 
                   pp. 198–207.

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_najibi_suite_one (float): Computed static Young's modulus magnitude 
                from the Najibi Suite 1 correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 0.014
        exponent = 1.96
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_najibi_suite_one_gpa = multiplier * yme_dyn_gpa**exponent
        yme_sta_najibi_suite_one = UnitConverter.convert_pressure(yme_sta_najibi_suite_one_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_najibi_suite_one)

    @staticmethod
    def dyn2sta_yme_sharifi(yme_dyn: float, unit: str) -> float:
        """(4) Convert dynamic to static Young's modulus using Sharifi Iranian Carbonates correlation.

        Equation type: Power law (y = a*x**b)
                x = Dynamic Young's modulus in GPa
                multiplier = 0.1098
                exponent = 1.2514

        Applicable for: Carbonate reservoirs (Sarvak and Ilam Formations).

        Reference: Sharifi, J., Mirzakhanian, M., Mondol, N. H., and Saberi, M. R., 
                   2017. Dynamic to Static Elastic Moduli Transformation of Carbonate 
                   Reservoir Core Plug Blocks from Southwest Iran. Journal of Geophysics 
                   and Engineering, Vol. 14, No. 6, pp. 1450–1462.

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_sharifi (float): Computed static Young's modulus magnitude 
                from the Sharifi correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 0.1098
        exponent = 1.2514
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_sharifi_gpa = multiplier * yme_dyn_gpa**exponent
        yme_sta_sharifi = UnitConverter.convert_pressure(yme_sta_sharifi_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_sharifi)

    @staticmethod
    def velocity2sta_yme_daoud(vp: float, unit: str) -> float:
        """(5) Convert compressional wave velocity to static Young's modulus using Daoud correlation.

        Equation type: Power law (y = a*(x/1000)**b)
                x = Compressional wave velocity in m/s
                multiplier = 0.295
                exponent = 3.4522

        Applicable for: Carbonate rocks from the Kurdistan Region.

        Reference: Daoud, H. S., Rashid, K. A., and Alshkane, Y. M., 2017. Static and 
                   Dynamic Young's Modulus Correlation for Carbonate Rocks from Kurdistan 
                   Region, Iraq. Journal of Engineering and Applied Sciences, Vol. 12, 
                   pp. 5600–5607.

        Args:
            vp (float): Compressional wave velocity. Unit: m/s
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_daoud (float): Computed static Young's modulus magnitude 
                from the Daoud correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 0.295
        exponent = 3.4522
        
        yme_sta_daoud_gpa = multiplier * (vp / 1000.0)**exponent
        yme_sta_daoud = UnitConverter.convert_pressure(yme_sta_daoud_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_daoud)

    @staticmethod
    def velocity2sta_yme_najibi(vp: float, unit: str) -> float:
        """(6) Convert compressional wave velocity to static Young's modulus using Najibi correlation.

        Equation type: Power law (y = a*(x/1000)**b)
                x = Compressional wave velocity in m/s
                multiplier = 0.169
                exponent = 3.324

        Applicable for: Carbonate core baselines.

        Reference: Najibi, A. et al., 2015. Velocity based transformation bypasses dynamic 
                   moduli parsing. Journal of Petroleum Science and Engineering.

        Args:
            vp (float): Compressional wave velocity. Unit: m/s
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_najibi (float): Computed static Young's modulus magnitude 
                from the Najibi correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 0.169
        exponent = 3.324
        
        yme_sta_najibi_gpa = multiplier * (vp / 1000.0)**exponent
        yme_sta_najibi = UnitConverter.convert_pressure(yme_sta_najibi_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_najibi)

        @staticmethod

    @staticmethod
    def dyn2sta_yme_ohen(yme_dyn: float, unit: str) -> float:
        """(7) Convert dynamic to static Young's modulus using Ohen Shale Inversion.

        Equation type: Power law (y = a*x**b) 
                x = Dynamic Young's modulus in GPa
                multiplier = 0.0158
                exponent = 2.74

        Applicable for: Shaly unconsolidated caprock integrities.

        Reference: Ohen, H. A., 2003. Dynamic-to-Static Moduli Transformations Mapped 
                   for Shaly Unconsolidated Caprock Integrities. Corporate Technical 
                   Monograph Series, H3 Tech Exploration Services / Chevron Nigeria JIP.

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_ohen (float): Computed static Young's modulus magnitude 
                from the Ohen correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 0.0158
        exponent = 2.74
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_ohen_gpa = multiplier * yme_dyn_gpa**exponent
        yme_sta_ohen = UnitConverter.convert_pressure(yme_sta_ohen_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_ohen)

    @staticmethod
    def dyn2sta_yme_vanheerden(yme_dyn: float, unit: str) -> float:
        """(8) Convert dynamic to static Young's modulus using Vanheerden Low Stress correlation.

        Equation type: Power law (y = a*x**b) 
                x = Dynamic Young's modulus in GPa
                multiplier = 0.097
                exponent = 1.485

        Applicable for: Low confinement baseline settings.

        Reference: Vanheerden, W. L., 1987. General Empirical Conversions Evaluated 
                   Under Low Confinement Baselines. Journal of the South African Institute 
                   of Mining and Metallurgy, Vol. 87, No. 10, pp. 325-331.

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_vanheerden (float): Computed static Young's modulus magnitude 
                from the Vanheerden correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 0.097
        exponent = 1.485
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_vanheerden_gpa = multiplier * yme_dyn_gpa**exponent
        yme_sta_vanheerden = UnitConverter.convert_pressure(yme_sta_vanheerden_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_vanheerden)

    @staticmethod
    def dyn2sta_yme_lacy_sandstone(yme_dyn: float, unit: str) -> float:
        """(9.1) Convert dynamic to static Young's modulus using Lacy Sandstone correlation.

        Equation type: Polynomial law (y = a*x**2 + b*x)
                x = Dynamic Young's modulus in GPa
                second_order_multiplier = 0.0026107
                first_order_multiplier = 0.4220

        Applicable for: Sandstone formations.

        Reference: Lacy, L. L., 1997. Dynamic and Static Modulus Relationships for 
                   Sands and Shales. SPE Drilling & Completion, Vol. 12, No. 2, 
                   pp. 116–122, SPE-29426-PA.

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_lacy_sandstone (float): Computed static Young's modulus magnitude 
                for sandstones. Unit: Pressure unit as defined in unit(str)
        """
        second_order_multiplier = 0.0026107
        first_order_multiplier = 0.4220
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_lacy_sandstone_gpa = (second_order_multiplier * yme_dyn_gpa**2.0) + (first_order_multiplier * yme_dyn_gpa)
        yme_sta_lacy_sandstone = UnitConverter.convert_pressure(yme_sta_lacy_sandstone_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_lacy_sandstone)

    @staticmethod
    def dyn2sta_yme_lacy_shale(yme_dyn: float, unit: str) -> float:
        """(9.2) Convert dynamic to static Young's modulus using Lacy Shale correlation.

        Equation type: Polynomial law (y = a*x**2 + b*x)
                x = Dynamic Young's modulus in GPa
                second_order_multiplier = 0.0042496
                first_order_multiplier = 0.4533

        Applicable for: Shale formations.

        Reference: Lacy, L. L., 1997. Dynamic and Static Modulus Relationships for 
                   Sands and Shales. SPE Drilling & Completion, Vol. 12, No. 2, 
                   pp. 116–122, SPE-29426-PA.

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_lacy_shale (float): Computed static Young's modulus magnitude 
                for shales. Unit: Pressure unit as defined in unit(str)
        """
        second_order_multiplier = 0.0042496
        first_order_multiplier = 0.4533
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_lacy_shale_gpa = (second_order_multiplier * yme_dyn_gpa**2.0) + (first_order_multiplier * yme_dyn_gpa)
        yme_sta_lacy_shale = UnitConverter.convert_pressure(yme_sta_lacy_shale_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_lacy_shale)

    @staticmethod
    def dyn2sta_yme_lacy_generic(yme_dyn: float, unit: str) -> float:
        """(9.3) Convert dynamic to static Young's modulus using Lacy Generic correlation.

        Equation type: Polynomial law (y = a*x**2 + b*x)
                x = Dynamic Young's modulus in GPa
                second_order_multiplier = 0.0061635
                first_order_multiplier = 0.2334

        Applicable for: Generic clastic rock classifications.

        Reference: Lacy, L. L., 1997. Dynamic and Static Modulus Relationships for 
                   Sands and Shales. SPE Drilling & Completion, Vol. 12, No. 2, 
                   pp. 116–122, SPE-29426-PA.

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_lacy_generic (float): Computed static Young's modulus magnitude 
                for generic lithologies. Unit: Pressure unit as defined in unit(str)
        """
        second_order_multiplier = 0.0061635
        first_order_multiplier = 0.2334
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_lacy_generic_gpa = (second_order_multiplier * yme_dyn_gpa**2.0) + (first_order_multiplier * yme_dyn_gpa)
        yme_sta_lacy_generic = UnitConverter.convert_pressure(yme_sta_lacy_generic_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_lacy_generic)

class StaticElasticPropertiesLinearLawConverter:
    """Convert dynamic elastic properties to static elastic properties using Linear Law correlations.

    ## Content:
    ### Linear Law Models (Es = a*X + b)
        (1) Ameen Arabian-Basin Matrix Model: Citation: Ameen, M. S., Smart, B. G. D., Somerville, J. M., Hammilton, S., and Naji, N. A., 2009; Static-Dynamic Relations for Stiffness and Strength Matrix Components in Carbonate Reservoirs (Saudi Arabia), Journal of Petroleum Science and Engineering, Vol. 68, No. 1, pp. 39–50.
        
        (2) Wang Empirical Shift: Citation: Wang, Z. N., 2000; Dynamic-to-Static Conversion Shift Matrix across Multiple Overburden Geographies, Core Laboratories Rock Mechanics Research Division Internal whitepaper.
        
        (3) Yasar and Erdogan Linear Velocity Conversion: Citation: Yasar, E., and Erdogan, Y., 2004; Correlations between Sound Velocity and Static Elastic Moduli of Carbonate Rocks, International Journal of Rock Mechanics and Mining Sciences, Vol. 41, No. 5, pp. 871-875.
        
        (4) Wang and Nur Mechanical Threshold Splits: Citation: Wang, Z., and Nur, A. A., 2000; Static and Dynamic Moduli of Sedimentary Rocks: General Core Cross-Plots, in Seismic and Acoustic Velocities in Reservoir Rocks, Volume 3, Society of Exploration Geophysicists (SEG).
        
        (5) Eissa and Kazi Global Sedimentary Line: Citation: Eissa, E. A., and Kazi, A., 1988; Relation between Static and Dynamic Moduli of Rocks, International Journal of Rock Mechanics and Mining Sciences & Geomechanics Abstracts, Vol. 25, No. 6, pp. 479 482.
        
        (6) Brotons Calcarenite Framework: Citation: Brotons, V., Tomás, R., Ivorra, S., and Grediaga, A., 2014; Empirical Correlation between Static and Dynamic Moduli of Porous Calcarenites, Rock Mechanics and Rock Engineering, Vol. 47, No. 6, pp. 2237 2247.
        
        (7) Sayers Schlumberger Private Stonefish Line: Citation: Sayers, C. M., 2005; Schlumberger Carbonate Geomechanics User Defined Function Reference Standard (Private Module Calibration).
    """

    @staticmethod
    def dyn2sta_yme_ameen(yme_dyn: float, unit: str) -> float:
        """Convert dynamic to static Young's modulus using Ameen Arabian-Basin Matrix Model.

        Equation type: Linear law (y = a*x + b)
                x = Dynamic Young's modulus in GPa
                multiplier = 0.541
                intercept = 12.852

        Applicable for: Carbonate reservoirs (Saudi Arabia).

        Reference: Ameen, M. S., Smart, B. G. D., Somerville, J. M., Hammilton, S., 
                   and Naji, N. A., 2009. Static-Dynamic Relations for Stiffness and 
                   Strength Matrix Components in Carbonate Reservoirs (Saudi Arabia). 
                   Journal of Petroleum Science and Engineering, Vol. 68, No. 1, pp. 39–50.

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_ameen (float): Computed static Young's modulus magnitude 
                from the Ameen correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 0.541
        intercept = 12.852
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_ameen_gpa = multiplier * yme_dyn_gpa + intercept
        yme_sta_ameen = UnitConverter.convert_pressure(yme_sta_ameen_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_ameen)

    @staticmethod
    def dyn2sta_yme_wang_shift(yme_dyn: float, unit: str) -> float:
        """Convert dynamic to static Young's modulus using Wang Empirical Shift.

        Equation type: Linear law (y = x - b)
                x = Dynamic Young's modulus in GPa
                intercept = 15.2

        Applicable for: Multiple overburden geographies.

        Reference: Wang, Z. N., 2000. Dynamic-to-Static Conversion Shift Matrix across 
                   Multiple Overburden Geographies. Core Laboratories Rock Mechanics 
                   Research Division Internal whitepaper.

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_wang_shift (float): Computed static Young's modulus magnitude 
                from the Wang shift correlation. Unit: Pressure unit as defined in unit(str)
        """
        intercept = 15.2
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_wang_shift_gpa = yme_dyn_gpa - intercept
        yme_sta_wang_shift = UnitConverter.convert_pressure(yme_sta_wang_shift_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_wang_shift)

    @staticmethod
    def velocity2sta_yme_yasar_erdogan(vp: float, unit: str) -> float:
        """Convert compressional wave velocity to static Young's modulus using Yasar and Erdogan correlation.

        Equation type: Linear law (y = a*(x/1000) + b)
                x = Compressional wave velocity in m/s
                multiplier = 10.67
                intercept = -18.71

        Applicable for: Carbonate rocks.

        Reference: Yasar, E., and Erdogan, Y., 2004. Correlations between Sound Velocity 
                   and Static Elastic Moduli of Carbonate Rocks. International Journal 
                   of Rock Mechanics and Mining Sciences, Vol. 41, No. 5, pp. 871–875.

        Args:
            vp (float): Compressional wave velocity. Unit: m/s
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_yasar_erdogan (float): Computed static Young's modulus magnitude 
                from the Yasar and Erdogan correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 10.67
        intercept = -18.71
        
        yme_sta_yasar_erdogan_gpa = multiplier * (vp / 1000.0) + intercept
        yme_sta_yasar_erdogan = UnitConverter.convert_pressure(yme_sta_yasar_erdogan_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_yasar_erdogan)

    @staticmethod
    def dyn2sta_yme_wang_nur_soft(yme_dyn: float, unit: str) -> float:
        """Convert dynamic to static Young's modulus using Wang and Nur Soft Sedimentary constraint.

        Equation type: Linear law (y = a*x + b)
                x = Dynamic Young's modulus in GPa
                multiplier = 0.41
                intercept = -1.06

        Applicable for: Soft Sedimentary Formations (Where Static Modulus < 15 GPa).

        Reference: Wang, Z., and Nur, A. A., 2000. Static and Dynamic Moduli of 
                   Sedimentary Rocks: General Core Cross-Plots. In Seismic and Acoustic 
                   Velocities in Reservoir Rocks, Volume 3, Society of Exploration Geophysicists (SEG).

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_wang_nur_soft (float): Computed static Young's modulus magnitude 
                for soft formations. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 0.41
        intercept = -1.06
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_wang_nur_soft_gpa = multiplier * yme_dyn_gpa + intercept
        yme_sta_wang_nur_soft = UnitConverter.convert_pressure(yme_sta_wang_nur_soft_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_wang_nur_soft)

    @staticmethod
    def dyn2sta_yme_wang_nur_hard(yme_dyn: float, unit: str) -> float:
        """Convert dynamic to static Young's modulus using Wang and Nur Hard Sedimentary constraint.

        Equation type: Linear law (y = a*x + b)
                x = Dynamic Young's modulus in GPa
                multiplier = 1.153
                intercept = -15.22

        Applicable for: Hard Sedimentary Formations (Where Static Modulus > 15 GPa).

        Reference: Wang, Z., and Nur, A. A., 2000. Static and Dynamic Moduli of 
                   Sedimentary Rocks: General Core Cross-Plots. In Seismic and Acoustic 
                   Velocities in Reservoir Rocks, Volume 3, Society of Exploration Geophysicists (SEG).

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_wang_nur_hard (float): Computed static Young's modulus magnitude 
                for hard formations. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 1.153
        intercept = -15.22
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_wang_nur_hard_gpa = multiplier * yme_dyn_gpa + intercept
        yme_sta_wang_nur_hard = UnitConverter.convert_pressure(yme_sta_wang_nur_hard_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_wang_nur_hard)

    @staticmethod
    def dyn2sta_yme_eissa_kazi(yme_dyn: float, unit: str) -> float:
        """Convert dynamic to static Young's modulus using Eissa and Kazi Global Sedimentary Line.

        Equation type: Linear law (y = a*x + b)
                x = Dynamic Young's modulus in GPa
                multiplier = 0.74
                intercept = -0.82

        Applicable for: Global sedimentary rock networks.

        Reference: Eissa, E. A., and Kazi, A., 1988. Relation between Static and Dynamic 
                   Moduli of Rocks. International Journal of Rock Mechanics and Mining 
                   Sciences & Geomechanics Abstracts, Vol. 25, No. 6, pp. 479–482.

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_eissa_kazi (float): Computed static Young's modulus magnitude 
                from the Eissa and Kazi correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 0.74
        intercept = -0.82
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_eissa_kazi_gpa = multiplier * yme_dyn_gpa + intercept
        yme_sta_eissa_kazi = UnitConverter.convert_pressure(yme_sta_eissa_kazi_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_eissa_kazi)

    @staticmethod
    def dyn2sta_yme_brotons(yme_dyn: float, unit: str) -> float:
        """Convert dynamic to static Young's modulus using Brotons Calcarenite Framework.

        Equation type: Linear law (y = a*x + b)
                x = Dynamic Young's modulus in GPa
                multiplier = 0.867
                intercept = -2.085

        Applicable for: Porous calcarenites structures.

        Reference: Brotons, V., Tomás, R., Ivorra, S., and Grediaga, A., 2014. Empirical 
                   Correlation between Static and Dynamic Moduli of Porous Calcarenites. 
                   Rock Mechanics and Rock Engineering, Vol. 47, No. 6, pp. 2237–2247.

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_brotons (float): Computed static Young's modulus magnitude 
                from the Brotons correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 0.867
        intercept = -2.085
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_brotons_gpa = multiplier * yme_dyn_gpa + intercept
        yme_sta_brotons = UnitConverter.convert_pressure(yme_sta_brotons_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_brotons)

    @staticmethod
    def dyn2sta_yme_sayers(yme_dyn: float, unit: str) -> float:
        """Convert dynamic to static Young's modulus using Sayers Schlumberger Private Stonefish Line.

        Equation type: Linear law (y = a*x)
                x = Dynamic Young's modulus in GPa
                multiplier = 0.7792

        Applicable for: Carbonate Geomechanics modules.

        Reference: Sayers, C. M., 2005. Schlumberger Carbonate Geomechanics User Defined 
                   Function Reference Reference Standard (Private Module Calibration).

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_sayers (float): Computed static Young's modulus magnitude 
                from the Sayers private module calibration. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 0.7792
        
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        yme_sta_sayers_gpa = multiplier * yme_dyn_gpa
        yme_sta_sayers = UnitConverter.convert_pressure(yme_sta_sayers_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_sayers)
   
class StaticElasticPropertiesExponentialLawConverter:
    """Convert dynamic elastic properties to static elastic properties.

    Reference:
        ### Exponential Law Porosity Models (Es = a * exp(b * por))
        (1) Ameen Ghawar Arab-D Dolomite: Citation: Ameen, M. S. et al., 2009; Ghawar Field Reservoirs - Mechanical Characterization and Dynamic to Static Stiffness Transforms, Journal of Petroleum Science and Engineering, Vol. 68, pp. 39-50.
        
        (2) Ameen Ghawar Arab-D Limestone: Citation: Ameen, M. S. et al., 2009; Core-calibrated limestone suite matrix from the Ghawar Field, Saudi Arabia.
        
        (3) Farquhar UK North Sea Carbonates: Citation: Farquhar, R. A., Somerville, J. M., and Smart, B. G., 1994; Porosity as a Direct Controlling Parameter for Static Young's Modulus in Central North Sea Chalk Core Plugs, Proceedings of the SPE European Petroleum Conference, SPE-28847-MS.
        
        (4) Faÿ-Gomord High-Porosity Chalk Matrix: Citation: Faÿ-Gomord, O., Descamps, F., Tshibangu, J. P., Vandycke, S., and Swennen, R., 2016; Stiffness, Strength, and Porosity Relationships in Fractured Chalk Reservoirs of Western Europe, Journal of Structural Geology, Vol. 84, pp. 122–136.
        
        (5) Hadi and Nygaard Basinal Sweeps: Citation: Hadi, F., and Nygaard, R., 2011; Empirical Modeling of Carbonate and Chalk Stiffness Inversion Functions for Drilling Optimization, Journal of Petroleum Exploration and Production Technology, Vol. 1, pp. 29–38.
        """

    @staticmethod
    def por2sta_yme_ameen_dolomite(por: float, unit: str) -> float:
        """Convert porosity directly to static Young's modulus using Ameen Ghawar Arab-D Dolomite Model.

        Equation type: Exponential law (y = a*exp(b*x)) 
                x = Total fractional porosity (por)
                multiplier = 92.612
                exponent_factor = -4.75

        Applicable for: Dolomite reservoir formations (Ghawar Field Arab-D).

        Reference: Ameen, M. S. et al., 2009. Ghawar Field Reservoirs - Mechanical 
                   Characterization and Dynamic to Static Stiffness Transforms. Journal 
                   of Petroleum Science and Engineering, Vol. 68, pp. 39-50.

        Args:
            por (float): Total fractional porosity magnitude. Unit: unitless
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_ameen_dolomite (float): Computed static Young's modulus magnitude 
                from the Ameen dolomite correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 92.612
        exponent_factor = -4.75
        
        yme_sta_gpa = multiplier * math.exp(exponent_factor * por)
        yme_sta_ameen_dolomite = UnitConverter.convert_pressure(yme_sta_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_ameen_dolomite)

    @staticmethod
    def por2sta_yme_ameen_limestone(por: float, unit: str) -> float:
        """Convert porosity directly to static Young's modulus using Ameen Ghawar Arab-D Limestone Model.

        Equation type: Exponential law (y = a*exp(b*x)) 
                x = Total fractional porosity (por)
                multiplier = 86.094
                exponent_factor = -5.34

        Applicable for: Limestone reservoir formations (Ghawar Field Arab-D).

        Reference: Ameen, M. S. et al., 2009. Core-calibrated limestone suite matrix 
                   from the Ghawar Field, Saudi Arabia. Journal of Petroleum Science 
                   and Engineering, Vol. 68, pp. 39-50.

        Args:
            por (float): Total fractional porosity magnitude. Unit: unitless
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_ameen_limestone (float): Computed static Young's modulus magnitude 
                from the Ameen limestone correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 86.094
        exponent_factor = -5.34
        
        yme_sta_gpa = multiplier * math.exp(exponent_factor * por)
        yme_sta_ameen_limestone = UnitConverter.convert_pressure(yme_sta_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_ameen_limestone)

    @staticmethod
    def por2sta_yme_farquhar(por: float, unit: str) -> float:
        """Convert porosity directly to static Young's modulus using Farquhar UK North Sea Carbonates Model.

        Equation type: Exponential law (y = a*exp(b*x)) 
                x = Total fractional porosity (por)
                multiplier = 69.05
                exponent_factor = -6.0

        Applicable for: Central North Sea Chalk Core Plugs.

        Reference: Farquhar, R. A., Somerville, J. M., and Smart, B. G., 1994. Porosity 
                   as a Direct Controlling Parameter for Static Young's Modulus in Central 
                   North Sea Chalk Core Plugs. Proceedings of the SPE European Petroleum 
                   Conference, SPE-28847-MS.

        Args:
            por (float): Total fractional porosity magnitude. Unit: unitless
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_farquhar (float): Computed static Young's modulus magnitude 
                from the Farquhar correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 69.05
        exponent_factor = -6.0
        
        yme_sta_gpa = multiplier * math.exp(exponent_factor * por)
        yme_sta_farquhar = UnitConverter.convert_pressure(yme_sta_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_farquhar)

    @staticmethod
    def por2sta_yme_fay_gomord(por: float, unit: str) -> float:
        """Convert porosity directly to static Young's modulus using Faÿ-Gomord High-Porosity Chalk Model.

        Equation type: Exponential law (y = a*exp(b*x)) 
                x = Total fractional porosity (por)
                multiplier = 249.968
                exponent_factor = -10.8

        Applicable for: Fractured chalk reservoirs of Western Europe.

        Reference: Faÿ-Gomord, O., Descamps, F., Tshibangu, J. P., Vandycke, S., and 
                   Swennen, R., 2016. Stiffness, Strength, and Porosity Relationships 
                   in Fractured Chalk Reservoirs of Western Europe. Journal of 
                   Structural Geology, Vol. 84, pp. 122–136.

        Args:
            por (float): Total fractional porosity magnitude. Unit: unitless
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_fay_gomord (float): Computed static Young's modulus magnitude 
                from the Faÿ-Gomord correlation. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 249.968
        exponent_factor = -10.8
        
        yme_sta_gpa = multiplier * math.exp(exponent_factor * por)
        yme_sta_fay_gomord = UnitConverter.convert_pressure(yme_sta_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_fay_gomord)

    @staticmethod
    def por2sta_yme_hadi_chalk(por: float, unit: str) -> float:
        """Convert porosity directly to static Young's modulus using Hadi and Nygaard Chalk Inversion.

        Equation type: Exponential law (y = a*exp(b*x)) 
                x = Total fractional porosity (por)
                multiplier = 52.6
                exponent_factor = -5.62

        Applicable for: Chalk reservoir structures.

        Reference: Hadi, F., and Nygaard, R., 2011. Empirical Modeling of Carbonate 
                   and Chalk Stiffness Inversion Functions for Drilling Optimization. 
                   Journal of Petroleum Exploration and Production Technology, Vol. 1, pp. 29–38.

        Args:
            por (float): Total fractional porosity magnitude. Unit: unitless
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_hadi_chalk (float): Computed static Young's modulus magnitude 
                for chalk formations. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 52.6
        exponent_factor = -5.62
        
        yme_sta_gpa = multiplier * math.exp(exponent_factor * por)
        yme_sta_hadi_chalk = UnitConverter.convert_pressure(yme_sta_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_hadi_chalk)

    @staticmethod
    def por2sta_yme_hadi_dolomite(por: float, unit: str) -> float:
        """Convert porosity directly to static Young's modulus using Hadi and Nygaard Dolomite Inversion.

        Equation type: Exponential law (y = a*exp(b*x)) 
                x = Total fractional porosity (por)
                multiplier = 95.0
                exponent_factor = -8.06

        Applicable for: Dolomite reservoir structures.

        Reference: Hadi, F., and Nygaard, R., 2011. Empirical Modeling of Carbonate 
                   and Chalk Stiffness Inversion Functions for Drilling Optimization. 
                   Journal of Petroleum Exploration and Production Technology, Vol. 1, pp. 29–38.

        Args:
            por (float): Total fractional porosity magnitude. Unit: unitless
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_hadi_dolomite (float): Computed static Young's modulus magnitude 
                for dolomite formations. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 95.0
        exponent_factor = -8.06
        
        yme_sta_gpa = multiplier * math.exp(exponent_factor * por)
        yme_sta_hadi_dolomite = UnitConverter.convert_pressure(yme_sta_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_hadi_dolomite)

    @staticmethod
    def por2sta_yme_hadi_limestone(por: float, unit: str) -> float:
        """Convert porosity directly to static Young's modulus using Hadi and Nygaard Limestone Inversion.

        Equation type: Exponential law (y = a*exp(b*x)) 
                x = Total fractional porosity (por)
                multiplier = 56.9
                exponent_factor = -4.46

        Applicable for: Limestone reservoir structures.

        Reference: Hadi, F., and Nygaard, R., 2011. Empirical Modeling of Carbonate 
                   and Chalk Stiffness Inversion Functions for Drilling Optimization. 
                   Journal of Petroleum Exploration and Production Technology, Vol. 1, pp. 29–38.

        Args:
            por (float): Total fractional porosity magnitude. Unit: unitless
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_hadi_limestone (float): Computed static Young's modulus magnitude 
                for limestone formations. Unit: Pressure unit as defined in unit(str)
        """
        multiplier = 56.9
        exponent_factor = -4.46
        
        yme_sta_gpa = multiplier * math.exp(exponent_factor * por)
        yme_sta_hadi_limestone = UnitConverter.convert_pressure(yme_sta_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_hadi_limestone)

class StaticElasticPropertiesOtherConverter:
    """Convert dynamic elastic properties to static elastic properties.

    ## Content:
    ### Logarithmic & Segmented Models (Morales Frameworks)
    (1) Morales Grain Support Modification: Citation: Morales, R. H., and Marcinew, R. P., 1993; Fracturing of High-Permeability Formations: Mechanical Properties Correlations, Proceedings of the SPE SPE/IADC Drilling Conference, SPE-25396-MS.
    
    (2) Morales Porosity Segmented Sandstones: Citation: Morales, R. H., and Marcinew, R. P., 1993; Raw equations mapped to psi log-space bounds.
    ### Generic Equations
    (3) Custom Power Law: Generic (y = a*x**b)
    
    (4) Custom Linear Law: Generic (y = a*x + b)
    ### Static Poisson's Ratio 
    (5) Linear Poisson Multiplier Shift: Citation: Technical Workflow Manual - Geomechanical Logging Transformations, Baker Hughes Inc., 2006."""



    @staticmethod
    def dyn2sta_yme_morales_grain(yme_dyn: float, por: float, unit: str) -> float:
        """Convert dynamic to static Young's modulus using Morales Grain Support Modification.

        Equation type: Segmented Law (y = (-2.21*phi_c + 0.965)*x)
                x = Dynamic Young's modulus in GPa
                phi_c = Critical porosity boundary condition (por capped at 0.35)
                por = Total fractional porosity

        Applicable for: High-permeability formations.

        Reference: Morales, R. H., and Marcinew, R. P., 1993. Fracturing of 
                   High-Permeability Formations: Mechanical Properties Correlations. 
                   Proceedings of the SPE SPE/IADC Drilling Conference, SPE-25396-MS.

        Args:
            yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Pressure unit as defined in unit(str)
            por (float): Total fractional porosity magnitude. Unit: unitless
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            yme_sta_morales_grain (float): Computed static Young's modulus magnitude 
                from the Morales correlation. Unit: Pressure unit as defined in unit(str)
        """
        yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn, unit_in=unit, unit_out='gpa')
        
        if por < 0.35:
            phi_c = por
        else:
            phi_c = 0.35
            
        yme_sta_morales_grain_gpa = (-2.21 * phi_c + 0.965) * yme_dyn_gpa
        yme_sta_morales_grain = UnitConverter.convert_pressure(yme_sta_morales_grain_gpa, unit_in='gpa', unit_out=unit)

        return float(yme_sta_morales_grain)

    @staticmethod
    def dyn2sta_yme_morales(yme_dyn: float, porosity: float, exclude_low_por: bool = True, unit: str) -> float:
        """Convert dynamic to static Young's modulus using Morales correlation.

        Equation type: Power law (y = a*x**b)
        Applicable for: Sandstones:
                          Consolidated (10%-15%)
                          Moderately Consolidated (15%-25%)
                          Weakly consolidated (>25%)

        Reference: R.H. Morales and R. P. Marcinew, Fracturing of High-Permeability Formations: Mechanical Properties Correlations, SPE 26561, October 1993.

        Args:
           yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Mpsi
           porosity (float): Formation porosity as a fraction (e.g. 0.20 for 20%). Unit: unitless
           exclude_low_por (bool): If True, returns -9999 for porosity below 0.10. Defaults to True.

        Returns:
           float: Static Young's modulus magnitude (Morales). Unit: Mpsi"""
        if porosity < 0.15:
            multiplier = 2.137
            exponent = 0.6612
        elif porosity < 0.25:
            multiplier = 1.829
            exponent = 0.6920
        elif porosity > 0.25:
            multiplier = -0.4575
            exponent = 0.9404 

        if porosity < 0.10 and exclude_low_por is True:
            yme_sta_morales = -9999
        else:
            yme_dyn_gpa = UnitConverter.convert_pressure(yme_dyn,unit_in = unit,unit_out='mpsi')
            yme_sta_morales_mpsi = multiplier * yme_dyn_gpa**exponent
            yme_sta_morales = UnitConverter.convert_pressure(yme_sta_morales_mpsi,unit_in='mpsi',unit_out=unit)
        return yme_sta_morales

    @staticmethod
    def convert_dyn2sta_yme_custom_power_law(yme_dyn: float, multiplier: float, exponent: float) -> float:
        """Convert dynamic to static Young's modulus using a custom power law.

        Equation type: Power law (y = a*x**b)
        Applicable for: Generic.

        Args:
           yme_dyn (float): Dynamic Young's modulus magnitude. Unit: Mpsi
           multiplier (float): Custom multiplier - from core data regression or alternative power law correlation Unit: unitless
           exponent (float): Custom exponent - from core data regression or alternative power law correlation  Unit: unitless

        Returns:
           yme_sta_morales: Static Young's modulus magnitude. Unit: Mpsi"""
        yme_sta_power_law = multiplier * yme_dyn**exponent

        return float(yme_sta_power_law)

    @staticmethod
    def dyn2sta_yme_custom_linear_law(yme_dyn: float, slope: float, intercept: float) -> float:
        """Convert dynamic to static Young's modulus using a custom linear law.

        Equation type: Linear law (y = a*x + b)
        Applicable for: Generic.

        Args:
           yme_dyn (float): Dynamic Young's modulus magnitude Unit: Mpsi
           slope (float):  custom slope - from core data regression or alternative linear law correlation Unit: unitless
           intercept (float):  custom intercept - from core data regression or alternative linear law correlation Unit: unitless

        Returns:
           float: Static Young's modulus magnitude. Unit: Mpsi"""
        yme_sta_linear_law = slope * yme_dyn + intercept

        return yme_sta_linear_law

    @staticmethod
    def dyn2sta_poissons_ratio(pr_dyn: float, multiplier: float) -> float:
        """Convert dynamic to static Poisson's ratio using a constant multiplier.

        Equation type: Linear scaling (y = a*b)
        Applicable for: Generic.

        Reference: -

        Args:
           pr_dyn (float): Dynamic Poisson's ratio. Unit: unitless
           multiplier (float): Constant multiplier. Unit: unitless. Default set to 1.

        Returns:
           float: Static Poisson's ratio. Unit: unitless"""
        pr_sta = pr_dyn * multiplier

        return pr_sta

class BiotCoefficientConverter:
    """Calculate the Poroelastic Parameter (Biot's Coefficient Alpha).

    Reference:
        Formulations capturing pore-fluid pressure volumetric coupling parameters for effective matrix stress calculations.
        (1) Klimentos Permeability Matrix Inversion: Citation: Klimentos, T., 2003; Direct Extraction of Biot's Coefficient utilizing Joint NMR/CMR Pore Geometry Identifiers, Transactions of the 44th SPWLA Annual Logging Symposium, Paper QQ.
        (2) Krief Porosity Volumetric Coupling: Citation: Krief, M., Garat, J., Stellingwerff, J., and Ventre, J., 1990; A Petrophysical Interpretation of Seismic Wave Velocities in Clay-Bearing, Shaly Sandstone Formations, The Log Analyst, Vol. 31, No. 5, pp. 355–369.
        (3) Mecpro Analytical Alpha: Citation: Mecpro Technical Reference Suite - Rock Intelligence Framework, CGG Veritas Geomechanic Modules, 2011.
        (4) Biot Coefficient constant law: Generic (y = c)"""
    
    @staticmethod
    def calculate_biot_coefficient_klimentos(permeability_md: float) -> float:
        """Calculate Biot's coefficient using Klimentos Permeability Matrix Inversion.

        Equation type: Logarithmic law (y = a * log(x) + b)
                x = Matrix permeability in mD (perm)
                multiplier = 0.045
                intercept = 0.7

        Applicable for: Complex pore geometries and shaly sandstones.

        Reference: Klimentos, T., 2003. Direct Extraction of Biot's Coefficient utilizing 
                   Joint NMR/CMR Pore Geometry Identifiers. Transactions of the 44th 
                   SPWLA Annual Logging Symposium, Paper QQ.

        Args:
            permeability_md (float): Matrix permeability magnitude. Unit: mD

        Returns:
            biot_coefficient_klimentos (float): Computed Biot's effective stress 
                coefficient. Unit: unitless
        """
        multiplier = 0.045
        intercept = 0.7

        biot_coefficient_klimentos = multiplier * math.log(permeability_md) + intercept

        return float(biot_coefficient_klimentos)

    @staticmethod
    def calculate_biot_coefficient_krief(por: float) -> float:
        """Calculate Biot's coefficient using Krief Porosity Volumetric Coupling.

        Equation type: Power law (y = 1 - (1 - x)**(3 / (1 - x)))
                x = Total fractional porosity (por)

        Applicable for: Clay-bearing, shaly sandstone formations.

        Reference: Krief, M., Garat, J., Stellingwerff, J., and Ventre, J., 1990. 
                   A Petrophysical Interpretation of Seismic Wave Velocities in 
                   Clay-Bearing, Shaly Sandstone Formations. The Log Analyst, 
                   Vol. 31, No. 5, pp. 355–369.

        Args:
            por (float): Total fractional porosity magnitude. Unit: unitless

        Returns:
            biot_coefficient_krief (float): Computed Biot's effective stress 
                coefficient. Unit: unitless
        """
        porosity_complement = 1.0 - por
        exponent = 3.0 / porosity_complement

        biot_coefficient_krief = 1.0 - (porosity_complement ** exponent)

        return float(biot_coefficient_krief)

    @staticmethod
    def calculate_biot_coefficient_mecpro(bmk_dyn: float, bmk_solid: float) -> float:
        """Calculate Biot's coefficient using Mecpro Analytical Alpha framework.

        Equation type: Rational compliance law (y = 1 - x_frame / x_solid)
                x_frame = Dynamic framework bulk modulus (bmk_dyn)
                x_solid = Grain matrix mineral bulk modulus (bmk_solid)

        Applicable for: Generic clean and shaly sandstone/carbonate reservoir blocks.

        Reference: Mecpro Technical Reference Suite - Rock Intelligence Framework, 
                   CGG Veritas Geomechanic Modules, 2011.

        Args:
            bmk_dyn (float): Dynamic bulk modulus of the dry rock frame. Unit: Pressure unit as defined in unit(str)
            bmk_solid (float): Bulk modulus of the clean solid grain matrix. Unit: Pressure unit as defined in unit(str)
            unit (str): Target physical unit for conversion boundary logic. Unit: unitless

        Returns:
            biot_coefficient_mecpro (float): Computed Biot's effective stress 
                coefficient. Unit: unitless
        """
        biot_coefficient_mecpro = 1.0 - (bulk_modulus_frame / bulk_modulus_solid)

        return float(biot_coefficient_mecpro)

    @staticmethod
    def biot_coefficient_constant_law(constant: float) -> float:
        """Assign a constant value for the Biot coefficient.

        Equation type: Constant law (y = c)
        Applicable for: Generic.

        Reference: Jaeger, John Conrad, Neville GW Cook, and Robert Zimmerman. Fundamentals of rock mechanics. John Wiley & Sons, 2009.

        Args:
           constant (float): Constant value for Biot's coefficient. Unit: unitless

        Returns:
           float: Biot coefficient. Unit: unitless"""
        biot_constant = constant

        return biot_constant


    
    # @staticmethod
    # def dyn2sta_yme_bradord_array(yme_dyn: list[float]) -> list[float]:
    #     """Convert an array of dynamic Young's modulus values to static using Bradford correlation.

    #     Args:
    #        yme_dyn (list[float]): Dynamic Young's modulus values. Unit: Mpsi

    #     Returns:
    #        list[float]: Static Young's modulus values (Bradford). Unit: Mpsi"""
    #     return [
    #         StaticElasticPropertiesConverter.dyn2sta_yme_bradord(yme_dyn=value)
    #         for value in yme_dyn
    #     ]

    # @staticmethod
    # def dyn2sta_yme_najib_array(yme_dyn: list[float]) -> list[float]:
    #     """Convert an array of dynamic Young's modulus values to static using Najib correlation.

    #     Args:
    #        yme_dyn (list[float]): Dynamic Young's modulus values. Unit: Mpsi

    #     Returns:
    #        list[float]: Static Young's modulus values (Najib). Unit: Mpsi"""
    #     return [
    #         StaticElasticPropertiesConverter.dyn2sta_yme_najib(yme_dyn=value)
    #         for value in yme_dyn
    #     ]

    # @staticmethod
    # def dyn2sta_yme_fuller_array(yme_dyn: list[float]) -> list[float]:
    #     """Convert an array of dynamic Young's modulus values to static using Fuller correlation.

    #     Args:
    #        yme_dyn (list[float]): Dynamic Young's modulus values. Unit: GPa

    #     Returns:
    #        list[float]: Static Young's modulus values (Fuller). Unit: GPa"""
    #     return [
    #         StaticElasticPropertiesConverter.dyn2sta_yme_fuller(yme_dyn=value)
    #         for value in yme_dyn
    #     ]

    # @staticmethod
    # def dyn2sta_yme_morales_array(yme_dyn: list[float], porosity: list[float], exclude_low_por: bool = True) -> list[float]:
    #     """Convert arrays of dynamic Young's modulus and porosity to static Young's modulus using Morales correlation.

    #     Args:
    #        yme_dyn (list[float]): Dynamic Young's modulus values. Unit: Mpsi
    #        porosity (list[float]): Formation porosity values as fractions. Unit: unitless
    #        exclude_low_por (bool): If True, returns -9999 for porosity below 0.10. Defaults to True.

    #     Returns:
    #        list[float]: Static Young's modulus values (Morales). Unit: Mpsi"""
    #     return [
    #         StaticElasticPropertiesConverter.dyn2sta_yme_morales(
    #             yme_dyn=yme_value,
    #             porosity=porosity_value,
    #             exclude_low_por=exclude_low_por,
    #         )
    #         for yme_value, porosity_value in zip(yme_dyn, porosity, strict=True)
    #     ]

    # @staticmethod
    # def convert_dyn2sta_yme_custom_power_law_array(yme_dyn: list[float], multiplier: float, exponent: float) -> list[float]:
    #     """Convert an array of dynamic Young's modulus values to static using a custom power law.

    #     Args:
    #        yme_dyn (list[float]): Dynamic Young's modulus values. Unit: Mpsi
    #        multiplier (float): Custom multiplier from regression. Unit: unitless
    #        exponent (float): Custom exponent from regression. Unit: unitless

    #     Returns:
    #        list[float]: Static Young's modulus values. Unit: Mpsi"""
    #     return [
    #         StaticElasticPropertiesConverter.convert_dyn2sta_yme_custom_power_law(
    #             yme_dyn=value,
    #             multiplier=multiplier,
    #             exponent=exponent,
    #         )
    #         for value in yme_dyn
    #     ]

    # @staticmethod
    # def dyn2sta_yme_custom_linear_law_array(yme_dyn: list[float], slope: float, intercept: float) -> list[float]:
    #     """Convert an array of dynamic Young's modulus values to static using a custom linear law.

    #     Args:
    #        yme_dyn (list[float]): Dynamic Young's modulus values. Unit: Mpsi
    #        slope (float): Custom slope from regression. Unit: unitless
    #        intercept (float): Custom intercept from regression. Unit: unitless

    #     Returns:
    #        list[float]: Static Young's modulus values. Unit: Mpsi"""
    #     return [
    #         StaticElasticPropertiesConverter.dyn2sta_yme_custom_linear_law(
    #             yme_dyn=value,
    #             slope=slope,
    #             intercept=intercept,
    #         )
    #         for value in yme_dyn
    #     ]

    # @staticmethod
    # def dyn2sta_poissons_ratio_array(pr_dyn: list[float], multiplier: float) -> list[float]:

    #     """Convert an array of dynamic Poisson's ratio values to static using a constant multiplier.

    #     Args:
    #        pr_dyn (list[float]): Dynamic Poisson's ratio values. Unit: unitless
    #        multiplier (float): Constant multiplier. Unit: unitless

    #     Returns:
    #        list[float]: Static Poisson's ratio values. Unit: unitless"""
    #     return [
    #         StaticElasticPropertiesConverter.dyn2sta_poissons_ratio(pr_dyn=value, multiplier=multiplier)
    #         for value in pr_dyn
    #     ]
