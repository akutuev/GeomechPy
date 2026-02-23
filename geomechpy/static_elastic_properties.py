import math
from dataclasses import dataclass


class StaticElasticPropertiesConverter:
    """Convert dynamic elastic properties to static elastic properties

    Reference:
       Zhang, Yuliang, et al. "Extracting static elastic moduli of rock through elastic wave velocities." Acta Geophysica 72.2 (2024): 915-931.

    """

    @staticmethod
    def convert_dyn2sta_yme_custom_power_law(measurement: float, multiplier: float, exponent: float) -> float:
        """
        Convert dynamic to static Youngs modulus using Fuller correlation.

        Equation type: Power law (y = a*x**b)

        Applicable for: Sandstones:
                          Consolidated (10%-15%)
                          Moderately Consolidated (15%-25%)
                          Weakly consolidated (>25%)

        Reference: R.H. Morales and R. P. Marcinew, Fracturing of High-Permeability Formations: Mechanical Properties Correlations, SPE 26561, October 1993.

        Args:
           yme_dyn (float): Dynamic Young's modulus magnitude Unit: GPa

        Returns:
           yme_sta_morales: Static Young's modulus magnitude Morales Unit: GPa

        Raises:
           ValueError: In case if we need to check values

        Example for input in:
           >>> StaticElasticProperties.dyn2sta_yme_najib(20)
           StaticElasticProperties(yme_sta_bradford =4.967602387961677)
        """

        yme_sta_power_law = multiplier * measurement**exponent

        return float(yme_sta_power_law)

    @staticmethod
    def dyn2sta_yme_bradord(yme_dyn: float) -> float:
        """
        Convert dynamic to static Youngs modulus using Bradord correlation

        Equation type: Power law (y = a*x**b)

        Applicable for: Turbiditic sandstones, Everest Field, North Sea.

        Reference: Bradford, I. D. R., Fuller, J., Thompson, P. J., & Walsgrove, T. R. (1998). Benefits of assessing the solids production risk in a North Sea reservoir using elastoplastic modelling (SPE/ISRM 47360). In SPE/ISRM. Eurorock '98 (pp. 261-269).

        Args:
           yme_dyn (float): Dynamic Young's modulus magnitude Unit: GPa

        Returns:
           yme_sta_bradord: Static Young's modulus magnitude Bradford Unit: GPa

        Raises:
           ValueError: In case if we need to check values

        Example for input in:
           >>> StaticElasticProperties.dyn2sta_yme_bradord(20)
           StaticElasticProperties(yme_sta_bradford =5.862103654131427)
        """
        multiplier = 0.0018
        exponent = 2.7
        yme_sta_bradford = multiplier * yme_dyn**exponent

        return float(yme_sta_bradford)

    @staticmethod
    def dyn2sta_yme_najib(yme_dyn: float) -> float:
        """
        Convert dynamic to static Youngs modulus using Najib correlation.

        Equation type: Power law (y = a*x**b)

        Applicable for: Carbonates from Iran (Asmari and Sarvak limestone)

        Reference: Najibi, A., Ghafoori, M., Lashkaripur, G. and Asef, M., 2015. Empirical relations between strength and static and dynamic elastic properties of Asmari and Sarvak limestones, two main oil reservoirs in Iran, Journal of Petroleum Science and Engineering 126 (2015) 78-82.

        Args:
           yme_dyn (float): Dynamic Young's modulus magnitude Unit: GPa

        Returns:
           yme_sta_najib: Static Young's modulus magnitude Bradford Unit: GPa

        Raises:
           ValueError: In case if we need to check values

        Example for input in:
           >>> StaticElasticProperties.dyn2sta_yme_najib(20)
           StaticElasticProperties(yme_sta_bradford =4.967602387961677)
        """
        multiplier = 0.014
        exponent = 1.96
        yme_sta_najib = multiplier * yme_dyn**exponent

        return float(yme_sta_najib)

    @staticmethod
    def dyn2sta_yme_fuller(yme_dyn: float) -> float:
        """
        Convert dynamic to static Youngs modulus using Fuller correlation.

        Equation type: Power law (y = a*x**b)

        Applicable for: Sandstone/Shale

        Reference: Techlog help file

        Args:
           yme_dyn (float): Dynamic Young's modulus magnitude Unit: GPa

        Returns:
           yme_sta_fuller: Static Young's modulus magnitude Fuller Unit: GPa

        Raises:
           ValueError: In case if we need to check values

        Example for input in:
           >>> StaticElasticProperties.dyn2sta_yme_najib(20)
           StaticElasticProperties(yme_sta_bradford =4.967602387961677)
        """
        multiplier = 0.032
        exponent = 1.632
        yme_sta_fuller = multiplier * yme_dyn**exponent

        return float(yme_sta_fuller)

    @staticmethod
    def dyn2sta_yme_morales(yme_dyn: float, porosity: float, exclude_low_por: bool = True) -> float:
        """
        Convert dynamic to static Youngs modulus using Fuller correlation.

        Equation type: Power law (y = a*x**b)

        Applicable for: Sandstones:
                          Consolidated (10%-15%)
                          Moderately Consolidated (15%-25%)
                          Weakly consolidated (>25%)

        Reference: R.H. Morales and R. P. Marcinew, Fracturing of High-Permeability Formations: Mechanical Properties Correlations, SPE 26561, October 1993.

        Args:
           yme_dyn (float): Dynamic Young's modulus magnitude Unit: GPa

        Returns:
           yme_sta_morales: Static Young's modulus magnitude Morales Unit: GPa

        Raises:
           ValueError: In case if we need to check values

        Example for input in:
           >>> StaticElasticProperties.dyn2sta_yme_najib(20)
           StaticElasticProperties(yme_sta_bradford =4.967602387961677)
        """

        if porosity < 0.15:
            multiplier = 1.271160
            exponent = 0.6612
        elif porosity < 0.25:
            multiplier = 0.957194
            exponent = 0.6920
        elif porosity > 0.25:
            multiplier = 0.153073
            exponent = 0.9404

        if porosity < 0.15 and exclude_low_por is True:
            yme_sta_morales = -9999
        else:
            yme_sta_morales = multiplier * yme_dyn**exponent

        return yme_sta_morales

    @staticmethod
    def dyn2sta_yme_modified_morales(yme_dyn: float, porosity: float, exclude_low_por: bool = True) -> float:
        """
        Convert dynamic to static Youngs modulus using Fuller correlation.

        Equation type: Power law (y = a*x**b)

        Applicable for: Sandstones:
                          Consolidated (10%-15%)
                          Moderately Consolidated (15%-25%)
                          Weakly consolidated (>25%)

        Reference: R.H. Morales and R. P. Marcinew, Fracturing of High-Permeability Formations: Mechanical Properties Correlations, SPE 26561, October 1993.

        Args:
           yme_dyn (float): Dynamic Young's modulus magnitude Unit: GPa

        Returns:
           yme_sta_morales: Static Young's modulus magnitude Morales Unit: GPa

        Raises:
           ValueError: In case if we need to check values

        Example for input in:
           >>> StaticElasticProperties.dyn2sta_yme_najib(20)
           StaticElasticProperties(yme_sta_bradford =4.967602387961677)
        """

        multiplier = 1.271160
        exponent = 0.6612

        yme_sta_morales = multiplier * yme_dyn**exponent

        return float(yme_sta_morales)

    @staticmethod
    def dyn2sta_yme_lacy(yme_dyn: float, porosity: float, exclude_low_por: bool = True) -> float:
        """
        Convert dynamic to static Youngs modulus using Fuller correlation.

        Equation type: Power law (y = a*x**b)

        Applicable for: Sandstones:
                          Consolidated (10%-15%)
                          Moderately Consolidated (15%-25%)
                          Weakly consolidated (>25%)

        Reference: R.H. Morales and R. P. Marcinew, Fracturing of High-Permeability Formations: Mechanical Properties Correlations, SPE 26561, October 1993.

        Args:
           yme_dyn (float): Dynamic Young's modulus magnitude Unit: GPa

        Returns:
           yme_sta_morales: Static Young's modulus magnitude Morales Unit: GPa

        Raises:
           ValueError: In case if we need to check values

        Example for input in:
           >>> StaticElasticProperties.dyn2sta_yme_najib(20)
           StaticElasticProperties(yme_sta_bradford =4.967602387961677)
        """

        multiplier = 1.271160
        exponent = 0.6612

        yme_sta_morales = multiplier * yme_dyn**exponent

        return float(yme_sta_morales)

    @staticmethod
    def dyn2sta_yme_custom_linear_law(measurement: float, slope: float, intercept: float) -> float:
        """
        Convert dynamic to static Youngs modulus using Fuller correlation.

        Equation type: Power law (y = a*x**b)

        Applicable for: Sandstones:
                          Consolidated (10%-15%)
                          Moderately Consolidated (15%-25%)
                          Weakly consolidated (>25%)

        Reference: R.H. Morales and R. P. Marcinew, Fracturing of High-Permeability Formations: Mechanical Properties Correlations, SPE 26561, October 1993.

        Args:
           yme_dyn (float): Dynamic Young's modulus magnitude Unit: GPa

        Returns:
           yme_sta_morales: Static Young's modulus magnitude Morales Unit: GPa

        Raises:
           ValueError: In case if we need to check values

        Example for input in:
           >>> StaticElasticProperties.dyn2sta_yme_najib(20)
           StaticElasticProperties(yme_sta_bradford =4.967602387961677)
        """
        yme_sta_linear_law = slope * measurement + intercept

        return yme_sta_linear_law

    @staticmethod
    def dyn2sta_yme_custom_exponential_law(measurement: float, multiplier: float, exponent: float) -> float:
        """
        Convert dynamic to static Youngs modulus using Fuller correlation.

        Equation type: Power law (y = a*x**b)

        Applicable for: Sandstones:
                          Consolidated (10%-15%)
                          Moderately Consolidated (15%-25%)
                          Weakly consolidated (>25%)

        Reference: R.H. Morales and R. P. Marcinew, Fracturing of High-Permeability Formations: Mechanical Properties Correlations, SPE 26561, October 1993.

        Args:
           yme_dyn (float): Dynamic Young's modulus magnitude Unit: GPa

        Returns:
           yme_sta_morales: Static Young's modulus magnitude Morales Unit: GPa

        Raises:
           ValueError: In case if we need to check values

        Example for input in:
           >>> StaticElasticProperties.dyn2sta_yme_najib(20)
           StaticElasticProperties(yme_sta_bradford =4.967602387961677)
        """

        yme_sta_morales = multiplier * yme_dyn**exponent

        return float(yme_sta_morales)

    @staticmethod
    def dyn2sta_yme_custom_constant_law(constant: float) -> float:
        """
        Custom constant law to assign a constant value for static Young's modulus

        Equation type: Constant law (y = c)

        Applicable for: generic

        Reference: -

        Args:
           constant (integer): Constant value in magnitude Unit: GPa

        Returns:
           yme_sta_constant: Static Young's modulus magnitude  Unit: GPa

        Raises:
           ValueError: In case if we need to check values

        Example for input in:
           >>> StaticElasticProperties.dyn2sta_yme_najib(20)
           StaticElasticProperties(yme_sta_bradford =4.967602387961677)
        """

        yme_sta_constant = constant

        return yme_sta_constant
