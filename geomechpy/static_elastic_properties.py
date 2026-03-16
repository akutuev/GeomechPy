import math
from dataclasses import dataclass
    

class StaticElasticPropertiesConverter:
   """Convert dynamic elastic properties to static elastic properties

   Reference:
       Zhang, Yuliang, et al. "Extracting static elastic moduli of rock through elastic wave velocities." Acta Geophysica 72.2 (2024): 915-931.
   """

   @staticmethod
   def dyn2sta_yme_bradord(yme_dyn: float) -> float:
      """
      Convert dynamic to static Youngs modulus using Bradord correlation
      Equation type: Power law (y = a*x**b)
      Applicable for: Turbiditic sandstones, Everest Field, North Sea.

      Reference: Bradford, I. D. R., Fuller, J., Thompson, P. J., & Walsgrove, T. R. (1998). Benefits of assessing the solids production risk in a North Sea reservoir using elastoplastic modelling (SPE/ISRM 47360). In SPE/ISRM. Eurorock '98 (pp. 261-269).

      Args:
         yme_dyn (float): Dynamic Young's modulus magnitude Unit: Mpsi

      Returns:
         yme_sta_bradord: Static Young's modulus magnitude Bradford Unit: Mpsi
      """
      
      multiplier = 0.04794626440600849
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
         yme_dyn (float): Dynamic Young's modulus magnitude Unit: Mpsi

      Returns:
         yme_sta_najib: Static Young's modulus magnitude Bradford Unit: Mpsi
      """
      multiplier = 0.07277314417314575
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
      """
      multiplier = 0.08143824177457351
      exponent = 1.632
      yme_sta_fuller = multiplier * yme_dyn**exponent

      return float(yme_sta_fuller)

   @staticmethod
   def dyn2sta_yme_morales(yme_dyn: float, porosity: float, exclude_low_por: bool = True) -> float:
      """
      Convert dynamic to static Youngs modulus using Morales correlation.
      Equation type: Power law (y = a*x**b)
      Applicable for: Sandstones:
                        Consolidated (10%-15%)
                        Moderately Consolidated (15%-25%)
                        Weakly consolidated (>25%)

      Reference: R.H. Morales and R. P. Marcinew, Fracturing of High-Permeability Formations: Mechanical Properties Correlations, SPE 26561, October 1993.

      Args:
         yme_dyn (float): Dynamic Young's modulus magnitude Unit: Mpsi

      Returns:
         yme_sta_morales: Static Young's modulus magnitude Morales Unit: Mpsi
      """

      if porosity < 0.15:
         multiplier = 2.562214651764409
         exponent = 0.6612
      elif porosity < 0.25:
         multiplier = 0.5281242638335621
         exponent = 0.6920
      elif porosity > 0.25:
         multiplier = 0.3467028522374105
         exponent = 0.9404

      if porosity < 0.10 and exclude_low_por is True:
         yme_sta_morales = -9999
      else:
         yme_sta_morales = multiplier * yme_dyn**exponent

      return yme_sta_morales

   @staticmethod
   def convert_dyn2sta_yme_custom_power_law(measurement: float, multiplier: float, exponent: float) -> float:
      """
      Convert dynamic to static Youngs modulus using custom power law.
      Equation type: Power law (y = a*x**b)
      Applicable for: Generic
      
      Args:
         yme_dyn (float): Dynamic Young's modulus magnitude Unit: Mpsi

      Returns:
         yme_sta_morales: Static Young's modulus magnitude Morales Unit: Mpsi
      """

      yme_sta_power_law = multiplier * measurement**exponent

      return float(yme_sta_power_law)

   @staticmethod
   def dyn2sta_yme_custom_linear_law(measurement: float, slope: float, intercept: float) -> float:
      """
      Convert dynamic to static Youngs modulus using custome linear law
      Equation type: Linear law (y = a*x + b)
      Applicable for: generic

      Args:
         yme_dyn (float): Dynamic Young's modulus magnitude Unit: Mpsi

      Returns:
         yme_sta_morales: Static Young's modulus magnitude Morales Unit: Mpsi
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
      """

      yme_sta_morales = multiplier * yme_dyn**exponent

      return float(yme_sta_morales)

   @staticmethod
   def dyn2sta_yme_custom_constant_law(constant: float) -> float:
      """
      Convert dynamic Poisson's ratio to static Poisson's ratio using a constant multiplier
      Equation type: Constant law (y = c)
      Applicable for: generic

      Reference: -

      Args:
         constant (integer): Constant value in magnitude Unit: Mpsi

      Returns:
         yme_sta_constant: Static Young's modulus magnitude  Unit: Mpsi
      """

      yme_sta_constant = constant

      return yme_sta_constant
    
   @staticmethod
   def dyn2sta_poissons_ratio(pr_dyn: float, multiplier: float) -> float:
      """
      Custom constant law to assign a constant value for static Young's modulus

      Equation type: Constant law (y = a*b)

      Applicable for: generic

      Reference: -

      Args:
         pr_dyn (float): Dynamic Poisson's ratio Unit: [unitless]
         b (float): constant multiplier Default set to 1.

      Returns:
         pr_sta (float): Static Poisson's ratio Unit: [unitless]
      """

      yme_sta_constant = constant

      return yme_sta_constant

   @staticmethod
   def biot_coefficient_constant_law(constant: float) -> float:
     
     """
        Biot coefficient

        Equation type: Constant law (y = c)

        Applicable for: generic

        Reference: Jaeger, John Conrad, Neville GW Cook, and Robert Zimmerman. Fundamentals of rock mechanics. John Wiley & Sons, 2009.

        Args:
           constant (integer): Constant value for Biot's coefficient: unitless

        Returns:
           alpha_const: Biot coefficient  Unit: unitless
        """

      biot_constant = constant

      return biot_constant