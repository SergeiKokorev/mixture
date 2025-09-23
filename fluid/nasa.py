"""
Provides emperical equations for fitting thermodynamic functions\n

Molar Heat Capacity at Constant Pressure dimensionless form\n
Cp0(T) / Ru = a1 * T ** -2 + a2 * T ** -1 + a3 + a4 * T + a5 * T ** 2 + a6 * T ** 3 + a7 * T ** 4\n\n

Enthalpy\n
H0(T) / (Ru * T) = -a1 * T ** -1 + a2 * ln(T) / T + a3 + a4 * T / 2 + a5 * T ** 2 / 3 + a6 * T ** 3 / 4 + a7 * T ** 4 / 5 + b1 / t\n\n

Entropy\n
S0(T) / Ru = -a1 * T ** -2 / 2 - a2 * T ** -1 + a3 * ln(T) + a4 * T + a5 * T ** 2 / 2 + a6 * T ** 3 / 3 + a7 * T ** 4 / 4 + b2\n

where   Ru is the universal gas constant, Ru = 8.3145 [J * mol -1 * K -1]
        a1 ... a7 are the interpolation coefficients
        b1 ... b2 are the integration coefficients
        T is the temperature, [K]

NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species\n
Bonnie J. McBride, Michael J. Zehe, and Sanford Gordon.\n
Glenn Research Center, Cleveland, Ohio, 2002
"""

from typing import List
from math import log


def molar_heat_capacity_at_constant_pressure(temperature: float, a: List[float]) -> float:

    """
    Computes Demensionless Molar Heat Capacity at constant Pressure (Cp0(T) / Ru)

    :param temperature: float, Temperature in Kelvin at which molar heat capacity to be calculated
    :param a: List[float], NASA approximation coefficientsto compute molar heat capacity, number of coeficients is always 7
    
    :return cp0(T) / Ru: float, Demensionless Molar Heat Capacity at Constant Pressure
    """

    if len(a) < 7:
            raise ValueError(f'Number of approximation coefficient must be 7. Given {len(a)}')
    
    return sum([ai * temperature ** (i - 2) for i, ai in enumerate(a)])


def entropy(temperature: float, a: List[float], b: float) -> float:
      
    """
    Computes Demensionless Entropy (S0(T) / Ru)

    :param temperature: float, Temperature in Kelvin at which entropy to be calculated
    :param a: List[float], NASA approximation coefficients to calculate entropy, number of coefficients is always 7
    :param b: float, First NASA integration coeffcient to calculate entropy
    
    :return S0(T) / Ru: float, Demensionless Entropy
    """

    if len(a) < 7:
        raise ValueError(f'Number of approxiamtion coefficients must be 7. Given {len(a)}')
    
    s0 = (-a[0] * temperature ** -2) / 2 - a[1] * temperature ** -1 + a[2] * log(temperature)
    s0 += sum([(ai * temperature ** i) / i for i, ai in enumerate(a[3:], 1)]) + b
    
    return s0    


def molar_enthalpy(temperature: float, a: List[float], b: float) -> float:
    
    """
    Computes Demensionless Enthalpy H0(T) / (Ru * T)

    :param temperature: float, Temperature in Kelvin at which Enthalpy to be calculated
    :param a: List[float], NASA approxiamtion coefficients to calculate enthalpy, number of coefficients is always 7
    :param b: float, Second NASA integration coefficient

    :return H(T) / (Ru * T): float, Demensionles Enthalpy
    """

    if len(a) < 7:
        raise ValueError(f'Number of approximation coefficients must be 7. Given {len(a)}')
    
    h0 = -a[0] * temperature ** -2 + a[1] * log(temperature) * temperature ** -1
    h0 += sum([ai * temperature ** i / (i + 1) for i, ai in enumerate(a[2:])]) + b / temperature
    
    return h0


def compute(temperature: float, a: List[float], b: List[float]) -> dict:
     
    """
    Computes demensionless parameters using NASA approximation (Molar Heat Capacity at Constant Pressure,
    Entropy and Enthalpy)

    :param temperature: float, Temperature at which parameters to be calculated
    :param a: List[float], NASA approximation coefficients, number of coefficients is always 7
    :param b: List[float], NASA integration coefficients, number of coeffcients is always 2

    :return parameters: dict, Demensionless Molar Heat Capacity at Constant Pressure cp0,
    Enthalpy h0, Entropy s0
    """

    if len(a) < 7:
        raise ValueError(f'Number of approximation coefficients must be 7. Given {len(a)}')

    if len(b) < 2:
        raise ValueError(f'Number of integration coefficients must be 2. Given {len(b)}')

    return {
        'cp0': molar_heat_capacity_at_constant_pressure(temperature, a),
        'h0': molar_enthalpy(temperature, a, b[0]),
        's0': entropy(temperature, a, b[1])
    }
