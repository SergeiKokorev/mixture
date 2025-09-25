import os
import sys
import numpy as np
from typing import List, Tuple

PYTHON_PATH = os.path.abspath(os.path.join(
    os.path.dirname(__file__), os.pardir
))
sys.path.append(PYTHON_PATH)


from fluid.const import R


def __check_args__(args):
    if isinstance(args, float | int):
        return [args]
    elif hasattr(args, "__iter__"):
        for arg in args:
            arg = [arg]
    else:
        raise ValueError(f"Arguments must be int, float, list or tuple type. Give {type(args)}")
    return args



def redlich_kwong(T: float, P: float, ys: List[float], Tcs: List[float], Pcs: List[float], MWs: List[float], ws: List[float], vcs: List[float]=None, /, k_ij:List[List[float]]=None, eos: str="standard") -> Tuple[float]:

    """
    Computes Density, Compressibility, Acentric Factor of the mixture
    using family of Redlich-Kwong EoS

    Args:
        T (float): reference Temperatureof the mixture in K
        P (float): referencr Pressure of the mixture in Pa
        ys list(float): mole fractions of all gases in the mixture
        Tcs (list(float)): critical temperature of all gases in the mixture in K
        Pcs (list(float)): critical pressure of all gases in the mixture in Pa
        vcs (list(float)): critical specific molar volume of all gases in the mixture in m^3 mol^-1. Optional. Default None. Needed only for Aungier and Soave EoS.
        MWs (list(float): molecular weights of all gases in the mixture
        ws (list(float)): acentric factors of all gases in Mixture
        k_ij (list(float)): binary interaction parameters, optional. Defaul None, will be set to zero array
        eos (str): EoS, optional. Default "standard". Available eos "aungier" for Redlich-Kwong-Aungier EoS, "soave" for Soave Redlich-Kwong EoS and "standard for Standard Redlich-Kwong EoS

    Returns:
        tuple: Mixture compressibility factor, density and acrentric factor
    """

    # Initilize data
    n = len(ys)
    Tcs, Pcs, MWs, ws = np.array(Tcs), np.array(Pcs), np.array(MWs), np.array(ws)
    if k_ij is None:
        k_ij = np.zeros((n, n))
    else:
        k_ij = np.array(k_ij)

    w_mix = sum(ys * ws)
    MW_mix = sum(ys * MWs)

    b_i = 0.08664 * R * Tcs / Pcs
    a_i =( 0.42747 * R**2 * Tcs**2) / Pcs

    if eos == "standard":
        kappa_i = 0.5
        c_i = np.zeros((n,))
        alpha_i = a_i * (T / Tcs)**(-n)
    elif eos == "aungier":
        kappa_i = 0.4986 + 1.1735 * ws + 0.4754 * ws**2
        c_i = R * Tcs / (Pcs + a_i / (vcs * (vcs + b_i))) + b_i - vcs
        alpha_i = a_i * (T / Tcs)**(-kappa_i)
    elif eos == "soave":
        kappa_i = 0.480 + 1.574 * ws  - 0.176 * ws**2
        c_i = R * Tcs / (Pcs + a_i / (vcs * (vcs + b_i))) + b_i - vcs
        alpha_i = a_i * (1 + kappa_i * (1 - np.sqrt(T / Tcs)))**2

    # Define paramters of the mixture
    b_mix = sum(ys * b_i)
    c_mix = sum(ys * c_i)
    a_mix = 0.0

    for i in range(n):
        for j in range(n):
            a_cross = (1 - k_ij[i,j]) * np.sqrt(a_i[i] * alpha_i[i] * a_i[j] * alpha_i[j])
            a_mix += ys[i] * ys[j] * a_cross
        
    # Compute coefficients for solving a cubic equation for Z
    A = (a_mix * P) / (R * T)**2
    B = (b_mix * P) / (R * T)
    C = (c_mix * P) / (R * T)
    coeffs = [1, C - 1, B * (-B - 1 + C), A * (C - B)]

    # Compute Compressibility Factor and Density
    Z_roots = np.roots(coeffs)
    Z = max([root.real for root in Z_roots if abs(root.imag) < 1e-10])
    rho = (MW_mix * P) / (R * T * Z)
    return float(Z), float(rho), float(w_mix)
    


def peng_robinson(T: float, P: float, ys:List[float], Tcs: List[float], Pcs: List[float], MWs: List[float], ws: List[float], /, k_ij:List[List[float]]=None) -> Tuple[float]:

    """
    Computes Density, Compressibility, Acentric Factor of the mixture
    using Peng-Robinson EoS

    Args:
        T (float): reference Temperatureof the mixture in K
        P (float): referencr Pressure of the mixture in Pa
        ys list(float): mole fractions of all gases in the mixture
        Tcs (list(float)): critical temperature of all gases in the mixture in K
        Pcs (list(float)): critical pressure of all gases in the mixture in Pa
        MWs (list(float): molecular weights of all gases in the mixture
        ws (list(float)): acentric factors of all gases in Mixture
        k_ij (list(float)): binary interaction parameters, optional. Defaul None, will be set to zero array

    Returns:
        tuple: Mixture compressibility factor, density and acrentric factor
    """

    # Initialize data
    n = len(ys)
    Tcs, Pcs, MWs, ws = np.array(Tcs), np.array(Pcs), np.array(MWs), np.array(ws)
    if k_ij is None:
        k_ij = np.zeros((n, n))
    else:
        k_ij = np.array(k_ij)
    
    # Define parameters using a simple mole-fraction-weighted average
    omega_mix = sum(ys * ws)
    Mw_mix = sum(ys * MWs)
    
    # Define paramters of the mixture
    b_i = 0.0778 * R * T / Pcs
    a_i = 0.45724 * (R**2 * Tcs**2) / Pcs
    kappa_i = 0.37464 + 1.54226 * ws - 0.26993 * ws ** 2
    alpha_i = (1 + kappa_i * (1 - np.sqrt(T / Tcs)))**2

    b_mix = sum(ys * b_i)
    a_mix = 0.0

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            a_cross = (1 - k_ij[i,j]) * np.sqrt(a_i[i] * alpha_i[i] * a_i[j] * alpha_i[j])
            a_mix += ys[i] * ys[j] * a_cross
    
    # Compute cubic equation coefficints
    A = (a_mix * P) / (R * T)**2
    B = (b_mix * P) / (R * T)
    coeffs = [1, -(1 - B), (A - 2 * B - 3 * B**2), -(A * B - B**2 - B**3)]

    # Compute Compressibility Factor and Density
    Z_roots = np.roots(coeffs)
    Z = max([root.real for root in Z_roots if abs(root.imag) < 1e-10])
    rho = (Mw_mix * P) / (R * T * Z)
    return float(Z), float(rho), float(omega_mix)
