import os
import sys
import numpy as np


PYTHON_PATH = os.path.abspath(os.path.join(
    os.path.dirname(__file__), os.pardir
))
sys.path.append(PYTHON_PATH)

from fluid.const import R


def wilke(ys, xs, MWs) -> float:

    """
    Computes property of the mixture using Wilke's mixing rule

    Args:
        ys (list(float)): Mole fractions of each component in the mixture
        xs (list(float)): Parameter being computed for each omponent in the mixture
        MWs (list(float)): Molecular weights for each component in the mixture

    Returns:
        float: X paramter of the mixture
    """

    n = len(ys)
    x_mix = 0.0

    ys = np.array(ys)
    xs = np.array(xs)
    MWs = np.array(MWs)

    for i in range(n):
        denom = 0.0
        for j in range(n):
            if i == j:
                continue
            phi_ij = (1 + (xs[i] / xs[j])**0.5 * (MWs[j] / MWs[i])**0.25)**2 / np.sqrt(8 * (1 + MWs[i] / MWs[j]))
            denom += ys[j] * phi_ij
        x_mix += ys[i] * xs[i] / denom

    return x_mix


def kay(ys, xs):

    """
    Computes property of the mixture using Kay's mixing rule (simple weighted averaging)

    Args:
        ys (list(float)): Mole fractions of each component in the mixture
        xs (list(float)): Parameter being computed for each omponent in the mixture

    Returns:
        float: X paramter of the mixture
    """

    ys, xs = np.array(ys), np.array(xs)
    return sum(ys * xs)


def prausnitz_gunn(ys, Tcs, Vcs, ws, MWs, Tc_method="kay", Vc_method="kay") -> dict:

    """
    Computes Pseudo critical properties for the mixture using Prausnitz-Gunn approach

    Args:
        ys list(float): mole fractions of each component in the mixture
        Tcs list(float): Critical temperatures in K for each component in mixture
        Vcs list(float): Molar specific volumes in m^3/mol for each component in the mixture
        ws list(float): Acentric factors for each componet in the mixture
        MWs list(float): Molecular weight in kg/mol for each component in the mixture
        Tc_method (str): Method to mix critical temperature. Optional. Default "kay". Available methods "kay" and "binary"
        Vc_method (str): Method to mix critical molar specific volume. Optional. Default "kay". Available methods "kay" and "binary"
    
    Returns:
        dict: Dictionary of the pseudocritical parameters (Tc_mix [K], Pc_mix [Pa], Vc_mix [m^3/mol], rhoc_mix [kg/m^3], Zc_mix) of the mixture and acentric factor ws_mix of the mixture
    """

    n = len(ys)
    ys = np.array(ys)
    Tcs = np.array(Tcs)
    Vcs = np.array(Vcs)
    ws = np.array(ws)
    MWs = np.array(MWs)

    omega_mix = sum(ys * ws)

    if Tc_method == "kay":
        Tc_mix = kay(ys,  Tcs)
    elif Tc_method == "binary":
        Tc_mix = 0.0
        for i in range(n):
            for j in range(n):
                Tc_mix += ys[i] * ys[j] * np.sqrt(Tcs[i] * Tcs[j])

    if Vc_method == "kay":
        Vc_mix = kay(ys, Vcs)
    elif Vc_method == "binary":
        Vc_mix = 0.0
        for i in range(n):
            for j in range(n):
                Vc_ij = ((Vcs[i]**(1/3) + Vcs[j]**(1/3)) / 2)**3
                Vc_mix += ys[i] * ys[j] * Vc_ij
    
    Zc_mix = (0.2905 - 0.085 * omega_mix)
    Pc_mix = Zc_mix * R * Tc_mix / Vc_mix
    WM_mix = sum(ys * MWs)
    rhoc_mix = WM_mix / (Vc_mix)

    return {
        "Tc_mix": float(Tc_mix),
        "Pc_mix": float(Pc_mix),
        "Vc_mix": float(Vc_mix),
        "rhoc_mix": float(rhoc_mix),
        "Zc_mix": float(Zc_mix),
        "ws_mix": float(omega_mix)
    }
