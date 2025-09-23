import os
import sys
import numpy as np
import numbers

from typing import Tuple


PYTHON_PATH = os.path.abspath(os.path.join(
    os.path.dirname(__file__), os.pardir
))
sys.path.append(PYTHON_PATH)


from fluid.const import R


def redlich_kwong(T: float, P: float, Mw:float, Tc: float, Pc: float, rhoc: float = None, omega: float = None, model: str="standard") -> float:

    a0 = 0.42747 * R**2 * Tc**2 / Pc
    b = 0.08664 * R * Tc / Pc

    if model == "standard":
        n = 0.5
        c = 0.0
        a = a0 * (T / Tc)**(-n)
    elif model == "aungier":
        
        if  not (rhoc is None and omega is None):
            RuntimeError("rho and/or omega is not provided")

        vc = 1 / rhoc
        n = 0.4986 + 1.1735 * omega + 0.4754 * omega**2
        c = R * Tc / (Pc + a0 / (vc * (vc + b))) + b - vc
        a = a0 * (T / Tc)**(-n)
        if T < Tc:
            coeffs = [-3.80666e+3, 6.59754e+1, -3.92603e-1, 6.11597e-4, 2.74395e-6, -1.18587e-8, 1.26942e-11]
            n_liquid = sum([ai * T ** i for i, ai in enumerate(coeffs)])
            a *= (Tc / T)**(n_liquid / n)
        vm_roots = np.roots([P, P * c - R * T, P * b - R * T * b + a, -a * b + a * c])
        rho = float(Mw / max(vm_roots.real))
        Z = P * Mw / (R * T * rho)
        return Z, rho
    elif model == "soave":
        
        if  not (rhoc is None and omega is None):
            RuntimeError("rho and/or omega is not provided")

        n = 0.48 + 1.574 * omega - 0.176 * omega ** 2
        a = a0 * (1 + n * (1 - (T / Tc)**0.5))**2
    else:
        raise RuntimeError(f"EOS model has not been determined")

    A = a * P / (R * T)**2
    B = b * P / (R * T)
    Z_roots = np.roots([1, -1, (A - B - B**2), -A * B])
    Z = float(max(Z_roots.real))
    rho = P * Mw / (R * T * Z)

    return Z, rho


def peng_robinson(T, P, Mw, Tc, Pc, omega) -> Tuple[float]:
    
    b = 0.0778 * R * Tc / Pc
    a0 = 0.45724 * (R * Tc)**2 / Pc
    n = 0.37464 + 1.54226 * omega - 0.26993 * omega**2
    a = a0 * (1 + n * (1 - np.sqrt(T / Tc))) ** 2

    A = float(a * P / (R * T)**2)
    B = float(b * P / (R * T))
    coeffs = [1, -(1 - B), (A - 2 * B - 3 * B ** 3), -(A * B - B ** 2 - B ** 3)]
    Z_roots = np.roots(coeffs)
    Z = float(max(Z_roots.real))
    rho = P * Mw /(Z * R * T)

    return Z, rho
