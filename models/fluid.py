from __future__ import annotations
import os
import sys
from dataclasses import dataclass




PYTHON_PATH = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(PYTHON_PATH)


from models.properties import *
from fluid.nasa import molar_heat_capacity_at_constant_pressure as compute_cp
from fluid.eos import *

EoS = [
    {"eos": "redlich_kwong_standard", "fun": redlich_kwong, "args": {"method": "standard"}},
    {"eos": "redlich_kwong_soave", "fun": redlich_kwong, "args": {"method": "soave"}},
    {"eos": "redlich_kwong_aungier", "fun": redlich_kwong, "args": {"method": "aungier"}},
    {"eos": "peng_robinson", "fun": peng_robinson, "args": ()}
]


@dataclass
class Chemical:
    __name: str
    __formula: str
    __ref: Reference
    __general: General
    __perfect_gas: PerfectGas
    __real_gas: RealGas

    def __get_CP_coeffs(self, Tref, coeffs):
        # Get temperature range in which the Cp is determined
        T_mix, T_max = self.NASA_Glenn_T_range

        if T_max < Tref < T_mix:
            raise RuntimeError(f"The Cp coefficients for chemical {self.name} ({self.formula}) are not determined at the temperature {Tref} [K]")

        # Get Cp coefficients
        for coeff in coeffs:
            if coeff[0][0] <= Tref <= coeff[0][1]:
                return coeff[1]

        raise RuntimeError(f"Coeffcients are not found")

    @property
    def name(self):
        return self.__name
    
    @property
    def formula(self):
        return self.__formula

    @property
    def T(self):
        return self.__ref.Tref
    
    @property
    def P(self):
        return self.__ref.Pref

    @property
    def MW(self):
        return self.__general.MW

    @property
    def y(self):
        return self.__general.y

    @property
    def Rg(self):
        return self.__general.Rg

    @property
    def Cpg(self):
        def get_Cpg(T:float=self.ref.Tref):
            a = self.__get_CP_coeffs(T, self.__perfect_gas.NASA_Glenn_coeffs)
            return compute_cp(T, a)
        return get_Cpg

    @property
    def NASA_Glenn_T_range(self):
        """
        Returns temperature range in which NASA Glenn coefficients can be used with appropriate accuracy

        Returns:
            tuple: The temperature range to use NASA Glenn coefficients
        """
        T_range = [_[0] for _ in self.__perfect_gas.NASA_Glenn_coeffs]
        return T_range[0][0], T_range[-1][-1]

    @property
    def NASA_Glenn_coefficients(self):
        def get_a(T):
            return self.__get_CP_coeffs(T, self.__perfect_gas.NASA_Glenn_coeffs)
        return get_a

    @property
    def Tc(self):
        return self.__real_gas.Tc
    
    @property
    def Pc(self):
        return self.__real_gas.Pc
    
    @property
    def vc(self):
        return self.__real_gas.vc
    
    @property
    def rhoc(self):
        return self.__general.MW / self.__real_gas.vc
    
    @property
    def w(self):
        return self.__real_gas.w

    @property
    def rhog(self):
        def compute(T=None, P=None, equation_of_state="redlich_kwong_aungier", **kwargs):
            args = []
            for e in EoS:
                if e["eos"] == equation_of_state:
                    fun = e["fun"]
                    if (a := e.get("args")):
                        args.append(a)
            T = [self.__ref.Tref] if T is None else [T]
            P = [self.__ref.Pref] if P is None else [P]
            MW, y = [self.__general.MW], [self.__general.y]
            w = [self.__real_gas.w]
            Tc, Pc, vc = [self.__real_gas.Tc], [self.__real_gas.Pc], [self.__real_gas.vc]
            if (k_ij := kwargs.get("k_ij")):
                args.append(k_ij)
            Z, rho, w = fun(T, P, y, Tc, Pc, MW, w, *args)
            return Z, rho, w
        return compute[1]
    
    @property
    def Zg(self):
        def compute(T=None, P=None, equation_of_state="redlich_kwong_aungier", **kwargs):
            args = []
            for e in EoS:
                if e["eos"] == equation_of_state:
                    fun = e["fun"]
                    if (a := e.get("args")):
                        args.append(a)
            T = [self.__ref.Tref] if T is None else [T]
            P = [self.__ref.Pref] if P is None else [P]
            MW, y = [self.__general.MW], [self.__general.y]
            w = [self.__real_gas.w]
            Tc, Pc, vc = [self.__real_gas.Tc], [self.__real_gas.Pc], [self.__real_gas.vc]
            if (k_ij := kwargs.get("k_ij")):
                args.append(k_ij)
            Z, rho, w = fun(T, P, y, Tc, Pc, MW, w, *args)
            return Z, rho, w
        return compute[0]

    @T.seeter
    def T(self, T:float):
        if not isinstance(T, float | int):
            raise ValueError(f"Temperature must be float type. Given {type(T)}")
        self.__ref.Tref = T

    @P.seeter
    def P(self, P:float):
        if not isinstance(P, float | int):
            raise ValueError(f"Pressure must be float type. Given {type(P)}")
        self.__ref.Tref = P

    @y.setter
    def y(self, y: float):
        if not isinstance(y, float | int):
            raise ValueError(f"Molecular weight must be float type. Given {type(y)}")
        self.__general.y = y
