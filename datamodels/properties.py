import os
import sys
from dataclasses import dataclass
from typing import List


PYTHON_PATH = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(PYTHON_PATH)


from fluid.const import R
from fluid.nasa import molar_heat_capacity_at_constant_pressure as compute_Cp


@dataclass
class PropertyForm:

    def form(self):
        form_dict = {}
        for name, val in self.__dict__.items():
            if not (name.startswith("__")) or name.startswith("_"):
                form_dict[name] = val

        return form_dict

    @property
    def data(self):
        return {}
    
    def __str__(self):
        return str([f"{_["name"]} {_["symbol"]} = {_["val"]} [{_["units"]}]" for _ in self.data])[1:-1].replace(", ","\n").replace("'", "")

@dataclass
class Reference(PropertyForm):
    Tref: float
    Pref: float
    __name: str = "Reference State"

    @property
    def name(self):
        return self.__name
    
    @property
    def data(self):
        return [
            {"name": "Reference Temperature", "symbol": "Tref", "units": "K", "val": self.Tref},
            {"name": "Reference Pressure", "symbol": "Pref", "units": "Pa", "val": self.Pref}
        ]


@dataclass
class General(PropertyForm):
    MW: float
    mole_fraction: float = 1.0
    __name: str = "General"

    @property
    def name(self):
        return self.__name

    @property
    def data(self):
        return [
            {"name": "Molecular Weight", "symbol": "MW", "units": "m^3/mol", "val": self.MW},
            {"name": "Mole Fraction", "symbol": "y", "units": "-", "val": self.mole_fraction},
            {"name": "Specific Gas Constant", "symbol": "Rg", "units": "J/kg K", "val": self.Rg}
        ]

    @property
    def Rg(self):
        return R / self.MW
    

@dataclass
class PerfectGas(PropertyForm):
    Z: float
    rhog: float = None
    Cpg: float = None
    __name: str = "Perfect Gas"

    @property
    def name(self):
        return self.__name

    @property
    def data(self):
        return [
            {"symbol": "rhog", "name": "Density (gas)", "units": "m^3/kg", "val": self.rhog},
            {"symbol": "Z", "name": "Compressibility Factor (gas)","units": "-", "val": self.Z},
            {"symbol": "Pressure Constant Specific Heat Capasity (gas)", "units": "J/kg K", "val": self.Cpg}
        ]
    
@dataclass
class RealGas(PropertyForm):
    Tc: float
    Pc: float
    vc: float
    w: float
    NASA_Glenn_coeffs: List[float]
    __name: str = "Real Gas"

    def __get_coefficients(self, Tref, coeffs):
        # Get temperature range in which the Cp is determined
        T_range = [c[0] for c in coeffs]
        T_mix, T_max = T_range[0][0], T_range[-1][-1]

        if T_max < Tref < T_mix:
            raise RuntimeError(f"The Cp coefficients are not determined at the temperature {Tref} [K]")

        # Get Cp coefficients
        for coeff in coeffs:
            if coeff[0][0] <= Tref <= coeff[0][1]:
                return coeff[1]

        raise RuntimeError(f"Coeffcients are not found")

    @property
    def name(self):
        return self.__name
    
    @property
    def NASA_Glenn_T_range(self):
        T_range = [_[0] for _ in self.NASA_Glenn_coeffs]
        return T_range[0][0], T_range[-1][-1]

    @property
    def Cpg(self):
        def get_Cpg(T: float=288.15):
            a = self.__get_coefficients(T, self.NASA_Glenn_coeffs)
            return compute_Cp(T, a)
        return get_Cpg

    @property
    def data(self):
        return [
            {"name": "Critical Temperature", "symbol": "Tc", "units": "K", "val": self.Tc},
            {"name": "Critical Pressure", "symbol": "Pc", "units": "Pa", "val": self.Pc},
            {"name": "Critical Specific Molar Volume", "symbol": "vc", "units": "m^3/mol", "val": self.vc},
            {"name": "Acentric Factor", "symbol": "", "units": "m^3/mol", "val": self.vc},
            {"name": "Cp Temperature Range", "symbol": "Tmax - Tmin", "units": "K", "val": self.NASA_Glenn_T_range}

        ]


@dataclass
class Transport(PropertyForm):
    mu: float
    k: float
    __name = "Transport"

    @property
    def name(self):
        return self.__name
    
    @property
    def nu(self):
        def get_nu(rho: float=1.0):
            return self.mu / rho
        return get_nu

    @property
    def data(self):
        return [
            {"name": "Dynamic Viscosity", "symbol": "mu", "units": "Pa s", "val": self.mu},
            {"name": "Thermal Conductivity", "symbol": "k", "units": "W/m K", "val": self.k}
        ]


@dataclass
class Phase(PropertyForm):
    Ttr: float
    Ptr: float
    __coeffs: List[float]
    __name: str = "Phase Transition"

    @property
    def name(self):
        return self.__name

    @property
    def antonie_coefficients(self):
        return self.__coeffs

    @property
    def data(self):
        return [
            {"name": "Triple Point Temperature", "symbol": "Ttr", "units": "K", "val": self.Ttr},
            {"name": "Triple Point Pressure", "symbol": "Ptr", "units": "Pa", "val": self.Ptr},
        ]
