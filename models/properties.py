from __future__ import annotations

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

    """
    A simple base class for descibing inherited class properties and method

    Base dataclass provides convient way to get data for displaying GUI from

    Attributes:
        data (dict): Dictionary of the main class properties

    """

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
    """
    Dataclass provides convient way to describe and change reference state

    Attributes:
        Tref (float): Reference temperature in K
        Pref (float): Reference pressure in Pa
        __name (str): Desribes how class will be represented in GUI from
    """
    Tref: float
    Pref: float
    __name: str = "Reference State"

    @property
    def name(self):
        """
        Returns incapsulated attribute self.__name of the class

        Returns:
            str: The name to display class in the GUI form
        """
        return self.__name
    
    @property
    def data(self) -> List[dict]:
        """
        Returns list of the class main attributes to display them in the GUI form

        Returns:
            list: List of dictionaries describing class main attributes
        """
        return [
            {"name": "Reference Temperature", "symbol": "Tref", "units": "K", "val": self.Tref},
            {"name": "Reference Pressure", "symbol": "Pref", "units": "Pa", "val": self.Pref}
        ]


@dataclass
class General(PropertyForm):
    """
    Provides convinient way to describe the general properties of the fluid

    Attributes:
        MW (float): Molecular weight in kg/mol
        y (float): Mole fraction of a fluid in the mixture. Optional. Default 1.0, a pure gas.
        __name (str): Desribes how class will be represented in GUI from
    """
    MW: float
    y: float = 1.0
    __name: str = "General"

    @property
    def name(self):
        """
        Returns incapsulated attribute self.__name of the class

        Returns:
            str: The name to display class in the GUI form
        """
        return self.__name

    @property
    def data(self):
        """
        Returns list of the class main attributes to display them in the GUI form

        Returns:
            list: List of dictionaries describing class main attributes
        """
        return [
            {"name": "Molecular Weight", "symbol": "MW", "units": "m^3/mol", "val": self.MW},
            {"name": "Mole Fraction", "symbol": "y", "units": "-", "val": self.y},
            {"name": "Specific Gas Constant", "symbol": "Rg", "units": "J/kg K", "val": self.Rg}
        ]

    @property
    def Rg(self) -> float:
        """
        Specific gas constant in J/(kg K). To be computed as Ru[J/(mol K)] / MW[kg/mol],
        where Ru is the universal gas constant [J/(mol K)] and 
        MW is the molecular weight [kg/mol]
        """
        return R / self.MW
    

@dataclass
class PerfectGas(PropertyForm):
    """
    Provides properties of perfect gases

    Attributes:
        NASA_Glenn_coeffs list: The NASA Glenn coeffcients and temperature range to use it to compute demensionless constant pressure specific heat capacity
        __name (str): Desribes how class will be represented in GUI from
    """    
    NASA_Glenn_coeffs: List[float]
    __name: str = "Perfect Gas"

    @property
    def name(self):
        """
        Returns incapsulated attribute self.__name of the class

        Returns:
            str: The name to display class in the GUI form
        """
        return self.__name

    @property
    def data(self):
        """
        Returns list of the class main attributes to display them in the GUI form

        Returns:
            list: List of dictionaries describing class main attributes
        """
        return [
            {"symbol": "Cp coeff","name": "Cp coefficients", "units": "-", "val": self.NASA_Glenn_coeffs},
        ]


@dataclass
class RealGas(PropertyForm):
    """
    Provides convinient way to describe real gas properties

    Attributes:
        Tc (float): Critical temperature in K
        Pc (float): Critical pressure in Pa
        vc (float): Specific molar volume in m^3/mol
        Zg (float): Compressibility factor of gase phase. Optional. Default 1.0, ideal gase case
        rhog (float): Density of a gas phase in kg/m^3. Optional. Need to be commputed in class Fluid (for pure substances) or Mixture (for gases compounds)
        w (float): Acentric factor
        __name (str): Desribes how class will be represented in GUI from
    """
    Tc: float
    Pc: float
    vc: float
    w: float
    rhog: float = None
    Zg: float = 1.0
    __name: str = "Real Gas"

    @property
    def name(self):
        """
        Returns incapsulated attribute self.__name of the class

        Returns:
            str: The name to display class in the GUI form
        """
        return self.__name
    
    @property
    def data(self):
        """
        Returns list of the class main attributes to display them in the GUI form

        Returns:
            list: List of dictionaries describing class main attributes
        """
        return [
            {"name": "Critical Temperature", "symbol": "Tc", "units": "K", "val": self.Tc},
            {"name": "Critical Pressure", "symbol": "Pc", "units": "Pa", "val": self.Pc},
            {"name": "Critical Specific Molar Volume", "symbol": "vc", "units": "m^3/mol", "val": self.vc},
            {"symbol": "rhog", "name": "Density (gas)", "units": "m^3/kg", "val": self.rhog},
            {"name": "Acentric Factor", "symbol": "w", "units": "-", "val": self.w},
            {"name": "Compressibility Factor (gas)", "symbol": "Zg", "units": "-", "val": self.Zg}

        ]


@dataclass
class Transport(PropertyForm):
    """
    Provides convinient way to describe transport properties of gases

    Attributes:
        mu (float): Dynamic viscosity in Pa s
        k (float): Thermanl conductivity in W/(m K)
        __name (str): Desribes how class will be represented in GUI from
    """
    mu: float = None
    k: float = None
    __name = "Transport"

    @property
    def name(self):
        """
        Returns incapsulated attribute self.__name of the class

        Returns:
            str: The name to display class in the GUI form
        """
        return self.__name
    
    @property
    def nu(self):
        """
        Returns kinematic viscosity in m^2/s.

        Returns:
            float: Kinamatic viscosity in m^2/s
        """
        def get_nu(rho: float=1.0):
            return self.mu / rho
        return get_nu

    @property
    def data(self):
        """
        Returns list of the class main attributes to display them in the GUI form

        Returns:
            list: List of dictionaries describing class main attributes
        """
        return [
            {"name": "Dynamic Viscosity", "symbol": "mu", "units": "Pa s", "val": self.mu},
            {"name": "Thermal Conductivity", "symbol": "k", "units": "W/m K", "val": self.k}
        ]


@dataclass
class Phase(PropertyForm):
    """
    Provides convinient way to describe triple point of gases

    Attributes:
        Ttr (float): Triple point temperature in K
        Ptr (float): Triple point pressure in Pa
        __coeffs (list): Antonie coefficients to compute the bubble (boiling) temperature
        __name (str): Desribes how class will be represented in GUI from
    """
    Ttr: float
    Ptr: float
    __coeffs: List[float]
    __name: str = "Phase Transition"

    @property
    def name(self):
        """
        Returns incapsulated attribute self.__name of the class

        Returns:
            str: The name to display class in the GUI form
        """
        return self.__name

    @property
    def antonie_coefficients(self):
        return self.__coeffs

    @property
    def data(self):
        """
        Returns list of the class main attributes to display them in the GUI form

        Returns:
            list: List of dictionaries describing class main attributes
        """
        return [
            {"name": "Triple Point Temperature", "symbol": "Ttr", "units": "K", "val": self.Ttr},
            {"name": "Triple Point Pressure", "symbol": "Ptr", "units": "Pa", "val": self.Ptr},
        ]
