import os
import sys
import json


PYTHON_PATH = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))

sys.path.append(PYTHON_PATH)


from models.properties import *
from models.fluid import Chemical


if __name__ == "__main__":

    datafile = os.path.abspath(os.path.join(PYTHON_PATH, "lib", "thermo.json"))

    with open(datafile, "r") as fp:
            database = json.load(fp)

    T, P = 326.15, 7344000
    e = {"name": "Methane", "symbol": "CH4", "mole_fraction": 1.64 / 100, "Tc": 190.6, "Pc": 4.64e6, "rhoc": 162.66, "Tb": 111.65,"omega": 0.011, "sigma": 3.758, "epsilon_kb": 148.6, "lam": 0.0332}

    ref = Reference(T, P)
    symbol = e["symbol"]
    name = e["name"]
    mw = database[symbol]["molecular_weight"]*1e-3
    gen = General(mw)
    perfect_gas = PerfectGas(database[symbol]["coefficients"])
    vc = mw / e["rhoc"]
    real_gas = RealGas(e["Tc"], e["Pc"], vc, e["omega"])
    transport = Transport()

    chemical = Chemical(name, symbol, ref, gen, perfect_gas, real_gas)
    print(chemical.Zg())

