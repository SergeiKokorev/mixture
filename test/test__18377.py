import os
import sys
import json
import numpy as np
import thermo

from scipy.optimize import curve_fit

PYTHON_PATH = os.path.abspath(os.path.join(
    os.path.dirname(__file__), os.pardir
))
sys.path.append(PYTHON_PATH)

import fluid.eos as eos
import fluid.mixing as mixing
import fluid.nasa as nasa

from fluid.const import R


def mole_to_mass_fraction(mixture):

    d = sum(fluid["mole_fraction"] * fluid["mw"] for fluid in mixture)
    fluids = []
    for fluid in mixture:
        mass_fraction = fluid["mole_fraction"] * fluid["mw"] / d
        fluids.append((mass_fraction, fluid["name"], fluid["symbol"]))
    return fluids


def get_Cp_coefficients(Tref, coeffs, name, symbol=None) -> list:
    # Get temperature range in which the Cp is determined
    T_range = [c[0] for c in coeffs]
    T_mix, T_max = T_range[0][0], T_range[-1][-1]

    if T_max < Tref < T_mix:
        raise RuntimeError(f"The Cp coefficients are not determined for {name} ({symbol}) at the temperature {Tref} [K]")

    # Get Cp coefficients
    for coeff in coeffs:
        if coeff[0][0] <= Tref <= coeff[0][1]:
            return coeff[1]

    raise RuntimeError(f"Component {name} ({symbol}). Coeffcients are not found")


def objective_function(T, a0, a1, a2, a3, a4, a5, a6, a7):
    return a0 + a1 * T + a2 * T**2 + a3 * T**3 + a4 * T**4 + a5 * T**5 + a6 * T**6 + a7 * T**7


def polynomial_objective(T, a1, a2, a3, a4, a5):
    return a1 + a2 * T + a3 * T**2 + a4 * T**3 + a5 * T**4


if __name__ == "__main__":

    nasa_glenn_data_file = os.path.abspath(os.path.join(PYTHON_PATH, "lib", "thermo.json"))
    csv_path = os.path.abspath(os.path.join(PYTHON_PATH, "outdata"))
    
    with open(nasa_glenn_data_file, "r") as fp:
        database = json.load(fp)

    fluid_list = [
        {"name": "Hydrogen", "symbol": "H2", "mole_fraction": 96.24 / 100, "Tc": 32.95, "Pc": 1.297e6, "rhoc": 10**-26, "Tb": 20.25, "omega": -0.218, "sigma": 2.827, "epsilon_kb": 59.7, "lam": 0.180},
        {"name": "Nitrogen", "symbol": "N2", "mole_fraction": 100 * 1e-6},
        {"name": "Methane", "symbol": "CH4", "mole_fraction": 1.64 / 100, "Tc": 190.6, "Pc": 4.64e6, "rhoc": 162.66, "Tb": 111.65,"omega": 0.011, "sigma": 3.758, "epsilon_kb": 148.6, "lam": 0.0332},
        {"name": "Ethane", "symbol": "C2H6", "mole_fraction": 0.3 / 100, "Tc": 305.322, "Pc": 4.872e6, "rhoc": 206.6, "Tb": 184.65,"omega": 0.099, "sigma": 4.443, "epsilon_kb": 215.7, "lam": 0.0185},
        {"name": "Propane", "symbol": "C3H8", "mole_fraction": 0.62 / 100, "Tc": 369.8, "Pc": 4.25e6, "rhoc": 259.75, "Tb": 231.05,"omega": 0.152, "sigma": 4.807, "epsilon_kb": 248.9, "lam": 0.0146},
        {"name": "Isobutane", "symbol": "C4H10,isobutane", "mole_fraction": 0.28 / 100, "Tc": 408.13, "Pc": 3.64e6, "rhoc": 221, "Tb": 261.45,"omega": 0.1835, "sigma": 5.278, "epsilon_kb": 330.1, "lam": 0.0135},
        {"name": "n-butane", "symbol": "C4H10,n-butane", "mole_fraction": 0.46 / 100, "Tc": 425.1, "Pc": 3.796e6, "rhoc": 227.93, "Tb": 273.1,"omega": 0.2008, "sigma": 4.687, "epsilon_kb": 531.4, "lam": 0.0152},
        {"name": "n-hexane", "symbol": "C6H14,n-hexane", "mole_fraction": 0.28 / 100, "Tc": 507.6, "Pc": 3.025e6, "rhoc": 665, "Tb": 342,"omega": 0.329, "sigma": 5.949, "epsilon_kb": 399.3, "lam": 0.0130},
        {"name": "water", "symbol": "H2O", "mole_fraction": 0.17 / 100, "Tc": 647, "Pc": 22.064e6, "rhoc": 322, "Tb": 373.15, "omega": 0.344, "sigma": 2.641, "epsilon_kb": 809.1, "lam": 0.0181},
    ]

    T, P = 321, 7443000
    fluids = [fluid for fluid in fluid_list]

    print(f"Reference {T = } [K]; {P = } [Pa]")
    my_mixture = []
    thermo_fluids= []
    d = 0.0
    for i, fluid in enumerate(fluids):
        name = fluid.get("name")
        symbol = fluid.get("symbol")

        nasa_coeffs = database[symbol]["coefficients"]
        mole_fraction = fluid.get("mole_fraction")

        e = thermo.Chemical(name)
        Tc, Pc, rhoc, omega = e.Tc, e.Pc, e.rhoc, e.omega
        mw = e.MW
        thermo_fluids.append((name, mole_fraction))
        d += mole_fraction * mw
        my_mixture.append(
            {
                "name": name, "symbol": symbol,
                "Tc": Tc, "Pc": Pc, "rhoc": rhoc, "omega": omega,
                "mole_fraction": mole_fraction, "mw": mw * 1e-3,
                "nasa_coefficients": nasa_coeffs
            }
        )
    
    mass_fractions = mole_to_mass_fraction(my_mixture)

    for fluid in mass_fractions:
        print(fluid)

    mixture = thermo.Mixture(IDs=[fluid[0] for fluid in thermo_fluids], zs=[fluid[1] for fluid in thermo_fluids])

    ys = np.array([fluid["mole_fraction"] for fluid in my_mixture])
    Tcs = np.array([fluid["Tc"] for fluid in my_mixture])
    Pcs = np.array([fluid["Pc"] for fluid in my_mixture])
    rhocs = np.array([fluid["rhoc"] for fluid in my_mixture])
    MWs = np.array([fluid["mw"] for fluid in my_mixture])
    ws = np.array([fluid["omega"] for fluid in my_mixture])
    Vcs = np.array([fluid["mw"] / fluid["rhoc"] for fluid in my_mixture])

    
    Tc, Pc, Vc, rhoc, Zc, omega = mixing.prausnitz_gunn(ys, Tcs, Vcs, ws, MWs, Tc_method="binary", Vc_method="binary").values()
    MW = sum(ys * MWs)
    Rg = R / MW

    print("Solving ctirical point")
    print(f"\t{mixture.Tc = :.2f} : {Tc = :.2f}")
    print(f"\t{mixture.Pc = :.2f} : {Pc = :.2f}")
    print(f"\t{mixture.Zc = :.4f} : {Zc = :.4f}")
    print(f"\t{mixture.Vc = :.6f} : {Vc = :.6f}")
    print(f"\t{mixture.rhoc = :.6f} : {rhoc = :.6f}")

    print("Solving EoS")
    Z_mix_PR, rho_mix_PR, omega_mix_PR = eos.peng_robinson(T, P, ys, Tcs, Pcs, MWs, ws)
    Z_mix_RK, rho_mix_RK, omega_mix_RK = eos.redlich_kwong(T, P, ys, Tcs, Pcs, MWs, ws, Vcs, eos="aungier")
    mixture.T, mixture.P = T, P
    print(f"\t{mixture.rhog = } : {rho_mix_PR = } : {rho_mix_RK = }")
    print(f"\t{mixture.Zg = } : {Z_mix_PR = } : {Z_mix_RK = }")
    print(f"\t{mixture.omega = } : {omega_mix_PR} : {omega_mix_RK = }")
    print(f"\t{mixture.R_specific = }")
    print(f"\t{mixture.mug = }")
    print(f"\t{mixture.Cpg / mixture.Cvg = }")
    print(f"\t{mixture.kg = }")

    Ts = np.linspace(200, 5000, 20)
    print("Solution for Cp")
    for Ti in Ts:
        cps = []
        print(f"{Ti = :.1f}")
        for fluid in my_mixture:
            coeff = fluid["nasa_coefficients"]
            mw = fluid["mw"]
            mf = fluid["mole_fraction"]
            name = fluid["name"]
            symbol = fluid["symbol"]
            cp_coeff = get_Cp_coefficients(Ti, coeff, name, symbol)
            Cp = nasa.molar_heat_capacity_at_constant_pressure(Ti, cp_coeff)
            cps.append((Cp, mf))
            e = thermo.Chemical(name)
            e.T, e.P = Ti, P
            print(f"\t{name} ({symbol}) : {mw = :.4f} : {e.Cpg = :.1f} : {Cp * R / mw = :.1f} : {(abs(e.Cpg - Cp * R / mw) / (Cp * R / mw)) * 100} %")
        cps = np.array(cps)
        Cp = sum(cps[:,0] * cps[:,1]) * R / MW
        mixture.T = Ti
        print(f"\tMixture\t{mixture.Cpg = :.1f} : {Cp = :.1f} : {(abs(Cp - mixture.Cpg) / Cp) * 100:.2f} %", end="\n\n")

    mixture.T = T
    ts1 = np.linspace(100, 1000, 100)
    ts2 = np.linspace(1000, 5000, 100)
    ts = np.hstack((ts1, ts2))

    cps1 = []
    cps2 = []
    cps3 = []
    cps4 = []

    for ti in ts1:
        mixture.T = ti
        cps1.append(mixture.Cpg)
        cps3.append(mixture.Cpg)

    for ti in ts2:
        mixture.T = ti
        cps2.append(mixture.Cpg)
        cps4.append(mixture.Cpg)

    coeffs1 = curve_fit(objective_function, ts1, cps1)[0]
    coeffs2 = curve_fit(objective_function, ts2, cps2)[0]
    coeffs3 = curve_fit(polynomial_objective, ts1, cps3)[0]
    coeffs4 = curve_fit(polynomial_objective, ts2, cps4)[0]

    print(f"Cp coefficients in the temperature range {ts1[0]} - {ts1[-1]}")
    print(coeffs1)

    print(f"Cp coefficients in the temperature range {ts2[0]} - {ts2[-1]}")
    print(coeffs2)

    print(f"Polynomil Cp coefficients in the temperature range {ts1[0]} - {ts1[-1]}")
    print(coeffs3)

    print(f"Polynomil Cp coefficients in the temperature range {ts2[0]} - {ts2[-1]}")
    print(coeffs4)

    # for ti in ts1:
    #     cp1 = objective_function(ti, *coeffs1)
    #     cp2 = polynomial_objective(ti, *coeffs3)
    #     print(f"{ti} : {100*(abs(cp1 - cp2)/cp1):.2f} %")

    # for ti in ts2:
    #     cp1 = objective_function(ti, *coeffs2)
    #     cp2 = polynomial_objective(ti, *coeffs4)
    #     print(f"{ti} : {100*(abs(cp1 - cp2)/cp1):.2f} %")

    ts = np.linspace(200, 1000, 100)
    csv_file = os.path.join(csv_path, "mixture_18377.csv")
    
    with open(csv_file, "w", newline="") as fp:
        fp.write("T,Cp,Cv,gamma\n")
        for ti in ts:
            cps = []
            for fluid in my_mixture:
                a = get_Cp_coefficients(ti, fluid["nasa_coefficients"], fluid["name"], fluid["symbol"])
                y = fluid["mole_fraction"]
                cp = (nasa.molar_heat_capacity_at_constant_pressure(ti, a))
                cps.append((cp, y))
            cps = np.array(cps)
            cp = sum(cps[:,0] * cps[:,1]) * R / MW
            cv = cp - Rg
            gamma= cp / cv
            fp.write(f"{ti},{cp},{cv},{gamma}\n")

    
    Z, rho, w = eos.redlich_kwong(331, 8038847, ys, Tcs, Pcs, MWs, ws, Vcs, eos="aungier")
    print(Z, rho, w)

