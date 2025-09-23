import os
import sys
import thermo
import numpy as np
import json


PYTHON_PATH = os.path.abspath(os.path.join(
    os.path.dirname(__file__), os.pardir
))
sys.path.append(PYTHON_PATH)

import fluid.mix_eos as eos
import fluid.mixing as mixing
import fluid.nasa as nasa


from fluid.const import R


def get_coeffs(Tref, coeffs) -> list:
    for coef in coeffs:
        if Tref < coef[0][0]:
            raise ValueError(f"Coefficients are not defined for fluid at the temperature {Tref}")
        elif coef[0][0] <= Tref <= coef[0][1]:
            return coef[1]
    
    return None


if __name__ == "__main__":

    nasa_glenn_data_file = os.path.abspath(os.path.join(PYTHON_PATH, "lib", "thermo.json"))
    
    with open(nasa_glenn_data_file, "r") as fp:
        database = json.load(fp)

    fluid_list = [
        {"thermo": "H2O", "name": "Water", "symbol": "H2O", "mw": 18.02e-3, "sigma": 2.641, "epsilon_kb": 809.1, "Tc": 647.1, "pc": 22.064e6, "rhoc": 322, "Tb": 373.15,"omega": 0.344, "mole_fraction": 4.65e-2, "Ttr": 273.2, "Ptr": 612, "antoine_coef": [8.07131,1730.63,233.426,(273.15,373.15)]},
        {"thermo": "H2S", "name": "Hydrogen Sulfide", "symbol": "H2S", "mw": 34.08e-3, "sigma": 3.623, "epsilon_kb": 301.1, "Tc": 373.5, "pc": 89.63e5, "rhoc": 347.3, "Tb": 212.85, "omega":0.1, "mole_fraction": 2.27e-2, "Ttr": 187.7, "Ptr": 23259, "antoine_coef": [7.05496,1315.92,237.731, (197,349)]},
        {"thermo": "H2", "name": "Hydrogen", "symbol": "H2", "mw": 2.02e-3, "sigma": 2.827, "epsilon_kb": 59.7, "Tc": 33.2, "pc": 13e5, "omega": -0.218, "rhoc": 31.26, "Tb": 20.25, "mole_fraction": 5.86e-2,"Ttr": 14, "Ptr": 7358, "antoine_coef": [5.64473,123.758,12.684,(18,26)]},
        {"thermo": "N2", "name": "Nitrogen", "symbol": "N2", "mw": 28.1e-3, "sigma": 3.798, "epsilon_kb": 71.4, "Tc": 126.2, "pc": 33.9e5, "rhoc": 313.3,"omega": 0.037, "mole_fraction": 0.0, "Tb": 77.35,"Ttr": 63.2, "Ptr": 12520,"antoine_coef": [6.49457,255.68,266.55,(63,93)]},
        {"thermo": "NH3", "name": "Ammonia", "symbol": "NH3", "mw": 17.03e-3, "sigma": 2.9, "epsilon_kb": 558.3, "Tc": 405.65, "pc": 113.33e5, "rhoc": 233.3, "omega": 0.256, "mole_fraction": 0.44e-2, "Tb": 239.85,"Ttr": 195.5, "Ptr": 6091,"antoine_coef": [7.55466,1002.711,247.885,(190,333)]},
        {"thermo": "CO2", "name": "Carbon Dioxide", "symbol": "CO2", "mw": 44.01e-3, "sigma": 3.941, "epsilon_kb": 195.2, "Tc": 304.13, "pc": 73.75e5, "rhoc": 467.6, "omega":0.225, "mole_fraction": 0.32e-2, "Tb": 194.65,"Ttr": 216.6, "Ptr": 5.1796e5,"antoine_coef": [9.81064,1347.788,273.00,(183,217)]},
        {"thermo": "CH4","name": "Methane", "symbol": "CH4", "mw": 16.04e-3, "sigma": 3.758, "epsilon_kb": 148.6, "Tc": 190.6, "pc": 46e5, "rhoc": 162.7, "omega": 0.011, "mole_fraction": 30.3e-2, "Tb": 111.65,"Ttr": 90.7, "Ptr": 11696,"antoine_coef": [6.61184,389.93,266.00,(92,121)]},
        {"thermo": "C2H4", "name": "Ethylene", "symbol": "C2H4", "mw": 28.05e-3, "sigma": 4.163, "epsilon_kb": 224.7, "Tc": 282.35, "pc": 50.41e5, "rhoc": 214.2, "omega": 0.086, "mole_fraction": 2.57e-2, "Tb": 169.45,"Ttr": 104, "Ptr": 122, "antoine_coef": [6.74756,585.00,255.00,(165,198)]},
        {"thermo": "C2H6","name": "Ethane", "symbol": "C2H6", "mw": 30.07e-3, "sigma": 4.443, "epsilon_kb": 215.7, "Tc": 305.32, "pc": 48.72e5, "rhoc": 206.2,"omega": 0.099, "mole_fraction": 16.87e-2, "Tb": 184.55,"Ttr": 90.4, "Ptr": 1,"antoine_coef": [6.80266,656.40,256.00, (170,248)]},
        {"thermo": "C3H6,propylene", "name": "Propene", "symbol": "C3H6", "mw": 42.08e-3, "sigma": 4.678, "epsilon_kb": 298.9, "Tc": 365.0, "pc": 46e5, "rhoc": 229.6,"omega": 0.142, "mole_fraction": 3.47e-2, "Tb": 225.55,"Ttr": 273, "Ptr": 3.427e5,"antoine_coef": [6.81960,785.00,247.00, (185,225)]},
        {"thermo": "C3H8", "name": "Propane", "symbol": "C3H8", "mw": 44.1e-3, "sigma": 5.118, "epsilon_kb": 237.1, "Tc": 369.8, "pc": 42.48e5, "rhoc": 220.5,"omega":0.152, "mole_fraction": 8.28e-2, "Tb": 231.05,"Ttr": 85.5, "Ptr": 0,"antoine_coef": [6.80338,803.810,246.990,(165,248)]},
        {"thermo": "C4H10,isobutane", "name": "Isobutane", "symbol": "i-C4H10", "mw": 58.12e-3, "sigma": 5.278, "epsilon_kb": 330.1, "Tc": 408.1, "pc": 36.48e5, "rhoc": 225.5,"omega": 0.186, "mole_fraction": 1.39e-2, "Tb": 261.45,"Ttr": 113.7, "Ptr": 0,"antoine_coef": [6.52587,787.160,248.870,(166,261)]},
        {"thermo": "C4H8,isobutene", "name": "Isobutene", "symbol": "i-C4H8", "mw": 56.11e-3, "sigma": 4.863, "epsilon_kb": 331.4, "Tc": 418.9, "pc": 40.0e5, "rhoc": 234, "omega": 0.194, "mole_fraction": 1.29e-2, "Tb": 266.25,"Ttr":132.4, "Ptr":1,"antoine_coef": [6.53990,811.658,236.260,(182,266)]},
        {"thermo": "C4H8,1-butene", "name": "1-Butene", "symbol": "1-C4H8", "mw": 56.11e-3, "sigma": 4.678, "epsilon_kb": 331.4, "Tc": 419.5, "pc": 40.2e5, "rhoc": 234, "omega": 0.191, "mole_fraction": 1.91e-2, "Tb": 266.65,"Ttr":87.8, "Ptr":0,"antoine_coef": [6.84290,926.10,240.00,(193,267)]},
        {"thermo": "C4H8,tr2-butene", "name": "trans-2-Butene", "symbol": "t-2-C4H8", "mw": 56.11e-3, "sigma": 4.678, "epsilon_kb": 331.4, "Tc": 428.6, "pc": 41.4e5, "rhoc": 234, "omega": 0.210, "mole_fraction": 0.7e-2, "Tb": 274.05,"Ttr":167.6, "Ptr": 75,"antoine_coef": [6.84570,1002.771,237.553,(241,274)]},
        {"thermo": "C4H8,cis2-buten", "name": "cis-2-Butene", "symbol": "c-2-C4H8", "mw": 56.11e-3, "sigma": 4.678, "epsilon_kb": 331.4, "Tc": 435.5, "pc": 42.4e5, "rhoc": 234, "omega": 0.205, "mole_fraction": 0.45e-2, "Tb": 276.85,"Ttr":134.3, "Ptr":0,"antoine_coef": [6.86970,1016.980,237.483,(260,277)]},
        {"thermo": "C4H6,butadiene", "name": "1,3-Butadiene", "symbol": "C4H6", "mw": 54.09e-3, "sigma": 4.840, "epsilon_kb": 293.0, "Tc": 425.0, "pc": 43.2e5, "rhoc": 237.9, "omega": 0.195, "mole_fraction": 0.07e-2, "Tb": 268.75,"Ttr": 273.2, "Ptr": 612,"antoine_coef": [6.85941,935.531,239.554,(260,290)]},
        {"thermo": "C4H10,n-butane", "name": "Butane", "symbol": "n-C4H10", "mw": 58.12e-3, "sigma": .678, "epsilon_kb": 531.4, "Tc": 425.12, "pc": 37.96e5, "rhoc": 228, "omega": 0.199, "mole_fraction": 4.83e-2, "Tb": 272.65,"Ttr":134.9, "Ptr":1,"antoine_coef": [6.80776,935.860,238.730,(195,272)]},
        {"thermo": "C5H12,i-pentane", "name": "Isopentane", "symbol": "i-C5H12", "mw": 72.15e-3, "sigma": 5.392, "epsilon_kb": 341.1, "Tc": 460.4, "pc": 33.84e5, "rhoc": 236, "omega":0.227, "mole_fraction": 1.04e-2, "Tb": 300.95,"Ttr":112.6, "Ptr":0,"antoine_coef": [6.78967,1020.012,233.097,(216,301)]},
        {"thermo": "C5H10,1-pentene", "name": "1-Pentene", "symbol": "1-C5H10", "mw": 70.14e-3, "sigma": 5.180, "epsilon_kb": 341.0, "Tc": 464.8, "pc": 35.6e5, "rhoc": 228, "omega": 0.235, "mole_fraction": 2.7e-2, "Tb": 303.15,"Ttr": 273.2, "Ptr": 612,"antoine_coef": [6.85296,1130.649,231.719,(223,303)]},
        {"thermo": "C5H12,n-pentane", "name": "Pentene", "symbol": "n-C5H12", "mw": 72.15e-3, "sigma": 5.397, "epsilon_kb": 341.1, "Tc": 469.7, "pc": 33.7e5, "rhoc": 236, "omega": 0.251, "mole_fraction": 2.8e-2, "Tb": 309.25,"Ttr":143.5, "Ptr":0,"antoine_coef": [6.87632,1075.780,233.205,(223,331)]},
        {"thermo": "C5H10,cyclo-", "name": "Cyclopentane", "symbol": "C5H10", "mw": 70.14e-3, "sigma": 5.270, "epsilon_kb": 345.0, "Tc": 511.7, "pc": 45.02e5, "rhoc": 228, "omega": 0.196, "mole_fraction": 0.46e-2, "Tb": 322.45,"Ttr":179.7, "Ptr":9,"antoine_coef": [6.88685,1124.162,231.361,(223,322)]},
        {"thermo": "C6H6", "name": "Benzene", "symbol": "C6H6", "mw": 78.12e-3, "sigma": 5.349, "epsilon_kb": 412.3, "Tc": 562.05, "pc": 48.95e5, "rhoc": 304.8, "omega": 0.212, "mole_fraction": 7.31e-2, "Tb": 353.25,"Ttr":278.7, "Ptr":4784,"antoine_coef": [6.90565,1211.033,220.790,8,103,(281,376)]},
    ]

    for fluid in fluid_list:
        component_symbol = fluid["name"]
        thermo_name = fluid["thermo"]
        thermo_component = thermo.Chemical(component_symbol)
        thermo_component.T = 1000
        my_component = database[thermo_name]

        thermo_Cpg = thermo_component.Cpg
        coeffs = get_coeffs(thermo_component.T, my_component["coefficients"])
        my_Cp = nasa.molar_heat_capacity_at_constant_pressure(thermo_component.T, coeffs)
        my_Cp *= R / (my_component["molecular_weight"] * 1e-3)

        print(f"{thermo_component.T = } K : {thermo_component.formula} : {thermo_component.MW} : {thermo_component.name} : {thermo_component.formula} : {thermo_Cpg = :.2f} : {my_Cp = :.2f} : Discrepancy {100 * ((thermo_Cpg - my_Cp) / thermo_Cpg):.2f} %")

    T, P = 326.15, 172000

    # Test pure component
    fluids = [fluid for fluid in fluid_list]#fluid_list[6], fluid_list[10]
    y1 = 0.753
    mole_fractions = [fluid["mole_fraction"] for fluid in fluid_list]#[y1, 1 - y1]
    thermo_fluids = []
    my_fluids = []

    print(f"Reference {T = } [K]; {P = } [Pa]")
    for i, fluid in enumerate(fluids):
        name = fluid.get("name")
        symbol = fluid.get("symbol")
        mw = fluid.get("mw")
        Tc = fluid.get("Tc")
        Pc = fluid.get("pc")
        rhoc = fluid.get("rhoc")
        omega = fluid.get("omega")
        thermo_name = fluid.get("thermo")
        nasa_coeffs = database[thermo_name]["coefficients"]
        mole_fraction = mole_fractions[i]
        thermo_fluids.append((name, mole_fraction))

        my_fluids.append(
            {
                "name": name, "symbol": symbol,
                "Tc": Tc, "Pc": Pc, "rhoc": rhoc, "omega": omega,
                "mole_fraction": mole_fraction, "mw": mw,
                "nasa_coefficients": nasa_coeffs
            }
        )

        cp_coeff = get_coeffs(T, nasa_coeffs)
        cp_nasa = nasa.molar_heat_capacity_at_constant_pressure(T, cp_coeff) * R / mw
        e = thermo.Chemical(name)
        e.T = T
        cp_thermo = e.Cpg
    
    mixture = thermo.Mixture(IDs=[fluid[0] for fluid in thermo_fluids], zs=[fluid[1] for fluid in thermo_fluids])
    mixture.T, mixture.P = T, P
    
    ys = np.array([fluid["mole_fraction"] for fluid in my_fluids])
    Tcs = np.array([fluid["Tc"] for fluid in my_fluids])
    Pcs = np.array([fluid["Pc"] for fluid in my_fluids])
    rhocs = np.array([fluid["rhoc"] for fluid in my_fluids])
    Mws = np.array([fluid["mw"] for fluid in my_fluids])
    ws = np.array([fluid["omega"] for fluid in my_fluids])
    Vcs = np.array([fluid["mw"] / fluid["rhoc"] for fluid in my_fluids])

    Tc, Pc, Vc, Zc = mixing.prausnitz_gunn(ys, Tcs, Vcs, ws, Mws).values()

    print("Solving ctirical point")
    print(f"\t{mixture.Tc = :.2f} : {Tc = :.2f}")
    print(f"\t{mixture.Zc = :.4f} : {Zc = :.4f}")
    print(f"\t{mixture.Vc = :.6f} : {Vc = :.6f}")
    print(f"\t{mixture.Pc = :.0f} : {Pc = :.0f}")

    print("Solving EoS")
    Z_mix, rho_mix, omega_mix = eos.peng_robinson(T, P, ys, Tcs, Pcs, Mws, ws)
    print(f"\t{mixture.rhog = :.4f} : {rho_mix}")
    print(f"\t{mixture.Zg = :.4f} : {Z_mix = :.4f}")
    print(f"\t{mixture.omega = :.4f} : {omega_mix = :.4f}")

    Ts = np.linspace(200, 1000, 20)
    Mw = sum(ys * Mws)
    print("Solution for Cp")
    for Ti in Ts:
        cps = []
        for fluid in my_fluids:
            coeff = fluid["nasa_coefficients"]
            cp_coeff = get_coeffs(Ti, coeff)
            mw = fluid["mw"]
            mf = fluid["mole_fraction"]
            name = fluid["name"]
            cps.append((nasa.molar_heat_capacity_at_constant_pressure(Ti, cp_coeff), mf))
        cps = np.array(cps)
        Cp = sum(cps[:,0] * cps[:,1]) * R / Mw
        mixture.T = Ti
        print(f"\t{Ti = :.2f} : thermo Cp {mixture.Cpg:.1f} : my Cp {Cp:.1f} : Deviation {100 * (abs(mixture.Cpg - Cp) / Cp)} %")
