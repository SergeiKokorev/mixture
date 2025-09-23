if __name__ == "__main__":

    fluid_list = [
            {"law": "ig", "thermo": "H2O", "name": "Water", "symbol": "H2O", "mw": 18.02e-3, "sigma": 2.641, "epsilon_kb": 809.1, "Tc": 647.1, "pc": 220550, "omega": 0.344, "mole_fraction": 4.65e-2},
            {"thermo": "H2S", "name": "Hydrogen Sulfide", "symbol": "H2S", "mw": 34.08e-3, "sigma": 3.623, "epsilon_kb": 301.1, "Tc": 373.5, "pc": 89.63e5, "omega":0.1, "mole_fraction": 2.27e-2},
            {"thermo": "H2", "name": "Hydrogen", "symbol": "H2", "mw": 2.02e-3, "sigma": 2.827, "epsilon_kb": 59.7, "Tc": 33.2, "pc": 13e5, "omega": -0.218+.3, "mole_fraction": 5.86e-2},
            {"thermo": "N2", "name": "Nitrogen", "symbol": "N2", "mw": 28.1e-3, "sigma": 3.798, "epsilon_kb": 71.4, "Tc": 126.2, "pc": 33.9e5, "omega": 0.037, "mole_fraction": 0.0},
            {"thermo": "NH3", "name": "Ammonia", "symbol": "NH3", "mw": 17.03e-3, "sigma": 2.9, "epsilon_kb": 558.3, "Tc": 405.65, "pc": 113.33e5, "omega": 0.256, "mole_fraction": 0.44e-2},
            {"thermo": "CO2", "name": "Carbon Dioxide", "symbol": "CO2", "mw": 44.01e-3, "sigma": 3.941, "epsilon_kb": 195.2, "Tc": 304.13, "pc": 73.75e5, "omega":0.225, "mole_fraction": 0.32e-2},
            {"thermo": "CH4","name": "Methane", "symbol": "CH4", "mw": 16.04e-3, "sigma": 3.758, "epsilon_kb": 148.6, "Tc": 190.6, "pc": 46e5, "omega": 0.11, "mole_fraction": 30.3e-2},
            {"thermo": "C2H4", "name": "Ethylene", "symbol": "C2H4", "mw": 28.05e-3, "sigma": 4.163, "epsilon_kb": 224.7, "Tc": 282.35, "pc": 50.41e5, "omega": 0.086, "mole_fraction": 2.57e-2},
            {"thermo": "C2H6","name": "Ethane", "symbol": "C2H6", "mw": 30.07e-3, "sigma": 4.443, "epsilon_kb": 215.7, "Tc": 305.32, "pc": 48.72e5, "omega": 0.099, "mole_fraction": 16.87e-2},
            {"thermo": "C3H6,propylene", "name": "Propene", "symbol": "C3H6", "mw": 42.08e-3, "sigma": 4.678, "epsilon_kb": 298.9, "Tc": 365.0, "pc": 46e5, "omega": 0.142, "mole_fraction": 3.47e-2},
            {"thermo": "C3H8", "name": "Propane", "symbol": "C3H8", "mw": 44.1e-3, "sigma": 5.118, "epsilon_kb": 237.1, "Tc": 369.8, "pc": 42.48e5, "omega":0.152, "mole_fraction": 8.28e-2},
            {"thermo": "C4H10,isobutane", "name": "Isobutane", "symbol": "i-C4H10", "mw": 58.12e-3, "sigma": 5.278, "epsilon_kb": 330.1, "Tc": 408.1, "pc": 36.48e5, "omega": 0.186, "mole_fraction": 1.39e-2},
            {"thermo": "C4H8,isobutene", "name": "Isobutene", "symbol": "i-C4H8", "mw": 56.11e-3, "sigma": 4.863, "epsilon_kb": 331.4, "Tc": 418.9, "pc": 40.0e5, "omega": 0.194, "mole_fraction": 1.29e-2},
            {"thermo": "C4H8,1-butene", "name": "1-Butene", "symbol": "1-C4H8", "mw": 56.11e-3, "sigma": 4.678, "epsilon_kb": 331.4, "Tc": 419.5, "pc": 40.2e5, "omega": 0.191, "mole_fraction": 1.91e-2},
            {"thermo": "C4H8,tr2-butene", "name": "trans-2-Butene", "symbol": "t-2-C4H8", "mw": 56.11e-3, "sigma": 4.678, "epsilon_kb": 331.4, "Tc": 428.6, "pc": 41.4e5, "omega": 0.210, "mole_fraction": 0.7e-2},
            {"thermo": "C4H8,cis2-buten", "name": "cis-2-Butene", "symbol": "c-2-C4H8", "mw": 56.11e-3, "sigma": 4.678, "epsilon_kb": 331.4, "Tc": 435.5, "pc": 42.4e5, "omega": 0.205, "mole_fraction": 0.45e-2},
            {"thermo": "C4H6,butadiene", "name": "1,3-Butadiene", "symbol": "C4H6", "mw": 54.09e-3, "sigma": 4.840, "epsilon_kb": 293.0, "Tc": 425.0, "pc": 43.2e5, "omega": 0.195, "mole_fraction": 0.07e-2},
            {"thermo": "C4H10,n-butane", "name": "Butane", "symbol": "n-C4H10", "mw": 58.12e-3, "sigma": .678, "epsilon_kb": 531.4, "Tc": 425.12, "pc": 37.96e5, "omega": 0.199, "mole_fraction": 4.83e-2},
            {"thermo": "C5H12,i-pentane", "name": "Isopentane", "symbol": "i-C5H12", "mw": 72.15e-3, "sigma": 5.392, "epsilon_kb": 341.1, "Tc": 460.4, "pc": 33.84e5, "omega":0.227, "mole_fraction": 1.04e-2},
            {"thermo": "C5H10,1-pentene", "name": "1-Pentene", "symbol": "1-C5H10", "mw": 70.14e-3, "sigma": 5.180, "epsilon_kb": 341.0, "Tc": 464.8, "pc": 35.6e5, "omega": 0.235, "mole_fraction": 2.7e-2},
            {"thermo": "C5H12,n-pentane", "name": "Pentene", "symbol": "n-C5H12", "mw": 72.15e-3, "sigma": 5.397, "epsilon_kb": 341.1, "Tc": 469.7, "pc": 33.7e5, "omega": 0.251, "mole_fraction": 2.8e-2},
            {"thermo": "C5H10,cyclo-", "name": "Cyclopentane", "symbol": "C5H10", "mw": 70.14e-3, "sigma": 5.270, "epsilon_kb": 345.0, "Tc": 511.7, "pc": 45.02, "omega": 0.196, "mole_fraction": 0.46e-2},
            {"thermo": "C6H6", "name": "Benzene", "symbol": "C6H6", "mw": 78.12e-3, "sigma": 5.349, "epsilon_kb": 412.3, "Tc": 562.05, "pc": 48.95e5, "omega": 0.212, "mole_fraction": 7.31e-2},
        ]
    
