import os
import sys


PYTHON_PATH = os.path.abspath(os.path.join(
    os.path.dirname(__file__), os.pardir
))
sys.path.append(PYTHON_PATH)


from fluid.const import NA


def chapman_enskog(T: float, MW: float, sigma: float, epsilon_kb: float):
    """
    Computes dynamic viscosity using Chapman-Enskog model

    Args:
        T (float): Temperature in K at which dynamic viscosity be calculated
        MW (float): Molercular weight in kg/mol of the component
        sigma (float): Characteristic diameter. A measure of the size of the fas molecules
        epsilon_kb
    """
