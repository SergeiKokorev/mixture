import os
import sys


PYTHON_PATH = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))

sys.path.append(PYTHON_PATH)


from datamodels.properties import *


if __name__ == "__main__":

    ref = Reference(326.15, 7344000)
    print(ref, end="\n\n")
    transport = Transport(1e-6, 0.125)
    print(transport)
    print(transport.nu(1.25))
