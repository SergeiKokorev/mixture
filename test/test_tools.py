import os
import sys

PYTHON_PATH = os.path.abspath(os.path.join(
    os.path.dirname(__file__), os.pardir
))
sys.path.append(PYTHON_PATH)


from fluid.eos import __check_args__

if __name__ == "__main__":

    args = (1, 2, 3)

    args = __check_args__(args)
    print(args)

    args = 1
    args = __check_args__(args)
    print(args)

