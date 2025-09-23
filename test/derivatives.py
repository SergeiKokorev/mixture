import sympy as sp


if __name__ == "__main__":

    v, R, t, a, b, n = sp.symbols('v, R, t, a, b, n')
    p = R * t / (v - b) - a / (v**2 + 2 * b * v - b**2)
    
    partial_p_v = sp.diff(p, v)
    print(partial_p_v)
    partial_2p_v2 = sp.diff(partial_p_v, v)
    print(partial_2p_v2)

    p = n * R * t / (v - n * b) - n**2 * a / v**2
    partial_p_v = sp.diff(p, v)
    partial_2p_v2 = sp.diff(partial_p_v, v)

    print(partial_p_v)
    print(partial_2p_v2)
