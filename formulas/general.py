from sympy import symbols, exp, sqrt, Rational, pi
from formulas.base import Formula1D

def linear(x="x", a="a", b="b"):
    x, a, b = symbols((x, a, b))
    formula = a*x+b
    return Formula1D(formula, x, [a, b])

def power_law(x="x", K="K", x_scale="x_scale", alpha="alpha"):
    """
    A power-law formula.

    Parameters
    ----------
    x : string
        The symbol for the independent variable.
    K : string
        The symbol for the normalization.
    x_scale : string
        The symbol for the scale value.
    alpha : string
        The symbol for the power-law index.
    """
    x, K, x_scale, alpha = symbols((x, K, x_scale, alpha))
    formula = K*(x/x_scale)**alpha
    return Formula1D(formula, x, [x_scale, K, alpha])

def exponential(x="x", A="A", x_0="x_0", x_s="x_s"):
    x, A, x_0, x_s = symbols((x, A, x_0, x_s))
    formula = A*exp((x-x_0)/x_s)
    return Formula1D(formula, x, [A, x_0, x_s])

def gaussian(x="x", A="A", mu="mu", sigma="sigma"):
    x, A, sigma = symbols((x, A, sigma))
    params = [A, sigma]
    if mu != 0:
        mu = symbols(mu)
        params.append(mu)
    formula = A*exp(-Rational(1,2)*((x-mu)/sigma)**2)/(sigma*sqrt(2*pi))
    return Formula1D(formula, x, params)
