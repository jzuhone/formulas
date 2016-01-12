from sympy import symbols, exp
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
    formula = A*exp(-(x-x_0)/x_s)
    return Formula1D(formula, x, [A, x_0, x_s])

def gaussian(x="x", A="A", mu="mu", sigma="sigma"):
    x, A, mu, sigma = symbols((x, A, mu, sigma))
    formula = A*exp(-0.5*((x-mu)/sigma)**2)
    return Formula1D(formula, x, [A, mu, sigma])
