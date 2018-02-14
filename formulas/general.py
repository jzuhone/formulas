from sympy import symbols, exp, sqrt, Rational, pi
from formulas.base import Formula1D

def linear(x="x", a="a", b="b"):
    """
    A linear formula.

    Parameters
    ----------
    x : string
        The symbol for the independent variable.
    a : string
        The symbol for the slope of the line.
    b : string
        The symbol for the intercept.
    """
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
    """
    An exponential formula.

    Parameters
    ----------
    x : string
        The symbol for the independent variable.
    A : string
        The symbol for the normalization.
    x_0 : string
        The symbol for the y-intercept.
    x_s : string
        The symbol for the scale parameter.
    """
    x, A, x_0, x_s = symbols((x, A, x_0, x_s))
    formula = A*exp((x-x_0)/x_s)
    return Formula1D(formula, x, [A, x_0, x_s])

def gaussian(x="x", A="A", mu="mu", sigma="sigma"):
    """
    A Gaussian formula.

    Parameters
    ----------
    x : string
        The symbol for the independent variable.
    A : string
        The symbol for the normalization.
    mu : string
        The symbol for the mean.
    sigma : string
        The symbol for the standard deviation.
    """
    x, A, mu, sigma = symbols((x, A, mu, sigma))
    params = [A, mu, sigma]
    formula = A*exp(-Rational(1,2)*((x-mu)/sigma)**2)/(sigma*sqrt(2*pi))
    return Formula1D(formula, x, params)
