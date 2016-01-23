from formulas.general import gaussian
from formulas.base import Formula1D
from sympy import pi, exp, Rational, sqrt, symbols

def maxwellian_velocity(v="v", A="A", sigma="sigma"):
    return gaussian(x=v, A=A, sigma=sigma, mu=0)

def maxwellian_speed(v="v", A="A", sigma="sigma"):
    v, A, sigma = symbols((v, A, sigma))
    params = [A, sigma]
    f = A*exp(-Rational(1,2)*((v/sigma)**2))/(sigma*sqrt(2*pi))**3
    return Formula1D(4*pi*v*v*f, v, params)
