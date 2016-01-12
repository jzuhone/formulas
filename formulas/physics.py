from formulas.base import variable
from formulas.general import gaussian
from formulas.constants import pi

def maxwellian_velocity(v="v", A="A", sigma="sigma"):
    return gaussian(x=v, A=A, sigma=sigma, mu=0)

def maxwellian_speed(v="v", A="A", sigma="sigma"):
    f = maxwellian_velocity(v=v, A=A, sigma=sigma)
    v_var = variable(v)
    return 4*pi*v_var**2*f