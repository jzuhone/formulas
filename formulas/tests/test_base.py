from formulas.base import Formula1D
import yt.units as u
from numpy.testing import assert_allclose

f_x = Formula1D("a*x**2+b","x",["a","b"])
a = 50.*u.km/u.s**2
b = 30.*u.km
f_x.set_param_values(a=a, b=b)
x = 60.*u.s

def test_unitless():
    assert_allclose(f_x(x).v, f_x.unitless()(x.v))

def test_clear():
    f_x.clear_param_values()
    for pv in f_x.param_values.values():
        assert pv is None
    f_x.set_param_values(a=a, b=b)