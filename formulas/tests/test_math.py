from formulas.base import Formula1D, Formula2D
from formulas.constants import FormulaConstant
from numpy.testing import assert_allclose
from sympy import Float
import yt.units as u

f_x = Formula1D("a*x**2+b","x",["a","b"])
g_x = Formula1D("p*cos(c*x)","x",["c","p"])
h_y = Formula1D("k*exp(d*y)","y",["k","d"])
j_o = Formula1D("2*v*o+z", "o", ["v","z"])
a = 50.*u.km/u.s**2
b = 30.*u.km
c = 10.*u.s**-1
d = 5.*u.s**-1
k = 10.*u.km
p = 20.*u.km
v = 3.0
z = 1.5
f_x.set_param_values(a=a, b=b)
g_x.set_param_values(p=p, c=c)
h_y.set_param_values(k=k, d=d)
j_o.set_param_values(v=v, z=z)
x = 60.*u.s
y = 15.*u.s
w = FormulaConstant("w", 2*u.km)
o = 3.4

def test_add():
    n_x = f_x+g_x
    assert n_x.ndim == 1
    assert isinstance(n_x, Formula1D)
    assert_allclose(n_x(x).v, (f_x(x)+g_x(x)).v)
    assert str(n_x(x).units) == "km"
    m_xy = n_x+h_y
    assert m_xy.ndim == 2
    assert isinstance(m_xy, Formula2D)
    assert_allclose(m_xy(x=x,y=y).v, (n_x(x)+h_y(y)).v)
    assert str(m_xy(x=x,y=y).units) == "km"
    assert (f_x + w)(x) == (w + f_x)(x)
    assert (j_o + 4)(o) == (4 + j_o)(o)
    assert (+f_x)(x) == +(f_x(x))

def test_subtract():
    n_x = f_x-g_x
    assert n_x.ndim == 1
    assert isinstance(n_x, Formula1D)
    assert_allclose(n_x(x).v, (f_x(x)-g_x(x)).v)
    assert str(n_x(x).units) == "km"
    m_xy = n_x-h_y
    assert m_xy.ndim == 2
    assert isinstance(m_xy, Formula2D)
    assert_allclose(m_xy(x=x,y=y).v, (n_x(x)-h_y(y)).v)
    assert str(m_xy(x=x,y=y).units) == "km"
    l_x = f_x-f_x
    assert len(l_x.params.values()) == 0
    assert len(l_x.param_values.values()) == 0
    assert len(l_x.var_symbols) == 0
    assert (-f_x)(x) == -(f_x(x))

def test_multiply():
    n_x = f_x*g_x
    assert n_x.ndim == 1
    assert isinstance(n_x, Formula1D)
    assert_allclose(n_x(x).v, (f_x(x)*g_x(x)).v)
    assert str(n_x(x).units) == "km**2"
    m_xy = n_x*h_y
    assert m_xy.ndim == 2
    assert isinstance(m_xy, Formula2D)
    assert_allclose(m_xy(x=x,y=y).v, (n_x(x)*h_y(y)).v)
    assert str(m_xy(x=x,y=y).units) == "km**3"
    assert (f_x * 2)(x) == (2 * f_x)(x)

def test_divide():
    n_x = f_x/g_x
    assert n_x.ndim == 1
    assert isinstance(n_x, Formula1D)
    assert_allclose(n_x(x).v, (f_x(x)/g_x(x)).v)
    assert str(n_x(x).units) == "dimensionless"
    m_xy = n_x/h_y
    assert m_xy.ndim == 2
    assert isinstance(m_xy, Formula2D)
    assert_allclose(m_xy(x=x,y=y).v, (n_x(x)/h_y(y)).v)
    assert str(m_xy(x=x,y=y).units) == "1/km"
    l_x = f_x/f_x
    assert len(l_x.params.values()) == 0
    assert len(l_x.param_values.values()) == 0
    assert len(l_x.var_symbols) == 0

def test_pow():
    n_x = f_x**3
    assert n_x.ndim == 1
    assert isinstance(n_x, Formula1D)
    assert_allclose(n_x(x).v, (f_x(x)**3).v)
    assert str(n_x(x).units) == "km**3"
    m_x = f_x**0.5
    assert n_x.ndim == 1
    assert isinstance(m_x, Formula1D)
    assert_allclose(m_x(x).v, (f_x(x)**0.5).v)
    assert str(m_x(x).units) == "sqrt(km)"
    l_x = f_x**Float(0.3)
    assert l_x.ndim == 1
    assert isinstance(l_x, Formula1D)
    assert_allclose(l_x(x).v, (f_x(x)**0.3).v)
    assert str(l_x(x).units) == "km**(3/10)"

def test_diff():
    n_x = f_x.diff("x")
    assert n_x.ndim == 1
    assert isinstance(n_x, Formula1D)
    assert_allclose(n_x(x).v, (2.*a*x).v)
    assert str(n_x(x).units) == "km/s"
    m_x = n_x.diff("x")
    assert m_x.ndim == 0
    assert_allclose(m_x().v, (2.*a).v)
    assert str(m_x().units) == "km/s**2"

