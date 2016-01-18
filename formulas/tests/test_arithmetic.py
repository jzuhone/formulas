from formulas.base import Formula1D, Formula2D
from numpy.testing import assert_allclose
import yt.units as u

f_x = Formula1D("a*x**2+b","x",["a","b"])
g_x = Formula1D("p*cos(c*x)","x",["c","p"])
h_y = Formula1D("k*exp(d*y)","y",["k","d"])
f_x.set_param_values(a=50.*u.km/u.s**2, b=30.*u.km)
g_x.set_param_values(p=20.*u.km, c=10.*u.s**-1)
h_y.set_param_values(k=10.*u.km, d=5.*u.s**-1)
x = 60.*u.s
y = 15.*u.s

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

def test_pow():
    n_x = f_x**3
    assert n_x.ndim == 1
    assert isinstance(n_x, Formula1D)
    assert_allclose(n_x(x).v, (f_x(x)**3).v)
    assert str(n_x(x).units) == "km**3"
