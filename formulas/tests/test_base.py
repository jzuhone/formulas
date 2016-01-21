from formulas.base import Formula1D, variable, Formula2D
from formulas.radial_profiles import NFW_density_profile
from formulas.general import gaussian
import yt.units as u
from numpy.testing import assert_allclose

f_x = Formula1D("a*x**2+b","x",["a","b"])
a = 50.*u.km/u.s**2
b = 30.*u.km
f_x.set_param_values(a=a, b=b)

g_x = gaussian(x="x", A="A_x", mu="mu_x", sigma="sigma_x")
g_y = gaussian(x="y", A="A_y", mu="mu_y", sigma="sigma_y")
g_z = gaussian(x="z", A="A_z", mu="mu_z", sigma="sigma_z")
A_x = A_y = A_z = 1.0
mu_x = mu_y = mu_z = 30*u.s
sigma_x = sigma_y = sigma_z = 40*u.s
g_x.set_param_values(A_x=A_x, mu_x=mu_x, sigma_x=sigma_x)
g_y.set_param_values(A_y=A_y, mu_y=mu_y, sigma_y=sigma_y)
g_z.set_param_values(A_z=A_z, mu_z=mu_z, sigma_z=sigma_z)

g_xy = g_x*g_y
g_xyz = g_xy*g_z

x = 60.*u.s
y = 25.*u.s
z = -10.*u.s

def test_unitless():
    assert_allclose(f_x(x).v, f_x.unitless()(x.v))
    assert_allclose(g_xy(x=x,y=y).v, g_xy.unitless()(x=x.v,y=y.v))
    assert_allclose(g_xyz(x=x,y=y,z=z).v, g_xyz.unitless()(x=x.v,y=y.v,z=z.v))

def test_clear():
    f_x.clear_param_values()
    for pv in f_x.param_values.values():
        assert pv is None
    f_x.set_param_values(a=a, b=b)

def test_variable():
    x = variable("x")
    y = variable("y")
    z1 = x*y
    z2 = Formula2D("x*y", "x", "y", [])
    assert z1.formula == z2.formula

def test_latex():
    pd = NFW_density_profile()
    lr = '\\frac{r_{s} \\rho_{s}}{r \\left(\\frac{r}{r_{s}} + 1\\right)^{2}}'
    assert pd.latex_representation() == lr