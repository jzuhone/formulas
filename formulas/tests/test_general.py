from formulas.general import \
    linear, exponential, gaussian, \
    power_law
import numpy as np
from numpy.testing import assert_allclose
import yt.units as ytu
import astropy.units as apu
from yt.units.yt_array import YTArray

def test_linear():
    p = linear(x="y", a="d", b="e")
    d = 1.0
    e = 5.0
    y = 0.5
    p.set_param_values(d=d, e=e)
    assert_allclose(p(y), d*y+e)

def test_power_law():
    p = power_law(x="E", K="F", x_scale="E_0")
    F = 1.0e-10*ytu.erg/ytu.s/ytu.cm**2/ytu.keV
    E_0 = 1.0*ytu.keV
    alpha = -1.1
    E = YTArray(np.linspace(0.1, 10.0, 1000), "keV")
    p.set_param_values(F=F, E_0=E_0, alpha=alpha)
    assert_allclose(p(E).v, (F*(E/E_0)**alpha).v)
    assert str(p(E).units) == str(F.units)

def test_exponential():
    p = exponential(x="z",A="P_0",x_0="z_0",x_s="z_s")
    P_0 = 100.0*apu.kPa
    z_0 = 0.0*apu.m
    z_s = -100.0*apu.m
    z = apu.Quantity(np.linspace(0,1000.,100),"m")
    p.set_param_values(P_0=P_0, z_0=z_0, z_s=z_s)
    assert_allclose(p(z).value, (P_0*np.exp((z-z_0)/z_s)).value)
    assert str(p(z).unit) == str(P_0.unit)

def test_gaussian():
    p = gaussian()
    A = 1.0
    mu = 3.0
    sigma = 2.0
    x = np.linspace(-10.,10.,100)
    p.set_param_values(A=A, mu=mu, sigma=sigma)
    assert_allclose(p(x), A*np.exp(-0.5*((x-mu)/sigma)**2)/(np.sqrt(2.*np.pi)*sigma))