from formulas.physics import \
    maxwellian_speed, \
    maxwellian_velocity
from numpy.testing import assert_allclose
import numpy as np
import yt.units as u
from yt.units.yt_array import YTArray

v_th = 1000*u.km/u.s
n_0 = 1.0*u.cm**-3*u.km**-3*u.s**3
v = YTArray(np.linspace(-1000.,1000.,100), "km/s")
A = 1.0

def test_maxwellian_velocity():
    p = maxwellian_velocity(sigma="v_th")
    p.set_param_values(A=A, v_th=v_th)
    assert_allclose(p(v).v, (np.exp(-0.5*(v/v_th)**2)/(np.sqrt(2.*np.pi)*v_th)).v)
    assert str(p(v).units) == str((1/v_th).units)

def test_maxwellian_speed():
    p = maxwellian_speed(sigma="v_th", A="n_0")
    p.set_param_values(n_0=n_0, v_th=v_th)
    assert_allclose(p(v).v, (4.*np.pi*v*v*n_0*np.exp(-0.5*(v/v_th)**2)/(np.sqrt(2.*np.pi)*v_th)**3).v)
    assert str(p(v).units) == str((n_0*v_th**-1).units)
