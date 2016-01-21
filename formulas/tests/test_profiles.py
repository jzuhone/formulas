from formulas.radial_profiles import \
    beta_model_profile, \
    vikhlinin_density_profile, \
    vikhlinin_temperature_profile, \
    baseline_entropy_profile, \
    AM06_density_profile, \
    AM06_temperature_profile, \
    NFW_density_profile, \
    hernquist_density_profile, \
    exponential_taper_profile, \
    rescale_profile_by_mass
from yt.units.yt_array import YTArray
import yt.units as u
import numpy as np
from numpy.testing import assert_allclose

r = YTArray(np.linspace(10., 1000., 10000), "kpc")

def test_nfw():
    p = NFW_density_profile()
    r_s = 100*u.kpc
    rho_s = 1.0e8*u.Msun/u.kpc**3
    p.set_param_values(r_s=r_s, rho_s=rho_s)
    x = r/r_s
    assert_allclose(p(r).v, (rho_s/(x*(1.+x)**2)).v)
    assert str(p(r).units) == str(rho_s.units)