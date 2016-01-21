from formulas.radial_profiles import \
    beta_model_profile, \
    vikhlinin_density_profile, \
    vikhlinin_temperature_profile, \
    baseline_entropy_profile, \
    AM06_density_profile, \
    AM06_temperature_profile, \
    NFW_density_profile, \
    NFW_mass_profile, \
    hernquist_density_profile, \
    hernquist_mass_profile, \
    exponential_taper_profile, \
    rescale_profile_by_mass
from yt.units.yt_array import YTArray
from astropy.units import Quantity
import yt.units as ytu
import astropy.units as apu
import numpy as np
from numpy.testing import assert_allclose

r_yt = YTArray(np.linspace(10., 1000., 10000), "kpc")
r_ap = Quantity(np.linspace(0.03, 2.0, 10000), "Mpc")

def test_nfw_density():
    p = NFW_density_profile()
    r_s = 100*ytu.kpc
    rho_s = 1.0e8*ytu.Msun/ytu.kpc**3
    p.set_param_values(r_s=r_s, rho_s=rho_s)
    x = r_yt/r_s
    assert_allclose(p(r_yt).v, (rho_s/(x*(1.+x)**2)).v)
    assert str(p(r_yt).units) == str(rho_s.units)

def test_hernquist_density():
    p = hernquist_density_profile()
    a = 200.*ytu.kpc
    M_0 = 1.0e14*ytu.Msun
    p.set_param_values(a=a, M_0=M_0)
    x = r_yt/a
    assert_allclose(p(r_yt).v, (M_0/(2*np.pi*a**3)/(x*(1.+x)**3)).v)
    assert str(p(r_yt).units) == str((M_0/a**3).units)

def test_hernquist_mass():
    p = hernquist_mass_profile()
    a = 350.*apu.kpc
    M_0 = 1.0e15*apu.Msun
    p.set_param_values(a=a, M_0=M_0)
    assert_allclose(p(r_ap).to("solMass").value, (M_0*r_ap**2/(r_ap+a)**2).value)

def test_beta_model():
    p = beta_model_profile()
    r_c = 100.*apu.kpc
    rho_c = 1.0e-25*apu.g/apu.cm**3
    beta = 1.0
    p.set_param_values(r_c=r_c, rho_c=rho_c, beta=beta)
    assert_allclose(p(r_ap).value, (rho_c*(1.+(r_ap/r_c)**2)**(-1.5*beta)).value)

def test_baseline_entropy():
    p = baseline_entropy_profile()
    K_0 = 10.*ytu.keV*ytu.cm**2
    K_200 = 10000.*ytu.keV*ytu.cm**2
    r_200 = 1.0*ytu.Mpc
    alpha = 1.1
    p.set_param_values(K_0=K_0, K_200=K_200, r_200=r_200, alpha=alpha)
    assert_allclose(p(r_yt).value, (K_0+K_200*(r_yt/r_200)**alpha).value)
    assert str(p(r_yt).units) == str(K_0.units)

def test_AM06_temperature():
    p = AM06_temperature_profile()
    a_c = 60*ytu.kpc
    a = 600.*ytu.kpc
    c = 0.17
    T_0 = 10.*ytu.keV
    p.set_param_values(T_0=T_0, a_c=a_c, c=c, a=a)
    assert_allclose(p(r_yt).value, (T_0/(1.+r_yt/a)*(c+r_yt/a_c)/(1.+r_yt/a_c)).value)
    assert str(p(r_yt).units) == str(T_0.units)

def test_mass_rescaling_yt():
    M = 6.0e14*ytu.Msun
    R = 1500.0*ytu.kpc
    pm = NFW_mass_profile()
    pd = NFW_density_profile()
    pd.set_param_values(r_s=350*ytu.kpc, rho_s=1.0*ytu.Msun/ytu.kpc**3)
    rescale_profile_by_mass(pd, ["rho_s"], M, R)
    pm.set_param_values(**pd.param_values)
    assert_allclose(pm(R).in_units("Msun").v, M.v)

def test_mass_rescaling_astropy():
    M = 6.0e14*apu.Msun
    R = 1500.0*apu.kpc
    pm = NFW_mass_profile()
    pd = NFW_density_profile()
    pd.set_param_values(r_s=350*apu.kpc, rho_s=1.0*apu.Msun/apu.kpc**3)
    rescale_profile_by_mass(pd, ["rho_s"], M, R)
    pm.set_param_values(**pd.param_values)
    assert_allclose(pm(R).to("Msun").value, M.value)
