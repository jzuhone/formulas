from formulas.base import Formula1D
from sympy import symbols, pi, log, exp
from mpmath import quad
from formulas.utils import in_cgs, get_units, in_units, check_type
import numpy as np


def beta_model_profile(r="r", rho_c="rho_c", r_c="r_c", beta="beta"):
    """
    A beta-model density profile (like for galaxy clusters,
    Cavaliere A., Fusco-Femiano R., 1976, A&A, 49, 137).

    Parameters
    ----------
    r : string
        The symbol for the radius variable.
    rho_c : string
        The symbol for the scale density.
    r_c : string
        The symbol for the core radius.
    beta : string
        The symbol for the beta parameter.
    """
    r, rho_c, r_c, beta = symbols((r, rho_c, r_c, beta))
    profile = rho_c*((1+(r/r_c)**2)**(-3*beta/2))
    return Formula1D(profile, r, [r_c, rho_c, beta])


def vikhlinin_density_profile(r="r", rho_0="rho_0", r_c="r_c",
                              r_s="r_s", alpha="alpha", beta="beta",
                              epsilon="epsilon", gamma=None):
    """
    A modified beta-model density profile for galaxy
    clusters from Vikhlinin, A., Kravtsov, A., Forman, W.,
    et al. 2006, ApJ, 640, 691.

    Parameters
    ----------
    r : string
        The symbol for the radius variable.
    rho_0 : string
        The symbol for the scale density of the profile.
    r_c : string
        The symbol for the core radius.
    r_s : string
        The symbol for the scale radius.
    alpha : string
        The symbol for the inner logarithmic slope parameter.
    beta : string
        The symbol for the middle logarithmic slope parameter.
    epsilon : string
        The symbol for the outer logarithmic slope parameter.
    gamma : string
        The symbol for controlling the width of the outer
        transition. If None, it will be gamma = 3 by default.
    """
    r, rho_0, r_c, r_s, alpha, beta, epsilon = \
        symbols((r, rho_0, r_c, r_s, alpha, beta, epsilon))
    params = [rho_0, r_c, r_s, alpha, beta, epsilon]
    if gamma is None:
        gamma = 3
    else:
        gamma = symbols(gamma)
        params += [gamma]
    profile = rho_0*(r/r_c)**(-alpha/2) * \
        (1+(r/r_c)**2)**(-3*beta/2+alpha/4) * \
        (1+(r/r_s)**gamma)**(-epsilon/gamma/2)
    return Formula1D(profile, r, params)


def hernquist_density_profile(r="r", M_0="M_0", a="a"):
    """
    A Hernquist density profile (Hernquist, L. 1990,
    ApJ, 356, 359).

    Parameters
    ----------
    r : string
        The symbol for the radius variable.
    M_0 : string
        The symbol for the total mass of the profile.
    a : string
        The symbol for the scale radius.
    """
    r, M_0, a = symbols((r, M_0, a))
    profile = M_0/(2*pi*a**3)/((r/a)*(1+r/a)**3)
    return Formula1D(profile, r, [M_0, a])


def hernquist_mass_profile(r="r", M_0="M_0", a="a"):
    """
    A Hernquist mass profile (Hernquist, L. 1990,
    ApJ, 356, 359).

    Parameters
    ----------
    r : string
        The symbol for the radius variable.
    M_0 : string
        The symbol for the total mass of the profile.
    a : string
        The symbol for the scale radius.
    """
    r, M_0, a = symbols((r, M_0, a))
    profile = M_0*r**2/(r+a)**2
    return Formula1D(profile, r, [M_0, a])



def NFW_density_profile(r="r", rho_s="rho_s", r_s="r_s"):
    """
    An NFW density profile (Navarro, J.F., Frenk, C.S.,
    & White, S.D.M. 1996, ApJ, 462, 563).

    Parameters
    ----------
    r : string
        The symbol for the radius variable.
    rho_s : string
        The symbol for the scale density of the profile.
    r_s : string
        The symbol for the scale radius.
    """
    r, rho_s, r_s = symbols((r, rho_s, r_s))
    profile = rho_s/((r/r_s)*(1+r/r_s)**2)
    return Formula1D(profile, r, [rho_s, r_s])


def NFW_mass_profile(r="r", rho_s="rho_s", r_s="r_s"):
    """
    An NFW mass profile (Navarro, J.F., Frenk, C.S.,
    & White, S.D.M. 1996, ApJ, 462, 563).

    Parameters
    ----------
    r : string
        The symbol for the radius variable.
    rho_s : string
        The symbol for the scale density of the profile.
    r_s : string
        The symbol for the scale radius.
    """
    r, rho_s, r_s = symbols((r, rho_s, r_s))
    x = r/r_s
    profile = 4*pi*rho_s*r_s**3*(log(1+x)-x/(1+x))
    return Formula1D(profile, r, [rho_s, r_s])


def sNFW_density_profile(r="r", M="M", a="a"):
    """
    A "super-NFW" density profile (Lilley, E. J.,
    Wyn Evans, N., & Sanders, J.L. 2018, MNRAS).

    Parameters
    ----------
    r : string
        The symbol for the radius variable.
    M : string
        The symbol for the total mass of the profile.
    a : string
        The symbol for the scale radius.
    """
    r, M, a = symbols((r, M, a))
    x = r/a
    profile = 3*M/(16*pi*a**3)/(x*(1+x)**(5/2))
    return Formula1D(profile, r, [M, a])


def sNFW_mass_profile(r="r", M="M", a="a"):
    """
    A "super-NFW" mass profile (Lilley, E. J.,
    Wyn Evans, N., & Sanders, J.L. 2018, MNRAS).

    Parameters
    ----------
    r : string
        The symbol for the radius variable.
    M : string
        The symbol for the total mass of the profile.
    a : string
        The symbol for the scale radius.
    """
    r, M, a = symbols((r, M, a))
    x = r/a
    profile = M*(1-(2+3*x)/(2*(1+x)**(3/2)))
    return Formula1D(profile, r, [M, a])


def einasto_density_profile(r="r", rho_0="rho_0", h="h", alpha="alpha"):
    """
    A density profile where the logarithmic slope is a 
    power-law. The form here is that given in Equation 5 of
    Retana-Montenegro et al. 2012, A&A, 540, A70.

    Parameters
    ----------
    r : string
        The symbol for the radius variable.
    rho_0 : string
        The symbol for the core density.
    h : string
        The symbol for the scale radius.
    alpha : string
        The symbol for the power-law index.
    """
    r, rho_0, h, alpha = symbols((r, rho_0, h, alpha))
    x = r/h
    profile = rho_0*exp(-x**alpha)
    return Formula1D(profile, r, [rho_0, h, alpha])


def AM06_density_profile(r="r", rho_0="rho_0", a="a", a_c="a_c", c="c",
                         alpha="alpha", beta="beta"):
    """
    The density profile for galaxy clusters suggested by
    Ascasibar, Y., & Markevitch, M. 2006, ApJ, 650, 102.
    Works best in concert with the ``AM06_temperature_profile``.

    Parameters
    ----------
    r : string
        The symbol for the radius variable.
    rho_0 : string
        The symbol for the scale density of the profile.
    a : string
        The symbol for the scale radius.
    a_c : string
        The symbol for the scale radius of the cool core.
    c : string
        The symbol for the scale of the temperature drop of the cool
        core.
    alpha : string
        The symbol for the first slope parameter.
    beta : string
        The symbol for the second slope parameter.
    """
    r, rho_0, a, a_c, c, alpha, beta = symbols((r, rho_0, a, a_c, c, 
                                                alpha, beta))
    profile = rho_0*(1+r/a_c)*(1+r/a_c/c)**alpha*(1+r/a)**beta
    return Formula1D(profile, r, [rho_0, a, a_c, c, alpha, beta])


def vikhlinin_temperature_profile(r="r", T_0="T_0", a="a", b="b", c="c",
                                  r_t="r_t", T_min="T_min",
                                  r_cool="r_cool", a_cool="a_cool"):
    """
    A temperature profile for galaxy clusters from
    Vikhlinin, A., Kravtsov, A., Forman, W., et al.
    2006, ApJ, 640, 691.

    Parameters
    ----------
    r : string
        The symbol for the radius variable.
    T_0 : string
        The symbol for the scale temperature of the profile.
    a : string
        The symbol for the inner logarithmic slope.
    b : string
        The symbol for the width of the transition region.
    c : string
        The symbol for the outer logarithmic slope.
    r_t : string
        The symbol for the scale radius.
    T_min : string
        The symbol for the minimum temperature.
    r_cool : string
        The symbol for the cooling radius.
    a_cool : string
        The symbol for the logarithmic slope in the cooling region.
    """
    r, T_0, a, b, c, r_t, T_min, r_cool, a_cool = \
        symbols((r, T_0, a, b, c, r_t, T_min, r_cool, a_cool))
    x = (r/r_cool)**a_cool
    t = (r/r_t)**(-a)/((1.+(r/r_t)**b)**(c/b))
    profile = T_0*t*(x+T_min/T_0)/(x+1)
    return Formula1D(profile, r, [T_0, a, b, c, r_t, T_min, r_cool, a_cool])


def AM06_temperature_profile(r="r", T_0="T_0", a="a", a_c="a_c", c="c"):
    """
    The temperature profile for galaxy clusters suggested by
    Ascasibar, Y., & Markevitch, M. 2006, ApJ, 650, 102.
    Works best in concert with the ``AM06_density_profile``.

    Parameters
    ----------
    r : string
        The symbol for the radius variable.
    T_0 : string
        The symbol for the scale temperature of the profile.
    a : string
        The symbol for the scale radius.
    a_c : string
        The symbol for the scale radius of the cool core.
    c : string
        The symbol for the scale of the temperature drop of the cool
        core.
    """
    r, T_0, a, a_c, c = symbols((r, T_0, a, a_c, c))
    profile = T_0/(1+r/a)*(c+r/a_c)/(1+r/a_c)
    return Formula1D(profile, r, [T_0, a, a_c, c])


def baseline_entropy_profile(r="r", K_0="K_0", K_200="K_200", r_200="r_200",
                             alpha="alpha"):
    """
    The baseline entropy profile for galaxy clusters (Voit, G.M.,
    Kay, S.T., & Bryan, G.L. 2005, MNRAS, 364, 909).

    Parameters
    ----------
    r : string
        The symbol for the radius variable.
    K_0 : string
        The symbol for the central entropy floor.
    K_200 : string
        The symbol for the entropy at the radius *r_200*.
    r_200 : string
        The symbol for the virial radius of the profile.
    alpha : string
        The symbol for the logarithmic slope of the profile.
    """
    r, K_0, K_200, r_200, alpha = symbols((r, K_0, K_200, r_200, alpha))
    profile = K_0 + K_200*(r/r_200)**alpha
    return Formula1D(profile, r, [K_0, K_200, r_200, alpha])


def rescale_profile_by_mass(profile, param, mass, radius):
    """
    Rescale a density profile by a total mass
    within some radius.

    Parameters
    ----------
    profile : Formula1D
        Formula that is a radial density profile.
    param : string
        The density-valued parameter that needs to be rescaled.
    mass : YTQuantity
        The mass of the object.
    radius : YTQuantity
        The radius that the *mass* corresponds to.

    Examples
    --------
    >>> import yt.units as u
    >>> import numpy as np
    >>> gas_density = AM06_density_profile()
    >>> a = 600.0*u.kpc
    >>> a_c = 60.0*u.kpc
    >>> c = 0.17
    >>> alpha = -2.0
    >>> beta = -3.0
    >>> M0 = 1.0e14*u.Msun
    >>> # Don't set the density parameter rho_0!
    >>> gas_density.set_param_values(a=a, a_c=a_c, c=c, 
    ...                              alpha=alpha, beta=beta)
    >>> rescale_profile_by_mass(gas_density, "rho_0", M0, np.inf*u.kpc)
    """
    R = float(in_cgs(radius).value)
    M = float(in_cgs(mass).value)
    rho = profile.copy()
    values = {}
    for n, p in profile.param_values.items():
        if n == param:
            # Set the density parameter to change to unity
            values[n] = 1.0
        elif isinstance(p, float):
            values[n] = p
        else:
            values[n] = float(in_cgs(p).value)
    rho.set_param_values(**values)
    mass_int = lambda r: rho(r)*r*r
    scale = float(M/(4.*np.pi*quad(mass_int, [0, R])))
    u = get_units(mass/radius**3)
    quan = check_type(mass)(scale, "g/cm**3")
    profile.param_values[param] = in_units(quan, u)


def _nfw_factor(conc):
    return 1.0/(np.log(conc+1.0)-conc/(1.0+conc))


def compute_NFW_scale_density(conc, z=0.0, delta=200.0, hubble=0.7):
    """
    Compute a scale density parameter for an NFW profile
    given a concentration parameter, and optionally
    a redshift, overdensity, and cosmology.

    Parameters
    ----------
    conc : float
        The concentration parameter for the halo, which should 
        correspond the selected overdensity (which has a default
        of 200). 
    z : float, optional
        The redshift of the halo formation. Default: 0.0
    delta : float, optional
        The overdensity parameter for which the concentration
        is defined. Default: 200.0
    hubble : float, optional
        The Hubble parameter at the current epoch, which 
        determines the value of the critical density.
    """
    import unyt as u
    H0 = hubble*100.0*u.km/u.s/u.Mpc
    rho_crit = (3.0*H0**2/(8.0*np.pi*u.G)).to("Msun/kpc**3")
    rho_crit *= (1.0+z)**3
    rho_s = delta*rho_crit*conc**3*_nfw_factor(conc)
    return rho_s


def convert_NFW_to_hernquist(M200, r200, conc):
    """
    Given M200, r200, and a concentration parameter for an
    NFW profile, return the Hernquist mass and scale radius
    parameters.

    Parameters
    ----------
    M200 : YTQuantity
        The mass of the halo at r200.
    r200 : YTQuantity
        The radius corresponding to the overdensity of 200 times the
        critical density of the universe.
    conc : float
        The concentration parameter r200/r_s for the NFW profile.
    """
    a = r200/(np.sqrt(0.5*conc*conc*_nfw_factor(conc))-1.0)
    M0 = M200*(r200+a)**2/r200**2
    return M0, a
