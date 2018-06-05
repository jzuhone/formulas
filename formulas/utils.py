import numpy as np

class NotAvailable(object):
    pass


try:
    from yt.units.yt_array import YTArray, YTQuantity
except ImportError:
    YTArray = NotAvailable
    YTQuantity = NotAvailable

try:
    from astropy.units import Quantity
except ImportError:
    Quantity = NotAvailable

try:
    from pint import UnitRegistry
    ureg = UnitRegistry(system='cgs')
    PintQuantity = ureg.Quantity
except ImportError:
    PintQuantity = NotAvailable

try:
    from unyt import unyt_array, unyt_quantity
except ImportError:
    unyt_quantity = NotAvailable
    unyt_array = NotAvailable


def check_type(x):
    if isinstance(x, (YTArray, YTQuantity)):
        return YTArray
    elif isinstance(x, (unyt_array, unyt_quantity)):
        return unyt_array
    elif isinstance(x, Quantity):
        return Quantity
    elif isinstance(x, PintQuantity):
        return PintQuantity
    else:
        return np.ndarray


def in_cgs(x):
    if isinstance(x, (YTArray, unyt_array)):
        return x.in_cgs()
    elif isinstance(x, Quantity):
        return x.cgs
    elif isinstance(x, PintQuantity):
        return x.to_base_units()
    else:
        return x


def in_units(x, units):
    return x.to(units)


def get_units(x):
    if hasattr(x, "units"):
        return x.units
    elif hasattr(x, "unit"):
        return x.unit
    else:
        raise RuntimeError("No units are attached to this object!!")


def latexify_units(x):
    if isinstance(x, (YTArray, unyt_array)):
        return "$"+x.units.latex_representation()+"$"
    elif isinstance(x, Quantity):
        return x.unit.to_string("latex")
    elif isinstance(x, PintQuantity):
        return '{:L}'.format(x.units)