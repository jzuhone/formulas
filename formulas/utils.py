import numpy as np

try:
    from yt.units.yt_array import YTArray, YTQuantity
except ImportError:
    YTArray = object
    YTQuantity = object

try:
    from astropy.units import Quantity
except ImportError:
    Quantity = object

try:
    from pint import UnitRegistry
    ureg = UnitRegistry(system='cgs')
    PintQuantity = ureg.Quantity
except ImportError:
    PintQuantity = object

def check_type(x):
    if isinstance(x, (YTArray, YTQuantity)):
        return YTArray
    elif isinstance(x, Quantity):
        return Quantity
    elif isinstance(x, PintQuantity):
        return PintQuantity
    else:
        return np.ndarray

def in_cgs(x):
    if isinstance(x, YTArray):
        return x.in_cgs()
    elif isinstance(x, Quantity):
        return x.cgs
    elif isinstance(x, PintQuantity):
        return x.to_base_units()
    else:
        return x

def in_units(x, units):
    if isinstance(x, YTArray):
        return x.in_units(units)
    elif isinstance(x, (Quantity, PintQuantity)):
        return x.to(units)

def get_units(x):
    if isinstance(x, (YTArray, PintQuantity)):
        return x.units
    elif isinstance(x, Quantity):
        return x.unit

def latexify_units(x):
    if isinstance(x, YTArray):
        return "$"+x.units.latex_representation()+"$"
    elif isinstance(x, Quantity):
        return x.unit.to_string("latex")
    elif isinstance(x, PintQuantity):
        return '{:L}'.format(x.units)