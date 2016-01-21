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

def check_type(x):
    if isinstance(x, (YTArray,YTQuantity)):
        return YTArray
    elif isinstance(x, Quantity):
        return Quantity
    else:
        return np.ndarray

def in_cgs(x):
    if isinstance(x, YTArray):
        return x.in_cgs()
    elif isinstance(x, Quantity):
        return x.cgs
    else:
        return x

def in_units(x, units):
    if isinstance(x, YTArray):
        return x.in_units(units)
    elif isinstance(x, Quantity):
        return x.to(units)

def get_units(x):
    if isinstance(x, YTArray):
        return x.units
    elif isinstance(x, Quantity):
        return x.unit

def latexify_units(x):
    if isinstance(x, YTArray):
        return "$"+x.units.latex_representation()+"$"
    elif isinstance(x, Quantity):
        return x.unit.to_string("latex")