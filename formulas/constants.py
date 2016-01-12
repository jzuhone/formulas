from sympy.core.numbers import Pi
from formulas.base import FormulaConstant
try:
    import yt.utilities.physical_constants as yt_pc
except ImportError:
    yt_pc = None
try:
    import astropy.constants as astropy_pc
except ImportError:
    astropy_pc = None

pi = Pi()

yt_map = {"m_e":"me","m_p":"mp","m_h":"mh","k_B":"kboltz"}
astropy_map = {}

class PhysicalConstants(object):
    def __init__(self, constants, map):
        self.constants = constants
        self.map = map

    def __getattr__(self, item):
        const = self.map.get(item, item)
        return FormulaConstant(item, getattr(self.constants, const))

if yt_pc is not None:
    yt_constants = PhysicalConstants(yt_pc, yt_map)
if astropy_pc is not None:
    astropy_constants = PhysicalConstants(astropy_pc, astropy_map)