from sympy.core.numbers import Pi
from formulas.base import FormulaConstant
import yt.utilities.physical_constants as pc

pi = Pi()

class PhysicalConstants(object):
    def __init__(self, constants):
        self.constants = constants

    def __getattr__(self, item):
        return FormulaConstant(item, getattr(self.constants,item))

yt_constants = PhysicalConstants(pc)