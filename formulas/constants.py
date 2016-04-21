from sympy import symbols
from sympy import pi as sym_pi
from formulas.base import Formula
try:
    import yt.utilities.physical_constants as yt_pc
except ImportError:
    yt_pc = None
try:
    import astropy.constants as astropy_pc
except ImportError:
    astropy_pc = None
try:
    from pint import UnitRegistry
    pint_pc = UnitRegistry(system='cgs')
except ImportError:
    pint_pc = None

yt_map = {"m_e": "me",
          "m_p": "mp",
          "m_h": "mh",
          "k_B": "kboltz"}
astropy_map = {}
pint_map = {"G": "newtonian_constant_of_gravitation",
            "k_B": "k"}

class FormulaConstant(Formula):
    def __init__(self, name, value):
        name = symbols(name)
        super(FormulaConstant, self).__init__(name, [], [name])
        self.param_values[str(name)] = value
        self._value = value

    def set_param_values(self, **kwargs):
        """
        Set the values of one or more parameters.
        """
        if self.num_params > 0:
            raise RuntimeError("Can't change the value of a constant!")

    def clear_param_values(self):
        """
        Set all of the parameter values to None.
        """
        if self.num_params > 0:
            raise RuntimeError("Can't change the value of a constant!")

    @property
    def value(self):
        return self._value

class FormulaPi(FormulaConstant):
    def __init__(self):
        Formula.__init__(self, sym_pi, [], [])

    @property
    def value(self):
        return self.formula.evalf()

class PhysicalConstants(object):
    def __init__(self, constants, map):
        self.constants = constants
        self.map = map

    def __getattr__(self, item):
        const = self.map.get(item, item)
        return FormulaConstant(item, 1.0*getattr(self.constants, const))

pi = FormulaPi()

if yt_pc is not None:
    yt_constants = PhysicalConstants(yt_pc, yt_map)
if astropy_pc is not None:
    astropy_constants = PhysicalConstants(astropy_pc, astropy_map)
if pint_pc is not None:
    pint_constants = PhysicalConstants(pint_pc, pint_map)