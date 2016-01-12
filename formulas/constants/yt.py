from formulas.base import FormulaConstant
import yt.utilities.physical_constants as pc

# Some physical constants

G = FormulaConstant("G", pc.G)
c = FormulaConstant("c", pc.clight)
m_e = FormulaConstant("m_e", pc.me)
m_p = FormulaConstant("m_p", pc.mp)
k_B = FormulaConstant("k_B", pc.kboltz)
h = FormulaConstant("h", pc.hcgs)
e = FormulaConstant("e", pc.qp)
