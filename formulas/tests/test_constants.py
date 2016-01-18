from formulas.constants import yt_constants, astropy_constants
import yt.utilities.physical_constants as yt_pc
import astropy.constants as astropy_pc

def test_constants():
    k_B = yt_constants.k_B
    assert k_B.value == yt_pc.kboltz
    c = astropy_constants.c
    assert c.value == astropy_pc.c