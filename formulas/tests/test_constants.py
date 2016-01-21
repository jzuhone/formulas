from formulas.constants import yt_constants, astropy_constants, pi
import yt.utilities.physical_constants as yt_pc
import astropy.constants as astropy_pc
from numpy.testing import assert_allclose
import numpy as np

def test_constants():
    k_B = yt_constants.k_B
    assert k_B.value == yt_pc.kboltz
    c = astropy_constants.c
    assert c.value == astropy_pc.c
    assert_allclose(float(pi.value), np.pi)