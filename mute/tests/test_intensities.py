import mute.constants as mtc
import mute.underground as mtu

import numpy as np

mtc.clear()


def test_u_intensities():

    u_intensities_calc = mtu.calc_u_intensities(force=True, output=False)
    u_intensities_read = [
        8.38550525e-07,
        3.99126333e-07,
        1.54879991e-07,
        9.35155921e-08,
        3.84605857e-08,
        2.59495635e-08,
        1.99959054e-08,
        9.16859070e-09,
        5.79662798e-09,
        4.12477377e-09,
        1.79575398e-09,
        9.19703745e-10,
        4.80652437e-10,
        5.83854128e-10,
        2.04629714e-10,
        1.06581030e-10,
        9.71642206e-11,
        3.11745792e-11,
        4.15591218e-11,
        1.99809107e-11,
        7.48659868e-12,
        3.07673449e-12,
        5.15677696e-12,
    ]

    assert np.allclose(u_intensities_calc, u_intensities_read)
