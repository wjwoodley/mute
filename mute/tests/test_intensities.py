import os

import numpy as np

import mute.constants as mtc
import mute.underground as mtu


def test_u_intensities():

    mtc.clear()

    u_intensities_calc = mtu.calc_u_intensities(method="sd", output=False, force=True)
    u_intensities_read = [
        5.44204714e-06,
        1.36106429e-06,
        5.38845074e-07,
        2.47298461e-07,
        1.26498482e-07,
        6.82698286e-08,
        3.78720866e-08,
        2.13392590e-08,
        1.22844489e-08,
        7.19956984e-09,
        4.23236899e-09,
        2.49549138e-09,
        1.45351792e-09,
        8.50501192e-10,
        4.99272534e-10,
        2.93738116e-10,
        1.74149878e-10,
        1.03076706e-10,
        6.10735102e-11,
        3.65166249e-11,
        2.14842462e-11,
        1.27262101e-11,
        7.58114103e-12,
        4.54062606e-12,
        2.67074717e-12,
        1.60773886e-12,
        9.41137354e-13,
        5.52113709e-13,
    ]

    assert np.allclose(u_intensities_calc, u_intensities_read)
