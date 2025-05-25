import os

import numpy as np

import mute.constants as mtc
import mute.underground as mtu


def test_u_intensities():

    mtc.clear()

    u_intensities_calc = mtu.calc_u_intensities(
        method="sd", model="daemonflux", output=False, force=True
    )
    u_intensities_read = [
        6.94919658e-06,
        1.80906517e-06,
        7.24444624e-07,
        3.32127337e-07,
        1.70705322e-07,
        9.22534290e-08,
        5.07858934e-08,
        2.85312282e-08,
        1.63857159e-08,
        9.51599387e-09,
        5.56572979e-09,
        3.27576784e-09,
        1.89893763e-09,
        1.10476908e-09,
        6.45823164e-10,
        3.79235329e-10,
        2.23574504e-10,
        1.31776529e-10,
        7.82436124e-11,
        4.62467535e-11,
        2.74526405e-11,
        1.62848533e-11,
        9.64117005e-12,
        5.71209783e-12,
        3.39018880e-12,
        2.00295030e-12,
        1.17807861e-12,
        6.95181479e-13,
    ]

    assert np.allclose(u_intensities_calc, u_intensities_read)
