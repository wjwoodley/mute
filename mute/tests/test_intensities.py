import os

import numpy as np

import mute.constants as mtc
import mute.underground as mtu


def test_u_intensities():

    mtc.clear()

    u_intensities_calc = mtu.calc_u_intensities(output=False, force=True)
    u_intensities_read = [
        9.002276027162749e-07,
        3.538269015684193e-07,
        1.6283577432654826e-07,
        8.406651142850705e-08,
        4.6122098928456166e-08,
        2.621442646505831e-08,
        1.5111826767744374e-08,
        8.887248723219226e-09,
        5.316230203572155e-09,
        3.1847620816949173e-09,
        1.910456525776394e-09,
        1.1570508762100986e-09,
        7.024617216480318e-10,
        4.2590492177512545e-10,
        2.578458039056593e-10,
        1.5681305256751851e-10,
        9.494820560621108e-11,
        5.741622006482442e-11,
        3.496678873927601e-11,
        2.091468781658136e-11,
        1.2574666445503009e-11,
        7.593747060608253e-12,
        4.604313931291195e-12,
    ]

    assert np.allclose(u_intensities_calc, u_intensities_read)
