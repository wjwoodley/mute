import os

import numpy as np

import mute.constants as mtc
import mute.underground as mtu


def test_u_tot_flux():

    mtc.clear()

    mtc.set_output(False)
    mtc.set_vertical_depth(3.5)

    u_tot_flux_calc = mtu.calc_u_tot_flux(E_th = 100, force=True)
    u_tot_flux_read = 1.123374482778583e-08

    assert np.allclose(u_tot_flux_calc, u_tot_flux_read)
