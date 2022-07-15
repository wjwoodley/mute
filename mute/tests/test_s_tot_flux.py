import os

import numpy as np

import mute.constants as mtc
import mute.surface as mts


def test_u_tot_flux():

    mtc.clear()
    
    s_fluxes = mts.load_s_fluxes_from_file(primary_model = "HG")
    
    s_tot_flux_calc = mts.calc_s_tot_flux(s_fluxes = s_fluxes)
    s_tot_flux_read = 0.011555347252101776

    assert np.allclose(s_tot_flux_calc, s_tot_flux_read)
