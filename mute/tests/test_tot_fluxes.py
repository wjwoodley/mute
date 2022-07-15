import os

import numpy as np

import mute.constants as mtc
import mute.surface as mts
import mute.underground as mtu


def test_s_tot_flux():

    mtc.clear()

    s_fluxes = mts.load_s_fluxes_from_file(primary_model="HG")

    s_tot_flux_calc = mts.calc_s_tot_flux(s_fluxes=s_fluxes)
    s_tot_flux_read = 0.011555347252101776

    assert np.allclose(s_tot_flux_calc, s_tot_flux_read)


def test_u_tot_flux_flat():

    mtc.clear()

    mtc.set_output(False)
    mtc.set_overburden("flat")
    mtc.set_vertical_depth(3.5)

    u_fluxes = mtu.calc_u_fluxes(full_tensor=True)

    u_tot_flux_calc = mtu.calc_u_tot_flux(u_fluxes=u_fluxes, E_th=100, force=True)
    u_tot_flux_read = 1.123374482778583e-08

    assert np.allclose(u_tot_flux_calc, u_tot_flux_read)


def test_u_tot_flux_mountain():

    mtc.clear()

    mtc.set_output(False)
    mtc.set_overburden("mountain")

    mtc.load_mountain("example_mountain_profile.txt")

    u_tot_flux_calc = mtu.calc_u_tot_flux(interaction_model="DDM", force=True)
    u_tot_flux_read = 1.4121546671702574e-06

    assert np.allclose(u_tot_flux_calc, u_tot_flux_read)
