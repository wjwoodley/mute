import os

import numpy as np

import mute.constants as mtc
import mute.surface as mts
import mute.underground as mtu


def test_s_tot_flux():

    mtc.clear()

    s_fluxes = mts.load_s_fluxes_from_file(
        model="mceq", interaction_model="sibyll23c", primary_model="hg", output=False
    )

    s_tot_flux_calc = mts.calc_s_tot_flux(s_fluxes=s_fluxes)
    s_tot_flux_read = 0.011553938671563626

    assert np.allclose(s_tot_flux_calc, s_tot_flux_read)


def test_u_tot_flux_flat():

    mtc.clear()

    mtc.set_output(False)
    mtc.set_overburden("flat")
    mtc.set_vertical_depth(3.5)

    u_fluxes = mtu.calc_u_fluxes(full_tensor=True, model="daemonflux", output=False)

    u_tot_flux_calc = mtu.calc_u_tot_flux(u_fluxes=u_fluxes, E_th=100, force=True)
    u_tot_flux_read = 1.4480871762241838e-08

    assert np.allclose(u_tot_flux_calc, u_tot_flux_read)


def test_u_tot_flux_mountain():

    mtc.clear()

    mtc.set_output(False)
    mtc.set_overburden("mountain")

    mountain_path = os.path.join(
        os.path.dirname(__file__), "example_mountain_profile.txt"
    )

    mtc.load_mountain(mountain_path)

    u_tot_flux_calc = mtu.calc_u_tot_flux(
        model="mceq", interaction_model="ddm", output=False, force=True
    )
    u_tot_flux_read = 1.400753073325132e-06

    assert np.allclose(u_tot_flux_calc, u_tot_flux_read)
