import numpy as np
import pytest
import sys

import mute.constants as mtc
import mute.surface as mts

try:

    import crflux.models as pm

except ImportError:

    pass


def test_s_fluxes():

    mtc.clear()

    s_fluxes_calc = mts.calc_s_fluxes(model="daemonflux", output=False, test=True)
    s_fluxes_read = mts.load_s_fluxes_from_file(model="daemonflux", test=True)

    assert np.allclose(s_fluxes_calc, s_fluxes_read)


def test_s_fluxes_mceq():

    mtc.clear()

    s_fluxes_calc = mts.calc_s_fluxes(
        model="mceq", interaction_model="sibyll23c", output=False, test=True
    )
    s_fluxes_read = mts.load_s_fluxes_from_file(
        model="mceq", interaction_model="sibyll23c", test=True
    )

    assert np.allclose(s_fluxes_calc, s_fluxes_read)


@pytest.mark.skipif("crflux.models" not in sys.modules, reason="Requires crflux.")
def test_s_fluxes_primary_model():

    mtc.clear()

    s_fluxes_calc = mts.calc_s_fluxes(
        model="mceq",
        primary_model=(pm.GaisserStanevTilav, "3-gen"),
        output=False,
        test=True,
    )
    s_fluxes_read = mts.load_s_fluxes_from_file(
        model="mceq", primary_model="gst3", test=True
    )

    assert np.allclose(s_fluxes_calc, s_fluxes_read)


def test_s_fluxes_location():

    mtc.clear()

    s_fluxes_calc = mts.calc_s_fluxes(
        model="mceq",
        interaction_model="sibyll23c",
        atmosphere="msis00",
        month="July",
        location=(46.472, -81.187),
        output=False,
        test=True,
    )
    s_fluxes_read = mts.load_s_fluxes_from_file(
        model="mceq",
        interaction_model="sibyll23c",
        atmosphere="msis00",
        month="July",
        location=(46.472, -81.187),
        test=True,
    )

    assert np.allclose(s_fluxes_calc, s_fluxes_read)
