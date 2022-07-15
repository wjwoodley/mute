import numpy as np

import mute.constants as mtc
import mute.surface as mts


def test_s_fluxes():

    mtc.clear()
    mtc.set_directory(".")

    s_fluxes_calc = mts.calc_s_fluxes(output=False, test=True)
    s_fluxes_read = mts.load_s_fluxes_from_file(test=True)

    assert np.allclose(s_fluxes_calc, s_fluxes_read)
