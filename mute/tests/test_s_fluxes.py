import mute.constants as mtc
import mute.surface as mts

import numpy as np

mtc.clear()


def test_s_fluxes():

    s_fluxes_calc = mts.calc_s_fluxes(force=True, output=False)
    s_fluxes_read = mts.load_s_fluxes_from_file()

    assert np.allclose(s_fluxes_calc, s_fluxes_read)
