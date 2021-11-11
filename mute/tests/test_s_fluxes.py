import numpy as np

import mute.constants as mtc
import mute.surface as mts

mtc.clear()


def test_s_fluxes():

    s_fluxes_calc = mts.calc_s_fluxes(output=False, force=True)
    s_fluxes_read = mts.load_s_fluxes_from_file()

    assert np.allclose(s_fluxes_calc, s_fluxes_read)
