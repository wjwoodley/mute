import mute.constants as mtc
import mute.propagation as mtp

import numpy as np
import pytest
import sys

try:

    import proposal as pp

except ImportError:

    pass

mtc.clear()


@pytest.mark.skipif("proposal" not in sys.modules, reason="Requires proposal.")
def test_propagation():

    mtc.set_n_muon(3)
    mtp._create_propagator()
    pp.RandomGenerator.get().set_seed(500)

    u_energy_calc = mtp._propagation_loop(mtc.energies[50], mtc.slant_depths[0])
    u_energy_read = [6935594.383751289, 4372686.094864153, 7017531.178970211]

    assert np.allclose(u_energy_calc, u_energy_read)
