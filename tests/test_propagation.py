import pytest
import sys

import numpy as np

import mute.constants as mtc
import mute.propagation as mtp

try:

    import proposal as pp

except ImportError:

    pass


@pytest.mark.skipif("proposal" not in sys.modules, reason="Requires proposal.")
def test_propagation():

    mtc.clear()

    mtc.set_n_muon(3)
    mtp._create_propagator(force=True)
    pp.RandomGenerator.get().set_seed(500)

    u_energy_calc = mtp._propagation_loop(
        mtc.ENERGIES[50], mtc.slant_depths[0], force=True
    )
    u_energy_read = [8088234.57870225, 7698468.554854496, 7568635.147518327]

    assert np.allclose(u_energy_calc, u_energy_read)
