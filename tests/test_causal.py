import numpy as np
from lambwaves.causal import initial_condition_residual, switched_harmonic


def test_zero_initial_conditions():
    u,v=initial_condition_residual()
    assert u < 1e-13 and v < 1e-13


def test_transient_decays():
    t=np.array([0.,1000.])
    _,_,tr=switched_harmonic(t,1.2,.01,1.0)
    assert abs(tr[-1]) < abs(tr[0])*1e-3
