import numpy as np
from lambwaves.dispersion import characteristic, low_frequency_errors, solve_branches
from lambwaves.material import Material


def test_roots_have_small_residual():
    m=Material(); K=np.linspace(.05,2,30)
    roots,res=solve_branches(K,"S",m,count=2,nscan=3000,wmax=6)
    assert np.nanmax(res) < 1e-7
    assert np.isfinite(roots[:,0]).all()


def test_low_frequency_limits():
    a,s=low_frequency_errors(Material())
    assert a < .005
    assert s < .001


def test_family_validation():
    try: characteristic(1.0,.5,"X",Material())
    except ValueError: pass
    else: raise AssertionError("invalid family accepted")
