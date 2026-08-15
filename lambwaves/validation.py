from __future__ import annotations

import numpy as np

from .causal import initial_condition_residual
from .dispersion import low_frequency_errors, solve_branches
from .material import Material


def validation_report(material: Material, K: np.ndarray):
    AW, Ar = solve_branches(K, "A", material)
    SW, Sr = solve_branches(K, "S", material)
    AWc, _ = solve_branches(K, "A", material, nscan=3500)
    SWc, _ = solve_branches(K, "S", material, nscan=3500)
    aerr, serr = low_frequency_errors(material)
    ic_u, ic_v = initial_condition_residual()
    return {
        "A": AW, "S": SW,
        "max_residual": float(np.nanmax([Ar, Sr])),
        "scan_convergence": float(np.nanmax(abs(np.r_[AW-AWc, SW-SWc]))),
        "A0_low_frequency_error": aerr,
        "S0_low_frequency_error": serr,
        "initial_displacement_residual": float(ic_u),
        "initial_velocity_residual": float(ic_v),
    }
