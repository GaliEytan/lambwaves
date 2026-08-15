from __future__ import annotations

import numpy as np


def switched_harmonic(t, omega_n, gamma, Omega, force=1.0):
    """Exact complex response with zero initial displacement and velocity."""
    t = np.asarray(t)
    H = force / (omega_n**2 - Omega**2 - 2j*gamma*Omega)
    wd = np.sqrt(omega_n**2 - gamma**2 + 0j)
    sp, sm = -gamma + 1j*wd, -gamma - 1j*wd
    Cp = H*(1j*Omega + sm)/(sp-sm)
    Cm = -H-Cp
    steady = H*np.exp(-1j*Omega*t)
    transient = Cp*np.exp(sp*t) + Cm*np.exp(sm*t)
    return steady + transient, steady, transient


def initial_condition_residual(omega_n=1.2, gamma=.01, Omega=1.0):
    H = 1/(omega_n**2-Omega**2-2j*gamma*Omega)
    wd = np.sqrt(omega_n**2-gamma**2+0j)
    sp, sm = -gamma+1j*wd, -gamma-1j*wd
    Cp = H*(1j*Omega+sm)/(sp-sm); Cm = -H-Cp
    return abs(H+Cp+Cm), abs(-1j*Omega*H+sp*Cp+sm*Cm)
