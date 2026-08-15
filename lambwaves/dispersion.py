from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from .material import Material


def characteristic(W, K: float, family: str, material: Material):
    """Scaled determinant in W=omega*h/cT and K=k*h coordinates."""
    W = np.asarray(W)
    r = material.cT / material.cL
    p = np.sqrt((r * W) ** 2 - K**2 + 0j)
    q = np.sqrt(W**2 - K**2 + 0j)
    if family.upper() == "S":
        det = (q**2 - K**2) ** 2 * np.cos(p) * np.sin(q)
        det += 4 * K**2 * p * q * np.sin(p) * np.cos(q)
    elif family.upper() == "A":
        det = 4 * K**2 * p * q * np.cos(p) * np.sin(q)
        det += (q**2 - K**2) ** 2 * np.sin(p) * np.cos(q)
    else:
        raise ValueError("family must be 'S' or 'A'")
    scale = 1 + abs((q**2 - K**2) ** 2) + abs(4 * K**2 * p * q)
    return np.where(W < K, det.imag, det.real) / scale


def _threshold(W: float, K: float, material: Material, tol=5e-3) -> bool:
    r = material.cT / material.cL
    return abs(W - K) < tol or abs(W - K / r) < tol


def roots_at_k(K: float, family: str, material: Material, *, wmax=12.0,
               nscan=7000, count=4) -> np.ndarray:
    scan = np.linspace(1e-5, wmax, nscan)
    values = characteristic(scan, K, family, material)
    crossings = np.flatnonzero(values[:-1] * values[1:] < 0)
    roots: list[float] = []
    for j in crossings:
        if _threshold(scan[j], K, material):
            continue
        try:
            root = brentq(lambda w: float(characteristic(w, K, family, material)),
                          scan[j], scan[j + 1], xtol=1e-12)
        except ValueError:
            continue
        if _threshold(root, K, material) or any(abs(root-r) < 2e-5 for r in roots):
            continue
        roots.append(root)
    return np.asarray(sorted(roots)[:count])


def solve_branches(K: np.ndarray, family: str, material: Material, *, count=3,
                   nscan=7000, wmax=12.0) -> tuple[np.ndarray, np.ndarray]:
    """Solve roots and continue labels by nearest frequency prediction."""
    branches = np.full((len(K), count), np.nan)
    residuals = np.full_like(branches, np.nan)
    for i, kval in enumerate(K):
        roots = roots_at_k(kval, family, material, wmax=wmax, nscan=nscan,
                           count=count + 2)
        if i == 0:
            chosen = roots[:count]
        else:
            chosen = np.full(count, np.nan)
            available = list(roots)
            for b in range(count):
                if not available:
                    break
                prediction = branches[i-1, b]
                if i > 1 and np.isfinite(branches[i-2, b]):
                    prediction += branches[i-1, b] - branches[i-2, b]
                j = int(np.argmin(abs(np.asarray(available) - prediction)))
                chosen[b] = available.pop(j)
        branches[i, :len(chosen)] = chosen
        for b, root in enumerate(chosen):
            if np.isfinite(root):
                residuals[i, b] = abs(characteristic(root, kval, family, material))
    return branches, residuals


def dimensionalize(W: np.ndarray, K: np.ndarray, material: Material):
    omega = W * material.cT / material.h
    k = K / material.h
    vg = np.full_like(omega, np.nan)
    for b in range(W.shape[1]):
        good = np.isfinite(omega[:, b])
        vg[good, b] = np.gradient(omega[good, b], k[good])
    return omega, k, vg


def low_frequency_errors(material: Material) -> tuple[float, float]:
    K = np.linspace(0.002, 0.08, 20)
    AW, _ = solve_branches(K, "A", material, count=1, nscan=4000, wmax=.5)
    SW, _ = solve_branches(K, "S", material, count=1, nscan=4000, wmax=.5)
    Ao, k, _ = dimensionalize(AW, K, material)
    So, _, _ = dimensionalize(SW, K, material)
    D = material.young * material.thickness**3 / (12*(1-material.poisson**2))
    expected_a = np.sqrt(D/(material.rho*material.thickness))*k**2
    expected_s = material.plate_speed*k
    return (float(np.nanmedian(abs(Ao[:, 0]-expected_a)/expected_a)),
            float(np.nanmedian(abs(So[:, 0]-expected_s)/expected_s)))
