"""Reproduce the first Lamb-wave dispersion and causal-buildup study."""
from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq


OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)


@dataclass(frozen=True)
class Material:
    rho: float = 2700.0
    cL: float = 6320.0
    cT: float = 3130.0
    d: float = 1.0e-3

    @property
    def h(self) -> float:
        return self.d / 2


MAT = Material()


def characteristic(W: float | np.ndarray, K: float, family: str) -> np.ndarray:
    """Scaled, real Rayleigh--Lamb characteristic function."""
    W = np.asarray(W)
    ratio = MAT.cT / MAT.cL
    p = np.sqrt((ratio * W) ** 2 - K**2 + 0j)
    q = np.sqrt(W**2 - K**2 + 0j)
    if family == "S":
        det = (q**2 - K**2) ** 2 * np.cos(p) * np.sin(q)
        det += 4 * K**2 * p * q * np.sin(p) * np.cos(q)
    else:
        det = 4 * K**2 * p * q * np.cos(p) * np.sin(q)
        det += (q**2 - K**2) ** 2 * np.sin(p) * np.cos(q)
    scale = 1 + np.abs((q**2 - K**2) ** 2) + np.abs(4 * K**2 * p * q)
    return np.where(W < K, det.imag, det.real) / scale


def near_threshold(W: float, K: float) -> bool:
    r = MAT.cT / MAT.cL
    # The determinant form contains removable p=0 and q=0 factors. A tolerance
    # wider than either scan spacing prevents those bulk-wave thresholds from
    # being misidentified as Lamb roots under grid refinement.
    return abs(W - K) < 5e-3 or abs(W - K / r) < 5e-3


def solve_family(family: str, K: np.ndarray, wmax=12.0, nscan=7000, nbranches=3):
    grid = np.linspace(1e-4, wmax, nscan)
    result = np.full((K.size, nbranches), np.nan)
    residual = np.full_like(result, np.nan)
    for i, kval in enumerate(K):
        values = characteristic(grid, kval, family)
        crossings = np.flatnonzero(values[:-1] * values[1:] < 0)
        roots: list[float] = []
        for j in crossings:
            if near_threshold(grid[j], kval):
                continue
            try:
                root = brentq(lambda w: float(characteristic(w, kval, family)),
                              grid[j], grid[j + 1], xtol=1e-11)
            except ValueError:
                continue
            if near_threshold(root, kval) or any(abs(root - x) < 2e-4 for x in roots):
                continue
            roots.append(root)
        roots.sort()
        for branch, root in enumerate(roots[:nbranches]):
            result[i, branch] = root
            residual[i, branch] = abs(characteristic(root, kval, family))
    return result, residual


def kinematics(W: np.ndarray, K: np.ndarray):
    omega = W * MAT.cT / MAT.h
    k = K / MAT.h
    vg = np.full_like(omega, np.nan)
    curvature = np.full_like(omega, np.nan)
    for b in range(W.shape[1]):
        good = np.isfinite(omega[:, b])
        vg[good, b] = np.gradient(omega[good, b], k[good])
        curvature[good, b] = np.gradient(vg[good, b], k[good])
    return omega, k, vg, curvature


def local_minima(omega: np.ndarray, K: np.ndarray, family: str, max_k=1.5):
    found = []
    for b in range(1, omega.shape[1]):
        y = omega[:, b]
        for i in range(2, K.size - 2):
            if K[i] <= max_k and np.all(np.isfinite(y[i - 1:i + 2])) and y[i] < y[i - 1] and y[i] < y[i + 1]:
                # Quadratic refinement using the three neighboring grid points.
                coeff = np.polyfit(K[i - 1:i + 2], y[i - 1:i + 2], 2)
                kr = -coeff[1] / (2 * coeff[0])
                wr = np.polyval(coeff, kr)
                found.append((f"{family}{b}", kr, wr / (2*np.pi)*MAT.d/1000, wr/(2*np.pi)/1e6))
    return found


def modal_integral(wn, weights, Omega, gamma, t):
    H = weights / (wn**2 - Omega**2 - 2j*gamma*Omega)
    # Complex square root handles both underdamped and overdamped spectral
    # components without collapsing the two characteristic roots.
    wd = np.sqrt(wn**2 - gamma**2 + 0j)
    sp, sm = -gamma + 1j*wd, -gamma - 1j*wd
    Cp = H * (1j*Omega + sm) / (sp - sm)
    Cm = -H - Cp
    ss = H.sum() * np.exp(-1j*Omega*t)
    transient = (np.exp(np.outer(t, sp)) @ Cp + np.exp(np.outer(t, sm)) @ Cm)
    return ss + transient, ss, transient


def break_plot_jumps(y: np.ndarray, threshold=0.35) -> np.ndarray:
    """Insert NaNs around root-index discontinuities for honest plotting."""
    clean = y.copy()
    jumps = np.flatnonzero(np.abs(np.diff(clean)) > threshold)
    clean[jumps] = np.nan
    clean[jumps + 1] = np.nan
    return clean


def causal_study(K, k, Aomega, Somega, Avg, x=20e-3, width=1e-3,
                 damping_ratio=0.01, ntime=2501):
    target = np.nanargmin(abs(K - 1.0))
    Omega = Aomega[target, 0]
    gamma = damping_ratio * Omega
    # Isolate the driven modal band. Very-low-frequency off-resonant components
    # represent a quasistatic turn-on field and should not be included in the
    # settling metric of the frequency-filtered Lamb-wave signal.
    good_a = np.isfinite(Aomega[:, 0]) & (Aomega[:, 0] > 0.5*Omega) & (Aomega[:, 0] < 1.5*Omega)
    good_s = np.isfinite(Somega[:, 0]) & (Somega[:, 0] > 0.5*Omega) & (Somega[:, 0] < 1.5*Omega)
    ka, ks = k[good_a], k[good_s]
    weights_a = np.exp(-(ka*width)**2/4) * np.cos(ka*x) * np.mean(np.diff(k))/np.pi
    weights_s = np.exp(-(ks*width)**2/4) * np.cos(ks*x) * np.mean(np.diff(k))/np.pi
    period = 2*np.pi/Omega
    t = np.linspace(0, 250*period, ntime)
    ua, ssa, tra = modal_integral(Aomega[good_a, 0], weights_a, Omega, gamma, t)
    us, _, _ = modal_integral(Somega[good_s, 0], weights_s, Omega, gamma, t)
    ratio = abs(tra) / max(abs(ssa[0]), np.finfo(float).eps)
    future_max = np.maximum.accumulate(ratio[::-1])[::-1]
    settled = np.flatnonzero(future_max < 0.05)
    tss = t[settled[0]] if settled.size else np.nan
    purity = abs(ua)**2/(abs(ua)**2 + abs(us)**2 + np.finfo(float).eps)
    return dict(t=t, period=period, frequency=Omega/(2*np.pi), Omega=Omega,
                ua=ua, us=us, ssa=ssa, tra=tra, ratio=ratio, purity=purity,
                tss=tss, arrival=x/abs(Avg[target, 0]), N95=np.ceil(tss/period),
                x=x, width=width, gamma=gamma)


def main():
    K = np.linspace(0.04, 8.0, 260)
    print("Solving antisymmetric branches...")
    AW, Ares = solve_family("A", K)
    print("Solving symmetric branches...")
    SW, Sres = solve_family("S", K)
    Aomega, k, Avg, Acurv = kinematics(AW, K)
    Somega, _, Svg, Scurv = kinematics(SW, K)

    # Thin-plate validation over the smallest sampled K values.
    ratio = MAT.cT / MAT.cL
    E = MAT.rho*MAT.cT**2*(3*MAT.cL**2 - 4*MAT.cT**2)/(MAT.cL**2 - MAT.cT**2)
    nu = (MAT.cL**2 - 2*MAT.cT**2)/(2*(MAT.cL**2 - MAT.cT**2))
    cplate = np.sqrt(E/(MAT.rho*(1-nu**2)))
    rigidity = E*MAT.d**3/(12*(1-nu**2))
    Klow = np.linspace(0.002, 0.08, 20)
    AWlow, _ = solve_family("A", Klow, wmax=0.5, nscan=4000, nbranches=1)
    SWlow, _ = solve_family("S", Klow, wmax=0.5, nscan=4000, nbranches=1)
    Aolow, klow, _, _ = kinematics(AWlow, Klow)
    Solow, _, _, _ = kinematics(SWlow, Klow)
    s_expected = cplate*klow
    a_expected = np.sqrt(rigidity/(MAT.rho*MAT.d))*klow**2
    s_err = np.nanmedian(abs(Solow[:, 0]-s_expected)/s_expected)
    a_err = np.nanmedian(abs(Aolow[:, 0]-a_expected)/a_expected)
    zgvs = local_minima(Aomega, K, "A") + local_minima(Somega, K, "S")
    transient = causal_study(K, k, Aomega, Somega, Avg)

    # Convergence checks: repeat the characteristic scan at half resolution,
    # and repeat the transient calculation with half the time step.
    AWcoarse, _ = solve_family("A", K, nscan=3500)
    SWcoarse, _ = solve_family("S", K, nscan=3500)
    dispersion_convergence = np.nanmax(np.abs(np.r_[AWcoarse-AW, SWcoarse-SW]))
    transient_fine = causal_study(K, k, Aomega, Somega, Avg, ntime=5001)
    time_convergence_cycles = abs(transient_fine["tss"]-transient["tss"])/transient["period"]

    # Verify the analytic oscillator constants enforce zero initial conditions.
    test_wn, test_Omega, test_gamma = 1.2, 1.0, 0.01
    H = 1/(test_wn**2-test_Omega**2-2j*test_gamma*test_Omega)
    wd = np.sqrt(test_wn**2-test_gamma**2+0j)
    sp, sm = -test_gamma+1j*wd, -test_gamma-1j*wd
    Cp = H*(1j*test_Omega+sm)/(sp-sm); Cm = -H-Cp
    initial_displacement_error = abs(H+Cp+Cm)
    initial_velocity_error = abs(-1j*test_Omega*H+sp*Cp+sm*Cm)

    # Parameter map for the design observables requested in the manuscript.
    sweep = []
    for distance in (10e-3, 20e-3, 40e-3):
        for width in (0.5e-3, 1e-3, 2e-3):
            for damping in (0.005, 0.01, 0.02):
                q = causal_study(K, k, Aomega, Somega, Avg, x=distance,
                                 width=width, damping_ratio=damping, ntime=2501)
                sweep.append((distance*1e3, width*1e3, damping,
                              q["arrival"]*1e6, q["tss"]*1e6, q["N95"],
                              np.mean(q["purity"][-100:])))

    # Machine-readable output.
    headers = ["K"]
    columns = [K]
    for family, omega, vg in (("A", Aomega, Avg), ("S", Somega, Svg)):
        for b in range(3):
            headers += [f"{family}{b}_fd_MHz_mm", f"{family}{b}_vg_m_s"]
            columns += [omega[:, b]/(2*np.pi)*MAT.d/1000, vg[:, b]]
    with (OUT/"dispersion.csv").open("w", newline="") as f:
        writer = csv.writer(f); writer.writerow(headers); writer.writerows(zip(*columns))
    with (OUT/"parameter_sweep.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["distance_mm", "gaussian_width_mm", "gamma_over_omega",
                         "group_delay_us", "settling_time_us", "N95_cycles",
                         "late_time_A0_purity"])
        writer.writerows(sweep)
    np.savez_compressed(OUT/"calculation_data.npz", K=K, k=k, AW=AW, SW=SW,
                        Aomega=Aomega, Somega=Somega, Avg=Avg, Svg=Svg,
                        Acurvature=Acurv, Scurvature=Scurv, **transient)

    # Figures.
    colors = plt.cm.tab10.colors
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    for b in range(3):
        ay = break_plot_jumps(Aomega[:, b]/(2*np.pi)*MAT.d/1000)
        sy = break_plot_jumps(Somega[:, b]/(2*np.pi)*MAT.d/1000)
        ax.plot(k*MAT.d, ay, "--", color=colors[b], label=f"A{b}")
        ax.plot(k*MAT.d, sy, "-", color=colors[b], label=f"S{b}")
    ax.set(xlabel=r"$kd$", ylabel=r"$fd$ (MHz mm)", title="Rayleigh–Lamb dispersion: 1 mm aluminum plate")
    ax.grid(True, alpha=.3); ax.legend(ncol=2); fig.tight_layout(); fig.savefig(OUT/"dispersion.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.plot(Aomega[:, 0]/(2*np.pi)*MAT.d/1000, Avg[:, 0], label="A0")
    ax.plot(Somega[:, 0]/(2*np.pi)*MAT.d/1000, Svg[:, 0], label="S0")
    ax.axhline(0, color="k", lw=.8, ls=":"); ax.grid(True, alpha=.3); ax.legend()
    ax.set(xlabel=r"$fd$ (MHz mm)", ylabel="group velocity (m/s)", title="Fundamental-mode group velocity")
    fig.tight_layout(); fig.savefig(OUT/"group_velocity.png", dpi=180); plt.close(fig)

    cycles = transient["t"]/transient["period"]
    fig, axs = plt.subplots(2, 1, figsize=(8.2, 6.5), sharex=True)
    axs[0].plot(cycles, abs(transient["ua"])/max(abs(transient["ssa"])))
    axs[0].axhline(1, color="k", lw=.8, ls=":"); axs[0].set_ylabel(r"$|u_{A0}|/|u_{ss}|$")
    axs[1].plot(cycles, transient["purity"]); axs[1].set(xlabel="time (drive cycles)", ylabel="A0 modal purity", ylim=(0,1))
    for ax in axs: ax.grid(True, alpha=.3)
    fig.suptitle("Normalized causal buildup (scalar modal coupling)"); fig.tight_layout(); fig.savefig(OUT/"causal_buildup.png", dpi=180); plt.close(fig)

    with (OUT/"validation_summary.txt").open("w", encoding="utf-8") as f:
        f.write("CALCULATION AND VALIDATION SUMMARY\n")
        f.write(f"rho={MAT.rho:g} kg/m^3, cL={MAT.cL:g} m/s, cT={MAT.cT:g} m/s, d={MAT.d:g} m\n")
        f.write(f"Grid: {K.size} K values; 7000 frequency samples; Wmax=12\n\n")
        f.write(f"Maximum normalized characteristic residual: A={np.nanmax(Ares):.3e}, S={np.nanmax(Sres):.3e}\n")
        f.write(f"Low-frequency median relative error: A0 thin-plate={100*a_err:.3f}%, S0 plate-speed={100*s_err:.3f}%\n")
        f.write(f"Expected low-frequency S0 plate speed={cplate:.3f} m/s; cT/cL={ratio:.6f}\n\n")
        f.write("Convergence and analytic checks:\n")
        f.write(f"  max |W(3500 scan)-W(7000 scan)|={dispersion_convergence:.3e}\n")
        f.write(f"  settling-time change on halving dt={time_convergence_cycles:.3f} drive cycles\n")
        f.write(f"  oscillator initial displacement residual={initial_displacement_error:.3e}\n")
        f.write(f"  oscillator initial velocity residual={initial_velocity_error:.3e}\n\n")
        f.write("Detected higher-mode local minima (quadratically refined):\n")
        for mode, kval, fd, freq in zgvs:
            f.write(f"  {mode}: K={kval:.5f}, fd={fd:.6f} MHz mm, f={freq:.6f} MHz\n")
        f.write("\nNormalized causal study (not an absolute thermoelastic amplitude):\n")
        f.write(f"  drive={transient['frequency']/1e6:.6f} MHz; x={transient['x']*1e3:.3f} mm; width={transient['width']*1e3:.3f} mm; gamma/Omega=0.01\n")
        f.write(f"  group delay x/|vg|={transient['arrival']*1e6:.6f} us\n")
        f.write(f"  5% settling time={transient['tss']*1e6:.6f} us ({transient['tss']/transient['period']:.1f} cycles)\n")
        f.write(f"  N95 proxy={transient['N95']:.0f} cycles\n\n")
        finite_n = np.array([row[5] for row in sweep], dtype=float)
        f.write(f"Parameter sweep (27 cases): finite N95 range={np.nanmin(finite_n):.0f}--{np.nanmax(finite_n):.0f} cycles; see parameter_sweep.csv\n\n")
        f.write("Dispersion uses the full isotropic Rayleigh–Lamb equations. The causal study uses unit scalar modal overlap and constant damping; it validates the framework, not absolute experimental amplitude.\n")
    print((OUT/"validation_summary.txt").read_text())


if __name__ == "__main__":
    main()
