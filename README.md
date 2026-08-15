# Causal buildup and laser excitation of Lamb waves

Thesis-style LaTeX manuscript and reproducible numerical benchmark for the causal
buildup and selective thermoelastic laser excitation of Lamb waves.

The repository distinguishes three scientific levels:

- **Verified:** Rayleigh--Lamb dispersion, analytic low-frequency limits,
  determinant residuals, scan convergence, and oscillator initial conditions.
- **Preliminary:** preserved scalar-unit-overlap regression calculations. These
  are not physical modal purity or absolute thermoelastic predictions.
- **Pending:** normalized eigenfields, physical thermoelastic overlaps,
  energy/group-velocity checks, calibrated damping/detector models, and causal
  experimental predictions.

See `SCIENTIFIC_STATUS.md` before interpreting any numerical output.

## Build

With TeX Live and `latexmk` installed:

```powershell
latexmk -pdf main.tex
```

The result is `main.pdf`. The source is organized under `frontmatter/`,
`chapters/`, and `appendices/`.

## Set up and test the calculations

Create a virtual environment, install the pinned packages, and run:

```powershell
python -m venv .venv
.\.venv\Scripts\python -m pip install -r requirements.txt
.\.venv\Scripts\python -m pytest -q
.\.venv\Scripts\python -m lambwaves.run
```

Verified outputs are written to `results/verified/`, including metadata in JSON
and machine-readable CSV data. The benchmark uses explicitly recorded
representative aluminum inputs; they are not fitted specimen properties.

## Legacy preliminary calculation

`calculations/run_all.py` and `calculations/run_all.m` preserve the original
unit-overlap study. Their outputs live under `results/preliminary/`; do not cite
their settling time, `N95`, or mode ratio as experimental predictions.
