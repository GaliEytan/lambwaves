# Causal buildup of laser-excited Lamb waves

LaTeX manuscript and, later, numerical code for studying the transient formation
and selective excitation of Lamb waves by thermoelastic laser sources.

## Build

With a LaTeX distribution installed:

```powershell
latexmk -pdf main.tex
```

## Reproduce the calculations

Create a virtual environment, install the pinned packages, and run:

```powershell
python -m venv .venv
.\.venv\Scripts\python -m pip install -r requirements.txt
.\.venv\Scripts\python calculations\run_all.py
```

The Python script solves the isotropic Rayleigh--Lamb equations for a representative
1 mm aluminum plate, evaluates dispersion derivatives and ZGV candidates, and
runs a normalized causal modal-integral study. It writes figures, a CSV table,
the complete workspace, and a validation summary to `results/`.

The dispersion calculation is dimensional and uses the full characteristic
equations. The initial transient calculation deliberately uses scalar unit modal
overlap; absolute thermoelastic amplitudes require an experiment-specific optical
absorption profile and normalized Lamb eigenfields.

`calculations/run_all.m` is an equivalent MATLAB reference implementation, but
the reported repository results are generated and validated with the Python code.
