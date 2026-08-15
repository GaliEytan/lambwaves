from __future__ import annotations

import csv, json, os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT/".matplotlib"))
import matplotlib.pyplot as plt
import numpy as np

from .dispersion import dimensionalize
from .material import Material
from .validation import validation_report


def main():
    out = ROOT/"results"/"verified"; out.mkdir(parents=True, exist_ok=True)
    mat = Material(); K = np.linspace(.04, 8, 260)
    report = validation_report(mat, K)
    Ao, k, Avg = dimensionalize(report["A"], K, mat)
    So, _, Svg = dimensionalize(report["S"], K, mat)
    metadata = {"rho_kg_m3":mat.rho,"cL_m_s":mat.cL,"cT_m_s":mat.cT,
                "thickness_m":mat.thickness,"K_min":float(K.min()),
                "K_max":float(K.max()),"K_count":len(K),"root_scan":7000}
    summary = {key:value for key,value in report.items() if not isinstance(value,np.ndarray)}
    (out/"validation.json").write_text(json.dumps({"metadata":metadata,"validation":summary},indent=2))
    with (out/"dispersion.csv").open("w",newline="") as f:
        w=csv.writer(f); w.writerow(["K","A0_fd_MHz_mm","S0_fd_MHz_mm","A0_vg_m_s","S0_vg_m_s"])
        w.writerows(zip(K,Ao[:,0]/(2*np.pi)*mat.thickness/1000,
                         So[:,0]/(2*np.pi)*mat.thickness/1000,Avg[:,0],Svg[:,0]))
    def clean(y, jump=.35):
        y=y.copy(); j=np.flatnonzero(abs(np.diff(y))>jump); y[j]=np.nan; y[j+1]=np.nan
        return y
    fig,ax=plt.subplots(figsize=(8.2,5.2))
    colors=plt.cm.tab10.colors
    ax.plot(k*mat.thickness,Ao[:,0]/(2*np.pi)*mat.thickness/1000,"--",color=colors[0],label="A0")
    ax.plot(k*mat.thickness,So[:,0]/(2*np.pi)*mat.thickness/1000,"-",color=colors[0],label="S0")
    for b in (1,2):
        ay=clean(Ao[:,b]/(2*np.pi)*mat.thickness/1000)
        sy=clean(So[:,b]/(2*np.pi)*mat.thickness/1000)
        ay[K>3.2]=np.nan; sy[K>3.2]=np.nan
        ax.plot(k*mat.thickness,ay,"--",color=colors[b],label=f"A{b} segments")
        ax.plot(k*mat.thickness,sy,"-",color=colors[b],label=f"S{b} segments")
    ax.set(xlabel=r"$kd$",ylabel=r"$fd$ (MHz mm)",title="Rayleigh–Lamb dispersion benchmark")
    ax.grid(alpha=.3);ax.legend(ncol=2);fig.tight_layout();fig.savefig(out/"dispersion.png",dpi=180);plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__ == "__main__": main()
