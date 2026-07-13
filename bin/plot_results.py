#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gera painel PNG com as principais series temporais de results/2I9T-daidzeina/10_analysis/.
Uso:
    mamba activate md-gromacs   # matplotlib/numpy
    python bin/plot_results.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS_DIR = ROOT / "results" / "2I9T-daidzeina" / "10_analysis"
OUT_PNG = ANALYSIS_DIR / "painel_resumo.png"


def read_xvg(path):
    if not path.exists():
        return [], []
    x, y = [], []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(("@", "#")):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            x.append(float(parts[0]))
            y.append(float(parts[1]))
        except ValueError:
            continue
    return x, y


PANELS = [
    ("rmsd_backbone.xvg", "RMSD backbone receptor", "Tempo (ns)", "RMSD (nm)"),
    ("rmsd_ligante.xvg", "RMSD ligante (UNL)", "Tempo (ns)", "RMSD (nm)"),
    ("gyrate.xvg", "Raio de giro (receptor)", "Tempo (ns)", "Rg (nm)"),
    ("numcont_receptor_ligante.xvg", "Contatos receptor-ligante (<0,4nm)", "Tempo (ns)", "N contatos"),
    ("hbond.xvg", "Pontes de H receptor-ligante", "Tempo (ns)", "N pontes H"),
    ("dist_arg30.xvg", "Dist. minima Ligante-Arg30 (docking: 4,7-4,8 A)", "Tempo (ns)", "Dist. (nm)"),
    ("dist_glu279.xvg", "Dist. minima Ligante-Glu279 (docking: 1,9 A)", "Tempo (ns)", "Dist. (nm)"),
    ("sasa_ligante.xvg", "SASA ligante", "Tempo (ns)", "SASA (nm2)"),
]


def main():
    fig, axes = plt.subplots(4, 2, figsize=(12, 14))
    axes = axes.flatten()

    any_data = False
    for ax, (fname, title, xlabel, ylabel) in zip(axes, PANELS):
        x, y = read_xvg(ANALYSIS_DIR / fname)
        if x:
            any_data = True
            ax.plot(x, y, linewidth=0.8)
        else:
            ax.text(0.5, 0.5, "sem dados\n(rodar analyze.sh)",
                     ha="center", va="center", transform=ax.transAxes, fontsize=9)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(labelsize=7)

    fig.suptitle("2I9T (NF-kB) + Daidzeina — resumo da dinamica molecular", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150)
    print(f"[OK] Painel salvo em {OUT_PNG}")
    if not any_data:
        print("[AVISO] Nenhum .xvg encontrado ainda — rode bin/run_md.sh + bin/analyze.sh primeiro.")


if __name__ == "__main__":
    main()
