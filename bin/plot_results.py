#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gera painel PNG com as principais series temporais de uma pasta de analise.
Segue a mesma convencao de CLI usada em MD-gromacs/bin/plot_results.py.

Uso:
    mamba activate md-gromacs   # matplotlib/numpy
    python bin/plot_results.py --analise-dir results/2I9T-daidzeina/analise \\
        --titulo "2I9T (NF-kB) + Daidzeina" --window-ns 5
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


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


def moving_average(x, y, window_ns):
    if not x or window_ns <= 0 or len(x) < 3:
        return x, y
    dt = (x[-1] - x[0]) / max(len(x) - 1, 1)
    if dt <= 0:
        return x, y
    win = max(1, int(round(window_ns / dt)))
    if win <= 1:
        return x, y
    smoothed = []
    for i in range(len(y)):
        lo = max(0, i - win // 2)
        hi = min(len(y), i + win // 2 + 1)
        smoothed.append(sum(y[lo:hi]) / (hi - lo))
    return x, smoothed


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
    ap = argparse.ArgumentParser()
    ap.add_argument("--analise-dir", required=True, help="pasta com os .xvg do gmx")
    ap.add_argument("--titulo", default="Dinamica molecular — resumo")
    ap.add_argument("--window-ns", type=float, default=5.0, help="janela da media movel (ns), 0 desativa")
    ap.add_argument("--out", default=None, help="PNG de saida (default: <analise-dir>/painel_resumo.png)")
    args = ap.parse_args()

    analise_dir = Path(args.analise_dir)
    out_png = Path(args.out) if args.out else analise_dir / "painel_resumo.png"

    fig, axes = plt.subplots(4, 2, figsize=(12, 14))
    axes = axes.flatten()

    any_data = False
    for ax, (fname, title, xlabel, ylabel) in zip(axes, PANELS):
        x, y = read_xvg(analise_dir / fname)
        if x:
            any_data = True
            ax.plot(x, y, linewidth=0.5, alpha=0.4, color="tab:blue")
            xs, ys = moving_average(x, y, args.window_ns)
            ax.plot(xs, ys, linewidth=1.3, color="tab:blue")
        else:
            ax.text(0.5, 0.5, "sem dados\n(rodar analyze.sh / ANALYSES)",
                     ha="center", va="center", transform=ax.transAxes, fontsize=9)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(labelsize=7)

    fig.suptitle(args.titulo, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    print(f"[OK] Painel salvo em {out_png}")
    if not any_data:
        print(f"[AVISO] Nenhum .xvg encontrado em {analise_dir} ainda.")


if __name__ == "__main__":
    main()
