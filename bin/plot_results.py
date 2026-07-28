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
    ("sasa_ligante.xvg", "SASA ligante", "Tempo (ns)", "SASA (nm2)"),
]

# Uma cor distinta por painel (ciclo), pra diferenciar as figuras separadas
# umas das outras e no painel combinado -- antes tudo saia tab:blue.
PALETTE = [
    "tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple",
    "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan",
]


def discover_residue_panels(analise_dir):
    """Um painel por dist_<Residuo>.xvg encontrado (numero de residuos-chave
    varia por sistema — ver meta.key_residues/ANALYSES_RESIDUES)."""
    panels = []
    for f in sorted(analise_dir.glob("dist_*.xvg")):
        label = f.stem.replace("dist_", "")
        panels.append((f.name, f"Dist. minima Ligante-{label}", "Tempo (ns)", "Dist. (nm)"))
    return panels


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analise-dir", required=True, help="pasta com os .xvg do gmx")
    ap.add_argument("--titulo", default="Dinamica molecular — resumo")
    ap.add_argument("--window-ns", type=float, default=5.0, help="janela da media movel (ns), 0 desativa")
    ap.add_argument("--out", default=None, help="PNG de saida (default: <analise-dir>/painel_resumo.png)")
    ap.add_argument("--figures-dir", default=None,
                     help="pasta pras figuras individuais (default: <out>/figuras)")
    args = ap.parse_args()

    analise_dir = Path(args.analise_dir)
    out_png = Path(args.out) if args.out else analise_dir / "painel_resumo.png"
    figures_dir = Path(args.figures_dir) if args.figures_dir else out_png.parent / "figuras"
    figures_dir.mkdir(parents=True, exist_ok=True)

    panels = PANELS + discover_residue_panels(analise_dir)
    ncols = 2
    nrows = -(-len(panels) // ncols)  # ceil
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 3.5 * nrows))
    axes = axes.flatten()
    for ax in axes[len(panels):]:
        ax.axis("off")

    any_data = False
    for i, (ax, (fname, title, xlabel, ylabel)) in enumerate(zip(axes, panels)):
        color = PALETTE[i % len(PALETTE)]
        x, y = read_xvg(analise_dir / fname)
        if x:
            any_data = True
            ax.plot(x, y, linewidth=0.5, alpha=0.4, color=color)
            xs, ys = moving_average(x, y, args.window_ns)
            ax.plot(xs, ys, linewidth=1.3, color=color)
        else:
            ax.text(0.5, 0.5, "sem dados\n(rodar analyze.sh / ANALYSES)",
                     ha="center", va="center", transform=ax.transAxes, fontsize=9)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(labelsize=7)

        # Figura individual do mesmo painel, mesma cor, arquivo proprio —
        # facilita usar uma metrica isolada em slides/artigo sem recortar
        # do painel combinado.
        slug = Path(fname).stem
        ind_fig, ind_ax = plt.subplots(figsize=(6, 4))
        if x:
            ind_ax.plot(x, y, linewidth=0.6, alpha=0.4, color=color)
            xs, ys = moving_average(x, y, args.window_ns)
            ind_ax.plot(xs, ys, linewidth=1.8, color=color)
        else:
            ind_ax.text(0.5, 0.5, "sem dados\n(rodar analyze.sh / ANALYSES)",
                         ha="center", va="center", transform=ind_ax.transAxes, fontsize=10)
        ind_ax.set_title(f"{args.titulo}\n{title}", fontsize=10)
        ind_ax.set_xlabel(xlabel, fontsize=9)
        ind_ax.set_ylabel(ylabel, fontsize=9)
        ind_fig.tight_layout()
        ind_fig.savefig(figures_dir / f"{slug}.png", dpi=150)
        plt.close(ind_fig)

    fig.suptitle(args.titulo, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    print(f"[OK] Painel salvo em {out_png}")
    print(f"[OK] Figuras individuais salvas em {figures_dir}/")
    if not any_data:
        print(f"[AVISO] Nenhum .xvg encontrado em {analise_dir} ainda.")


if __name__ == "__main__":
    main()
