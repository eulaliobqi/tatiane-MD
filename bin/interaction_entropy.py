#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Estimativa de energia livre de ligação via Interaction Entropy (Duan, Liu,
Zhang, JACS 2016) sobre a energia de interação Coulomb-SR + LJ-SR extraída
de um rerun do GROMACS com energygrps Receptor/Ligante.

APROXIMAÇÃO DELIBERADA: energia de interação em vácuo (sem termo de
solvatação implícita PB/GB) + correção entrópica IE. Não é equivalente a um
ΔG de MM-GBSA/MM-PBSA rigoroso -- serve para comparação relativa entre
isoformas da mesma série calculadas pelo mesmo método, não como valor
absoluto de afinidade.

Uso:
  interaction_entropy.py --xvg interaction_energy.xvg --temperature 300 \
      --out-dir . [--titulo T]
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R_KJ = 0.0083144621  # kJ/mol/K


def logsumexp(x):
    """ln(sum(exp(x))), estavel numericamente (sem dependencia de scipy)."""
    x_max = np.max(x)
    return x_max + np.log(np.sum(np.exp(x - x_max)))


def load_xvg(path):
    data = []
    for ln in open(path):
        if ln.startswith(('#', '@')):
            continue
        parts = ln.split()
        if len(parts) < 3:
            continue
        try:
            data.append([float(x) for x in parts])
        except ValueError:
            pass
    return np.array(data)


def rolling_stats(values, window):
    n = len(values)
    w = max(1, min(window, n))
    v = values.astype(float)
    pad = w // 2
    vp = np.pad(v, (pad, w - pad - 1), mode='edge')
    cs = np.cumsum(np.insert(vp, 0, 0.0))
    mean = (cs[w:] - cs[:-w]) / w
    cs2 = np.cumsum(np.insert(vp ** 2, 0, 0.0))
    var = (cs2[w:] - cs2[:-w]) / w - mean ** 2
    std = np.sqrt(np.maximum(var, 0.0))
    return mean[:n], std[:n]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xvg", required=True,
                     help="xvg do gmx energy com colunas: tempo, Coul-SR, LJ-SR")
    ap.add_argument("--temperature", type=float, default=300.0)
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--titulo", default="Interaction Energy")
    ap.add_argument("--window", type=int, default=10,
                     help="Janela da média móvel, em frames (default 10)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    report_path = os.path.join(args.out_dir, "free_energy_estimate.txt")
    png_path = os.path.join(args.out_dir, "interaction_energy.png")

    arr = load_xvg(args.xvg)
    if arr.size == 0 or arr.shape[1] < 3:
        with open(report_path, "w") as fh:
            fh.write("=== Estimativa de Energia Livre (Interaction Entropy) ===\n")
            fh.write(f"ERRO: {args.xvg} vazio ou sem colunas suficientes "
                     "(esperado: tempo, Coul-SR, LJ-SR)\n")
        print(f"ERRO: xvg invalido: {args.xvg}", file=sys.stderr)
        return

    time_ps = arr[:, 0]
    coul_sr = arr[:, 1]
    lj_sr = arr[:, 2]
    dE_int = coul_sr + lj_sr
    n = len(dE_int)

    T = args.temperature
    kT = R_KJ * T
    mean_dE = dE_int.mean()
    sd_dE = dE_int.std()

    x = (dE_int - mean_dE) / kT
    ln_avg_exp = logsumexp(x) - np.log(n)
    minus_T_dS = kT * ln_avg_exp

    dG_estimate = mean_dE + minus_T_dS
    Kd_estimate_M = np.exp(dG_estimate / kT)  # mesma kT (R*T), Kd em unidade "molar relativa"

    unstable_warning = ""
    if sd_dE / kT > 10:
        unstable_warning = (
            "\nAVISO: desvio padrao da energia de interacao e' > 10 kT "
            f"({sd_dE/kT:.1f} kT) -- o metodo Interaction Entropy e' conhecido "
            "por ficar numericamente instavel/superestimar a entropia nesse "
            "regime (ver ressalva metodologica em docs/frontier-tech-roadmap.md, "
            "secao 3). Tratar o -TdS como qualitativo, nao quantitativo, neste caso.\n"
        )

    with open(report_path, "w") as fh:
        fh.write("=== Estimativa de Energia Livre de Ligacao (Interaction Entropy) ===\n\n")
        fh.write("*** APROXIMACAO -- LER ANTES DE USAR ***\n")
        fh.write("Energia de interacao Coulomb-SR + LJ-SR em vacuo (rerun do GROMACS\n")
        fh.write("com energygrps Receptor/Ligante) + correcao entropica via Interaction\n")
        fh.write("Entropy (Duan, Liu, Zhang, JACS 2016). NAO inclui termo de solvatacao\n")
        fh.write("implicita (PB/GB) -- nao e equivalente a um resultado de MM-GBSA/\n")
        fh.write("MM-PBSA rigoroso. Valido para comparacao RELATIVA entre isoformas da\n")
        fh.write("mesma serie calculadas pelo mesmo metodo, nao como valor absoluto de\n")
        fh.write("afinidade de ligacao.\n\n")
        fh.write("ESPERADO (nao e' bug): sem a blindagem do solvente, a eletrostatica em\n")
        fh.write("vacuo entre duas cadeias de proteina/peptideo tende a dar Delta_G/Kd em\n")
        fh.write("magnitude fisicamente exagerada (centenas de kJ/mol, Kd muito abaixo de\n")
        fh.write("qualquer faixa experimental) -- e exatamente por isso que MM-GBSA/MM-PBSA\n")
        fh.write("de verdade somam um termo de solvatacao implicita para cancelar a maior\n")
        fh.write("parte dessa eletrostatica. Sem esse termo aqui, so o RANKING relativo\n")
        fh.write("entre sistemas desta serie (calculados pelo mesmo metodo) e informativo\n")
        fh.write("-- o valor absoluto de Delta_G/Kd abaixo nao deve ser citado isoladamente.\n\n")
        fh.write(f"N. frames                    : {n}\n")
        fh.write(f"Temperatura                  : {T:.1f} K\n\n")
        fh.write(f"<Coul-SR>                    : {coul_sr.mean():.2f} +/- {coul_sr.std():.2f} kJ/mol\n")
        fh.write(f"<LJ-SR>                      : {lj_sr.mean():.2f} +/- {lj_sr.std():.2f} kJ/mol\n")
        fh.write(f"<Delta_E_MM> (Coul+LJ)        : {mean_dE:.2f} +/- {sd_dE:.2f} kJ/mol\n")
        fh.write(f"-T*Delta_S (Interaction Entropy) : {minus_T_dS:.2f} kJ/mol\n")
        fh.write(f"Delta_G estimado              : {dG_estimate:.2f} kJ/mol "
                 f"({dG_estimate/4.184:.2f} kcal/mol)\n")
        fh.write(f"Kd estimado                   : {Kd_estimate_M:.3e} (unidade molar relativa, ver nota)\n")
        fh.write(unstable_warning)
        fh.write("\nNota sobre Kd: derivado algebricamente de Delta_G = R*T*ln(Kd) usando\n")
        fh.write("o mesmo R*T da correcao IE acima -- herda toda a aproximacao do\n")
        fh.write("Delta_G (sem solvatacao explicita), nao deve ser reportado como Kd\n")
        fh.write("experimental ou como substituto de um ensaio de afinidade real.\n")
    print(f"[OK] {report_path}", file=sys.stderr)

    # ── interaction_energy.png ─────────────────────────────────────────────
    time_ns = time_ps / 1000.0
    rm, rs = rolling_stats(dE_int, args.window)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(time_ns, dE_int, lw=0.5, color="darkred", alpha=0.30)
    ax.fill_between(time_ns, rm - rs, rm + rs, alpha=0.18, color="darkred")
    ax.plot(time_ns, rm, lw=1.6, color="darkred")
    ax.axhline(mean_dE, ls='--', color='black', lw=1.0,
               label=f"Mean: {mean_dE:.1f} kJ/mol")
    ax.set_xlabel("Frame (ns equivalente)")
    ax.set_ylabel("Coul-SR + LJ-SR (kJ/mol)")
    ax.set_title(f"{args.titulo}\nΔG estimado (IE): {dG_estimate:.1f} kJ/mol "
                 f"({dG_estimate/4.184:.1f} kcal/mol)")
    ax.legend()
    ax.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[OK] {png_path}", file=sys.stderr)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
