#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fingerprint de interacao (ProLIF) receptor-ligante ao longo da trajetoria
(por fase — bound/relocated, ver PHASE_SPLIT).

Adaptado de Milena-MD/bin/prolif_fingerprint.py: ligante aqui e' molecula
pequena (daidzeina, resname UNL), nao peptideo — selecao via `resname UNL`
em vez de `moltype Protein_chain_B` (ver contact_map.py para o mesmo
raciocinio). ProLIF foi desenhado primariamente para esse caso de uso
(ligante pequeno), entao essa e' na verdade uma aplicacao mais natural do
que a original em Milena-MD.

Modulo de maior risco de dependencia do conjunto (prolif+rdkit+MDAnalysis)
-- escreve log de diagnostico ANTES de tentar rodar o ProLIF e cai num
fallback com CSV vazio em caso de erro, em vez de estourar traceback cru
(mesmo padrao defensivo ja usado em MMGBSA).

Uso:
  prolif_fingerprint.py --complexo-pdb complexo.pdb --tpr md.tpr \
      --xtc bound.xtc --out-dir . [--stride 1]
"""
import argparse
import os
import re
import sys
import traceback


def read_chain_residues_ordered(complexo_pdb, chain):
    """Lista (resid, resname) da cadeia, na ordem do arquivo (numeracao real do PDB)."""
    seen = set()
    residues = []
    with open(complexo_pdb) as fh:
        for line in fh:
            if line.startswith(('ATOM', 'HETATM')) and line[21:22].strip() == chain:
                resid = int(line[22:26])
                if resid not in seen:
                    seen.add(resid)
                    residues.append((resid, line[17:20].strip()))
    if not residues:
        raise SystemExit(f"ERRO: cadeia {chain} nao encontrada em {complexo_pdb}")
    return residues


def build_resid_map(atomgroup, pdb_residues, label):
    internal_resids = sorted(set(int(r) for r in atomgroup.resids))
    if len(internal_resids) != len(pdb_residues):
        raise SystemExit(
            f"ERRO: contagem de residuos do {label} nao bate "
            f"(MDAnalysis={len(internal_resids)}, PDB={len(pdb_residues)})"
        )
    return dict(zip(internal_resids, pdb_residues))


_LABEL_RE = re.compile(r"^([A-Za-z]+)(\d+)(\..*)?$")


def remap_label(label, resid_map):
    """'UNL1.A' (numeracao interna do GROMACS) -> 'UNL301.A' (numeracao PDB real)."""
    m = _LABEL_RE.match(str(label))
    if not m:
        return label
    resname, resid_str, suffix = m.groups()
    resid = int(resid_str)
    if resid not in resid_map:
        return label
    pdb_resid, pdb_resname = resid_map[resid]
    return f"{pdb_resname}{pdb_resid}{suffix or ''}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complexo-pdb", required=True)
    ap.add_argument("--tpr", required=True)
    ap.add_argument("--xtc", required=True)
    ap.add_argument("--stride", type=int, default=1,
                     help="Usa 1 a cada N frames (default 1 = todos)")
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    log_path = os.path.join(args.out_dir, "prolif_log.txt")
    csv_path = os.path.join(args.out_dir, "prolif_fingerprint.csv")
    png_path = os.path.join(args.out_dir, "prolif_heatmap.png")

    def log(msg):
        print(msg, file=sys.stderr)
        with open(log_path, "a") as fh:
            fh.write(msg + "\n")

    log("=== PROLIF_FINGERPRINT ===")

    try:
        import MDAnalysis as mda
        import prolif as plf
        import pandas as pd
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        rec_pdb = read_chain_residues_ordered(args.complexo_pdb, "A")
        lig_pdb = read_chain_residues_ordered(args.complexo_pdb, "B")

        u = mda.Universe(args.tpr, args.xtc)
        rec_sel = u.select_atoms("protein")
        lig_sel = u.select_atoms("resname UNL")
        log(f"[OK] Selecoes: ligante={len(lig_sel)} atomos, receptor={len(rec_sel)} atomos")

        rec_map = build_resid_map(rec_sel, rec_pdb, "receptor")
        lig_map = build_resid_map(lig_sel, lig_pdb, "ligante")

        fp = plf.Fingerprint([
            "HBAcceptor", "HBDonor", "Hydrophobic", "PiStacking",
            "Anionic", "Cationic", "VdWContact",
        ])
        traj = u.trajectory[::args.stride]
        log(f"[INFO] Rodando ProLIF sobre {len(traj)} frames (stride={args.stride})...")
        fp.run(traj, lig_sel, rec_sel)

        df = fp.to_dataframe()
        df.columns = df.columns.set_levels(
            [remap_label(l, lig_map) for l in df.columns.levels[0]], level=0
        ).set_levels(
            [remap_label(l, rec_map) for l in df.columns.levels[1]], level=1
        )
        df.to_csv(csv_path)
        log(f"[OK] {csv_path} ({df.shape[0]} frames x {df.shape[1]} interacoes)")

        freq = df.mean(axis=0)
        freq_df = freq.reset_index()
        freq_df.columns = ["ligand_residue", "protein_residue", "interaction", "frequency"]
        pivot = freq_df.pivot_table(index="protein_residue", columns="interaction",
                                     values="frequency", aggfunc="sum", fill_value=0.0)

        if pivot.empty:
            raise ValueError("Nenhuma interacao detectada pelo ProLIF (pivot vazio)")

        fig, ax = plt.subplots(figsize=(8, max(4, len(pivot) * 0.25)))
        im = ax.imshow(pivot.values, aspect='auto', cmap='magma', vmin=0, vmax=1)
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels(pivot.columns, rotation=45, ha='right', fontsize=8)
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(pivot.index, fontsize=7)
        ax.set_title("ProLIF: Frequencia de Interacao por Residuo do Receptor\n(numeracao PDB)")
        fig.colorbar(im, ax=ax, label="Frequencia")
        plt.tight_layout()
        plt.savefig(png_path, dpi=150, bbox_inches='tight')
        plt.close()
        log(f"[OK] {png_path}")
        log("[OK] ProLIF concluido com sucesso")

    except Exception:
        log("ERRO: ProLIF falhou -- ver traceback abaixo")
        log(traceback.format_exc())
        if not os.path.exists(csv_path):
            with open(csv_path, "w") as fh:
                fh.write("ligand_residue,protein_residue,interaction,frequency\n")
        log("[FALLBACK] CSV vazio gravado -- ProLIF nao disponivel/falhou nesta rodada")


if __name__ == "__main__":
    main()
