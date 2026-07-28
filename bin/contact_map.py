#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mapa de contatos receptor-ligante por residuo, a partir da trajetoria
GROMACS ja processada (por fase — bound/relocated, ver PHASE_SPLIT).

Adaptado de Milena-MD/bin/contact_map.py: la o ligante era um peptideo
(segunda cadeia proteica, selecionado via `moltype Protein_chain_B`); aqui
o ligante e' molecula pequena (daidzeina, resname UNL, CGenFF) — nao tem
moltype "Protein_chain_B", entao a selecao usa `resname UNL` diretamente
(mesma identificacao ja usada em modules/local/analyses, testada em
producao). Receptor selecionado via `protein` (MDAnalysis reconhece HSD/
HSE/HSP do CHARMM como residuos proteicos padrao desde mdanalysis>=2.0).

O receptor ainda usa a leitura de cadeia A do complexo.pdb + remapeamento
por posicao para numeracao PDB real, mesma tecnica do Milena-MD (pdb2gmx
pode alterar numeracao interna). Como o ligante e' um unico residuo (UNL),
o "mapa" colapsa numa unica coluna — ainda assim util: da a frequencia de
contato POR RESIDUO DO RECEPTOR ao longo da fase, revelando pra onde a
daidzeina foi apos a transicao de ~65-78ns.

Uso:
  contact_map.py --complexo-pdb complexo.pdb --tpr md.tpr --xtc bound.xtc \
      --out-dir . [--cutoff 4.0]
"""
import argparse
import os
import sys

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


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
    """Mapa {resid_interno_MDAnalysis: (resid_pdb, resname_pdb)} por posicao."""
    internal_resids = sorted(set(int(r) for r in atomgroup.resids))
    if len(internal_resids) != len(pdb_residues):
        raise SystemExit(
            f"ERRO: contagem de residuos do {label} nao bate "
            f"(MDAnalysis={len(internal_resids)}, PDB={len(pdb_residues)})"
        )
    return dict(zip(internal_resids, pdb_residues))


def residue_onehot(atomgroup, internal_resids):
    """Retorna matriz esparsa n_atoms x n_res one-hot, na ordem de internal_resids."""
    idx = {r: i for i, r in enumerate(internal_resids)}
    rows = np.arange(len(atomgroup))
    cols = np.array([idx[int(r)] for r in atomgroup.resids])
    data = np.ones(len(atomgroup), dtype=np.int32)
    return sparse.csr_matrix((data, (rows, cols)), shape=(len(atomgroup), len(internal_resids)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complexo-pdb", required=True)
    ap.add_argument("--tpr", required=True)
    ap.add_argument("--xtc", required=True)
    ap.add_argument("--cutoff", type=float, default=4.0,
                     help="Angstrom (default 4.0 = 0.4 nm, mesmo cutoff de ANALYSES)")
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    rec_pdb = read_chain_residues_ordered(args.complexo_pdb, "A")
    lig_pdb = read_chain_residues_ordered(args.complexo_pdb, "B")
    print(f"[CONTACT_MAP] Receptor (PDB): {len(rec_pdb)} residuos, "
          f"Ligante (PDB): {len(lig_pdb)} residuos", file=sys.stderr)

    u = mda.Universe(args.tpr, args.xtc)
    rec = u.select_atoms("protein and not name H*")
    lig = u.select_atoms("resname UNL and not name H*")

    if len(lig) == 0 or len(rec) == 0:
        raise SystemExit(f"ERRO: selecao vazia (lig={len(lig)} atomos, rec={len(rec)} atomos)")

    rec_map = build_resid_map(rec, rec_pdb, "receptor")
    lig_map = build_resid_map(lig, lig_pdb, "ligante")
    rec_internal = sorted(rec_map)
    lig_internal = sorted(lig_map)

    rec_onehot = residue_onehot(rec, rec_internal)
    lig_onehot = residue_onehot(lig, lig_internal)

    n_res_rec, n_res_lig = len(rec_internal), len(lig_internal)
    contact_counts = np.zeros((n_res_rec, n_res_lig), dtype=np.int64)

    n_frames = len(u.trajectory)
    for ts in u.trajectory:
        d = distances.distance_array(rec.positions, lig.positions, box=ts.dimensions)
        atom_contact = sparse.csr_matrix(d < args.cutoff)
        res_contact = (rec_onehot.T @ atom_contact @ lig_onehot) > 0
        contact_counts += res_contact.toarray()

    freq = contact_counts / max(n_frames, 1)

    os.makedirs(args.out_dir, exist_ok=True)

    rec_labels = [f"{rec_map[r][1]}{rec_map[r][0]}" for r in rec_internal]
    lig_labels = [f"{lig_map[r][1]}{lig_map[r][0]}" for r in lig_internal]

    # ── contact_map.csv ────────────────────────────────────────────────────
    csv_path = os.path.join(args.out_dir, "contact_map.csv")
    with open(csv_path, "w") as fh:
        fh.write("residue," + ",".join(lig_labels) + "\n")
        for i, rl in enumerate(rec_labels):
            fh.write(rl + "," + ",".join(f"{v:.4f}" for v in freq[i]) + "\n")
    print(f"[OK] {csv_path}", file=sys.stderr)

    # ── contact_map.png ────────────────────────────────────────────────────
    fig_h = max(6, n_res_rec * 0.12)
    fig_w = max(8, n_res_lig * 0.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(freq, aspect='auto', cmap='viridis', vmin=0, vmax=1,
                    interpolation='nearest')
    ax.set_xlabel("Ligand residue (numeracao PDB)")
    ax.set_ylabel("Receptor residue (numeracao PDB)")
    ax.set_title("Receptor-Ligand Contact Frequency Map\n(fraction of frames, cutoff "
                  f"{args.cutoff/10:.1f} nm)")
    step_x = max(1, n_res_lig // 25)
    step_y = max(1, n_res_rec // 40)
    ax.set_xticks(range(0, n_res_lig, step_x))
    ax.set_xticklabels([lig_labels[i] for i in range(0, n_res_lig, step_x)],
                        rotation=90, fontsize=6)
    ax.set_yticks(range(0, n_res_rec, step_y))
    ax.set_yticklabels([rec_labels[i] for i in range(0, n_res_rec, step_y)], fontsize=6)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Contact frequency")
    plt.tight_layout()
    png_path = os.path.join(args.out_dir, "contact_map.png")
    plt.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[OK] {png_path}", file=sys.stderr)

    # ── interface_residues.csv ─────────────────────────────────────────────
    rec_max = freq.max(axis=1) if n_res_lig else np.zeros(n_res_rec)
    lig_max = freq.max(axis=0) if n_res_rec else np.zeros(n_res_lig)
    iface_path = os.path.join(args.out_dir, "interface_residues.csv")
    with open(iface_path, "w") as fh:
        fh.write("chain,resid,resname,max_contact_freq\n")
        for r, m in zip(rec_internal, rec_max):
            pdb_resid, pdb_resname = rec_map[r]
            fh.write(f"receptor,{pdb_resid},{pdb_resname},{m:.4f}\n")
        for r, m in zip(lig_internal, lig_max):
            pdb_resid, pdb_resname = lig_map[r]
            fh.write(f"ligand,{pdb_resid},{pdb_resname},{m:.4f}\n")
    print(f"[OK] {iface_path}", file=sys.stderr)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
