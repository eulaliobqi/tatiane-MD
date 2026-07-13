#!/usr/bin/env python3
"""
Mescla topologia de receptor (pdb2gmx, CHARMM36) com ligante pequeno (CGenFF/cgenff_charmm2gmx.py).

Adaptado de MD-gromacs/bin/merge_small_molecule_topology.py (pipeline BEN, AMBER+ACPYPE/GAFF2)
para o par CHARMM36 (receptor) + CGenFF (ligante) — mesma logica de mescla de .gro/.top,
so troca a origem da topologia do ligante.

Operacoes:
1. Mescla receptor.gro + ligand.gro -> complexo.gro (com renumeracao de atomos)
2. Renumera residuo do ligante para sequencia apos o receptor
3. Patcheia receptor.top para incluir ligand.itp e adicionar ligante em [molecules]

Uso:
    python3 bin/merge_small_molecule_topology.py \\
        --protein-gro receptor.gro \\
        --ligand-gro  unl_ini.gro \\
        --protein-top receptor.top \\
        --ligand-itp  unl.itp \\
        --ligand-mol  UNL \\
        --out-gro     complexo.gro \\
        --out-top     topol.top
"""
import argparse
import re
import sys
from pathlib import Path


def read_gro(gro_path):
    """Le arquivo GRO -> (titulo, natoms, linhas_atomos, linha_box)."""
    text = Path(gro_path).read_text()
    lines = text.splitlines()
    title  = lines[0]
    natoms = int(lines[1].strip())
    atoms  = lines[2 : 2 + natoms]
    box    = lines[2 + natoms]
    return title, natoms, atoms, box


def last_resnum(atom_lines):
    """Retorna numero do ultimo residuo de uma lista de linhas GRO."""
    last = 1
    for line in atom_lines:
        try:
            last = int(line[:5].strip())
        except ValueError:
            pass
    return last


def merge_gro(protein_gro, ligand_gro, out_gro, ligand_mol):
    """Mescla dois GROs num unico arquivo com atomos renumerados."""
    _, n_prot, prot_atoms, box = read_gro(protein_gro)
    _, n_lig, lig_atoms, _    = read_gro(ligand_gro)

    last_res = last_resnum(prot_atoms)
    lig_res  = last_res + 1

    # Atualiza numero de residuo do ligante (col 0-4) para seguir receptor
    lig_renumbered = []
    for line in lig_atoms:
        lig_renumbered.append(f"{lig_res:5d}" + line[5:])

    combined = prot_atoms + lig_renumbered

    # Renumera atomos sequencialmente (col 15-19)
    final = []
    for i, line in enumerate(combined, 1):
        final.append(line[:15] + f"{i % 100000:5d}" + line[20:])

    total = n_prot + n_lig
    with open(out_gro, "w") as f:
        f.write(f"Complex protein + {ligand_mol} (CHARMM36/CGenFF)\n")
        f.write(f"{total:5d}\n")
        for line in final:
            f.write(line + "\n")
        f.write(box + "\n")

    print(f"GRO mesclado: {n_prot} prot + {n_lig} {ligand_mol} = {total} atomos | "
          f"{ligand_mol} = residuo {lig_res}")
    return lig_res


def patch_topology(protein_top, ligand_itp, ligand_prm, ligand_mol, out_top):
    """
    Insere #include do PRM (parametros extra) e do ITP do ligante no topol.top do receptor.
    Adiciona entrada do ligante em [molecules].
    """
    content = Path(protein_top).read_text()
    itp_name = Path(ligand_itp).name
    prm_name = Path(ligand_prm).name if ligand_prm else None

    # Insere #include apos a linha de include do forcefield (prm ANTES do itp — o itp
    # referencia os parametros extras gerados pelo cgenff_charmm2gmx.py)
    includes = f'#include "{prm_name}"\n' if prm_name else ""
    includes += f'#include "{itp_name}"'

    ff_pat = r'(#include\s+"[^"]+\.ff/forcefield\.itp")'
    if re.search(ff_pat, content):
        content = re.sub(ff_pat, rf'\1\n{includes}', content, count=1)
    else:
        # Fallback: antes do primeiro [ moleculetype ]
        content = re.sub(
            r'(\[\s*moleculetype\s*\])',
            rf'{includes}\n\n\1',
            content, count=1
        )

    # Adiciona ligante em [molecules] (ao final do arquivo)
    content = content.rstrip("\n") + f"\n{ligand_mol:<20} 1\n"

    Path(out_top).write_text(content)
    print(f"TOP patcheado: #include {itp_name!r}"
          + (f" + {prm_name!r}" if prm_name else "")
          + f" + {ligand_mol} 1 -> {out_top}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein-gro", required=True)
    ap.add_argument("--ligand-gro",  required=True)
    ap.add_argument("--protein-top", required=True)
    ap.add_argument("--ligand-itp",  required=True)
    ap.add_argument("--ligand-prm",  default=None,
                     help="unl.prm gerado por cgenff_charmm2gmx.py (parametros bonded extra)")
    ap.add_argument("--ligand-mol",  default="UNL")
    ap.add_argument("--out-gro",     required=True)
    ap.add_argument("--out-top",     required=True)
    args = ap.parse_args()

    merge_gro(args.protein_gro, args.ligand_gro, args.out_gro, args.ligand_mol)
    patch_topology(args.protein_top, args.ligand_itp, args.ligand_prm, args.ligand_mol, args.out_top)
    return 0


if __name__ == "__main__":
    sys.exit(main())
