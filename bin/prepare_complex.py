#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Monta o complexo receptor+ligante para o pipeline CHARMM36/CGenFF:
  - Cadeia A = receptor (numeracao original do PDB preservada, ex. 17-291)
  - Cadeia B = ligante UNL (renumerado para logo apos o ultimo residuo do receptor)

Forca cadeia/numeracao explicitamente (nao confia em pdb2pqr preservar chain ID).

Uso:
    prepare_complex.py --receptor receptor_ph.pdb --ligand unl_ini.pdb \\
        --out complex.pdb
"""
import argparse


def read_atom_lines(path):
    lines = []
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                lines.append(line.rstrip("\n"))
    return lines


def last_resseq(atom_lines):
    last = 0
    for line in atom_lines:
        try:
            resseq = int(line[22:26])
            last = max(last, resseq)
        except ValueError:
            pass
    return last


def renumber_atoms(lines, start=1):
    out = []
    for i, line in enumerate(lines, start):
        serial = i % 100000
        out.append(f"{line[:6]}{serial:5d}{line[11:]}")
    return out


def set_record_resname_chain_resseq(line, record, resname, chain, resseq):
    """
    Reconstroi uma linha ATOM/HETATM com record/resName/chainID/resSeq novos,
    preservando os demais campos de largura fixa do PDB (serial, atom name,
    altLoc, iCode, coordenadas, etc.) pelas colunas padrao:
      1-6 record | 7-11 serial | 13-16 name | 17 altLoc |
      18-20 resName | 21 blank | 22 chainID | 23-26 resSeq | 27+ resto
    """
    return (
        f"{record:<6s}"      # cols  1-6
        f"{line[6:16]}"      # cols  7-16 (serial, blank, atom name) — preservados
        f"{line[16]}"        # col  17 (altLoc) — preservado
        f"{resname:<3s}"     # cols 18-20
        f" "                 # col  21 (blank)
        f"{chain:1s}"        # col  22
        f"{resseq:4d}"       # cols 23-26
        f"{line[26:]}"       # col  27+ (iCode, coordenadas, etc.) — preservado
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--ligand", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--receptor-chain", default="A")
    ap.add_argument("--ligand-chain", default="B")
    ap.add_argument("--ligand-resname", default="UNL")
    args = ap.parse_args()

    receptor_lines = read_atom_lines(args.receptor)
    ligand_lines = read_atom_lines(args.ligand)

    if not receptor_lines:
        raise SystemExit(f"ERRO: nenhum atomo ATOM/HETATM lido de {args.receptor}")
    if not ligand_lines:
        raise SystemExit(f"ERRO: nenhum atomo ATOM/HETATM lido de {args.ligand}")

    # Forca cadeia A + record ATOM em todo o receptor (pdb2pqr pode ter alterado/limpo a coluna de cadeia)
    fixed_receptor = []
    for line in receptor_lines:
        resname = line[17:20].strip()
        resseq = int(line[22:26])
        fixed_receptor.append(
            set_record_resname_chain_resseq(line, "ATOM", resname, args.receptor_chain, resseq)
        )

    rec_last_res = last_resseq(fixed_receptor)
    lig_resseq = rec_last_res + 1

    # Forca cadeia B, HETATM, novo numero de residuo e resname UNL em todo o ligante
    fixed_ligand = []
    for line in ligand_lines:
        fixed_ligand.append(
            set_record_resname_chain_resseq(
                line, "HETATM", args.ligand_resname, args.ligand_chain, lig_resseq
            )
        )

    combined = fixed_receptor + fixed_ligand
    combined = renumber_atoms(combined, start=1)

    with open(args.out, "w") as f:
        for line in combined:
            f.write(line + "\n")
        f.write("TER\n")
        f.write("END\n")

    print(f"[OK] Receptor: {len(fixed_receptor)} atomos, cadeia {args.receptor_chain}, "
          f"residuos ate {rec_last_res}")
    print(f"[OK] Ligante ({args.ligand_resname}): {len(fixed_ligand)} atomos, cadeia "
          f"{args.ligand_chain}, residuo {lig_resseq}")
    print(f"[OK] Complexo escrito em {args.out} ({len(combined)} atomos totais)")


if __name__ == "__main__":
    main()
