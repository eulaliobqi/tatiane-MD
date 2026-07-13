#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Limpa output do PDB2PQR (--ff CHARMM --ffout CHARMM) para uso com pdb2gmx + charmm36.ff:
  - Remove hidrogenios adicionados pelo PDB2PQR (pdb2gmx os readiciona com nomenclatura
    CHARMM correta a partir do nome de residuo, que ja codifica o estado de protonacao)
  - Normaliza variantes de nomenclatura de titulacao (algumas versoes do PDB2PQR usam
    tokens no estilo AMBER mesmo com --ff CHARMM) para os nomes esperados pelo
    charmm36-mar2019.ff/merged.rtp: HSD/HSE/HSP (His) e ASPP/GLUP (Asp/Glu protonados)

Uso: pdb2pqr_process_charmm.py input_pdb2pqr.pdb output_gromacs.pdb
"""
import sys
import os

# Normalizacao defensiva: cobre nomenclatura CHARMM nativa (passthrough) e
# possiveis variantes AMBER-style que alguma versao do PDB2PQR possa emitir
# mesmo com --ff CHARMM.
RENAME = {
    'HSD': 'HSD', 'HSE': 'HSE', 'HSP': 'HSP',   # CHARMM nativo (passthrough, documentado)
    'HID': 'HSD', 'HIE': 'HSE', 'HIP': 'HSP',   # fallback estilo AMBER -> CHARMM
    'HISD': 'HSD', 'HISE': 'HSE', 'HISH': 'HSP',
    'ASPP': 'ASPP', 'GLUP': 'GLUP',             # CHARMM nativo (passthrough)
    'ASH': 'ASPP', 'GLH': 'GLUP',               # fallback estilo AMBER -> CHARMM
}


def process(infile, outfile):
    # ── Passo 1: le e filtra (remove nao-ATOM/HETATM e hidrogenios) ──────────
    lines = []
    with open(infile) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                # Descarta REMARK/CRYST1/TER/etc do PDB2PQR — so nos interessam
                # atomos; qualquer TER de fragmento interno seria mal-interpretado
                # como residuo por ferramentas a jusante se deixado passar
                continue

            atom = line[12:16].strip()
            # Remove H adicionados pelo PDB2PQR (pdb2gmx readiciona via charmm36.ff/*.hdb)
            if atom.startswith('H') or (len(atom) > 1 and atom[0].isdigit() and atom[1] == 'H'):
                continue

            lines.append(line)

    # ── Passo 2: normaliza residuo terminal com nome de patch inconsistente ──
    # PDB2QR (--ff CHARMM) as vezes rotula so os atomos do patch de terminal-N
    # (N, CA) com o nome do patch ("TER") em vez do nome real do residuo,
    # deixando o MESMO residuo (mesma cadeia+resSeq) com dois resnames
    # distintos nos seus atomos (ex.: N/CA="TER", C/O/CB/...="HIS"). Isso
    # gera um erro fatal no pdb2gmx ("chain ... do not have a consistent
    # type") porque "TER" nao e um nome de residuo valido em residuetypes.dat.
    # Descoberto rodando de verdade no servidor (residuo 17 do 2I9T).
    key_names = {}
    for line in lines:
        key = (line[21], line[22:26])  # (chainID, resSeq)
        resname = line[17:20].strip()
        key_names.setdefault(key, set()).add(resname)

    fixes = {}
    for key, names in key_names.items():
        if len(names) > 1 and 'TER' in names:
            real_names = names - {'TER'}
            if len(real_names) == 1:
                fixes[key] = real_names.pop()
            else:
                print(f"  AVISO: residuo {key} tem 'TER' + multiplos outros "
                      f"nomes {real_names} — nao foi possivel normalizar automaticamente")

    if fixes:
        for key, real_name in fixes.items():
            print(f"  Normalizado: residuo cadeia={key[0]!r} resSeq={key[1].strip()} "
                  f"'TER' -> {real_name!r} (patch de terminal do PDB2PQR)")

    # ── Passo 3: aplica normalizacao de patch + tabela de titulacao, escreve ─
    kept = 0
    renamed = 0
    with open(outfile, 'w') as out:
        for line in lines:
            key = (line[21], line[22:26])
            resname = line[17:20].strip()

            if resname == 'TER' and key in fixes:
                resname = fixes[key]
                line = line[:17] + f'{resname:<3s}' + line[20:]

            new_name = RENAME.get(resname, resname)
            if new_name != resname:
                line = line[:17] + f'{new_name:<3s}' + line[20:]
                renamed += 1
                print(f"  Renomeado: {resname} -> {new_name} (res {line[22:26].strip()})")
            kept += 1
            out.write(line)
    print(f"  {kept} atomos escritos em {outfile} ({renamed} residuos de titulacao renomeados, "
          f"{len(fixes)} residuos com patch 'TER' normalizados)")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(f"Uso: {sys.argv[0]} input.pdb output.pdb")
    if not os.path.exists(sys.argv[1]):
        sys.exit(f"Arquivo nao encontrado: {sys.argv[1]}")
    process(sys.argv[1], sys.argv[2])
