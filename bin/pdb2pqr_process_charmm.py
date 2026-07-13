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
    kept = 0
    renamed = 0
    with open(infile) as f, open(outfile, 'w') as out:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom = line[12:16].strip()
                # Remove H adicionados pelo PDB2PQR (pdb2gmx readiciona via charmm36.ff/*.hdb)
                if atom.startswith('H') or (len(atom) > 1 and atom[0].isdigit() and atom[1] == 'H'):
                    continue
                resname = line[17:20].strip()
                new_name = RENAME.get(resname, resname)
                if new_name != resname:
                    line = line[:17] + f'{new_name:<3s}' + line[20:]
                    renamed += 1
                    print(f"  Renomeado: {resname} -> {new_name} (res {line[22:26].strip()})")
                kept += 1
            out.write(line)
    print(f"  {kept} atomos escritos em {outfile} ({renamed} residuos de titulacao renomeados)")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(f"Uso: {sys.argv[0]} input.pdb output.pdb")
    if not os.path.exists(sys.argv[1]):
        sys.exit(f"Arquivo nao encontrado: {sys.argv[1]}")
    process(sys.argv[1], sys.argv[2])
