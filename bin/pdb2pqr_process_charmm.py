#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Limpa output do PDB2PQR (--ff CHARMM --ffout CHARMM) para uso com pdb2gmx + charmm36.ff:
  - Remove hidrogenios adicionados pelo PDB2PQR (pdb2gmx os readiciona com nomenclatura
    CHARMM correta a partir do nome de residuo/estado de titulacao)
  - Normaliza HSD/HSE/HSP (His protonado, 3 letras) — cabe na coluna padrao do PDB,
    gmx pdb2gmx le direto sem interatividade.
  - Para Asp/Glu protonados (ASPP/GLUP, 4 letras): gmx pdb2gmx NAO consegue ler um
    resname de 4 letras em NENHUMA posicao de coluna do PDB (confirmado rodando de
    verdade no servidor pro 6FWC — sempre le so as 3 colunas fixas padrao, cortando
    a 1a letra: "GLUP"->"LUP", erro fatal "residue LUP437 is of type 'Other'"). A
    forma suportada pelo GROMACS pra protonacao de Asp/Glu e o par de flags
    interativos `pdb2gmx -asp -glu` (um prompt 0/1 por residuo Asp/Glu, na ordem
    em que aparecem no arquivo) — NAO embutir a titulacao no nome do residuo do
    PDB. Por isso: reverte ASPP/GLUP pro nome padrao de 3 letras (ASP/GLU, cabe
    na coluna normal) e grava quais residuos devem responder "1" em
    <output>.protonated.txt, que o modulo TOPOLOGY usa via
    bin/build_asp_glu_answers.py pra montar as respostas do prompt do pdb2gmx.

Uso: pdb2pqr_process_charmm.py input_pdb2pqr.pdb output_gromacs.pdb
Gera tambem: output_gromacs.pdb.protonated.txt (linhas "chain resSeq resname_base")
"""
import sys
import os

# His protonado: 3 letras, cabe na coluna padrao do PDB, gmx le direto sem
# interatividade — passthrough/normalizacao de variantes AMBER-style.
RENAME = {
    'HSD': 'HSD', 'HSE': 'HSE', 'HSP': 'HSP',   # CHARMM nativo (passthrough, documentado)
    'HID': 'HSD', 'HIE': 'HSE', 'HIP': 'HSP',   # fallback estilo AMBER -> CHARMM
    'HISD': 'HSD', 'HISE': 'HSE', 'HISH': 'HSP',
}
# Titulacao de Asp/Glu (nomes de 4 letras, NUNCA escritos no PDB de saida —
# ver docstring): nome-detectado -> nome-base-3-letras.
TITRATION = {
    'ASPP': 'ASP', 'ASH': 'ASP',
    'GLUP': 'GLU', 'GLH': 'GLU',
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
        # Janela de 4 colunas [16:20) -- no PDB2PQR bruto, nomes de titulacao
        # de 4 letras (ASPP/GLUP) "emprestam" a coluna de altLoc (17),
        # alinhados a direita ate a col.20; strip() cobre o caso de 3 letras
        # sem quebrar nada (so remove o espaco extra a esquerda).
        resname = line[16:20].strip()
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

    # ── Passo 3: resolve titulacao/patch, escreve sempre no campo padrao de
    #     3 colunas (18-20) — ASPP/GLUP nunca chegam a ser escritos como tal,
    #     so seus residuos-base (ASP/GLU) + entrada na lista protonated. ──
    kept = 0
    renamed = 0
    protonated = []
    with open(outfile, 'w') as out:
        for line in lines:
            key = (line[21], line[22:26])
            resname = line[16:20].strip()

            if resname == 'TER' and key in fixes:
                resname = fixes[key]

            if resname in TITRATION:
                write_name = TITRATION[resname]
                protonated.append((key[0], key[1].strip(), write_name))
            else:
                write_name = RENAME.get(resname, resname)

            if write_name != resname:
                renamed += 1

            # Sobrescreve as 4 colunas [16:20) por inteiro (nao so [17:20)) --
            # quando o nome original era de titulacao de 4 letras (ASPP/GLUP),
            # a coluna 16 (altLoc) tinha a 1a letra do nome ("G" de "GLUP");
            # escrever so em [17:20) deixava essa letra sobrando na saida
            # ("GGLU" em vez de "GLU"). Pego rodando de verdade pro 6FWC.
            line = line[:16] + f' {write_name:<3s}' + line[20:]
            kept += 1
            out.write(line)

    # Sempre grava o arquivo (mesmo vazio, 0 residuos protonados) -- o modulo
    # TOPOLOGY declara isso como input fixo (nao opcional), mais simples que
    # lidar com arquivo ausente no Nextflow.
    protonated_unique = sorted(set(protonated))
    prot_file = outfile + '.protonated.txt'
    with open(prot_file, 'w') as pf:
        for chain, resid, base in protonated_unique:
            pf.write(f"{chain} {resid} {base}\n")
    if protonated_unique:
        listagem = ", ".join(f"{base}{resid}" for _, resid, base in protonated_unique)
        print(f"  {len(protonated_unique)} residuo(s) Asp/Glu protonados pelo PROPKA "
              f"(revertidos pro nome padrao no PDB, lista em {prot_file} p/ "
              f"pdb2gmx -asp -glu): {listagem}")
    else:
        print(f"  Nenhum residuo Asp/Glu protonado pelo PROPKA ({prot_file} vazio)")

    print(f"  {kept} atomos escritos em {outfile} ({renamed} residuos de titulacao/patch "
          f"normalizados, {len(fixes)} residuos com patch 'TER')")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(f"Uso: {sys.argv[0]} input.pdb output.pdb")
    if not os.path.exists(sys.argv[1]):
        sys.exit(f"Arquivo nao encontrado: {sys.argv[1]}")
    process(sys.argv[1], sys.argv[2])
