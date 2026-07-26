#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gera a sequencia de respostas para os prompts interativos do
`gmx pdb2gmx -asp -glu` (0=nao protonado/padrao, 1=protonado).

IMPORTANTE (confirmado rodando de verdade no servidor): o pdb2gmx NAO
pergunta na ordem sequencial do arquivo misturando ASP/GLU -- ele agrupa
por TIPO, perguntando primeiro TODOS os residuos ASP (na ordem em que
aparecem no PDB) e só depois TODOS os residuos GLU (idem). Alimentar as
respostas na ordem errada (sequencial misturada) faz a resposta cair no
residuo errado silenciosamente -- sem erro, mas com o estado de
protonacao trocado (confirmado: GLU437 marcado "1" saiu como GLU padrao
carga -1 em vez de GLUP carga 0 na primeira tentativa deste fix).

Ver bin/pdb2pqr_process_charmm.py para de onde vem protonated.txt (lista
de residuos que o PROPKA determinou como protonados, ASPP/GLUP).

Uso: build_asp_glu_answers.py receptor.pdb protonated.txt
Imprime uma resposta (0 ou 1) por linha no stdout, pronta pra virar stdin
do pdb2gmx (concatenar com a resposta final de terminal N/C se necessario).
Se protonated.txt nao existir (nenhum residuo protonado), ainda imprime
"0" pra cada Asp/Glu -- pdb2gmx sempre pergunta quando -asp/-glu estao
ativos, independente de haver algum protonado ou nao.
"""
import sys
import os


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Uso: {sys.argv[0]} receptor.pdb protonated.txt")
    receptor_pdb, protonated_file = sys.argv[1], sys.argv[2]

    protonated = set()
    if os.path.exists(protonated_file):
        with open(protonated_file) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    protonated.add((parts[0], parts[1]))  # (chain, resSeq)

    seen = set()
    asp_answers = []
    glu_answers = []
    with open(receptor_pdb) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            resname = line[17:20].strip()
            if resname not in ('ASP', 'GLU'):
                continue
            key = (line[21], line[22:26].strip())
            if key in seen:
                continue
            seen.add(key)
            answer = '1' if key in protonated else '0'
            (asp_answers if resname == 'ASP' else glu_answers).append(answer)

    # Todas as respostas de ASP primeiro, depois todas as de GLU -- ver nota
    # acima sobre o agrupamento por tipo do pdb2gmx.
    print('\n'.join(asp_answers + glu_answers))


if __name__ == '__main__':
    main()
