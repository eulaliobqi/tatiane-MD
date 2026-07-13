# Tatiana-MD

Dinamica molecular do dominio de ligacao a DNA de **NF-kB** (PDB **2I9T**, p65/p50
heterodimero — cadeia A usada isoladamente, residuos 17-291) em complexo com a
isoflavona **daidzeina** (resname `UNL`), pose obtida por docking com AutoDock Vina.

Campo de forca: **CHARMM36m** (proteina) + **CGenFF 5.0** via ParamChem (ligante) —
diferente dos demais pipelines deste laboratorio (MD-gromacs, Milena-MD), que usam
AMBER99SB-ILDN + GAFF2/ACPYPE. A escolha segue o que ja havia sido preparado
manualmente (submissao ao ParamChem ja feita antes desta sessao) e o protocolo
CHARMM36+CGenFF+`cgenff_charmm2gmx.py` e o padrao-ouro documentado para proteina +
pequena molecula em GROMACS (tutorial Lemkul, T4-lisozima+JZ4).

## Status

- [x] Topologia do ligante convertida e verificada localmente (`bin/cgenff_charmm2gmx.py`)
- [x] Scripts de preparo do complexo, topologia, box/solvatacao/ions, EM/NVT/NPT/producao
- [x] Script de analises (RMSD/RMSF/Rg/contatos/H-bonds/SASA + Arg30/Glu279)
- [x] Gerador de relatorio (`docs/artigo_md.md`)
- [ ] **Simulacao ainda NAO rodou** — nada foi executado no servidor, so preparado localmente
- [ ] Push para GitHub (a fazer por voce)

## Antes de rodar — decisoes ja tomadas, revisar se necessario

| Decisao | Valor usado | Por que |
|---|---|---|
| pH | **7,4** (nao 8,2) | Alvo humano (NF-kB), nao midgut de inseto — ver `bioinformatics.md` |
| Cation | **Na+** (nao K+) | "usar NA para mamiferos" — comentario ja existente em `MD-gromacs/nextflow.config` |
| Concentracao ionica | **0,15 M** (nao 0,10 M) | Fisiologico humano padrao |
| Box | cubica, **1,2 nm** (nao 2,0 nm) | Complexo globular+molecula pequena, nao proteina-peptideo alongada; 1,2nm e o padrao do tutorial CHARMM36 |
| Nao-bonded | `vdw-modifier = Force-switch`, `rvdw-switch=1,0`, `DispCorr=no` | Recomendacao oficial CHARMM36 p/ GROMACS — NAO copiar mdp AMBER de outros projetos |
| Protonacao do receptor | **pdb2pqr + PROPKA re-executado** (nao reusar `receptor-2I9T-fixed.pdb` como estava) | O PDBFixer previamente aplicado deixou **as 9 histidinas em HSP** (dupla protonacao) uniformemente — padrao suspeito de default ingenuo, nao calculo de pKa real. `bin/run_md.sh` roda PROPKA do zero a partir de `inputs/receptor-2I9T-original.pdb` |
| Tempo de producao | 100 ns (padrao), ajustavel via `TIME_NS=200 bash bin/run_md.sh` | Convencao do laboratorio |

**⚠️ Qualidade da topologia do ligante:** o CGenFF retornou *param penalty = 53*
(acima do limiar de 50 que a propria ferramenta define como "requer validacao
extensa" — analogia pobre para parte do anel cromona da daidzeina). A simulacao e
cientificamente valida para uma primeira rodada exploratoria, mas resultados de
energia livre/ligacao devem ser tratados como preliminares ate validacao adicional
(ex. reotimizacao QM dos dihedros de maior penalidade). Ver `inputs/ligand-UNL.str`
(comentario `param penalty= 53.000 ; charge penalty= 23.263`).

## Como rodar amanha

```bash
# 1. Local (Windows) — revisar o diff, depois push
cd ~/.claude/Tatiana-MD   # ou o caminho local equivalente
git remote add origin git@github.com:eulaliobqi/Tatiana-MD.git   # se ainda nao existir
git push -u origin main

# 2. Servidor
ssh eulalio@200.235.143.10
screen -S tatiana-2i9t-daidzeina
mamba activate md-gromacs
cd ~/gromacs && git clone git@github.com:eulaliobqi/Tatiana-MD.git   # ou git pull se ja existir
cd Tatiana-MD
df -h /home   # checar espaco antes (regra do laboratorio)

# 3. Dinamica completa (PREPARE_PH -> ... -> PRODUCTION -> POSTPROCESS), retomavel
bash bin/run_md.sh 2>&1 | tee run_md.log
# para producao mais longa: TIME_NS=200 bash bin/run_md.sh

# 4. Analises
bash bin/analyze.sh 2>&1 | tee analyze.log
python bin/plot_results.py

# 5. Relatorio (preenche docs/artigo_md.md com os numeros reais)
python bin/gerar_artigo_md.py
```

Cada etapa de `run_md.sh` e retomavel: se cair a conexao/screen, rodar o mesmo
comando de novo pula as etapas ja concluidas (checa se o arquivo de output final
de cada etapa ja existe).

## Depois que a simulacao terminar (pedidos feitos durante o preparo)

1. **`python bin/gerar_artigo_md.py`** — gera `docs/artigo_md.md` no formato
   padrao dos outros projetos (Resumo/Introducao/Metodologia/Resultados),
   com os numeros reais da simulacao. Passar por `/humanizer` antes de usar
   em qualquer documento final.
2. **Checklist de convergencia** (secao 3.3 do `docs/artigo_md.md`, gerada como
   TODO) — comparar robustez com:
   - outros pipelines ja validados do laboratorio (mesmo protocolo de
     equilibracao/cutoffs, ver `~/.claude/.claude/agents/bioinformatics.md`)
   - literatura publicada sobre daidzeina/isoflavonas e NF-kB (buscar e citar
     explicitamente — **nao fabricar numeros de terceiros**, ver skill
     `auditing-academic-sources`)
   - persistencia dos contatos Arg30/Glu279 previstos no docking original
     (`results/2I9T-daidzeina/10_analysis/dist_arg30.xvg`,
     `dist_glu279.xvg` — comparar com os 4,7-4,8 Å e 1,9 Å do docking)

## Estrutura do repositorio

```
inputs/     - PDBs/mol2/str originais (receptor bruto, receptor "fixed" por
              PDBFixer [nao usado diretamente — ver nota de protonacao acima],
              ligante docado, saida do ParamChem/CGenFF)
ff/         - CHARMM36 port (charmm36-mar2019.ff, vendorizado do GitHub
              intbio/gromacs_ff — servidor bloqueia HTTPS externo exceto github.com)
bin/        - todos os scripts (conversao de topologia, preparo de complexo,
              driver da MD, analises, plot, gerador de relatorio)
mdp/        - parametros de simulacao (em/nvt/npt/production + ions.mdp)
docs/       - artigo_md.md (gerado)
results/    - saida do pipeline (gitignored, regeneravel)
work/       - scratch de testes locais (gitignored)
```

## Arquivos-fonte relevantes

- `inputs/receptor-2I9T-original.pdb` — cadeia A crua do PDB 2I9T (SEM hidrogenios,
  usada como entrada real do `bin/run_md.sh`, via pdb2pqr+PROPKA)
- `inputs/receptor-2I9T-fixed.pdb` — saida do PDBFixer de uma sessao anterior;
  **mantida so como referencia**, NAO usada pelo pipeline (motivo: todas as 9
  histidinas saíram HSP, ver tabela acima)
- `inputs/ligand-UNL.pdb` — pose de docking original (21 atomos, H parciais)
- `inputs/ligand-UNL.cgenff.mol2` / `ligand-UNL.str` — saida do ParamChem/CGenFF
  (29 atomos com H completo; `RESI`/nome da molecula corrigidos de `ligand-U`
  para `UNL` nesta sessao — ver nota tecnica abaixo)
- `inputs/complex-2I9T-daidzeina-docked.pdb` — complexo docado original, referencia

## Nota tecnica: por que o resname foi renomeado de `ligand-U` para `UNL`

O ParamChem/CGenFF truncou o nome do arquivo mol2 submetido (`ligand-UNL`) em
`RESI ligand-U` (8 caracteres) no `.str`, e o mesmo texto aparecia como nome da
molecula no cabecalho `@<TRIPOS>MOLECULE` do mol2. Isso quebra dois pontos do
pipeline se nao corrigido:
1. Nome de residuo PDB tem largura fixa de 3-4 colunas — `ligand-U` (8 chars,
   com hifen) corrompe o alinhamento de colunas do PDB.
2. `cgenff_charmm2gmx.py` le o nome da molecula do mol2 (`@<TRIPOS>MOLECULE`) e
   **sobrescreve** o nome do residuo lido do `.str` — descoberto rodando a
   conversao de verdade e inspecionando o `.itp` gerado (o `[ moleculetype ]`
   saiu `ligand-U` mesmo depois de eu corrigir soh o `.str`; so ficou
   consistente apos corrigir os dois arquivos).

Ambos os arquivos em `inputs/` ja estao corrigidos (`RESI UNL` no `.str`, `UNL`
no mol2) — a conversao foi re-executada e verificada localmente (29/29 atomos,
sem tokens de erro/placeholder no `.itp`/`.prm` gerados).
