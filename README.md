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
- [x] Pipeline **Nextflow DSL2** completo (`main.nf` + `modules/local/*`), mirror da
      arquitetura do MD-gromacs/Milena-MD, 15 processos: PREPARE_PH →
      LIGAND_TOPOLOGY → PREPARE_COMPLEX → TOPOLOGY → BOX_SOLVATE_IONS →
      MINIMIZATION → NVT → NPT → PRODUCTION → POSTPROCESS → ANALYSES +
      ANALYSES_RESIDUES → MMGBSA → PLOT + REPORT
- [x] Scripts bash equivalentes mantidos como fallback/referencia (`bin/run_md.sh`,
      `bin/analyze.sh`) — mesma logica, uteis se Nextflow nao estiver disponivel
- [x] MM-GBSA (`gmx_MMPBSA`) — reescrito do zero evitando o bug real que
      derrubou essa mesma etapa no projeto irmao Milena-MD (ver secao dedicada)
- [x] Gerador de relatorio (`docs/artigo_md.md`, com secao de MM-GBSA)
- [ ] **Simulacao ainda NAO rodou** — nada foi executado no servidor, so preparado localmente
- [ ] Push para `github.com/eulaliobqi/tatiane-MD` (a fazer por voce — remoto ja
      configurado localmente, so falta autenticar e confirmar)

## Arquitetura Nextflow

Pipeline DSL2 completo em `main.nf` + `modules/local/*/main.nf`, samplesheet-driven
(`assets/samplesheet.csv`, hoje com uma unica linha — `2I9T-daidzeina` — mas a
arquitetura ja suporta adicionar outros pares receptor-ligante sem reescrever nada).
Perfis `local`/`slurm`/`conda` em `conf/*.config`, mesma convencao do MD-gromacs.

```bash
mamba activate md-gromacs
cd ~/gromacs/Tatiana-MD
nextflow run main.nf --outdir ~/gromacs/results-tatiana -profile local,conda
```

**Nao roda `gh`/Nextflow localmente nesta sessao** (Windows sem Nextflow/Java 17+
instalado, sem GROMACS) — o pipeline foi escrito e revisado por inspecao estatica
contra os mesmos padroes ja validados em producao no MD-gromacs/Milena-MD, mas
**precisa ser testado de verdade no servidor antes de confiar 100%** — ver
"Riscos conhecidos" abaixo.

### Riscos conhecidos (nao testados em ambiente real)

- Todo o wiring de canais Groovy (`main.nf`) foi escrito seguindo os padroes
  exatos ja validados no MD-gromacs/Milena-MD (join por `meta.id`, nao por
  `meta` Map inteiro — bug real documentado no Milena-MD que fazia o join
  nunca emitir nada, silenciosamente), mas nunca foi executado. Primeira
  rodada real deve ser acompanhada de perto (`nextflow run ... -with-report
  -with-trace`).
- Modulo MMGBSA foi escrito pelo agente `bioinformatics` (a pedido explicito),
  aplicando a causa-raiz mais provavel do bug que derrubou essa mesma etapa 3x
  no Milena-MD — mas nao pode ser testado sem `gmx_MMPBSA` real. `errorStrategy
  'ignore'` garante que uma falha aqui nao derruba o resto do pipeline.

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
# 1. Local (Windows) — revisar o diff, depois push (repo ja existe, vazio)
cd ~/.claude/Tatiana-MD
git push -u origin main

# 2. Servidor
ssh eulalio@200.235.143.10
screen -S tatiana-2i9t-daidzeina
mamba activate md-gromacs
cd ~/gromacs && git clone https://github.com/eulaliobqi/tatiane-MD.git   # ou git pull se ja existir
cd tatiane-MD
df -h /home   # checar espaco antes (regra do laboratorio)

# 2b. Envs conda dedicados (so na PRIMEIRA vez — nao existem ainda neste
#     projeto; nao reusar os do Milena-MD, isolamento entre projetos-execucao)
mamba create -n plot-env-tatiana python=3.11 numpy matplotlib -y
mamba create -n mmgbsa-env -c conda-forge -c bioconda gmx_mmpbsa python=3.11 ambertools -y
# gmx_MMPBSA tambem precisa de um binario GROMACS no PATH desse env (usa
# editconf/trjconv internamente) — o modulo MMGBSA avisa no log se faltar:
mamba install -n mmgbsa-env -c conda-forge gromacs -y   # build CPU, sem CUDA

# 3a. Via Nextflow (recomendado — arquitetura completa, retomavel via -resume)
nextflow run main.nf --outdir ~/gromacs/results-tatiana -profile local,conda \
    -with-report -with-trace 2>&1 | tee nextflow.log
# producao mais longa: --time_ns 200

# 3b. OU via bash (fallback/referencia, mesma logica sem o overhead do Nextflow)
bash bin/run_md.sh 2>&1 | tee run_md.log   # TIME_NS=200 bash bin/run_md.sh p/ mais longo
bash bin/analyze.sh 2>&1 | tee analyze.log
python bin/plot_results.py --analise-dir results/2I9T-daidzeina/analise --titulo "2I9T + Daidzeina"
python bin/gerar_artigo_md.py --analise-dir results/2I9T-daidzeina/analise \
    --mmgbsa-dir results/2I9T-daidzeina/mmgbsa --out docs/artigo_md.md
```

Cada etapa e retomavel: no Nextflow, `-resume` reaproveita o cache de tasks ja
concluidas; no `bin/run_md.sh`, o script checa se o arquivo de output final de
cada etapa ja existe antes de rodar de novo.

Se `gh auth login` nao estiver configurado, o `git push` acima usa o Git
Credential Manager do Windows (deve abrir o navegador pra autenticar) — o
repositorio `github.com/eulaliobqi/tatiane-MD` ja existe e esta vazio,
confirmado via `git ls-remote` antes de eu configurar o remote localmente.

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
main.nf     - pipeline Nextflow DSL2 (orquestracao dos 15 processos)
modules/local/*/main.nf - um processo por etapa (mirror MD-gromacs/Milena-MD)
conf/       - local.config / slurm.config / base.config (perfis de recursos)
assets/     - samplesheet.csv (sample_id,receptor,ligand_mol2,ligand_str)
inputs/     - PDBs/mol2/str originais (receptor bruto, receptor "fixed" por
              PDBFixer [nao usado diretamente — ver nota de protonacao acima],
              ligante docado, saida do ParamChem/CGenFF)
ff/         - CHARMM36 port (charmm36-mar2019.ff, vendorizado do GitHub
              intbio/gromacs_ff — servidor bloqueia HTTPS externo exceto github.com)
bin/        - scripts reutilizados pelos processos Nextflow (conversao de
              topologia, preparo de complexo, merge, plot, gerador de relatorio)
              + bin/run_md.sh e bin/analyze.sh (fallback bash standalone)
mdp/        - parametros de simulacao de referencia p/ o fallback bash
              (os processos Nextflow geram os mdp inline, mesmo conteudo)
docs/       - artigo_md.md (gerado)
results/    - saida do pipeline (gitignored, regeneravel)
work/       - work dir do Nextflow + scratch de testes locais (gitignored)
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

## MM-GBSA — historico e o que mudou

`modules/local/mmgbsa/main.nf` calcula energia livre de ligacao com `gmx_MMPBSA`
(protocolo de trajetoria unica, `igb=2`, decomposicao por residuo). **Esta mesma
ferramenta ja falhou de forma irreconciliavel no projeto irmao Milena-MD**
(serie trypsin×GORE12T, `modules/local/mmgbsa_robust/main.nf`), 3 tentativas de
fix, oficialmente abandonada (ver `bioinformatics.md`). Eu investiguei os logs
reais dessa falha (`Milena-MD/*/mmgbsa/mmgbsa.log`) e encontrei o sintoma exato:
a chamada final ao `gmx_MMPBSA` saiu com `-cs -ct -ci` **vazios** (`Several args
are duplicated`) mesmo depois de uma tentativa de fix que passava os caminhos
como argumentos posicionais de um script bash intermediario via `mamba run`.

O modulo aqui foi **reescrito do zero** (nao e um patch do antigo) evitando essa
camada intermediaria de script/argumentos posicionais inteiramente — a chamada
ao `gmx_MMPBSA` interpola os caminhos de arquivo diretamente no bloco `script:`
do Nextflow (que ja faz a substituicao de texto antes de qualquer shell rodar),
sem nenhum `run_mmgbsa.sh` nem `$1`/`$2`/`$3` no meio. Ambiente conda
`mmgbsa-env` continua separado do `md-gromacs` (AmberTools conflita com o build
CUDA do GROMACS — mesma razao documentada no laboratorio). Precisa existir no
servidor antes de rodar (ver comentario no topo do modulo se precisar criar).

`errorStrategy 'ignore'` no processo: se `gmx_MMPBSA` falhar de novo por
qualquer motivo, **o resto do pipeline continua normalmente** — RMSD/RMSF/Rg/
contatos/H-bonds/SASA/Arg30/Glu279 nao dependem do MM-GBSA (ver comentario em
`modules/local/report/main.nf` sobre por que REPORT nao usa o canal Nextflow do
MMGBSA — outro bug real documentado no Milena-MD, join contra canal de processo
com `errorStrategy 'ignore'` travou o PLOT inteiro em producao).

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
