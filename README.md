# Tatiana-MD

Dinamica molecular de 4 pares receptor-ligante (flavonoides/isoflavonas x alvos
humanos), triados por docking com AutoDock Vina: **2I9T** (NF-kB) + Daidzeina,
**4R3C** (p38α MAPK) + Liquiritigenina, **6FWC** (MAO-B) + Formononetina,
**7K2S** (Keap1) + Biochanina A. Pipeline samplesheet-driven, roda os 4 sistemas
em loop sequencial na GPU unica do servidor (ver "Status em 2026-07-26" abaixo).

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
- [x] Push para `github.com/eulaliobqi/tatiane-MD`
- [x] **2I9T-daidzeina: pipeline completo** (2026-07-13/15), incl. analises extra
      por fase (achado: pose de docking dissocia em ~65-78ns, migra pra outro
      sitio — ver memoria do projeto)
- [x] **Generalizado pra 4 sistemas** (2026-07-26): samplesheet com colunas
      por amostra (`key_residues`, `target_name`, `organism`, `ligand_name`,
      `pdb_id`, `run_extra_analyses`); `ANALYSES_RESIDUES`/`PLOT`/
      `gerar_artigo_md.py` parametrizados em vez de hardcoded pro 2I9T;
      `maxForks=1` nos labels de GPU serializa o loop na GPU unica
- [x] Ligantes novos (4R3C/6FWC/7K2S) parametrizados via CGenFF/ParamChem
      (`cgenff.com`) — estrutura quimica da PubChem (InChI) alinhada a pose
      real do docking via correspondencia de grafo + Kabsch (`.pdbqt` do
      AutoDock nao tem H apolares nem ordem de ligacao explicita)
- [x] **Rodando de verdade no servidor** (2026-07-26): 3 sistemas novos em
      loop sequencial, screen `tatiana-loop-3sistemas`. 2 bugs reais nao
      vistos no 2I9T encontrados/corrigidos rodando pro 6FWC (ver "Bugs reais"):
      coluna de resname CHARMM de 4 letras (ASPP/GLUP) e ordem dos prompts
      `pdb2gmx -asp -glu`. Sem erros desde entao.
- [ ] `PRODUCTION` (100ns cada, ~11h de GPU por sistema) → `POSTPROCESS` →
      analises → `MMGBSA`/`PLOT`/`REPORT` dos 3 sistemas novos — retomar
      verificando `nextflow-loop.log` no screen `tatiana-loop-3sistemas`

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

**Nao rodei `gh`/Nextflow localmente** (Windows sem Nextflow/Java 17+ instalado,
sem GROMACS) — o pipeline foi escrito e revisado por inspecao estatica antes da
primeira execucao real, mas so a execucao real no servidor pega certas classes
de bug (ver abaixo). Rode sempre com `-with-report -with-trace` e acompanhe de
perto as primeiras etapas.

### Bugs reais encontrados na primeira execucao no servidor (2026-07-13)

Nenhum destes apareceu na revisao estatica — todos so surgiram rodando de
verdade. Todos ja corrigidos e commitados:

1. **`pdb2pqr_process_charmm.py`**: bug de indentacao deixava passar qualquer
   linha do PDB2PQR (REMARK/CRYST1/TER), nao so ATOM/HETATM.
2. **Residuo N-terminal com dois nomes**: o PDB2PQR (`--ff CHARMM`) rotula os
   atomos do patch de terminal-N (N, CA) com o nome do patch (`TER`) em vez do
   nome real do residuo — o mesmo residuo (HIS17) saia com `N`/`CA`="TER" e o
   resto="HIS", fazendo o `pdb2gmx` abortar ("chain ... do not have a
   consistent type"). Corrigido com normalizacao em duas passadas no mesmo
   script (agrupa por cadeia+resSeq, substitui "TER" pelo nome real).
3. **`cp` redundante em `TOPOLOGY`**: `posre_UNL.itp` ja chegava staged no cwd
   como input declarado; copiar pra si mesmo falhava ("are the same file").
4. **Exclusao de input do glob de output** (`modules/local/topology/main.nf`):
   o Nextflow exclui por padrao arquivos de input do casamento de
   `path("*.itp")`/`path("*.prm")` — como `unl.itp`/`unl.prm` estavam staged
   direto no cwd (mesmo nome do input), o glob ficava vazio e a task falhava
   com exitcode 0 mas "Missing output file(s)". Corrigido com o mesmo padrao
   de staging em subdiretorio + cópia explicita já usado em
   `box_solvate_ions`/`minimization`/`nvt`/`npt`/`production`.
5. **Descompasso de versao do CGenFF** (o mais sério): `ff/charmm36-mar2019.ff`
   (port original, vendorizado no inicio do projeto) bundla CGenFF **4.1**,
   mas o `.str` do ParamChem para a daidzeina e CGenFF **5.0** — causando 8
   erros reais de `grompp` ("No default Proper Dih./U-B types") concentrados
   no oxigenio do anel pirano da cromona (a mesma regiao que o CGenFF ja
   sinalizava com `param penalty=53`). **Trocado o FF port** para
   `charmm36-feb2026_cgenff-5.0.ff` (MacKerell Lab, ferramenta `charmm2gmx`,
   Wacha & Lemkul JCIM 2023), que bundla CGenFF 5.0 e cobre os 8 termos
   (confirmado por grep contra a conversao real, nao suposicao). Considerado
   e descartado: trocar para AMBER99SB-ILDN+GAFF2/ACPYPE (como
   MD-gromacs/Milena-MD) resolveria tambem, mas exigiria refazer protonacao,
   todos os `.mdp` e o modulo MM-GBSA — mudanca bem maior sem garantia de
   nao ter sua propria lacuna de parametrizacao nessa mesma regiao do anel.

**O que ainda nao foi validado de verdade em 2026-07-13**: nada além de `TOPOLOGY`/
`BOX_SOLVATE_IONS` tinha rodado ate entao (nenhum GPU/mdrun real). Pipeline
completo desde entao pro 2I9T — ver "Status em 2026-07-26" abaixo pros
sistemas novos.

### Bugs reais encontrados generalizando pro 6FWC (2026-07-26)

Nenhum destes apareceu no 2I9T porque nenhum residuo daquele sistema exigiu
protonacao de Asp/Glu (so His, que nao tem esse problema — ver abaixo). Os
dois foram descobertos rodando de verdade contra o receptor do 6FWC (MAO-B),
onde o PROPKA protonou Glu437:

6. **`gmx pdb2gmx` nao le resname CHARMM de 4 letras** (`ASPP`/`GLUP`, usados
   pelo PDB2PQR pra Asp/Glu protonados) **em nenhuma posicao de coluna do
   PDB** — confirmado testando `pdb2gmx` diretamente no servidor com o nome
   em varias posicoes; ele sempre le a janela padrao de 3 colunas e corta a
   1a letra ("GLUP"->"LUP"), erro fatal `residue LUP437 is of type 'Other'`.
   A forma suportada pelo GROMACS e o par de flags interativos
   `pdb2gmx -asp -glu` (prompt 0/1 por residuo), nao embutir a titulacao no
   nome do residuo do PDB. `bin/pdb2pqr_process_charmm.py` agora reverte
   ASPP/GLUP pro nome padrao (ASP/GLU) e grava a lista de residuos
   protonados em `<receptor_ph.pdb>.protonated.txt`.
7. **Ordem dos prompts do `pdb2gmx -asp -glu` nao e sequencial pelo
   arquivo** — ele pergunta **todos os residuos ASP primeiro** (na ordem do
   arquivo), **depois todos os GLU** (idem), nao intercalado por posicao
   como seria natural assumir. A 1a tentativa de `bin/build_asp_glu_answers.py`
   assumiu ordem sequencial misturada e a resposta "protonado" caiu no
   residuo errado **silenciosamente** (sem erro fatal — GLU437 saiu com
   carga -1/nao-protonado em vez de 0/protonado). Corrigido agrupando as
   respostas por tipo (ASP primeiro, GLU depois) — confirmado no topol.top
   final: `residue 437 GLUH rtp GLUP q 0.0`.

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

## Status em 2026-07-13 (fim da sessao) — como continuar

Pipeline pushado, rodando de verdade no servidor desde hoje. `PREPARE_PH` →
`LIGAND_TOPOLOGY` → `PREPARE_COMPLEX` → `TOPOLOGY` → `BOX_SOLVATE_IONS` →
`MINIMIZATION` completos com sucesso (5 bugs reais encontrados e corrigidos
no processo, ver secao dedicada abaixo). `NVT` em andamento quando a sessao
terminou — CUDA confirmado ativo (`nvidia-smi` + `gmx_mpi -version` checados
antes do inicio).

**Pra continuar amanha**, no servidor, dentro do mesmo screen:

```bash
screen -r tatiana-2i9t-daidzeina   # ou screen -d -r se já estiver attached em outro lugar
# conferir onde parou:
tail -30 nextflow.log
nvidia-smi   # confirma se ainda esta rodando ou se terminou/caiu

# se caiu e precisar retomar:
cd ~/gromacs/tatiane-MD
git pull   # so por seguranca, caso eu tenha ajustado mais alguma coisa
nextflow run main.nf --outdir ~/gromacs/results-tatiana -profile local,conda \
    -with-report -with-trace -resume 2>&1 | tee -a nextflow.log
```

**Depois que `PRODUCTION`/`POSTPROCESS`/`ANALYSES` terminarem** (a etapa longa,
100 ns), os passos finais rodam automaticamente dentro do proprio Nextflow
(`MMGBSA`, `PLOT`, `REPORT`) — mas os 2 envs conda dedicados deste projeto
(`plot-env-tatiana`, `mmgbsa-env`) **ainda nao foram criados** no servidor
(nao confirmamos isso na sessao). Rodar antes do Nextflow chegar em `PLOT`:

```bash
mamba create -n plot-env-tatiana python=3.11 numpy matplotlib -y
mamba create -n mmgbsa-env -c conda-forge -c bioconda gmx_mmpbsa python=3.11 ambertools -y
mamba install -n mmgbsa-env -c conda-forge gromacs -y   # gmx_MMPBSA precisa de um gmx no PATH desse env
```

Fallback bash equivalente (se preferir nao usar Nextflow numa etapa especifica):

```bash
bash bin/analyze.sh 2>&1 | tee analyze.log
python bin/plot_results.py --analise-dir results/2I9T-daidzeina/analise --titulo "2I9T + Daidzeina"
python bin/gerar_artigo_md.py --analise-dir results/2I9T-daidzeina/analise \
    --mmgbsa-dir results/2I9T-daidzeina/mmgbsa --out docs/artigo_md.md
```

## Status em 2026-07-26 (fim da sessao) — como continuar

Pipeline generalizado de 1 pra 4 sistemas (2I9T ja concluido + 4R3C/6FWC/7K2S
novos), 2 bugs reais de protonacao Asp/Glu encontrados e corrigidos rodando
de verdade (ver secao "Bugs reais... 2026-07-26" acima). Loop rodando sem
erros desde a correcao: `6FWC` completou toda a equilibracao (NPT) e entrou
em `PRODUCTION`; `4R3C` em `NPT`; `7K2S` na fila de `PRODUCTION` (esperando a
GPU liberar — `maxForks=1` serializa corretamente).

**Pra continuar amanha**, no servidor:

```bash
ssh eulalio@200.235.143.10
cd ~/gromacs/tatiane-MD
screen -r tatiana-loop-3sistemas   # ou screen -d -r se ja estiver attached em outro lugar
tail -60 nextflow-loop.log
nvidia-smi   # confirma se ainda esta rodando ou se terminou/caiu

# se caiu e precisar retomar:
git pull
nextflow run main.nf -profile local,conda -resume -with-report -with-trace 2>&1 | tee -a nextflow-loop.log
```

Estimativa: ~11h de GPU de `PRODUCTION` por sistema novo (2 restantes a
comecar do zero + 1 em andamento), total ainda em aberto quando a sessao
terminou. Depois que os 3 sistemas novos terminarem `PRODUCTION`→`POSTPROCESS`
→`ANALYSES`, conferir:

- `docs/<sample_id>/artigo_md.md` — relatorio parametrizado por sistema
  (alvo/PDB/organismo/ligante/residuos-chave, ver `bin/gerar_artigo_md.py`)
- `~/gromacs/results-tatiana/<sample_id>/analise/painel_resumo.png` — RMSD/
  RMSF/Rg/contatos/H-bonds/SASA + distancias dos residuos-chave por sistema
- Os 6 modulos extra (`PHASE_SPLIT` etc.) **nao rodam** pros 3 sistemas
  novos nesta rodada (`run_extra_analyses=false` no samplesheet) — decidir
  caso a caso, depois de ver RMSD/distancias, se algum justifica uma 2a
  passada com esses modulos ligados

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
ff/         - CHARMM36 port (charmm36-feb2026_cgenff-5.0.ff, MacKerell Lab,
              ferramenta charmm2gmx — Wacha & Lemkul, JCIM 2023 — bundla
              CGenFF 5.0, casando com a versao do .str do ParamChem; baixado
              de mackerell.umaryland.edu e vendorizado aqui porque o servidor
              bloqueia HTTPS externo exceto github.com)
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
