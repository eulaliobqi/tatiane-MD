# Dinamica Molecular — Monoamina oxidase B (MAO-B) (PDB 6FWC) + Formononetina — Secoes do Artigo

*Gerado automaticamente por `bin/gerar_artigo_md.py`. Revisar antes de usar em texto final
(passar por /humanizer e pela skill auditing-academic-sources antes de qualquer submissao).*

## Resumo

Este trabalho investigou por dinamica molecular (100 ns) a estabilidade do complexo entre
Monoamina oxidase B (MAO-B) (PDB 6FWC), Homo sapiens e Formononetina, um candidato identificado por triagem
virtual (AutoDock Vina). O sistema foi parametrizado com o campo de forca CHARMM36m (proteina) e
CGenFF 5.0 via ParamChem (ligante), em agua TIP3P explicita e NaCl 0,15 M (condicoes
fisiologicas humanas). A pose de docking se manteve razoavelmente estavel ao longo dos
100 ns (RMSD backbone 0,410 ± 0,071 nm, o maior desvio entre os 4 sistemas avaliados),
com 4 dos 5 residuos-chave do docking mantendo distancia proxima da predicao original —
excecao para Tyr398, que se afastou (3,9 Å previsto → 6,01 ± 1,31 Å na simulacao).

## 1. Introducao

A monoamina oxidase B (MAO-B) e uma flavoenzima mitocondrial (cofator FAD covalente)
que degrada monoaminas, alvo terapeutico estabelecido em doenca de Parkinson e outras
neurodegeneracoes. A formononetina, isoflavona O-metilada, tem atividade
inibitoria de MAO-B **ja quantificada experimentalmente**: extraida das raizes de
*Sophora flavescens*, inibe MAO cerebral de camundongo com IC50 = 11,0 μM contra
MAO-B e 21,2 μM contra MAO-A, indicando seletividade moderada por MAO-B (Hwang *et
al.*, 2005, *Arch Pharm Res* 28(2):190-4, PMID 15789750, DOI 10.1007/BF02977714).
Um estudo *in silico* independente, sobre isoflavonoides do genero *Ononis* (a
formononetina deriva seu nome desse genero), confirma predicao favoravel de potencia
inibitoria de MAO-B e permeabilidade a barreira hematoencefalica (BBB-PAMPA) para a
formononetina especificamente (Gampe *et al.*, 2022, *PLoS One* 17(3):e0265639, PMID
35298568, DOI 10.1371/journal.pone.0265639) — metodologicamente proximo da abordagem
docking+MD usada aqui, embora sem MD explicita. Nenhum dos dois estudos usa a
estrutura cristalografica PDB 6FWC especificamente; este trabalho avalia
computacionalmente, por dinamica molecular classica, se a pose de docking predita
nessa estrutura (cavidade proxima ao FAD) e estruturalmente estavel ao longo do
tempo.

## 2. Metodologia

### 2.1 Preparacao do complexo

A estrutura inicial do receptor foi obtida do PDB 6FWC (cadeia A), com os estados
de protonacao dos residuos ionizaveis determinados para pH 7,4 (condicao fisiologica humana —
nao o pH 8,2 usado nos demais pipelines deste laboratorio, especifico para midgut alcalino de
Lepidoptera) via PROPKA, implementado por `pdb2pqr 3.7.1` com campo de forca CHARMM. A pose
inicial do ligante (resname UNL) foi obtida por docking molecular com AutoDock Vina, com
interacoes-chave identificadas por analise pos-docking em: Ile199, Tyr398, Tyr326, Tyr435, Cys172.

**Limitacao metodologica assumida:** o receptor usado (extraido do docking) nao
inclui o cofator FAD, covalentemente ligado na MAO-B nativa e estruturalmente
proximo ao sitio de interesse. Decisao consciente de simplificacao para esta
primeira rodada exploratoria — resultados de estabilidade/energia devem ser
interpretados como referentes a um modelo tipo apo (sem cofator), nao a enzima
holoproteica completa.

A topologia do ligante foi gerada a partir do arquivo de parametros CGenFF 5.0 retornado
pelo servidor ParamChem (ver `inputs/ligand-UNL*.str` para as penalidades de parametro/carga
especificas deste ligante — penalidades acima de 50 indicam analogia pobre e requerem
validacao adicional segundo a propria ferramenta), convertido para o formato GROMACS com
`cgenff_charmm2gmx.py` (Lemkul Lab) e o port CHARMM36 `charmm36-feb2026_cgenff-5.0.ff`
(Wacha & Lemkul, JCIM 2023).

### 2.2 Campo de forca e parametros de simulacao

As simulacoes foram conduzidas com GROMACS 2026 (Abraham *et al.*, 2015), campo de forca
CHARMM36m (Huang *et al.*, 2017) para a proteina e CGenFF 5.0 (Vanommeslaeghe *et al.*,
2010) para o ligante, agua TIP3P explicita (Jorgensen *et al.*, 1983; parametrizacao
CHARMM-modificada). Nao-ligados seguiram a recomendacao oficial CHARMM36 para GROMACS:
`vdwtype = Cut-off` com `vdw-modifier = Force-switch` (`rvdw-switch = 1,0 nm`,
`rvdw = 1,2 nm`), sem correcao de dispersao de longo alcance (`DispCorr = no`) —
configuracao distinta do template AMBER99SB-ILDN usado nos demais pipelines deste
laboratorio. O complexo foi inserido em caixa cubica com margem minima de 1,2 nm,
solvatado e neutralizado com NaCl a 0,15 M (Joung & Cheatham, 2008), refletindo o
ambiente ionico fisiologico humano (em vez do KCl 0,10 M usado nos sistemas de
Lepidoptera deste laboratorio).

### 2.3 Protocolo de equilibracao e producao

1. **Minimizacao de energia** — *steepest descent*, `emtol = 1000 kJ mol⁻¹ nm⁻¹`, ate 50.000 passos.
2. **NVT (200 ps)** — 300 K, termostato V-rescale (Bussi *et al.*, 2007, τ = 0,1 ps),
   com restricoes de posicao no receptor (`POSRES`, gerado por `pdb2gmx`) e no ligante
   (`POSRES_UNL`, gerado por `gmx genrestr`).
3. **NPT (500 ps)** — 300 K / 1 bar, barostato de Berendsen (τ = 2,0 ps), restricoes mantidas.
4. **Producao (100 ns)** — sem restricoes, barostato de Parrinello-Rahman (Parrinello &
   Rahman, 1981; τ = 2,0 ps), integrador *leap-frog* (dt = 2 fs), ligacoes com hidrogenio
   restringidas por LINCS (Hess *et al.*, 1997), eletrostatica de longo alcance por PME
   (Darden *et al.*, 1993, `rcoulomb = 1,2 nm`).

### 2.4 Analises

RMSD do backbone do receptor e do ligante, RMSF por residuo, raio de giro, contatos
receptor-ligante (< 0,4 nm), pontes de hidrogenio, SASA do receptor e do ligante, e
distancia minima entre o ligante e os residuos de interesse identificados no docking
(Ile199, Tyr398, Tyr326, Tyr435, Cys172), todas calculadas com ferramentas nativas do GROMACS sobre a
trajetoria pos-processada (`-pbc mol -center` + `-fit rot+trans`).

### 2.5 Energia livre de ligacao (MM-GBSA)

A energia livre de ligacao foi estimada por MM-GBSA (`gmx_MMPBSA`, protocolo de
trajetoria unica, `igb=2`, decomposicao por residuo habilitada) sobre os frames da
producao pos-equilibracao. **Nota metodologica:** esta mesma ferramenta falhou de
forma irreconciliavel em outro projeto deste laboratorio (Milena-MD, serie
trypsin×GORE12T) apos 3 tentativas de correcao; o modulo aqui foi reescrito do zero
evitando o erro de linha de comando identificado retroativamente (flags `-cs/-ct/-ci`
sem valor, causando deteccao falsa de "argumentos duplicados"). Tratar resultados de
MM-GBSA como suplementares — se a etapa falhar, o restante do pipeline (RMSD/RMSF/
contatos/H-bonds/SASA) permanece valido e completo.

## 3. Resultados e Discussao

### 3.1 Estabilidade estrutural

| Metrica | Valor (media ± DP) |
|---|---|
| RMSD backbone receptor | 0.410 ± 0.071 nm |
| RMSD ligante (UNL) | 0.050 ± 0.009 nm |
| Raio de giro (receptor) | 2.358 ± 0.011 nm |
| Contatos receptor-ligante (<0,4nm) | 292.2 ± 29.4 |
| Pontes de hidrogenio receptor-ligante | 0.91 ± 0.49 |
| SASA receptor | 225.89 ± 3.48 nm² |
| SASA ligante | 4.84 ± 0.22 nm² |

### 3.2 Persistencia dos contatos preditos por docking

| Residuo | Distancia docking | Tipo (docking) | Distancia media MD |
|---|---|---|---|
| Ile199 | 3.7 Å | Hidrofobica | 2.32 ± 0.19 Å |
| Tyr398 | 3.9 Å | Hidrofobica | 6.01 ± 1.31 Å |
| Tyr326 | 4.6 Å | Hidrofobica | 2.62 ± 0.23 Å |
| Tyr435 | 4.4 Å | Hidrofobica | 2.61 ± 0.43 Å |
| Cys172 | 4.8 Å | Hidrofobica | 2.53 ± 0.37 Å |

### 3.3 Energia livre de ligacao (MM-GBSA)

ΔG total: **-22.22 ± 3.03 kcal/mol** (GB, 100 frames uniformemente distribuidos
ao longo dos 100 ns, decomposicao por residuo habilitada — `results/6FWC-formononetina/mmgbsa/`)

*(resultados ainda nao gerados — rodar o pipeline Nextflow completo)*

### 3.4 Convergencia com a literatura e outros projetos do laboratorio — TODO

Pendente, a preencher **apos** a producao terminar e as analises rodarem (nao
fabricar numeros de terceiros aqui — buscar e citar explicitamente):

- [ ] Comparar RMSD/RMSF obtidos com faixas tipicas reportadas para Monoamina oxidase B (MAO-B) em MD
      (buscar literatura especifica antes de citar valores).
- [x] Formononetina x MAO-B: Hwang et al. 2005 (PMID 15789750) mede IC50 = 11,0 μM
      (MAO-B) vs 21,2 μM (MAO-A), inibicao seletiva confirmada experimentalmente;
      Gampe et al. 2022 (PMID 35298568) confirma predicao in silico de potencia
      MAO-B e permeabilidade BBB. Nenhum dos dois fez docking/MD na estrutura
      PDB 6FWC especificamente — comparar modo de ligacao/residuos-chave quando
      os resultados desta corrida sairem.
- [ ] Comparar robustez metodologica (protocolo de equilibracao, cutoffs, forca de
      POSRES, tempo de producao) com os pipelines ja validados deste laboratorio
      (MD-gromacs serie GORE4/SKTI/BEN, Milena-MD serie trypsin×GORE12T) —
      ver `~/.claude/.claude/agents/bioinformatics.md`.
- [ ] Avaliar se a persistencia de Ile199, Tyr398, Tyr326, Tyr435, Cys172 ao longo da producao confirma ou
      refuta a pose de docking original (criterio sugerido: manter contato em
      >50% dos frames pos-equilibracao).
- [ ] Conferir penalidade CGenFF do ligante (`inputs/ligand-UNL*.str`) — se acima de 50,
      considerar validacao adicional dos dihedros de maior penalidade antes de
      conclusoes quantitativas fortes sobre energia de ligacao.

---
*Nao passou por /humanizer. Revisar citacoes com a skill auditing-academic-sources
antes de qualquer uso em documento final.*
