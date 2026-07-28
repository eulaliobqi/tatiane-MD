# Dinamica Molecular — NF-kB (dominio de ligacao a DNA p65/p50) (PDB 2I9T) + Daidzeina — Secoes do Artigo

*Gerado automaticamente por `bin/gerar_artigo_md.py`. Revisar antes de usar em texto final
(passar por /humanizer e pela skill auditing-academic-sources antes de qualquer submissao).*

## Resumo

Este trabalho investigou por dinamica molecular (100 ns) a estabilidade do complexo entre
NF-kB (dominio de ligacao a DNA p65/p50) (PDB 2I9T), Mus musculus e Daidzeina, um candidato identificado por triagem
virtual (AutoDock Vina). O sistema foi parametrizado com o campo de forca CHARMM36m (proteina) e
CGenFF 5.0 via ParamChem (ligante), em agua TIP3P explicita e NaCl 0,15 M (condicoes
fisiologicas humanas). Resultados preliminares indicam RMSD do backbone de 0.316 ± 0.104 nm e 186 contatos receptor-ligante em media.

## 1. Introducao

NF-kB (dominio de ligacao a DNA p65/p50) (PDB 2I9T), Mus musculus foi selecionado como alvo de triagem virtual para
Daidzeina, com a pose de docking (AutoDock Vina) avaliada quanto a estabilidade temporal
por dinamica molecular classica. *(Secao a expandir com contexto biologico especifico do alvo
e revisao da literatura sobre o ligante — ver checklist na secao 3.4; nenhuma afirmacao sobre
relevancia biologica ou precedente na literatura deve ser incluida aqui sem verificacao
explicita, ver skill auditing-academic-sources.)*

## 2. Metodologia

### 2.1 Preparacao do complexo

A estrutura inicial do receptor foi obtida do PDB 2I9T (cadeia A), com os estados
de protonacao dos residuos ionizaveis determinados para pH 7,4 (condicao fisiologica humana —
nao o pH 8,2 usado nos demais pipelines deste laboratorio, especifico para midgut alcalino de
Lepidoptera) via PROPKA, implementado por `pdb2pqr 3.7.1` com campo de forca CHARMM. A pose
inicial do ligante (resname UNL) foi obtida por docking molecular com AutoDock Vina, com
interacoes-chave identificadas por analise pos-docking em: Arg30, Glu279.

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
(Arg30, Glu279), todas calculadas com ferramentas nativas do GROMACS sobre a
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
| RMSD backbone receptor | 0.316 ± 0.104 nm |
| RMSD ligante (UNL) | 0.052 ± 0.021 nm |
| Raio de giro (receptor) | 2.337 ± 0.021 nm |
| Contatos receptor-ligante (<0,4nm) | 185.9 ± 42.4 |
| Pontes de hidrogenio receptor-ligante | 1.50 ± 0.85 |
| SASA receptor | 159.123 ± 2.661 nm² |
| SASA ligante | 4.525 ± 0.200 nm² |

### 3.2 Persistencia dos contatos preditos por docking

| Residuo | Distancia docking | Tipo (docking) | Distancia media MD |
|---|---|---|---|
| Arg30 | 4.7-4.8 Å | Hidrofobica | 2.61 Å |
| Glu279 | 1.9 Å | Hidrogenio | 4.02 Å |

### 3.3 Energia livre de ligacao (MM-GBSA)

ΔG total: **-15.19 ± 5.23 kcal/mol**



### 3.4 Convergencia com a literatura e outros projetos do laboratorio — TODO

Pendente, a preencher **apos** a producao terminar e as analises rodarem (nao
fabricar numeros de terceiros aqui — buscar e citar explicitamente):

- [ ] Comparar RMSD/RMSF obtidos com faixas tipicas reportadas para NF-kB (dominio de ligacao a DNA p65/p50) em MD
      (buscar literatura especifica antes de citar valores).
- [ ] Buscar na literatura estudos computacionais ou experimentais de Daidzeina
      ligando NF-kB (dominio de ligacao a DNA p65/p50) (ou alvos homologos) e comparar modo de ligacao /
      residuos-chave / valores de ΔG de ligacao.
- [ ] Comparar robustez metodologica (protocolo de equilibracao, cutoffs, forca de
      POSRES, tempo de producao) com os pipelines ja validados deste laboratorio
      (MD-gromacs serie GORE4/SKTI/BEN, Milena-MD serie trypsin×GORE12T) —
      ver `~/.claude/.claude/agents/bioinformatics.md`.
- [ ] Avaliar se a persistencia de Arg30, Glu279 ao longo da producao confirma ou
      refuta a pose de docking original (criterio sugerido: manter contato em
      >50% dos frames pos-equilibracao).
- [ ] Conferir penalidade CGenFF do ligante (`inputs/ligand-UNL*.str`) — se acima de 50,
      considerar validacao adicional dos dihedros de maior penalidade antes de
      conclusoes quantitativas fortes sobre energia de ligacao.

---
*Nao passou por /humanizer. Revisar citacoes com a skill auditing-academic-sources
antes de qualquer uso em documento final.*
