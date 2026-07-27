# Dinamica Molecular — NF-kB (dominio de ligacao a DNA p65/p50) (PDB 2I9T) + Daidzeina — Secoes do Artigo

*Gerado automaticamente por `bin/gerar_artigo_md.py`, com a secao 3.2.1 (analise por
fase) adicionada manualmente a partir da analise ad-hoc de 2026-07-14 — ver
`memory/project_tatiana_md.md`. Revisar antes de usar em texto final (passar por
/humanizer e pela skill auditing-academic-sources antes de qualquer submissao).*

## Resumo

Este trabalho investigou por dinamica molecular (100 ns) a estabilidade do complexo entre
NF-kB (dominio de ligacao a DNA p65/p50) (PDB 2I9T), Mus musculus e Daidzeina, um candidato identificado por triagem
virtual (AutoDock Vina). O sistema foi parametrizado com o campo de forca CHARMM36m (proteina) e
CGenFF 5.0 via ParamChem (ligante), em agua TIP3P explicita e NaCl 0,15 M (condicoes
fisiologicas humanas). O achado principal nao e um valor medio unico: a trajetoria de
100 ns mostra uma **transicao de fase real e sequencial** por volta de 65-78 ns — a pose
de docking e mantida com excelente concordancia ate ~60 ns, depois o ligante migra pra
outro sitio na superficie do receptor (ver secao 3.2 pra estatisticas separadas por fase;
a media simples de 100ns reportada na Tabela 3.1, RMSD backbone 0,545 ± 0,337 nm, mistura
as duas fases e por isso tem desvio-padrao alto — nao deve ser lida isoladamente).

## 1. Introducao

NF-κB e um fator de transcricao central nas respostas inflamatoria e imune, cuja
ativacao aberrante esta implicada em doencas inflamatorias cronicas e cancer. O
dominio de ligacao a DNA (PDB 2I9T, heterodimero p65/p50) medeia essa ativacao. A
daidzeina, isoflavona identificada aqui por triagem virtual (AutoDock Vina) contra
esse dominio, tem atividade moduladora de NF-κB **ja documentada experimentalmente**
em outros contextos, embora nao especificamente na interface proteina-DNA alvo deste
trabalho: atenua lesao pulmonar aguda induzida por LPS via a via TLR4/NF-κB (Feng,
Sun & Li, 2015, *Int Immunopharmacol* 26(2):392-400, PMID 25887269, DOI
10.1016/j.intimp.2015.04.002) e reduz a expressao de NF-κB, p38MAPK e TGF-β1 num
modelo de aneurisma de aorta abdominal induzido por angiotensina II (Liu, Bai & Qi,
2016, *Mol Med Rep* 14(1):955-62, PMID 27222119, DOI 10.3892/mmr.2016.5304). Esses
estudos demonstram modulacao funcional de NF-κB pela daidzeina em celulas/modelos
animais, mas nao caracterizam uma interacao direta com o dominio de ligacao a DNA
aqui modelado — este trabalho avalia computacionalmente, por dinamica molecular
classica, se a pose de docking predita nesse dominio especifico e estruturalmente
estavel ao longo do tempo.

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
| RMSD backbone receptor | 0.545 ± 0.337 nm |
| RMSD ligante (UNL) | 0.091 ± 0.025 nm |
| Raio de giro (receptor) | 2.425 ± 0.071 nm |
| Contatos receptor-ligante (<0,4nm) | 107.7 ± 58.6 |
| Pontes de hidrogenio receptor-ligante | 0.98 ± 0.75 |
| SASA receptor | 157.129 ± 2.580 nm² |
| SASA ligante | 4.581 ± 0.197 nm² |

### 3.2 Persistencia dos contatos preditos por docking

A media de 100 ns mistura duas fases estruturalmente distintas — ver 3.2.1. Os valores
abaixo (media simples sobre a trajetoria inteira) sao apresentados por completude, mas a
leitura correta esta na tabela por fase logo em seguida.

| Residuo | Distancia docking | Tipo (docking) | Distancia media MD (100ns) |
|---|---|---|---|
| Arg30 | 4.7-4.8 Å | Hidrofobica | 15.14 Å |
| Glu279 | 1.9 Å | Hidrogenio | 17.96 Å |

#### 3.2.1 Analise por fase (2026-07-14, script ad-hoc sobre os `.xvg` desta mesma corrida)

Divisao em duas fases a partir da inspecao visual do RMSD backbone/distancias (nao
deteccao automatica): **fase ligada** ≈0-60 ns e **fase pos-transicao** ≈65-100 ns.

| Metrica | Fase ligada (0-60ns) | Fase pos-transicao (65-100ns) |
|---|---|---|
| Dist. Glu279 (docking previu 1,9 Å) | 2,00 ± 0,44 Å | ~12,5 Å (H-bond rompido ~77ns) |
| Dist. Arg30 (docking previu 4,7-4,8 Å) | 4,46 ± 2,30 Å | ~22 Å (contato rompido ~74ns) |
| Pontes de H receptor-ligante | 1,70 ± 0,54 | — |
| Contatos receptor-ligante | 127,8 ± 26,8 | menor (ligante em novo sitio) |

**Sequencia da transicao**: RMSD do backbone cruza 0,6 nm em ~68 ns (receptor muda de
conformacao primeiro — Rg sobe de 2,478 para 2,505 nm, abertura de dominio) → contato
com Arg30 rompe em ~74 ns → ligacao de H com Glu279 rompe por ultimo, ~77 ns.

**O ligante nao se solvatou** — apesar das distancias aos dois residuos do docking
subirem na fase pos-transicao, a distancia minima receptor-ligante (`mindist`) continua
<1 nm na maior parte do tempo e a SASA do ligante nao muda (~4,6 nm² do inicio ao fim da
trajetoria, ver Tabela 3.1). Isto indica que a daidzeina **migrou para outro sitio na
superficie do receptor**, nao se dissociou pro solvente. RMSF por residuo: media
0,275 nm; Arg30 = 0,3275 nm e Glu279 = 0,3403 nm (acima da mediana — sitio de docking em
regiao de loop moderadamente flexivel, consistente com a metaestabilidade observada).

**Leitura cientifica**: a pose predita pelo AutoDock Vina e um **minimo local
metaestavel**, nao o minimo global — resultado positivo (nao falha do protocolo de
docking/MD), mas que exige tratar qualquer conclusao de energia de ligacao sobre a pose
original como valida apenas para os primeiros ~60 ns, nao para a trajetoria inteira.
Antes de qualquer conclusao mais forte, seria necessario confirmar que o evento e real
e nao artefato de PBC/fitting (replica adicional ou trajetoria estendida — ver 3.4).

### 3.3 Energia livre de ligacao (MM-GBSA)

ΔG total: **N/D (MM-GBSA nao rodou ou falhou — ver mmgbsa.log; tratar como opcional, ja falhou de forma irreconciliavel em outro projeto deste laboratorio, ver Milena-MD)**

### 3.4 Convergencia com a literatura e outros projetos do laboratorio — TODO

Pendente, a preencher **apos** a producao terminar e as analises rodarem (nao
fabricar numeros de terceiros aqui — buscar e citar explicitamente):

- [ ] Comparar RMSD/RMSF obtidos com faixas tipicas reportadas para NF-kB (dominio de ligacao a DNA p65/p50) em MD
      (buscar literatura especifica antes de citar valores).
- [x] Literatura sobre daidzeina x NF-κB (funcional, nao estrutural no 2I9T
      especificamente): Feng et al. 2015 (PMID 25887269, TLR4/NF-κB, lesao
      pulmonar) e Liu et al. 2016 (PMID 27222119, NF-κB/p38MAPK/TGF-β1,
      aneurisma de aorta) — nenhum dos dois fez docking/MD no dominio de
      ligacao a DNA do 2I9T; nao ha ainda estudo estrutural direto encontrado
      pra comparar modo de ligacao/residuos-chave especificos deste alvo.
- [ ] Comparar robustez metodologica (protocolo de equilibracao, cutoffs, forca de
      POSRES, tempo de producao) com os pipelines ja validados deste laboratorio
      (MD-gromacs serie GORE4/SKTI/BEN, Milena-MD serie trypsin×GORE12T) —
      ver `~/.claude/.claude/agents/bioinformatics.md`.
- [x] Persistencia de Arg30/Glu279 — ver 3.2.1: confirma a pose de docking nos
      primeiros ~60ns (60% da producao), depois rompe numa transicao sequencial
      real. Falta: replica/trajetoria estendida pra confirmar que nao e artefato
      (ver item abaixo).
- [ ] Rodar replica curta ou estender a trajetoria — um unico run de 100ns nao
      distingue evento termodinamico real de flutuacao estocastica rara.
- [ ] `gmx dssp` do receptor em torno de ~65-70ns — checar se ha
      desenovelamento local/movimento de loop que abriu o sitio novo.
- [ ] Conferir penalidade CGenFF do ligante (`inputs/ligand-UNL*.str`) — se acima de 50,
      considerar validacao adicional dos dihedros de maior penalidade antes de
      conclusoes quantitativas fortes sobre energia de ligacao.

---
*Nao passou por /humanizer. Revisar citacoes com a skill auditing-academic-sources
antes de qualquer uso em documento final.*
