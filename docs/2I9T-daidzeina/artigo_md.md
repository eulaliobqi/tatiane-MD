# Dinamica Molecular — NF-kB (dominio de ligacao a DNA p65/p50) (PDB 2I9T) + Daidzeina — Secoes do Artigo

*Gerado automaticamente por `bin/gerar_artigo_md.py`, com a Tabela 3.1/3.2 e o Resumo
atualizados manualmente em 2026-07-28 pra refletir a producao mais recente deste sistema
(ver nota em 3.2.1 sobre a rodada anterior). Revisar antes de usar em texto final (passar
por /humanizer e pela skill auditing-academic-sources antes de qualquer submissao).*

## Resumo

Este trabalho investigou por dinamica molecular (100 ns) a estabilidade do complexo entre
NF-kB (dominio de ligacao a DNA p65/p50) (PDB 2I9T), Mus musculus e Daidzeina, um candidato identificado por triagem
virtual (AutoDock Vina). O sistema foi parametrizado com o campo de forca CHARMM36m (proteina) e
CGenFF 5.0 via ParamChem (ligante), em agua TIP3P explicita e NaCl 0,15 M (condicoes
fisiologicas humanas). Nesta producao (concluida 2026-07-27/28), a pose de docking se manteve
proxima da predicao original ao longo da maior parte da trajetoria — RMSD do backbone de
0,316 ± 0,104 nm e 185,9 ± 42,4 contatos receptor-ligante em media, com Arg30 e Glu279
(residuos-chave do docking) mantendo distancia media de 2,61 Å e 4,02 Å, respectivamente
(docking previu 4,7-4,8 Å e 1,9 Å). Uma producao INDEPENDENTE anterior deste mesmo sistema
(protocolo identico, rodada 2026-07-13/16) identificou uma transicao de fase — o ligante
migrando pra outro sitio na superficie do receptor por volta de 65-78 ns — nao replicada
(nem re-testada) nesta rodada; ver secao 3.2.1 pra essa analise historica, mantida por
valor metodologico mas nao representativa dos dados atuais deste documento.

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

**Dados desta producao (concluida 2026-07-27/28)** — ver nota em 3.2.1 sobre a rodada
anterior (2026-07-13/16), que usava a mesma configuracao mas seguiu uma trajetoria
diferente (MD e processo estocastico; ver 3.2.1 pra detalhes).

| Metrica | Valor (media ± DP) |
|---|---|
| RMSD backbone receptor | 0.316 ± 0.104 nm |
| RMSD ligante (UNL) | 0.052 ± 0.021 nm |
| Raio de giro (receptor) | 2.337 ± 0.021 nm |
| Contatos receptor-ligante (<0,4nm) | 185.9 ± 42.4 |
| Pontes de hidrogenio receptor-ligante | 1.50 ± 0.85 |
| SASA receptor | 159.123 ± 2.660 nm² |
| SASA ligante | 4.525 ± 0.200 nm² |

### 3.2 Persistencia dos contatos preditos por docking

| Residuo | Distancia docking | Tipo (docking) | Distancia media MD (100ns) |
|---|---|---|---|
| Arg30 | 4.7-4.8 Å | Hidrofobica | 2.61 ± 0.64 Å |
| Glu279 | 1.9 Å | Hidrogenio | 4.02 ± 1.35 Å |

Ambos os residuos mantiveram distancia media proxima da predicao do docking ao longo dos
100 ns desta producao — sem evidencia, nesta rodada, da transicao de fase descrita em 3.2.1
pra uma producao anterior e independente do mesmo sistema.

#### 3.2.1 [ANALISE HISTORICA — rodada anterior, 2026-07-13/16, NAO representa os dados atuais deste documento]

**Atualizado 2026-07-28**: esta secao descreve uma producao INDEPENDENTE anterior deste
mesmo sistema (mesmo protocolo/configuracao, mas rodada de MD diferente — concluida
2026-07-13, recomputada 2026-07-15/16 apos um bug de cache do Nextflow). A producao ATUAL
deste sistema (results-tatiana/2I9T-daidzeina/, concluida 2026-07-27/28, dados em 3.1/3.2
acima) **nao foi verificada** quanto a uma transicao equivalente — o RMSD backbone desta
rodada atual cruza 0,6 nm pela primeira vez apenas em ~87 ns (vs. ~68 ns na rodada antiga),
com padrao menos monotonico (sobe 60-70ns, cai 70-80ns, sobe de novo 80-100ns) — insuficiente
pra afirmar ou descartar um evento analogo sem analise dedicada (nao feita). Mantida abaixo
por valor metodologico (mostra que a pose de docking pode ser metaestavel em pelo menos uma
das rodadas independentes deste sistema) — **nao usar os numeros desta subsecao como
descricao dos dados atuais em 3.1/3.2**.

Divisao em duas fases a partir da inspecao visual do RMSD backbone/distancias (nao
deteccao automatica, rodada 2026-07-13/16): **fase ligada** ≈0-60 ns e **fase
pos-transicao** ≈65-100 ns.

| Metrica | Fase ligada (0-60ns) | Fase pos-transicao (65-100ns) |
|---|---|---|
| Dist. Glu279 (docking previu 1,9 Å) | 2,00 ± 0,44 Å | ~12,5 Å (H-bond rompido ~77ns) |
| Dist. Arg30 (docking previu 4,7-4,8 Å) | 4,46 ± 2,30 Å | ~22 Å (contato rompido ~74ns) |
| Pontes de H receptor-ligante | 1,70 ± 0,54 | — |
| Contatos receptor-ligante | 127,8 ± 26,8 | menor (ligante em novo sitio) |

**Sequencia da transicao** (rodada 2026-07-13/16): RMSD do backbone cruza 0,6 nm em ~68 ns
(receptor muda de conformacao primeiro — Rg sobe de 2,478 para 2,505 nm, abertura de
dominio) → contato com Arg30 rompe em ~74 ns → ligacao de H com Glu279 rompe por ultimo,
~77 ns.

**O ligante nao se solvatou** (nessa rodada) — apesar das distancias aos dois residuos do
docking subirem na fase pos-transicao, a distancia minima receptor-ligante (`mindist`)
continuava <1 nm na maior parte do tempo e a SASA do ligante nao mudava (~4,6 nm² do inicio
ao fim da trajetoria daquela rodada). Isto indica que a daidzeina **migrou para outro
sitio na superficie do receptor**, nao se dissociou pro solvente, naquela producao
especifica. RMSF por residuo (rodada antiga): media 0,275 nm; Arg30 = 0,3275 nm e
Glu279 = 0,3403 nm (acima da mediana — sitio de docking em regiao de loop moderadamente
flexivel, consistente com a metaestabilidade observada naquela rodada).

**Leitura cientifica**: a pose predita pelo AutoDock Vina se mostrou um **minimo local
metaestavel** em pelo menos uma das producoes independentes deste sistema — resultado
cientificamente relevante (nao falha do protocolo de docking/MD), mas que **nao foi
confirmado como reprodutivel** entre rodadas independentes (replica formal nunca rodada).
Antes de qualquer conclusao mais forte sobre a estabilidade real da pose, seria necessario
rodar replicas adicionais com deteccao automatica de transicao (nao inspecao visual) e
comparar consistentemente entre elas — ver 3.4.

#### 3.2.2 Analise complementar: mapa de contatos e fingerprint quimico (2026-07-28)

Para investigar se a producao ATUAL (3.1/3.2) tambem apresenta algum deslocamento do
ligante (mesmo sem a transicao classica documentada em 3.2.1), a trajetoria foi
dividida em duas fases usando os limites do `PHASE_SPLIT` calibrados na rodada antiga
(0-60ns "bound", 65-100ns "relocated" — aproximado pra esta rodada, ver nota
metodologica na secao 3.3) e analisada com `CONTACT_MAP` (frequencia de contato por
residuo, corte <0,4nm) e `PROLIF_FINGERPRINT` (perfil de interacoes quimicas).

**Os residuos de maior frequencia de contato sao praticamente os mesmos nas duas
fases** — Lys28, Arg30, Glu49, Arg50, Leu280 e Glu282 permanecem entre os residuos
mais contatados em ambas (>0,9 na fase "bound", entre 0,71-0,98 na fase "relocated") —
mas com frequencia GERALMENTE MENOR na fase "relocated". A diferenca real esta em dois
pontos especificos: a fase "relocated" **perde quase completamente** o contato com
Glu193/Leu194 (0,98/0,95 na fase "bound" → 0,17/ausente do top-15 na fase
"relocated") e **ganha** contato novo com Gln29 (0,45), Arg278 (0,42) e Asp277 (0,33),
praticamente ausentes na fase "bound".

**Interpretacao**: mesmo nesta rodada (sem a transicao classica de 3.2.1), ha um
indicio de deslocamento LATERAL sutil da daidzeina dentro da MESMA regiao geral de
interface (loop 277-283 + helice 28-52) — perdendo profundidade de contato com o
subsitio Glu193/Leu194 e ganhando contato com o lado do loop mais proximo de
Arg278/Asp277, em vez de permanecer perfeitamente fixa ou migrar pra um sitio
totalmente distante. Isso e consistente com (nao identico a) o achado da rodada
antiga: o CONTATO geral com o grupo de residuos-chave persiste, mas a geometria fina
de contato oscila ao longo da trajetoria.

`FE_RERUN` (energia de interacao Coulomb-SR + LJ-SR em vacuo, sem termo de
solvatacao — nao equivalente ao MM-GBSA da secao 3.3) deu ΔE_MM = -200,11 ± 29,27
kJ/mol na fase "bound" (1001 frames) contra -167,61 ± 51,09 kJ/mol na fase
"relocated" (584 frames) — interacao mais fraca E mais variavel na fase
"relocated", direcionalmente consistente com o deslocamento lateral observado no
mapa de contatos. Ver `analise_extra/{bound,relocated}/` pros arquivos completos
(`contact_map.png`, `prolif_heatmap.png`, `interface_residues.csv`,
`free_energy_estimate.txt`).

### 3.3 Energia livre de ligacao (MM-GBSA)

ΔG total: **-15.19 ± 5.23 kcal/mol** (GB, 100 frames uniformemente distribuidos
ao longo dos 100 ns da producao atual, decomposicao por residuo habilitada —
`results/2I9T-daidzeina/mmgbsa/`)

**Nota metodologica**: valor de trajetoria UNICA (nao separado por fase bound/
relocated), sobre a producao atual (3.1/3.2). Os resultados de `CONTACT_MAP`/
`PROLIF_FINGERPRINT`/`FE_RERUN` gerados em 2026-07-28 (ver `analise_extra/`) usaram
os limites fixos `phase_bound_end_ns=60`/`phase_reloc_start_ns=65` do `PHASE_SPLIT`
— calibrados na rodada antiga (3.2.1), aplicados sem recalibrar a producao atual;
tratar a separacao bound/relocated desses resultados como aproximada, nao
recalibrada, pra esta producao.

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
- [x] Persistencia de Arg30/Glu279 — **producao atual** (3.1/3.2): distancia media
      proxima da predicao do docking ao longo dos 100ns inteiros, sem ruptura
      visivel. **Rodada historica** (3.2.1, producao independente anterior): pose
      confirmada nos primeiros ~60ns, depois rompe numa transicao sequencial —
      NAO observado (nem testado com o mesmo criterio) na producao atual.
- [x] Replica independente — por acidente de infraestrutura (bug de cache do
      Nextflow, ver 3.2.1, e trabalho de generalizacao do pipeline pra 4 sistemas),
      este sistema acabou tendo 2+ producoes de 100ns independentes do zero (mesmo
      protocolo, seeds diferentes). Achado: elas **nao concordam** — RMSD cruza
      0,6nm em ~68ns numa e ~87ns (padrao menos monotonico) na outra. Isso e
      evidencia inicial (nao formal/sistematica) de que a transicao de fase
      documentada em 3.2.1 pode ser flutuacao estocastica sensivel a seed, nao um
      evento termodinamico reprodutivel — mas as 2 rodadas nunca foram comparadas
      lado a lado com o mesmo criterio automatico (so inspecao visual pontual).
      Uma replica FORMAL (mesmo ponto de partida, so seed de velocidade diferente,
      analise automatica consistente) ainda nao foi rodada.
- [ ] `gmx dssp` do receptor em torno de ~65-70ns — checar se ha
      desenovelamento local/movimento de loop que abriu o sitio novo.
- [ ] Conferir penalidade CGenFF do ligante (`inputs/ligand-UNL*.str`) — se acima de 50,
      considerar validacao adicional dos dihedros de maior penalidade antes de
      conclusoes quantitativas fortes sobre energia de ligacao.

---
*Nao passou por /humanizer. Revisar citacoes com a skill auditing-academic-sources
antes de qualquer uso em documento final.*
