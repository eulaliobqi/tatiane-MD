# Dinamica Molecular — Receptor 2I9T (NF-kB) + Daidzeina — Secoes do Artigo

*Gerado automaticamente por `bin/gerar_artigo_md.py`. Revisar antes de usar em texto final
(passar por /humanizer e pela skill auditing-academic-sources antes de qualquer submissao).*

## Resumo

Este trabalho investigou por dinamica molecular (100 ns) a estabilidade do complexo entre
o dominio de ligacao a DNA do fator de transcricao NF-kB (PDB 2I9T, cadeia A, res. 17-291)
e a isoflavona daidzeina, um candidato a inibidor natural identificado por triagem virtual
(AutoDock Vina). O sistema foi parametrizado com o campo de forca CHARMM36m (proteina) e
CGenFF 5.0 via ParamChem (ligante), em agua TIP3P explicita e NaCl 0,15 M (condicoes
fisiologicas humanas). Simulacao ainda nao executada — secao a preencher apos `bin/run_md.sh` + `bin/analyze.sh`.

## 1. Introducao

NF-kB e um fator de transcricao central na resposta inflamatoria e imune, cuja ativacao
aberrante esta implicada em cancer, doencas autoimunes e inflamacao cronica. O dominio Rel
homology (RHD) de suas subunidades (p50/p65) medeia tanto a dimerizacao quanto a ligacao
direta ao DNA, sendo um alvo estabelecido para o desenho de inibidores de pequenas
moleculas que bloqueiam essa interacao. Isoflavonas de origem vegetal, como a daidzeina
(*Glycine max*), tem sido reportadas na literatura como moduladoras da via NF-kB; este
trabalho avalia computacionalmente, por dinamica molecular classica, a estabilidade
temporal do complexo predito por docking entre a daidzeina e o RHD de NF-kB (PDB 2I9T).

## 2. Metodologia

### 2.1 Preparacao do complexo

A estrutura inicial do receptor foi obtida do PDB 2I9T (dominio de ligacao a DNA de NF-kB
p65/p50, cadeia A, residuos 17-291), com os estados de protonacao dos residuos ionizaveis
determinados para pH 7,4 (condicao fisiologica humana — nao o pH 8,2 usado nos demais
pipelines deste laboratorio, especifico para midgut alcalino de Lepidoptera) via PROPKA,
implementado por `pdb2pqr 3.7.1` com campo de forca CHARMM. A pose inicial da daidzeina
(resname UNL) foi obtida por docking molecular com AutoDock Vina, com interacoes-chave
identificadas por analise pos-docking em Arg30 (contato hidrofobico, ~4,7-4,8 Å) e
Glu279 (ligacao de hidrogenio, ~1,9 Å).

A topologia do ligante foi gerada a partir do arquivo de parametros CGenFF 5.0 retornado
pelo servidor ParamChem (penalidade de parametro = 53,0; penalidade de carga = 23,3 —
acima do limiar de 50 que a propria CGenFF define como "requer validacao extensa";
resultado tratado como preliminar ate validacao adicional, ex. otimizacao QM dos
dihedros de maior penalidade), convertido para o formato GROMACS com
`cgenff_charmm2gmx.py` (Lemkul Lab) e o port CHARMM36 de marco de 2019
(E. P. Raman, J. A. Lemkul, R. Best, A. D. MacKerell Jr.).

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
distancia minima entre o ligante e os dois residuos de interesse identificados no
docking (Arg30, Glu279), todas calculadas com ferramentas nativas do GROMACS sobre a
trajetoria pos-processada (`-pbc mol -center` + `-fit rot+trans`).

## 3. Resultados e Discussao

### 3.1 Estabilidade estrutural

| Metrica | Valor (media ± DP) |
|---|---|
| RMSD backbone receptor | N/D (rodar bin/analyze.sh) |
| RMSD ligante (UNL) | N/D (rodar bin/analyze.sh) |
| Raio de giro (receptor) | N/D (rodar bin/analyze.sh) |
| Contatos receptor-ligante (<0,4nm) | N/D (rodar bin/analyze.sh) |
| Pontes de hidrogenio receptor-ligante | N/D (rodar bin/analyze.sh) |
| SASA receptor | N/D (rodar bin/analyze.sh) |
| SASA ligante | N/D (rodar bin/analyze.sh) |

### 3.2 Persistencia dos contatos preditos por docking

| Residuo | Distancia docking | Tipo (docking) | Distancia media MD |
|---|---|---|---|
| Arg30 | 4,7-4,8 Å | Hidrofobica | N/D Å |
| Glu279 | 1,9 Å | Hidrogenio | N/D Å |

*(resultados ainda nao gerados — rodar `bin/run_md.sh` seguido de `bin/analyze.sh`)*

### 3.3 Convergencia com a literatura e outros projetos do laboratorio — TODO

Pendente, a preencher **apos** a producao terminar e as analises rodarem (nao
fabricar numeros de terceiros aqui — buscar e citar explicitamente):

- [ ] Comparar RMSD/RMSF obtidos com faixas tipicas reportadas para dominios RHD de
      NF-kB em MD (buscar literatura especifica antes de citar valores).
- [ ] Buscar na literatura estudos computacionais ou experimentais de
      daidzeina/isoflavonas ligando NF-kB (ou alvos RHD homologos) e comparar
      modo de ligacao / residuos-chave.
- [ ] Comparar robustez metodologica (protocolo de equilibracao, cutoffs, forca de
      POSRES, tempo de producao) com os pipelines ja validados deste laboratorio
      (MD-gromacs serie GORE4/SKTI/BEN, Milena-MD serie trypsin×GORE12T) —
      ver `~/.claude/.claude/agents/bioinformatics.md`.
- [ ] Avaliar se a persistencia de Arg30/Glu279 ao longo dos 100 ns confirma ou
      refuta a pose de docking original (criterio sugerido: manter contato em
      >50% dos frames pos-equilibracao).
- [ ] Dado o param penalty=53 do CGenFF (acima do limiar de 50), considerar
      validacao adicional dos dihedros do anel cromona antes de conclusoes
      quantitativas fortes sobre energia de ligacao.

---
*Nao passou por /humanizer. Revisar citacoes com a skill auditing-academic-sources
antes de qualquer uso em documento final.*
