#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gera docs/artigo_md.md a partir dos resultados de results/2I9T-daidzeina/10_analysis/*.xvg
Mirror do formato usado em MD-gromacs/artigo_md.md (Resumo / Introducao / Metodologia /
Resultados e Discussao), adaptado para um unico par receptor-ligante (nao uma serie).

IMPORTANTE: a secao "Convergencia com a literatura e outros projetos" e deixada como
checklist TODO — nao preenche comparacoes com a literatura automaticamente (ver skill
auditing-academic-sources: nenhum numero de artigo de terceiros deve ser citado sem
verificacao explicita via busca).

Uso (depois que bin/run_md.sh e bin/analyze.sh tiverem rodado):
    python bin/gerar_artigo_md.py
"""
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS_DIR = ROOT / "results" / "2I9T-daidzeina" / "10_analysis"
OUT_MD = ROOT / "docs" / "artigo_md.md"


def read_xvg(path):
    """Le um .xvg do GROMACS, ignorando cabecalhos @/# -> lista de (x, y[, y2...])."""
    if not path.exists():
        return []
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(("@", "#")):
            continue
        parts = line.split()
        try:
            rows.append(tuple(float(p) for p in parts))
        except ValueError:
            continue
    return rows


def mean_sd(values):
    if not values:
        return None, None
    if len(values) == 1:
        return values[0], 0.0
    return statistics.mean(values), statistics.stdev(values)


def fmt(mean, sd, unit="", nd=3):
    if mean is None:
        return "N/D (rodar bin/analyze.sh)"
    return f"{mean:.{nd}f} ± {sd:.{nd}f} {unit}".strip()


def summarize(filename, col=1):
    rows = read_xvg(ANALYSIS_DIR / filename)
    if not rows:
        return None, None
    values = [r[col] for r in rows if len(r) > col]
    return mean_sd(values)


def last_value(filename, col=1):
    rows = read_xvg(ANALYSIS_DIR / filename)
    if not rows:
        return None
    return rows[-1][col] if len(rows[-1]) > col else None


def count_events_above(filename, threshold, col=1):
    """Conta quantos frames tem valor > threshold (ex.: quebras de contato)."""
    rows = read_xvg(ANALYSIS_DIR / filename)
    if not rows:
        return None
    return sum(1 for r in rows if len(r) > col and r[col] > threshold)


def build_report():
    rmsd_rec_mean, rmsd_rec_sd = summarize("rmsd_backbone.xvg")
    rmsd_lig_mean, rmsd_lig_sd = summarize("rmsd_ligante.xvg")
    rg_mean, rg_sd = summarize("gyrate.xvg")
    contacts_mean, contacts_sd = summarize("numcont_receptor_ligante.xvg")
    hbond_mean, hbond_sd = summarize("hbond.xvg")
    sasa_rec_mean, sasa_rec_sd = summarize("sasa_receptor.xvg")
    sasa_lig_mean, sasa_lig_sd = summarize("sasa_ligante.xvg")
    arg30_mean, arg30_sd = summarize("dist_arg30.xvg")
    glu279_mean, glu279_sd = summarize("dist_glu279.xvg")

    has_results = rmsd_rec_mean is not None

    arg30_nm = f"{arg30_mean*10:.2f}" if arg30_mean is not None else "N/D"
    glu279_nm = f"{glu279_mean*10:.2f}" if glu279_mean is not None else "N/D"

    md = f"""# Dinamica Molecular — Receptor 2I9T (NF-kB) + Daidzeina — Secoes do Artigo

*Gerado automaticamente por `bin/gerar_artigo_md.py`. Revisar antes de usar em texto final
(passar por /humanizer e pela skill auditing-academic-sources antes de qualquer submissao).*

## Resumo

Este trabalho investigou por dinamica molecular (100 ns) a estabilidade do complexo entre
o dominio de ligacao a DNA do fator de transcricao NF-kB (PDB 2I9T, cadeia A, res. 17-291)
e a isoflavona daidzeina, um candidato a inibidor natural identificado por triagem virtual
(AutoDock Vina). O sistema foi parametrizado com o campo de forca CHARMM36m (proteina) e
CGenFF 5.0 via ParamChem (ligante), em agua TIP3P explicita e NaCl 0,15 M (condicoes
fisiologicas humanas). {"Resultados preliminares indicam RMSD do backbone de " + fmt(rmsd_rec_mean, rmsd_rec_sd, "nm") + " e " + f"{contacts_mean:.0f}" + " contatos receptor-ligante em media." if has_results else "Simulacao ainda nao executada — secao a preencher apos `bin/run_md.sh` + `bin/analyze.sh`."}

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
| RMSD backbone receptor | {fmt(rmsd_rec_mean, rmsd_rec_sd, "nm")} |
| RMSD ligante (UNL) | {fmt(rmsd_lig_mean, rmsd_lig_sd, "nm")} |
| Raio de giro (receptor) | {fmt(rg_mean, rg_sd, "nm")} |
| Contatos receptor-ligante (<0,4nm) | {fmt(contacts_mean, contacts_sd, "", 1)} |
| Pontes de hidrogenio receptor-ligante | {fmt(hbond_mean, hbond_sd, "", 2)} |
| SASA receptor | {fmt(sasa_rec_mean, sasa_rec_sd, "nm²")} |
| SASA ligante | {fmt(sasa_lig_mean, sasa_lig_sd, "nm²")} |

### 3.2 Persistencia dos contatos preditos por docking

| Residuo | Distancia docking | Tipo (docking) | Distancia media MD |
|---|---|---|---|
| Arg30 | 4,7-4,8 Å | Hidrofobica | {arg30_nm} Å |
| Glu279 | 1,9 Å | Hidrogenio | {glu279_nm} Å |

{"*(resultados ainda nao gerados — rodar `bin/run_md.sh` seguido de `bin/analyze.sh`)*" if not has_results else ""}

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
"""
    return md


def main():
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text(build_report(), encoding="utf-8")
    print(f"[OK] Relatorio gerado em {OUT_MD}")
    if not (ANALYSIS_DIR / "rmsd_backbone.xvg").exists():
        print("[AVISO] Nenhum resultado de analise encontrado ainda — "
              "rode bin/run_md.sh e bin/analyze.sh antes para preencher os numeros.")


if __name__ == "__main__":
    main()
