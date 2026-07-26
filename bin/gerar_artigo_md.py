#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gera artigo_md.md a partir dos resultados de uma pasta de analise (.xvg) e,
opcionalmente, de uma pasta de MM-GBSA (FINAL_RESULTS_MMGBSA.dat).
Mirror do formato usado em MD-gromacs/artigo_md.md (Resumo / Introducao /
Metodologia / Resultados e Discussao), adaptado para um unico par
receptor-ligante (nao uma serie). Generico por sistema — nome do alvo, PDB
ID, organismo, ligante e residuos-chave sao parametros de linha de comando,
nao fabricados/assumidos aqui.

IMPORTANTE: a secao "Convergencia com a literatura e outros projetos" e deixada
como checklist TODO — nao preenche comparacoes com a literatura automaticamente
(ver skill auditing-academic-sources: nenhum numero de artigo de terceiros deve
ser citado sem verificacao explicita via busca).

Uso:
    python bin/gerar_artigo_md.py \\
        --analise-dir results/2I9T-daidzeina/analise \\
        --mmgbsa-dir  results/2I9T-daidzeina/mmgbsa \\
        --target-name "NF-kB (dominio de ligacao a DNA)" --pdb-id 2I9T \\
        --organism "Mus musculus" --ligand-name Daidzeina \\
        --key-residues "30:ARG:4.7-4.8:Hidrofobica;279:GLU:1.9:Hidrogenio" \\
        --out docs/2I9T-daidzeina/artigo_md.md
"""
import argparse
import re
import statistics
from pathlib import Path


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
        return "N/D (rodar as analises)"
    return f"{mean:.{nd}f} ± {sd:.{nd}f} {unit}".strip()


def summarize(analise_dir, filename, col=1):
    rows = read_xvg(analise_dir / filename)
    if not rows:
        return None, None
    values = [r[col] for r in rows if len(r) > col]
    return mean_sd(values)


def parse_key_residues(spec):
    """'resid:resname:dockdist:doctype;...' -> lista de dicts com label (ex. Arg30)."""
    entries = []
    for chunk in (spec or "").split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        resid, resname, dockdist, doctype = chunk.split(":")
        label = f"{resname.lower().capitalize()}{resid}"
        entries.append({"resid": resid, "resname": resname, "dockdist": dockdist,
                         "doctype": doctype, "label": label})
    return entries


def read_mmgbsa_total(mmgbsa_dir):
    """Extrai DELTA TOTAL (kcal/mol) do FINAL_RESULTS_MMGBSA.dat, se existir e valido."""
    if mmgbsa_dir is None:
        return None
    dat = mmgbsa_dir / "FINAL_RESULTS_MMGBSA.dat"
    if not dat.exists():
        return None
    text = dat.read_text(errors="ignore")
    if "No results" in text or "failed" in text.lower():
        return None
    # Formato tipico do gmx_MMPBSA: bloco "DELTA TOTAL" com media/desvio
    m = re.search(r"DELTA TOTAL\s+(-?\d+\.?\d*)\s+(\d+\.?\d*)", text)
    if m:
        return float(m.group(1)), float(m.group(2))
    return None


def build_report(analise_dir, mmgbsa_dir, time_ns, target_name, pdb_id, organism,
                  ligand_name, key_residues):
    rmsd_rec_mean, rmsd_rec_sd = summarize(analise_dir, "rmsd_backbone.xvg")
    rmsd_lig_mean, rmsd_lig_sd = summarize(analise_dir, "rmsd_ligante.xvg")
    rg_mean, rg_sd = summarize(analise_dir, "gyrate.xvg")
    contacts_mean, contacts_sd = summarize(analise_dir, "numcont_receptor_ligante.xvg")
    hbond_mean, hbond_sd = summarize(analise_dir, "hbond.xvg")
    sasa_rec_mean, sasa_rec_sd = summarize(analise_dir, "sasa_receptor.xvg")
    sasa_lig_mean, sasa_lig_sd = summarize(analise_dir, "sasa_ligante.xvg")
    mmgbsa = read_mmgbsa_total(mmgbsa_dir)

    has_results = rmsd_rec_mean is not None

    residue_entries = parse_key_residues(key_residues)
    residue_rows = []
    for r in residue_entries:
        dmean, dsd = summarize(analise_dir, f"dist_{r['label']}.xvg")
        dist_nm = f"{dmean*10:.2f}" if dmean is not None else "N/D"
        residue_rows.append(f"| {r['label']} | {r['dockdist']} Å | {r['doctype']} | {dist_nm} Å |")
    residue_table = "\n".join(residue_rows) if residue_rows else "| — | — | — | — |"
    residue_names_txt = ", ".join(r["label"] for r in residue_entries) or "N/D"

    pdb_line = f" (PDB {pdb_id})" if pdb_id else ""
    organism_line = f", {organism}" if organism else ""

    mmgbsa_line = (f"{mmgbsa[0]:.2f} ± {mmgbsa[1]:.2f} kcal/mol"
                   if mmgbsa is not None else
                   "N/D (MM-GBSA nao rodou ou falhou — ver mmgbsa.log; tratar "
                   "como opcional, ja falhou de forma irreconciliavel em outro "
                   "projeto deste laboratorio, ver Milena-MD)")

    md = f"""# Dinamica Molecular — {target_name}{pdb_line} + {ligand_name} — Secoes do Artigo

*Gerado automaticamente por `bin/gerar_artigo_md.py`. Revisar antes de usar em texto final
(passar por /humanizer e pela skill auditing-academic-sources antes de qualquer submissao).*

## Resumo

Este trabalho investigou por dinamica molecular ({time_ns} ns) a estabilidade do complexo entre
{target_name}{pdb_line}{organism_line} e {ligand_name}, um candidato identificado por triagem
virtual (AutoDock Vina). O sistema foi parametrizado com o campo de forca CHARMM36m (proteina) e
CGenFF 5.0 via ParamChem (ligante), em agua TIP3P explicita e NaCl 0,15 M (condicoes
fisiologicas humanas). {"Resultados preliminares indicam RMSD do backbone de " + fmt(rmsd_rec_mean, rmsd_rec_sd, "nm") + " e " + (f"{contacts_mean:.0f}" if contacts_mean is not None else "N/D") + " contatos receptor-ligante em media." if has_results else "Simulacao ainda nao executada — secao a preencher apos o pipeline Nextflow rodar."}

## 1. Introducao

{target_name}{pdb_line}{organism_line} foi selecionado como alvo de triagem virtual para
{ligand_name}, com a pose de docking (AutoDock Vina) avaliada quanto a estabilidade temporal
por dinamica molecular classica. *(Secao a expandir com contexto biologico especifico do alvo
e revisao da literatura sobre o ligante — ver checklist na secao 3.4; nenhuma afirmacao sobre
relevancia biologica ou precedente na literatura deve ser incluida aqui sem verificacao
explicita, ver skill auditing-academic-sources.)*

## 2. Metodologia

### 2.1 Preparacao do complexo

A estrutura inicial do receptor foi obtida do PDB {pdb_id or 'N/D'} (cadeia A), com os estados
de protonacao dos residuos ionizaveis determinados para pH 7,4 (condicao fisiologica humana —
nao o pH 8,2 usado nos demais pipelines deste laboratorio, especifico para midgut alcalino de
Lepidoptera) via PROPKA, implementado por `pdb2pqr 3.7.1` com campo de forca CHARMM. A pose
inicial do ligante (resname UNL) foi obtida por docking molecular com AutoDock Vina, com
interacoes-chave identificadas por analise pos-docking em: {residue_names_txt}.

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
4. **Producao ({time_ns} ns)** — sem restricoes, barostato de Parrinello-Rahman (Parrinello &
   Rahman, 1981; τ = 2,0 ps), integrador *leap-frog* (dt = 2 fs), ligacoes com hidrogenio
   restringidas por LINCS (Hess *et al.*, 1997), eletrostatica de longo alcance por PME
   (Darden *et al.*, 1993, `rcoulomb = 1,2 nm`).

### 2.4 Analises

RMSD do backbone do receptor e do ligante, RMSF por residuo, raio de giro, contatos
receptor-ligante (< 0,4 nm), pontes de hidrogenio, SASA do receptor e do ligante, e
distancia minima entre o ligante e os residuos de interesse identificados no docking
({residue_names_txt}), todas calculadas com ferramentas nativas do GROMACS sobre a
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
{residue_table}

### 3.3 Energia livre de ligacao (MM-GBSA)

ΔG total: **{mmgbsa_line}**

{"*(resultados ainda nao gerados — rodar o pipeline Nextflow completo)*" if not has_results else ""}

### 3.4 Convergencia com a literatura e outros projetos do laboratorio — TODO

Pendente, a preencher **apos** a producao terminar e as analises rodarem (nao
fabricar numeros de terceiros aqui — buscar e citar explicitamente):

- [ ] Comparar RMSD/RMSF obtidos com faixas tipicas reportadas para {target_name} em MD
      (buscar literatura especifica antes de citar valores).
- [ ] Buscar na literatura estudos computacionais ou experimentais de {ligand_name}
      ligando {target_name} (ou alvos homologos) e comparar modo de ligacao /
      residuos-chave / valores de ΔG de ligacao.
- [ ] Comparar robustez metodologica (protocolo de equilibracao, cutoffs, forca de
      POSRES, tempo de producao) com os pipelines ja validados deste laboratorio
      (MD-gromacs serie GORE4/SKTI/BEN, Milena-MD serie trypsin×GORE12T) —
      ver `~/.claude/.claude/agents/bioinformatics.md`.
- [ ] Avaliar se a persistencia de {residue_names_txt} ao longo da producao confirma ou
      refuta a pose de docking original (criterio sugerido: manter contato em
      >50% dos frames pos-equilibracao).
- [ ] Conferir penalidade CGenFF do ligante (`inputs/ligand-UNL*.str`) — se acima de 50,
      considerar validacao adicional dos dihedros de maior penalidade antes de
      conclusoes quantitativas fortes sobre energia de ligacao.

---
*Nao passou por /humanizer. Revisar citacoes com a skill auditing-academic-sources
antes de qualquer uso em documento final.*
"""
    return md


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analise-dir", required=True)
    ap.add_argument("--mmgbsa-dir", default=None)
    ap.add_argument("--time-ns", type=int, default=100)
    ap.add_argument("--target-name", required=True)
    ap.add_argument("--pdb-id", default="")
    ap.add_argument("--organism", default="")
    ap.add_argument("--ligand-name", required=True)
    ap.add_argument("--key-residues", default="",
                     help="'resid:resname:dockdist:doctype;...' — mesmo formato do samplesheet")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    analise_dir = Path(args.analise_dir)
    mmgbsa_dir = Path(args.mmgbsa_dir) if args.mmgbsa_dir else None
    out_path = Path(args.out)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        build_report(analise_dir, mmgbsa_dir, args.time_ns, args.target_name, args.pdb_id,
                     args.organism, args.ligand_name, args.key_residues),
        encoding="utf-8")
    print(f"[OK] Relatorio gerado em {out_path}")
    if not (analise_dir / "rmsd_backbone.xvg").exists():
        print(f"[AVISO] Nenhum resultado de analise encontrado em {analise_dir} ainda.")


if __name__ == "__main__":
    main()
