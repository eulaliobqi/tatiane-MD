#!/usr/bin/env bash
# ============================================================================
# Analises pos-MD: receptor 2I9T + daidzeina (UNL) — mirror do padrao
# ANALYSES_BEN + ANALYSES_TRIAD do MD-gromacs, adaptado para 2 residuos de
# interesse especificos (nao a triade catalitica de uma protease):
#   Arg30  — contato hidrofobico previsto no docking (Vina/PLIP: 4.7-4.8 A)
#   Glu279 — ligacao de hidrogenio prevista no docking (Vina/PLIP: 1.9 A)
# Roda depois de bin/run_md.sh terminar (usa md.tpr + md_fit.xtc + complex.pdb)
#
#   mamba activate md-gromacs
#   cd ~/gromacs/Tatiana-MD
#   bash bin/analyze.sh 2>&1 | tee analyze.log
# ============================================================================
set -euo pipefail

SAMPLE_ID="2I9T-daidzeina"
GMX="gmx_mpi"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${ROOT_DIR}/results/${SAMPLE_ID}"
COMPLEX_PDB="${OUT}/02_complex/complex.pdb"
MD_TPR="${OUT}/09_postprocess/md.tpr"
MD_FIT_XTC="${OUT}/09_postprocess/md_fit.xtc"
ANALYSIS_DIR="${OUT}/10_analysis"

# Residuos de interesse (numeracao original do receptor, cadeia A — preservada
# pelo pipeline: ver bin/prepare_complex.py)
RES_ARG30=30
RES_GLU279=279

mkdir -p "${ANALYSIS_DIR}"; cd "${ANALYSIS_DIR}"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

[ -s "${MD_TPR}" ] || { log "ERRO: ${MD_TPR} nao existe — rode bin/run_md.sh primeiro"; exit 1; }
[ -s "${MD_FIT_XTC}" ] || { log "ERRO: ${MD_FIT_XTC} nao existe — rode bin/run_md.sh primeiro"; exit 1; }

# --- Detecta o ligante (UNL, cadeia B, HETATM) no complexo ------------------
LIG_FIRST=$(awk '(/^ATOM/ || /^HETATM/) && substr($0,22,1)=="B" {print substr($0,23,4)+0; exit}' "${COMPLEX_PDB}")
LIG_LAST=$(awk '(/^ATOM/ || /^HETATM/) && substr($0,22,1)=="B" {r=substr($0,23,4)+0} END{print r}' "${COMPLEX_PDB}")
[ -n "${LIG_FIRST}" ] && [ -n "${LIG_LAST}" ] || { log "ERRO: cadeia B (UNL) nao detectada em ${COMPLEX_PDB}"; exit 1; }
log "Ligante UNL: residuo ${LIG_FIRST}-${LIG_LAST} (cadeia B)"

# --- Indice: Ligante, Receptor, Arg30, Glu279 --------------------------------
N_DEFAULT=$(echo q | "${GMX}" make_ndx -f "${MD_TPR}" -o _default.ndx 2>&1 \
    | grep -cE "^ *[0-9]+ +[A-Za-z]")
rm -f _default.ndx
LIG_IDX=${N_DEFAULT}
REC_IDX=$((N_DEFAULT + 1))
ARG30_IDX=$((N_DEFAULT + 2))
GLU279_IDX=$((N_DEFAULT + 3))

"${GMX}" make_ndx -f "${MD_TPR}" -o analysis.ndx << EOF
r ${LIG_FIRST}-${LIG_LAST}
name ${LIG_IDX} Ligante
1 & ! ${LIG_IDX}
name ${REC_IDX} Receptor
r ${RES_ARG30}
name ${ARG30_IDX} Arg30
r ${RES_GLU279}
name ${GLU279_IDX} Glu279
q
EOF

# --- RMSD backbone do receptor -----------------------------------------------
log "RMSD backbone (receptor)"
printf 'Backbone\nBackbone\n' | "${GMX}" rms \
    -s "${MD_TPR}" -f "${MD_FIT_XTC}" -n analysis.ndx -o rmsd_backbone.xvg -tu ns

# --- RMSD ligante -------------------------------------------------------------
log "RMSD ligante (UNL)"
printf 'Ligante\nLigante\n' | "${GMX}" rms \
    -s "${MD_TPR}" -f "${MD_FIT_XTC}" -n analysis.ndx -o rmsd_ligante.xvg -tu ns

# --- RMSF por residuo (backbone receptor) ------------------------------------
log "RMSF por residuo"
printf 'Backbone\n' | "${GMX}" rmsf \
    -s "${MD_TPR}" -f "${MD_FIT_XTC}" -n analysis.ndx -o rmsf_residuos.xvg -res -fit

# --- Raio de giro (proteina) --------------------------------------------------
log "Raio de giro"
printf 'Protein\n' | "${GMX}" gyrate \
    -s "${MD_TPR}" -f "${MD_FIT_XTC}" -n analysis.ndx -o gyrate.xvg -tu ns

# --- Contatos receptor-ligante < 0.4 nm (mesmo cutoff do pipeline BEN) -------
log "Contatos receptor-ligante (<0.4nm)"
printf 'Receptor\nLigante\n' | "${GMX}" mindist \
    -s "${MD_TPR}" -f "${MD_FIT_XTC}" -n analysis.ndx \
    -od mindist_receptor_ligante.xvg -on numcont_receptor_ligante.xvg -d 0.4 -tu ns

# --- Pontes de hidrogenio receptor-ligante -----------------------------------
# CGenFF usa nomes de atomo CHARMM padrao (N/O + H) — gmx hbond costuma
# reconhecer normalmente, mas mantem fallback defensivo (mesmo padrao do
# pipeline BEN/GAFF2, onde isso de fato falhou em producao)
log "Pontes de hidrogenio receptor-ligante"
if ! printf 'Receptor\nLigante\n' | "${GMX}" hbond \
        -s "${MD_TPR}" -f "${MD_FIT_XTC}" -n analysis.ndx \
        -num hbond.xvg -tu ns 2>&1 | tee hbond.log; then
    log "AVISO: gmx hbond falhou — gravando placeholder vazio"
    printf '# gmx hbond falhou — ver hbond.log\n@ title "Number of Hydrogen Bonds"\n@ xaxis label "Time (ns)"\n@ yaxis label "Number"\n0.000 0\n' > hbond.xvg
fi

# --- SASA receptor e ligante ---------------------------------------------------
log "SASA receptor"
"${GMX}" sasa -s "${MD_TPR}" -f "${MD_FIT_XTC}" -n analysis.ndx \
    -surface 'Protein' -output 'Protein' -o sasa_receptor.xvg -tu ns

log "SASA ligante (valores baixos = UNL enterrado na interface)"
"${GMX}" sasa -s "${MD_TPR}" -f "${MD_FIT_XTC}" -n analysis.ndx \
    -surface 'Ligante' -output 'Ligante' -o sasa_ligante.xvg -tu ns

# --- Distancia minima Ligante-Arg30 (contato hidrofobico, docking: 4.7-4.8A) -
log "Distancia Ligante-Arg30 (docking previu contato hidrofobico ~4.7-4.8 A)"
printf 'Ligante\nArg30\n' | "${GMX}" mindist \
    -s "${MD_TPR}" -f "${MD_FIT_XTC}" -n analysis.ndx -od dist_arg30.xvg -tu ns

# --- Distancia minima Ligante-Glu279 (H-bond, docking: 1.9A) -----------------
log "Distancia Ligante-Glu279 (docking previu ligacao H ~1.9 A)"
printf 'Ligante\nGlu279\n' | "${GMX}" mindist \
    -s "${MD_TPR}" -f "${MD_FIT_XTC}" -n analysis.ndx -od dist_glu279.xvg -tu ns

cat > residues_of_interest.txt << EOF
Arg30	4.7-4.8	Hidrofobica	docking_AutoDock_Vina
Glu279	1.9	Hidrogenio	docking_AutoDock_Vina
EOF

log "=== Analises concluidas -> ${ANALYSIS_DIR} ==="
log "Proximo passo: python bin/gerar_artigo_md.py"
