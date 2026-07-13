#!/usr/bin/env bash
# ============================================================================
# Dinamica Molecular: receptor 2I9T (dominio de ligacao a DNA de NF-kB,
# cadeia A, res 17-291) + ligante daidzeina (isoflavona, resname UNL,
# pose de docking AutoDock Vina) — CHARMM36m + CGenFF 5.0 (ParamChem)
#
# Executar no servidor (eulalio@200.235.143.10), dentro de `screen`, com o
# ambiente conda `md-gromacs` ativo (GROMACS 2026.0 CUDA, pdb2pqr 3.7.1):
#
#   screen -S tatiana-2i9t-daidzeina
#   mamba activate md-gromacs
#   cd ~/gromacs/Tatiana-MD   # apos git pull
#   bash bin/run_md.sh 2>&1 | tee run_md.log
#
# Retomavel: cada etapa so roda se o output final dela ainda nao existir.
# Para forcar uma etapa a rodar de novo, apague o arquivo de output dela
# (ou o diretorio da etapa) e rode o script de novo.
# ============================================================================
set -euo pipefail

# --- Config -----------------------------------------------------------------
SAMPLE_ID="2I9T-daidzeina"
PH=7.4                       # fisiologico, alvo humano (NF-kB) — nao usar 8.2
                              # (esse era especifico p/ midgut de inseto, ver bioinformatics.md)
TIME_NS="${TIME_NS:-100}"    # producao, sobrescrevivel: TIME_NS=200 bash bin/run_md.sh
TEMPERATURE=300
NACL_CONC=0.15                # fisiologico p/ alvo humano (lab usa KCl 0.10M so p/ inseto)
CATION="NA"
BOX_DIST=1.2                  # nm; padrao do tutorial CHARMM36 p/ complexo globular+lig pequeno
                               # (nao usar os 2.0nm do lab — aquele valor eh p/ complexos
                               # proteina-peptideo alongados, nao se aplica aqui)
BOX_TYPE="cubic"

GMX="gmx_mpi"                 # regra nao-negociavel: gmx puro nao existe no servidor
MPI="mpirun -np 1"
NTOMP=8
USE_GPU=1
GPU_ID=0
MAXWARN=2

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUTS="${ROOT_DIR}/inputs"
BIN="${ROOT_DIR}/bin"
MDP="${ROOT_DIR}/mdp"
FF_SRC="${ROOT_DIR}/ff/charmm36-mar2019.ff"
OUT="${ROOT_DIR}/results/${SAMPLE_ID}"

mkdir -p "${OUT}"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

gpu_flags_nb_only() { [ "${USE_GPU}" = "1" ] && echo "-nb gpu -gpu_id ${GPU_ID}" || echo ""; }
gpu_flags_full()    { [ "${USE_GPU}" = "1" ] && echo "-nb gpu -pme gpu -bonded gpu -gpu_id ${GPU_ID}" || echo ""; }

# ============================================================================
# STAGE 1 — PREPARE_PH: pdb2pqr + PROPKA (pH 7.4, saida CHARMM) no receptor
# ============================================================================
stage_prepare_ph() {
    local d="${OUT}/00_prep_ph"
    local final="${d}/receptor_ph.pdb"
    [ -s "${final}" ] && { log "PREPARE_PH: ja concluido, pulando"; return; }
    mkdir -p "${d}"; cd "${d}"
    log "PREPARE_PH: pdb2pqr --ff CHARMM --ffout CHARMM --with-ph ${PH} (PROPKA)"

    pdb2pqr --ff CHARMM --ffout CHARMM \
        --titration-state-method propka --with-ph "${PH}" \
        --pdb-output receptor_raw.pdb \
        --nodebump \
        "${INPUTS}/receptor-2I9T-original.pdb" receptor.pqr

    python "${BIN}/pdb2pqr_process_charmm.py" receptor_raw.pdb "${final}"
    log "PREPARE_PH: OK -> ${final}"
}

# ============================================================================
# STAGE 2 — LIGAND_TOPOLOGY: cgenff_charmm2gmx.py (str do ParamChem -> itp/prm)
# ============================================================================
stage_ligand_topology() {
    local d="${OUT}/01_ligand_topology"
    local final="${d}/unl.itp"
    [ -s "${final}" ] && { log "LIGAND_TOPOLOGY: ja concluido, pulando"; return; }
    mkdir -p "${d}"; cd "${d}"
    log "LIGAND_TOPOLOGY: cgenff_charmm2gmx.py UNL (daidzeina, CGenFF 5.0, penalty=53/23)"

    python "${BIN}/cgenff_charmm2gmx.py" UNL \
        "${INPUTS}/ligand-UNL.cgenff.mol2" \
        "${INPUTS}/ligand-UNL.str" \
        "${FF_SRC}" 2>&1 | tee cgenff_charmm2gmx.log

    [ -s unl.itp ] || { log "ERRO: unl.itp nao foi gerado — ver cgenff_charmm2gmx.log"; exit 1; }

    # Converte a pose (todos os atomos, com H) para .gro via editconf, e gera
    # a restricao posicional do ligante (usada em NVT/NPT sob -DPOSRES_UNL)
    "${GMX}" editconf -f unl_ini.pdb -o unl_ini.gro
    "${GMX}" genrestr -f unl_ini.pdb -o posre_UNL.itp -fc 1000 1000 1000 <<< "0"

    cat >> unl.itp << 'EOF'

#ifdef POSRES_UNL
#include "posre_UNL.itp"
#endif
EOF
    log "LIGAND_TOPOLOGY: OK -> ${final} (+ posre_UNL.itp, unl_ini.gro)"
}

# ============================================================================
# STAGE 3 — PREPARE_COMPLEX: mescla receptor_ph.pdb (cadeia A) + unl_ini.pdb
#           (cadeia B, renumerado para logo apos o receptor)
# ============================================================================
stage_prepare_complex() {
    local d="${OUT}/02_complex"
    local final="${d}/complex.pdb"
    [ -s "${final}" ] && { log "PREPARE_COMPLEX: ja concluido, pulando"; return; }
    mkdir -p "${d}"; cd "${d}"
    log "PREPARE_COMPLEX: mesclando receptor (cadeia A) + UNL (cadeia B)"

    python "${BIN}/prepare_complex.py" \
        --receptor "${OUT}/00_prep_ph/receptor_ph.pdb" \
        --ligand   "${OUT}/01_ligand_topology/unl_ini.pdb" \
        --out      "${final}"
    log "PREPARE_COMPLEX: OK -> ${final}"
}

# ============================================================================
# STAGE 4 — TOPOLOGY: pdb2gmx (receptor, CHARMM36) + merge com topologia UNL
# ============================================================================
stage_topology() {
    local d="${OUT}/03_topology"
    local final="${d}/topol.top"
    [ -s "${final}" ] && { log "TOPOLOGY: ja concluido, pulando"; return; }
    mkdir -p "${d}"; cd "${d}"
    log "TOPOLOGY: pdb2gmx (charmm36-mar2019, tip3p) + merge com ligante CGenFF"

    # gmx pdb2gmx localiza NAME.ff por nome apenas se a pasta existir no cwd
    # (ou via GMXLIB) — copia local do FF vendorizado
    [ -d charmm36-mar2019.ff ] || cp -r "${FF_SRC}" charmm36-mar2019.ff

    awk '/^ATOM/ && substr($0,22,1)=="A" {print}' "${OUT}/02_complex/complex.pdb" > receptor.pdb
    echo "TER" >> receptor.pdb
    echo "END" >> receptor.pdb

    NATOM_REC=$(grep -c "^ATOM" receptor.pdb || echo 0)
    log "  Receptor: ${NATOM_REC} atomos ATOM (cadeia A)"

    # Sem -ter: usa terminos padrao (NH3+/COO-) sem prompt interativo —
    # os indices numericos do prompt -ter DIFEREM entre AMBER e CHARMM
    # (nao reusar o "printf '0\n0\n'" do pipeline BEN/AMBER aqui).
    # -ignh e redundante (receptor ja sem H apos pdb2pqr) mas mantido por seguranca.
    "${GMX}" pdb2gmx \
        -f receptor.pdb \
        -o receptor.gro \
        -p receptor.top \
        -i posre.itp \
        -ff charmm36-mar2019 \
        -water tip3p \
        -ignh \
        2>&1 | tee pdb2gmx.log

    [ -s receptor.gro ] || { log "ERRO: pdb2gmx falhou — ver pdb2gmx.log"; exit 1; }

    cp "${OUT}/01_ligand_topology/unl.itp" .
    cp "${OUT}/01_ligand_topology/unl.prm" .
    cp "${OUT}/01_ligand_topology/posre_UNL.itp" .
    cp "${OUT}/01_ligand_topology/unl_ini.gro" .

    python "${BIN}/merge_small_molecule_topology.py" \
        --protein-gro receptor.gro \
        --ligand-gro  unl_ini.gro \
        --protein-top receptor.top \
        --ligand-itp  unl.itp \
        --ligand-prm  unl.prm \
        --ligand-mol  UNL \
        --out-gro     complexo.gro \
        --out-top     "${final}"

    [ -s complexo.gro ] || { log "ERRO: merge_small_molecule_topology.py falhou"; exit 1; }
    NTOTAL=$(awk 'NR==2{print $1}' complexo.gro)
    log "TOPOLOGY: OK -> ${final} (complexo.gro: ${NTOTAL} atomos)"
}

# ============================================================================
# STAGE 5 — BOX_SOLVATE_IONS
# ============================================================================
stage_box_solvate_ions() {
    local d="${OUT}/04_box"
    local final="${d}/ions.gro"
    [ -s "${final}" ] && { log "BOX_SOLVATE_IONS: ja concluido, pulando"; return; }
    mkdir -p "${d}"; cd "${d}"
    log "BOX_SOLVATE_IONS: caixa ${BOX_TYPE} ${BOX_DIST}nm + TIP3P + ${CATION}Cl ${NACL_CONC}M"

    local topo_dir="${OUT}/03_topology"
    cp "${topo_dir}/topol.top" topol.top
    cp "${topo_dir}"/*.itp . 2>/dev/null || true
    cp "${topo_dir}"/*.prm . 2>/dev/null || true
    [ -d charmm36-mar2019.ff ] || cp -r "${FF_SRC}" charmm36-mar2019.ff

    "${GMX}" editconf -f "${topo_dir}/complexo.gro" -o box.gro \
        -c -d "${BOX_DIST}" -bt "${BOX_TYPE}"

    # spc216.gro (bundled com o GROMACS) e so o template de coordenadas da
    # caixa de agua equilibrada — os parametros reais vem do #include
    # "charmm36-mar2019.ff/tip3p.itp" ja presente no topol.top (via pdb2gmx
    # -water tip3p). Mesma pratica usada no pipeline AMBER do laboratorio.
    "${GMX}" solvate -cp box.gro -cs spc216.gro \
        -p topol.top -o solv.gro

    "${GMX}" grompp -f "${MDP}/ions.mdp" -c solv.gro \
        -p topol.top -o ions.tpr -maxwarn "${MAXWARN}"

    echo "SOL" | "${GMX}" genion \
        -s ions.tpr -o "${final}" \
        -p topol.top -pname "${CATION}" -nname CL \
        -neutral -conc "${NACL_CONC}"

    log "BOX_SOLVATE_IONS: OK -> ${final}"
}

# ============================================================================
# STAGE 6 — MINIMIZATION
# ============================================================================
stage_minimization() {
    local d="${OUT}/05_em"
    local final="${d}/em.gro"
    [ -s "${final}" ] && { log "MINIMIZATION: ja concluido, pulando"; return; }
    mkdir -p "${d}"; cd "${d}"
    log "MINIMIZATION: steepest descent, emtol=1000"

    # topol.top referencia "charmm36-mar2019.ff/forcefield.itp" como caminho
    # relativo ao cwd do grompp — precisa da pasta do FF presente em CADA
    # etapa, nao so na que rodou o pdb2gmx
    [ -d charmm36-mar2019.ff ] || cp -r "${FF_SRC}" charmm36-mar2019.ff

    cp "${OUT}/04_box/topol.top" topol.top
    cp "${OUT}/04_box"/*.itp . 2>/dev/null || true
    cp "${OUT}/04_box"/*.prm . 2>/dev/null || true

    "${GMX}" grompp -f "${MDP}/em.mdp" -c "${OUT}/04_box/ions.gro" \
        -p topol.top -o em.tpr -maxwarn "${MAXWARN}"

    ${MPI} "${GMX}" mdrun -v -deffnm em -ntomp "${NTOMP}" -pin on $(gpu_flags_nb_only)

    log "MINIMIZATION: OK -> ${final}"
}

# ============================================================================
# STAGE 7 — NVT (200 ps, V-rescale, POSRES+POSRES_UNL)
# ============================================================================
stage_nvt() {
    local d="${OUT}/06_nvt"
    local final="${d}/nvt.gro"
    [ -s "${final}" ] && { log "NVT: ja concluido, pulando"; return; }
    mkdir -p "${d}"; cd "${d}"
    log "NVT: 200 ps, 300 K, V-rescale, posicoes restritas"

    [ -d charmm36-mar2019.ff ] || cp -r "${FF_SRC}" charmm36-mar2019.ff

    cp "${OUT}/05_em/topol.top" topol.top
    cp "${OUT}/05_em"/*.itp . 2>/dev/null || true
    cp "${OUT}/05_em"/*.prm . 2>/dev/null || true

    "${GMX}" grompp -f "${MDP}/nvt.mdp" \
        -c "${OUT}/05_em/em.gro" -r "${OUT}/05_em/em.gro" \
        -p topol.top -o nvt.tpr -maxwarn "${MAXWARN}"

    ${MPI} "${GMX}" mdrun -v -deffnm nvt -ntomp "${NTOMP}" -pin on $(gpu_flags_full)

    log "NVT: OK -> ${final}"
}

# ============================================================================
# STAGE 8 — NPT (500 ps, Berendsen, POSRES+POSRES_UNL)
# ============================================================================
stage_npt() {
    local d="${OUT}/07_npt"
    local final="${d}/npt.gro"
    [ -s "${final}" ] && { log "NPT: ja concluido, pulando"; return; }
    mkdir -p "${d}"; cd "${d}"
    log "NPT: 500 ps, 1 bar, Berendsen, posicoes restritas"

    [ -d charmm36-mar2019.ff ] || cp -r "${FF_SRC}" charmm36-mar2019.ff

    cp "${OUT}/06_nvt/topol.top" topol.top
    cp "${OUT}/06_nvt"/*.itp . 2>/dev/null || true
    cp "${OUT}/06_nvt"/*.prm . 2>/dev/null || true

    "${GMX}" grompp -f "${MDP}/npt.mdp" \
        -c "${OUT}/06_nvt/nvt.gro" -r "${OUT}/06_nvt/nvt.gro" -t "${OUT}/06_nvt/nvt.cpt" \
        -p topol.top -o npt.tpr -maxwarn "${MAXWARN}"

    ${MPI} "${GMX}" mdrun -v -deffnm npt -ntomp "${NTOMP}" -pin on $(gpu_flags_full)

    log "NPT: OK -> ${final}"
}

# ============================================================================
# STAGE 9 — PRODUCTION (padrao 100 ns, Parrinello-Rahman, sem restricoes)
# ============================================================================
stage_production() {
    local d="${OUT}/08_prod"
    local final="${d}/md.gro"
    [ -s "${final}" ] && { log "PRODUCTION: ja concluido, pulando"; return; }
    mkdir -p "${d}"; cd "${d}"
    log "PRODUCTION: ${TIME_NS} ns, 300 K, 1 bar, Parrinello-Rahman"

    [ -d charmm36-mar2019.ff ] || cp -r "${FF_SRC}" charmm36-mar2019.ff

    cp "${OUT}/07_npt/topol.top" topol.top
    cp "${OUT}/07_npt"/*.itp . 2>/dev/null || true
    cp "${OUT}/07_npt"/*.prm . 2>/dev/null || true

    local prod_steps=$(( TIME_NS * 500000 ))
    sed "s/^nsteps.*=.*/nsteps               = ${prod_steps}/" "${MDP}/production.mdp" > md.mdp

    if [ ! -s md.tpr ]; then
        "${GMX}" grompp -f md.mdp \
            -c "${OUT}/07_npt/npt.gro" -t "${OUT}/07_npt/npt.cpt" \
            -p topol.top -o md.tpr -maxwarn "${MAXWARN}"
    fi

    local cpi=""
    [ -s md.cpt ] && cpi="-cpi md.cpt"

    ${MPI} "${GMX}" mdrun -v -deffnm md ${cpi} -ntomp "${NTOMP}" -pin on $(gpu_flags_full)

    log "PRODUCTION: OK -> ${final}"
}

# ============================================================================
# STAGE 10 — POSTPROCESS: remove PBC, centraliza, ajusta (fit) no receptor
# ============================================================================
stage_postprocess() {
    local d="${OUT}/09_postprocess"
    local final="${d}/md_fit.xtc"
    [ -s "${final}" ] && { log "POSTPROCESS: ja concluido, pulando"; return; }
    mkdir -p "${d}"; cd "${d}"
    log "POSTPROCESS: trjconv -pbc mol -center + fit rot+trans no backbone do receptor"

    local prod="${OUT}/08_prod"
    echo "Protein System" | "${GMX}" trjconv -s "${prod}/md.tpr" -f "${prod}/md.xtc" \
        -pbc mol -center -o md_center.xtc

    echo "Backbone System" | "${GMX}" trjconv -s "${prod}/md.tpr" -f md_center.xtc \
        -fit rot+trans -o "${final}"

    cp "${prod}/md.tpr" .
    log "POSTPROCESS: OK -> ${final}"
}

main() {
    log "=== run_md.sh: ${SAMPLE_ID} | pH=${PH} | ${TIME_NS} ns | ${CATION}Cl ${NACL_CONC}M ==="
    stage_prepare_ph
    stage_ligand_topology
    stage_prepare_complex
    stage_topology
    stage_box_solvate_ions
    stage_minimization
    stage_nvt
    stage_npt
    stage_production
    stage_postprocess
    log "=== MD concluida. Proximo passo: bash bin/analyze.sh ==="
}

main "$@"
