// Portado de Milena-MD/modules/local/fe_rerun. Roda uma vez POR FASE
// (bound/relocated) sobre a subtrajetoria reduzida de CLUSTERING —
// comparacao direta de energia de interacao entre a pose de docking
// original e o sitio novo pos-transicao (~65-78ns). NAO substitui um
// MM-GBSA rigoroso (sem termo de solvatacao implicita) — ver ressalvas em
// bin/interaction_entropy.py. Usa BOX_SOLVATE_IONS.out.system (topologia
// JA solvatada/com ions), NAO TOPOLOGY.out.topology — mesma regra
// nao-negociavel do MMGBSA deste projeto (numero de coordenadas precisa
// bater com a topologia).
process FE_RERUN {
    tag "${meta.id}:${meta.phase}"
    label 'process_medium'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/analise_extra/${meta.phase}/fe_estimate" }, mode: 'copy'

    input:
    tuple val(meta), path(complexo_gro), path(top, stageAs: 'input.top'),
                     path(itps, stageAs: 'itp_in/*'), path(fe_xtc), path(lig_ndx)

    output:
    tuple val(meta), path("interaction_energy.xvg"), emit: energy
    tuple val(meta), path("fe_rerun.log"),            emit: log

    script:
    """
    echo "=== FE_RERUN: ${meta.id} (${meta.phase}) ===" >&2
    cp ${top} topol.top
    cp itp_in/*.itp .
    cp -r ${projectDir}/ff/charmm36-feb2026_cgenff-5.0.ff .

    {
        echo "=== FE_RERUN: rerun com energygrps Receptor/Ligante ==="
        echo "Sistema: ${meta.id}  Fase: ${meta.phase}"
        echo ""
    } > fe_rerun.log

    cat > rerun.mdp << MDP_EOF
integrator           = md
dt                   = 0.002
nsteps               = 0
nstenergy            = 1
cutoff-scheme        = Verlet
nstlist              = 20
coulombtype          = PME
rcoulomb             = 1.2
vdwtype              = Cut-off
vdw-modifier         = Force-switch
rvdw-switch          = 1.0
rvdw                 = 1.2
DispCorr             = no
constraints          = h-bonds
constraint-algorithm = LINCS
pbc                  = xyz
energygrps           = Receptor Ligante
MDP_EOF

    # NOTA (mesma licao de MMGBSA/CLUSTERING neste projeto): Nextflow roda com
    # "set -e" implicito -- cada comando arriscado precisa de "|| true" para
    # nao abortar o script inteiro antes do fallback abaixo rodar.
    ${params.gmx_cmd} grompp \\
        -f rerun.mdp \\
        -c ${complexo_gro} \\
        -p topol.top \\
        -n ${lig_ndx} \\
        -o rerun.tpr \\
        -maxwarn ${params.maxwarn} \\
        2>&1 | tee -a fe_rerun.log || true

    if [ -f rerun.tpr ]; then
        ${params.mpi_cmd} ${params.gmx_cmd} mdrun \\
            -s rerun.tpr -rerun ${fe_xtc} \\
            -deffnm rerun -ntomp ${params.ntomp} \\
            2>&1 | tee -a fe_rerun.log || true
    else
        echo "ERRO: rerun.tpr nao gerado (grompp falhou)" | tee -a fe_rerun.log
    fi

    if [ -f rerun.edr ]; then
        printf 'Coul-SR:Receptor-Ligante\\nLJ-SR:Receptor-Ligante\\n0\\n' | \\
            ${params.gmx_cmd} energy -f rerun.edr -o interaction_energy.xvg \\
            2>&1 | tee -a fe_rerun.log || true
    fi

    if [ ! -f interaction_energy.xvg ]; then
        echo "ERRO: interaction_energy.xvg nao gerado" | tee -a fe_rerun.log
        echo "# rerun falhou -- ver fe_rerun.log" > interaction_energy.xvg
    else
        echo "[OK] interaction_energy.xvg gerado" | tee -a fe_rerun.log
    fi
    """
}
