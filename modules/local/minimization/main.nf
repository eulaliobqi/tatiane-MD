process MINIMIZATION {
    tag "${meta.id}"
    label 'process_gpu'

    publishDir { "${params.outdir}/${meta.id}/em" }, mode: 'copy'

    input:
    tuple val(meta), path(ions_gro), path(top, stageAs: 'input.top'), path(itps, stageAs: 'itp_in/*'),
                     path(prms, stageAs: 'prm_in/*')

    output:
    tuple val(meta), path("em.gro"), path("topol.top"), path("*.itp"), path("*.prm"), emit: system

    script:
    // steep NAO suporta PME GPU -> usa apenas -nb gpu (sem -pme gpu)
    def gpu_flags = params.use_gpu ? "-nb gpu -gpu_id ${params.gpu_id}" : ""
    def mpi       = params.mpi_cmd  ?: ""
    """
    cp ${top} topol.top
    cp itp_in/*.itp .
    cp prm_in/*.prm . 2>/dev/null || true
    # topol.top referencia charmm36-mar2019.ff/forcefield.itp por caminho
    # relativo ao cwd do grompp — precisa da pasta do FF em CADA etapa
    cp -r ${projectDir}/ff/charmm36-mar2019.ff .

    cat > em.mdp << 'MDP_EOF'
; Minimizacao — CHARMM36 + CGenFF: Force-switch + DispCorr=no (nao copiar
; mdp AMBER de outros pipelines deste laboratorio para uso aqui)
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
cutoff-scheme   = Verlet
nstlist         = 20
coulombtype     = PME
rcoulomb        = 1.2
vdwtype         = Cut-off
vdw-modifier    = Force-switch
rvdw-switch     = 1.0
rvdw            = 1.2
DispCorr        = no
pbc             = xyz
MDP_EOF

    ${params.gmx_cmd} grompp \\
        -f em.mdp -c ${ions_gro} \\
        -p topol.top -o em.tpr \\
        -maxwarn ${params.maxwarn}

    ${mpi} ${params.gmx_cmd} mdrun \\
        -v -deffnm em \\
        -ntomp ${params.ntomp} \\
        -pin on ${gpu_flags}
    """
}
