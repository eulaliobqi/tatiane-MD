process PRODUCTION {
    tag "${meta.id}"
    label 'process_gpu_long'

    publishDir { "${params.outdir}/${meta.id}/prod" }, mode: 'copy'

    input:
    tuple val(meta), path(npt_gro), path(npt_cpt), path(top, stageAs: 'input.top'),
                     path(itps, stageAs: 'itp_in/*'), path(prms, stageAs: 'prm_in/*')

    output:
    tuple val(meta), path("md.tpr"), path("md.xtc"), emit: traj
    tuple val(meta), path("md.gro"), path("md.cpt"), path("md.edr"), emit: checkpoint

    script:
    def gpu_flags  = params.use_gpu ? "-nb gpu -pme gpu -bonded gpu -gpu_id ${params.gpu_id}" : ""
    def mpi        = params.mpi_cmd ?: ""
    def prod_steps = (params.time_ns as long) * 500000L
    def temp       = params.temperature
    """
    cp ${top} topol.top
    cp itp_in/*.itp .
    cp prm_in/*.prm . 2>/dev/null || true
    cp -r ${projectDir}/ff/charmm36-feb2026_cgenff-5.0.ff .

    cat > md.mdp << MDP_EOF
integrator           = md
dt                   = 0.002
nsteps               = ${prod_steps}
nstxout              = 0
nstvout              = 0
nstenergy            = 5000
nstlog               = 5000
nstxout-compressed   = 5000
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
continuation         = yes
tcoupl               = V-rescale
tc-grps              = Protein Non-Protein
tau-t                = 0.1 0.1
ref-t                = ${temp} ${temp}
pcoupl               = Parrinello-Rahman
pcoupltype           = isotropic
tau-p                = 2.0
ref-p                = 1.0
compressibility      = 4.5e-5
pbc                  = xyz
comm-mode            = Linear
nstcomm              = 100
MDP_EOF

    ${params.gmx_cmd} grompp \\
        -f md.mdp \\
        -c ${npt_gro} -t ${npt_cpt} \\
        -p topol.top -o md.tpr \\
        -maxwarn ${params.maxwarn}

    # Retoma checkpoint se existir (recover de execucao interrompida)
    CPI=""
    [ -s md.cpt ] && CPI="-cpi md.cpt"

    ${mpi} ${params.gmx_cmd} mdrun \\
        -v -deffnm md \${CPI} \\
        -ntomp ${params.ntomp} \\
        -pin on ${gpu_flags}
    """
}
