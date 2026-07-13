process NPT {
    tag "${meta.id}"
    label 'process_gpu'

    publishDir { "${params.outdir}/${meta.id}/npt" }, mode: 'copy'

    input:
    tuple val(meta), path(nvt_gro), path(nvt_cpt), path(top, stageAs: 'input.top'),
                     path(itps, stageAs: 'itp_in/*'), path(prms, stageAs: 'prm_in/*')

    output:
    tuple val(meta), path("npt.gro"), path("npt.cpt"), path("topol.top"), path("*.itp"), path("*.prm"),
                     emit: system

    script:
    def gpu_flags = params.use_gpu ? "-nb gpu -pme gpu -bonded gpu -gpu_id ${params.gpu_id}" : ""
    def mpi       = params.mpi_cmd ?: ""
    def temp      = params.temperature
    """
    cp ${top} topol.top
    cp itp_in/*.itp .
    cp prm_in/*.prm . 2>/dev/null || true
    cp -r ${projectDir}/ff/charmm36-mar2019.ff .

    cat > npt.mdp << MDP_EOF
define          = -DPOSRES -DPOSRES_UNL
integrator      = md
dt              = 0.002
nsteps          = 250000
nstxout         = 2500
nstvout         = 2500
nstenergy       = 500
cutoff-scheme   = Verlet
nstlist         = 20
coulombtype     = PME
rcoulomb        = 1.2
vdwtype         = Cut-off
vdw-modifier    = Force-switch
rvdw-switch     = 1.0
rvdw            = 1.2
DispCorr        = no
constraints     = h-bonds
constraint-algorithm = LINCS
continuation    = yes
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau-t           = 0.1 0.1
ref-t           = ${temp} ${temp}
pcoupl          = Berendsen
pcoupltype      = isotropic
tau-p           = 2.0
ref-p           = 1.0
compressibility = 4.5e-5
refcoord_scaling = com
pbc             = xyz
MDP_EOF

    ${params.gmx_cmd} grompp \\
        -f npt.mdp \\
        -c ${nvt_gro} -r ${nvt_gro} -t ${nvt_cpt} \\
        -p topol.top -o npt.tpr \\
        -maxwarn ${params.maxwarn}

    ${mpi} ${params.gmx_cmd} mdrun \\
        -v -deffnm npt \\
        -ntomp ${params.ntomp} \\
        -pin on ${gpu_flags}
    """
}
