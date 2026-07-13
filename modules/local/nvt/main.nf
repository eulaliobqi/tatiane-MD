process NVT {
    tag "${meta.id}"
    label 'process_gpu'

    publishDir { "${params.outdir}/${meta.id}/nvt" }, mode: 'copy'

    input:
    tuple val(meta), path(em_gro), path(top, stageAs: 'input.top'), path(itps, stageAs: 'itp_in/*'),
                     path(prms, stageAs: 'prm_in/*')

    output:
    tuple val(meta), path("nvt.gro"), path("nvt.cpt"), path("topol.top"), path("*.itp"), path("*.prm"),
                     emit: system

    script:
    def gpu_flags = params.use_gpu ? "-nb gpu -pme gpu -bonded gpu -gpu_id ${params.gpu_id}" : ""
    def mpi       = params.mpi_cmd ?: ""
    def temp      = params.temperature
    """
    cp ${top} topol.top
    cp itp_in/*.itp .
    cp prm_in/*.prm . 2>/dev/null || true
    cp -r ${projectDir}/ff/charmm36-feb2026_cgenff-5.0.ff .

    cat > nvt.mdp << MDP_EOF
; POSRES restringe o receptor (posre.itp do pdb2gmx), POSRES_UNL restringe
; o ligante (posre_UNL.itp de gmx genrestr)
define          = -DPOSRES -DPOSRES_UNL
integrator      = md
dt              = 0.002
nsteps          = 100000
nstxout         = 1000
nstvout         = 1000
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
continuation    = no
gen-vel         = yes
gen-temp        = ${temp}
gen-seed        = -1
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau-t           = 0.1 0.1
ref-t           = ${temp} ${temp}
pcoupl          = no
pbc             = xyz
MDP_EOF

    ${params.gmx_cmd} grompp \\
        -f nvt.mdp \\
        -c ${em_gro} -r ${em_gro} \\
        -p topol.top -o nvt.tpr \\
        -maxwarn ${params.maxwarn}

    ${mpi} ${params.gmx_cmd} mdrun \\
        -v -deffnm nvt \\
        -ntomp ${params.ntomp} \\
        -pin on ${gpu_flags}
    """
}
