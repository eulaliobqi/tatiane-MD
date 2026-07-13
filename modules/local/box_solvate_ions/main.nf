process BOX_SOLVATE_IONS {
    tag "${meta.id}"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/box" }, mode: 'copy'

    input:
    tuple val(meta), path(gro), path(top, stageAs: 'input.top'), path(itps, stageAs: 'itp_in/*'),
                     path(prms, stageAs: 'prm_in/*')

    output:
    tuple val(meta), path("ions.gro"), path("topol.top"), path("*.itp"), path("*.prm"), emit: system

    script:
    """
    cp ${top} topol.top
    cp itp_in/*.itp .
    cp prm_in/*.prm . 2>/dev/null || true
    cp -r ${projectDir}/ff/charmm36-feb2026_cgenff-5.0.ff .

    # Caixa cubica com margem de ${params.box_dist} nm — 1,2nm (nao os 2,0nm
    # usados nos pipelines proteina-peptideo deste laboratorio): complexo
    # globular receptor + molecula pequena, sem a elongacao de um peptideo
    ${params.gmx_cmd} editconf \\
        -f ${gro} -o box.gro \\
        -c -d ${params.box_dist} -bt ${params.box_type}

    # spc216.gro (bundled com GROMACS) e so o template de coordenadas da caixa
    # de agua equilibrada — os parametros reais vem do #include
    # charmm36-feb2026_cgenff-5.0.ff/tip3p.itp ja presente no topol.top
    ${params.gmx_cmd} solvate \\
        -cp box.gro -cs spc216.gro \\
        -p topol.top -o solv.gro

    cat > ions.mdp << 'MDP_EOF'
integrator    = steep
nsteps        = 0
cutoff-scheme = Verlet
MDP_EOF

    ${params.gmx_cmd} grompp \\
        -f ions.mdp -c solv.gro \\
        -p topol.top -o ions.tpr \\
        -maxwarn ${params.maxwarn}

    # Na+/Cl- fisiologico (alvo humano) — NAO o K+ usado nos sistemas de
    # Lepidoptera deste laboratorio ("usar NA para mamiferos")
    echo "SOL" | ${params.gmx_cmd} genion \\
        -s ions.tpr -o ions.gro \\
        -p topol.top -pname ${params.cation} -nname CL \\
        -neutral -conc ${params.nacl_conc}
    """
}
