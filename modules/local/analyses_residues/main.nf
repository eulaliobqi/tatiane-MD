// Mirror de MD-gromacs/modules/local/analyses_triad, mas para os 2 residuos
// de interesse identificados no docking (nao a triade catalitica de uma
// protease): Arg30 (contato hidrofobico, Vina/PLIP ~4,7-4,8 A) e Glu279
// (ligacao de hidrogenio, Vina/PLIP ~1,9 A). Numeracao original do receptor
// (cadeia A), preservada pelo pipeline — ver bin/prepare_complex.py.
process ANALYSES_RESIDUES {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/analise" }, mode: 'copy'

    input:
    tuple val(meta), path(md_tpr), path(md_fit_xtc), path(lig_ndx)

    output:
    tuple val(meta), path("dist_arg30.xvg"), path("dist_glu279.xvg"),
                     path("sasa_arg30.xvg"), path("sasa_glu279.xvg"), emit: residues
    tuple val(meta), path("residues_info.txt"), emit: info

    script:
    def res_arg30  = params.res_arg30  ?: 30
    def res_glu279 = params.res_glu279 ?: 279
    """
    echo "=== ANALYSES_RESIDUES: ${meta.id} ===" >&2
    echo "Residuos de interesse: Arg${res_arg30} (hidrofobica, docking 4.7-4.8A) / Glu${res_glu279} (H-bond, docking 1.9A)" >&2

    cat > residues_info.txt << EOF
Arg30	4.7-4.8	Hidrofobica	docking_AutoDock_Vina
Glu279	1.9	Hidrogenio	docking_AutoDock_Vina
EOF

    N_CURR=\$(echo q | ${params.gmx_cmd} make_ndx \\
        -f ${md_tpr} -n ${lig_ndx} -o _tmp_count.ndx 2>&1 \\
        | grep -cE "^ *[0-9]+ +[A-Za-z]")
    rm -f _tmp_count.ndx
    ARG30_IDX=\${N_CURR}
    GLU279_IDX=\$((N_CURR + 1))

    ${params.gmx_cmd} make_ndx -f ${md_tpr} -n ${lig_ndx} -o residues.ndx << MNDX
r ${res_arg30}
name \${ARG30_IDX} Arg30
r ${res_glu279}
name \${GLU279_IDX} Glu279
q
MNDX

    printf 'Ligante\\nArg30\\n' | ${params.gmx_cmd} mindist \\
        -s ${md_tpr} -f ${md_fit_xtc} -n residues.ndx -od dist_arg30.xvg -tu ns

    printf 'Ligante\\nGlu279\\n' | ${params.gmx_cmd} mindist \\
        -s ${md_tpr} -f ${md_fit_xtc} -n residues.ndx -od dist_glu279.xvg -tu ns

    # SASA por residuo (surface = proteina completa; output = residuo
    # individual) — valores baixos indicam residuo enterrado/em contato
    printf 'Protein\\nArg30\\n' | ${params.gmx_cmd} sasa \\
        -s ${md_tpr} -f ${md_fit_xtc} -n residues.ndx -o sasa_arg30.xvg -tu ns

    printf 'Protein\\nGlu279\\n' | ${params.gmx_cmd} sasa \\
        -s ${md_tpr} -f ${md_fit_xtc} -n residues.ndx -o sasa_glu279.xvg -tu ns

    echo "[OK] Distancias e SASA de Arg30/Glu279 concluidos para ${meta.id}" >&2
    """
}
