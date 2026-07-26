// Distancia minima + SASA entre o ligante e os residuos-chave identificados
// no docking (numero variavel por sistema — nao mais um par fixo). Lista
// vem de meta.key_residues (samplesheet, coluna key_residues), formato
// "resid:resname:distancia_docking:tipo" separado por ";", ex.:
// "30:ARG:4.7-4.8:Hidrofobica;279:GLU:1.9:Hidrogenio". Numeracao original
// do receptor (cadeia A), preservada pelo pipeline — ver bin/prepare_complex.py.
process ANALYSES_RESIDUES {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/analise" }, mode: 'copy'

    input:
    tuple val(meta), path(md_tpr), path(md_fit_xtc), path(lig_ndx)

    output:
    tuple val(meta), path("dist_*.xvg"), path("sasa_*.xvg"), emit: residues
    tuple val(meta), path("residues_info.txt"), emit: info

    script:
    if (!meta.key_residues) {
        error "meta.key_residues vazio para ${meta.id} — preencha a coluna key_residues no samplesheet"
    }
    def entries = meta.key_residues.split(';').collect { e ->
        def (resid, resname, dockdist, doctype) = e.split(':')
        [resid: resid, resname: resname, dockdist: dockdist, doctype: doctype,
         label: "${resname.toLowerCase().capitalize()}${resid}"]
    }
    def info_lines = entries.collect { "${it.label}\t${it.dockdist}\t${it.doctype}\tdocking_AutoDock_Vina" }.join('\n')
    def mindist_cmds = entries.collect { e ->
        "printf 'Ligante\\n${e.label}\\n' | ${params.gmx_cmd} mindist -s ${md_tpr} -f ${md_fit_xtc} -n residues.ndx -od dist_${e.label}.xvg -tu ns"
    }.join('\n\n    ')
    def sasa_cmds = entries.collect { e ->
        "printf 'Protein\\n${e.label}\\n' | ${params.gmx_cmd} sasa -s ${md_tpr} -f ${md_fit_xtc} -n residues.ndx -o sasa_${e.label}.xvg -tu ns"
    }.join('\n\n    ')
    """
    echo "=== ANALYSES_RESIDUES: ${meta.id} (${entries.size()} residuo(s)-chave) ===" >&2

    cat > residues_info.txt << EOF
${info_lines}
EOF

    N_CURR=\$(echo q | ${params.gmx_cmd} make_ndx \\
        -f ${md_tpr} -n ${lig_ndx} -o _tmp_count.ndx 2>&1 \\
        | grep -cE "^ *[0-9]+ +[A-Za-z]")
    rm -f _tmp_count.ndx

    {
    ${entries.withIndex().collect { e, i -> "echo \"r ${e.resid}\"; echo \"name \$((N_CURR + ${i})) ${e.label}\"" }.join('\n    ')}
    echo q
    } | ${params.gmx_cmd} make_ndx -f ${md_tpr} -n ${lig_ndx} -o residues.ndx

    ${mindist_cmds}

    # SASA por residuo (surface = proteina completa; output = residuo
    # individual) — valores baixos indicam residuo enterrado/em contato
    ${sasa_cmds}

    echo "[OK] Distancias e SASA de ${entries.size()} residuo(s)-chave concluidos para ${meta.id}" >&2
    """
}
