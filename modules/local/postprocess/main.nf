process POSTPROCESS {
    tag "${meta.id}"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/postprocess" }, mode: 'copy'

    input:
    tuple val(meta), path(md_tpr), path(md_xtc)

    output:
    tuple val(meta), path(md_tpr), path("md_fit.xtc"), emit: fit

    script:
    """
    echo "Protein System" | ${params.gmx_cmd} trjconv -s ${md_tpr} -f ${md_xtc} \\
        -pbc mol -center -o md_center.xtc

    echo "Backbone System" | ${params.gmx_cmd} trjconv -s ${md_tpr} -f md_center.xtc \\
        -fit rot+trans -o md_fit.xtc
    """
}
