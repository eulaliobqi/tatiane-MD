process PLOT {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/analise" }, mode: 'copy'

    input:
    tuple val(meta), path(xvgs), path(residue_xvgs)

    output:
    tuple val(meta), path("painel_resumo.png"), emit: png
    tuple val(meta), path("figuras/*.png"), emit: png_individual

    script:
    """
    mkdir -p analise_dir
    cp *.xvg analise_dir/ 2>/dev/null || true

    python3 ${projectDir}/bin/plot_results.py \\
        --analise-dir analise_dir \\
        --titulo "${meta.target_name ?: meta.id} + ${meta.ligand_name ?: 'ligante'} — resumo da dinamica molecular" \\
        --window-ns ${params.window_ns} \\
        --out painel_resumo.png \\
        --figures-dir figuras
    """
}
