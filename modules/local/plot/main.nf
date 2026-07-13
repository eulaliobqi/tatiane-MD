process PLOT {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/analise" }, mode: 'copy'

    input:
    tuple val(meta), path(xvgs), path(residue_xvgs)

    output:
    tuple val(meta), path("painel_resumo.png"), emit: png

    script:
    """
    mkdir -p analise_dir
    cp *.xvg analise_dir/ 2>/dev/null || true

    python3 ${projectDir}/bin/plot_results.py \\
        --analise-dir analise_dir \\
        --titulo "2I9T (NF-kB) + Daidzeina — resumo da dinamica molecular" \\
        --window-ns ${params.window_ns} \\
        --out painel_resumo.png
    """
}
