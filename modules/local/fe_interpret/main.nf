// Portado de Milena-MD/modules/local/fe_interpret sem alteracao de logica —
// so o titulo do grafico ganha a fase (bound/relocated) para diferenciar as
// duas estimativas no relatorio final.
process FE_INTERPRET {
    tag "${meta.id}:${meta.phase}"
    label 'process_low'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/analise_extra/${meta.phase}/fe_estimate" }, mode: 'copy'

    input:
    tuple val(meta), path(interaction_xvg)

    output:
    tuple val(meta), path("free_energy_estimate.txt"),
                     path("interaction_energy.png", optional: true), emit: results

    script:
    def titulo = "${meta.id} (${meta.phase}) — Interaction Energy (Coul-SR + LJ-SR)"
    """
    echo "=== FE_INTERPRET: ${meta.id} (${meta.phase}) ===" >&2

    interaction_entropy.py \\
        --xvg ${interaction_xvg} \\
        --temperature ${params.temperature} \\
        --titulo "${titulo}" \\
        --out-dir .

    echo "[OK] Estimativa de energia livre concluida para ${meta.id} (${meta.phase})" >&2
    """
}
