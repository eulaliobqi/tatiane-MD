// Portado de Milena-MD/modules/local/prolif_fingerprint. ProLIF foi
// desenhado primariamente para ligante molecula-pequena (o caso de uso mais
// comum em drug discovery) — encaixa melhor aqui do que no par
// proteina-proteina do Milena-MD. Roda uma vez POR FASE: compara o padrao
// de interacao qualitativo (H-bond/hidrofobica/pi-stacking) da pose de
// docking original contra o do sitio novo pos-transicao.
process PROLIF_FINGERPRINT {
    tag "${meta.id}:${meta.phase}"
    label 'process_medium'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/analise_extra/${meta.phase}" }, mode: 'copy'

    input:
    tuple val(meta), path(complexo_pdb), path(md_tpr), path(phase_xtc)

    output:
    tuple val(meta), path("prolif_fingerprint.csv"), path("prolif_log.txt"),
                     path("prolif_heatmap.png", optional: true), emit: results

    script:
    """
    echo "=== PROLIF_FINGERPRINT: ${meta.id} (${meta.phase}) ===" >&2

    prolif_fingerprint.py \\
        --complexo-pdb ${complexo_pdb} \\
        --tpr ${md_tpr} \\
        --xtc ${phase_xtc} \\
        --out-dir .

    echo "[OK] ProLIF (ou fallback) concluido para ${meta.id} (${meta.phase})" >&2
    """
}
