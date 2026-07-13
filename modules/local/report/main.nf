// NAO depende do canal de saida do MMGBSA (errorStrategy 'ignore' pode nunca
// emitir, e um join contra esse canal ja travou PLOT inteiro em producao no
// projeto irmao Milena-MD — ver comentario em main.nf). Le o diretorio do
// MM-GBSA como caminho de disco (publishDir), tratando ausencia com
// graciosidade (gerar_artigo_md.py ja retorna N/D se o arquivo nao existir).
// Reexecutar manualmente depois que o MMGBSA_ROBUST terminar atualiza a secao.
process REPORT {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}" }, mode: 'copy'
    publishDir "${projectDir}/docs", mode: 'copy', pattern: 'artigo_md.md'

    input:
    tuple val(meta), path(xvgs), path(residue_xvgs)

    output:
    tuple val(meta), path("artigo_md.md"), emit: report

    script:
    """
    mkdir -p analise_dir
    cp *.xvg analise_dir/ 2>/dev/null || true

    python3 ${projectDir}/bin/gerar_artigo_md.py \\
        --analise-dir analise_dir \\
        --mmgbsa-dir  ${params.outdir}/${meta.id}/mmgbsa \\
        --time-ns     ${params.time_ns} \\
        --out         artigo_md.md
    """
}
