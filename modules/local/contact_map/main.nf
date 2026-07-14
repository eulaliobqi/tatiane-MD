// Portado de Milena-MD/modules/local/contact_map. Roda uma vez POR FASE
// (bound/relocated, ver PHASE_SPLIT) — o objetivo real de portar isso aqui
// e' descobrir para ONDE a daidzeina migrou apos a transicao de ~65-78ns:
// a fase "bound" deve mostrar Arg30/Glu279 como picos de contato (confirma
// o docking), a fase "relocated" revela o(s) residuo(s) do NOVO sitio.
process CONTACT_MAP {
    tag "${meta.id}:${meta.phase}"
    label 'process_medium'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/analise_extra/${meta.phase}" }, mode: 'copy'

    input:
    tuple val(meta), path(complexo_pdb), path(md_tpr), path(phase_xtc)

    output:
    tuple val(meta), path("contact_map.csv"), path("contact_map.png"),
                     path("interface_residues.csv"), emit: results

    script:
    """
    echo "=== CONTACT_MAP: ${meta.id} (${meta.phase}) ===" >&2

    contact_map.py \\
        --complexo-pdb ${complexo_pdb} \\
        --tpr ${md_tpr} \\
        --xtc ${phase_xtc} \\
        --out-dir .

    echo "[OK] Mapa de contatos concluido para ${meta.id} (${meta.phase})" >&2
    """
}
