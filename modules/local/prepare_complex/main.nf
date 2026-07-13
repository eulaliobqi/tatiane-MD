process PREPARE_COMPLEX {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/complex" }, mode: 'copy'

    input:
    tuple val(meta), path(receptor_ph), path(ligand_ini_pdb)

    output:
    tuple val(meta), path("complex.pdb"), emit: complexo

    script:
    """
    echo "=== PREPARE_COMPLEX: ${meta.id} (receptor cadeia A + UNL cadeia B) ===" >&2

    python3 ${projectDir}/bin/prepare_complex.py \\
        --receptor ${receptor_ph} \\
        --ligand   ${ligand_ini_pdb} \\
        --out      complex.pdb

    if [ ! -s complex.pdb ]; then
        echo "ERRO: complex.pdb nao foi gerado"; exit 1
    fi
    echo "[OK] PREPARE_COMPLEX concluido para ${meta.id}" >&2
    """
}
