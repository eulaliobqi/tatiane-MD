process PREPARE_PH {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/prep_ph" }, mode: 'copy'

    input:
    tuple val(meta), path(receptor)

    output:
    tuple val(meta), path("receptor_ph.pdb"), emit: protonated
    tuple val(meta), path("*.propka"), optional: true, emit: propka

    script:
    def ph = params.pH
    """
    echo "=== PREPARE_PH: ${meta.id} (pH ${ph}, CHARMM) ===" >&2

    # Alvo humano (NF-kB) -> pH fisiologico, NAO o pH 8.2 (midgut de inseto)
    # usado nos demais pipelines deste laboratorio. --ff/--ffout CHARMM produz
    # nomenclatura de titulacao ja compativel com charmm36.ff (HSD/HSE/HSP,
    # ASPP/GLUP), sem precisar de tabela de renomeacao tipo AMBER (HID/HIE/HIP).
    pdb2pqr --ff CHARMM --ffout CHARMM \\
        --titration-state-method propka --with-ph ${ph} \\
        --pdb-output receptor_raw.pdb \\
        --nodebump \\
        ${receptor} receptor.pqr

    python3 ${projectDir}/bin/pdb2pqr_process_charmm.py receptor_raw.pdb receptor_ph.pdb

    if [ ! -s receptor_ph.pdb ]; then
        echo "ERRO: receptor_ph.pdb nao foi gerado"; exit 1
    fi
    echo "[OK] PREPARE_PH concluido para ${meta.id}" >&2
    """
}
