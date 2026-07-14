// Portado de Milena-MD/modules/local/clustering (GROMOS via gmx cluster
// sobre o RMSD do ligante). Adaptado: la o output "for_mmgbsa" alimentava
// MMGBSA_ROBUST; aqui alimenta FE_RERUN (Interaction Entropy), ja que este
// projeto nao roda um segundo MM-GBSA sobre subtrajetoria reduzida — so o
// MM-GBSA de trajetoria unica em modules/local/mmgbsa (que ja falhou, ver
// memoria do projeto). Roda uma vez POR FASE (bound/relocated, tag em
// meta.phase, ver PHASE_SPLIT) -- clusters da pose original vs. clusters do
// sitio novo sao o ponto de comparacao, nao um clustering unico da
// trajetoria inteira (que misturaria as duas populacoes).
process CLUSTERING {
    tag "${meta.id}:${meta.phase}"
    label 'process_medium'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/analise_extra/${meta.phase}/clustering" }, mode: 'copy'

    input:
    tuple val(meta), path(md_tpr), path(phase_xtc), path(lig_ndx)

    output:
    tuple val(meta), path(md_tpr), path("fe_input.xtc"), path(lig_ndx), emit: for_fe
    tuple val(meta), path("clusterid.xvg"), path("cluster_centers.pdb"), emit: clusters
    tuple val(meta), path("clustering_report.txt"), emit: report

    script:
    def cutoff = params.cluster_cutoff       ?: 0.2
    def n_cl   = params.fe_n_clusters        ?: 3
    def fpc    = params.fe_frames_per_cluster ?: 50
    def target_frames = n_cl * fpc
    """
    echo "=== CLUSTERING: ${meta.id} (${meta.phase}) ===" >&2

    {
        echo "=== Relatorio de Clustering ==="
        echo "Sistema  : ${meta.id}"
        echo "Fase     : ${meta.phase}"
        echo "Metodo   : GROMOS"
        echo "Grupo    : Ligante (pose de ligacao)"
        echo "Cutoff   : ${cutoff} nm"
        echo ""
    } > clustering_report.txt

    printf 'Ligante\\nLigante\\n' | ${params.gmx_cmd} cluster \\
        -s ${md_tpr} \\
        -f ${phase_xtc} \\
        -n ${lig_ndx} \\
        -method gromos \\
        -cutoff ${cutoff} \\
        -o clusters.xpm \\
        -g cluster.log \\
        -clid clusterid.xvg \\
        -cl cluster_centers.pdb \\
        -tu ns \\
        2>&1 | tee -a clustering_report.txt

    echo "" >> clustering_report.txt
    echo "--- Clusters detectados ---" >> clustering_report.txt
    grep -E "^(cl\\.|Total|Middle|Found)" cluster.log 2>/dev/null \\
        | head -50 >> clustering_report.txt || true

    CHECK_OUT=\$(${params.gmx_cmd} check -f ${phase_xtc} 2>&1)
    N_FRAMES=\$(echo "\${CHECK_OUT}" | grep -E "^Last frame" | awk '{print \$3}')
    if [ -z "\${N_FRAMES}" ]; then
        echo "WARN: nao foi possivel extrair 'Last frame' de gmx check -- usando fallback 1000" >&2
        N_FRAMES=1000
    fi

    SKIP=\$(python3 -c "print(max(1, int(\${N_FRAMES} // ${target_frames})))")

    echo "" >> clustering_report.txt
    echo "--- Subsampling para FE_RERUN ---" >> clustering_report.txt
    echo "Frames na fase   : \${N_FRAMES}" >> clustering_report.txt
    echo "Target frames    : ${target_frames}  (${n_cl} clusters x ${fpc} frames)" >> clustering_report.txt
    echo "Stride usado     : \${SKIP}" >> clustering_report.txt

    echo "System" | ${params.gmx_cmd} trjconv \\
        -s ${md_tpr} \\
        -f ${phase_xtc} \\
        -o fe_input.xtc \\
        -skip \${SKIP} \\
        2>&1 | tail -5

    N_OUT=\$(${params.gmx_cmd} check -f fe_input.xtc 2>&1 \\
        | grep -E "^Last frame" | awk '{print \$3}' || echo "?")
    echo "Frames para FE_RERUN : \${N_OUT}" | tee -a clustering_report.txt
    echo "[OK] fe_input.xtc: \${N_OUT} frames" >&2
    """
}
