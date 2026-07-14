// Substitui STABILITY_FILTER (usado em Milena-MD/MD-gromacs) por extracao de
// fase EXPLICITA, nao deteccao automatica de plateau. stability_filter.py
// assume o padrao "instavel no inicio -> estabiliza e fica estavel ate o
// fim" (procura o PRIMEIRO platô de RMSD e usa tudo dali pra frente) -- o
// oposto do que essa trajetoria mostra: 2I9T+daidzeina fica estavel de
// 0-60ns (RMSD backbone plateau ~0.45nm, contato Arg30/Glu279 batendo com o
// docking) e DESESTABILIZA de forma real e sequencial em ~65-78ns (RMSD
// backbone cruza 0.6nm em ~68ns -> Arg30 rompe em ~74ns -> Glu279 rompe por
// ultimo em ~77ns), com o ligante migrando pra outro sitio na superficie do
// receptor (mindist continua <1nm, SASA do ligante nao muda). Rodar
// stability_filter.py aqui detectaria o platô ERRADO (20-100ns, incluindo a
// fase ja desestabilizada) porque olha só pro primeiro trecho de SD baixo.
//
// Os limites de fase abaixo vieram da analise real da trajetoria completa
// (ver memoria do projeto, 2026-07-14) -- nao sao um algoritmo generico,
// sao os limites EXATOS desta corrida. Se re-rodar com parametros diferentes
// (outra semente, outro tempo de producao), os limites precisam ser
// reavaliados a partir dos gráficos antes de reusar este modulo.
process PHASE_SPLIT {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/analise_extra" }, mode: 'copy'

    input:
    tuple val(meta), path(md_tpr), path(md_fit_xtc), path(lig_ndx)

    output:
    tuple val(meta), path(md_tpr), path("bound.xtc"),     path(lig_ndx), emit: bound
    tuple val(meta), path(md_tpr), path("relocated.xtc"), path(lig_ndx), emit: relocated
    tuple val(meta), path("phase_split_report.txt"),                     emit: report

    script:
    def bound_end    = params.phase_bound_end_ns    ?: 60
    def reloc_start  = params.phase_reloc_start_ns  ?: 65
    """
    echo "=== PHASE_SPLIT: ${meta.id} ===" >&2

    {
        echo "=== Relatorio de Divisao de Fase ==="
        echo "Sistema           : ${meta.id}"
        echo "Fase 'bound'      : 0 - ${bound_end} ns (pose de docking mantida, Arg30/Glu279 concordam com docking)"
        echo "Fase 'relocated'  : ${reloc_start} - fim (ligante migrou pra outro sitio da superficie)"
        echo "Janela descartada : ${bound_end} - ${reloc_start} ns (transicao, evita contaminar as duas fases)"
        echo "Metodo            : limites explicitos da inspecao da trajetoria real, NAO deteccao automatica"
        echo ""
    } > phase_split_report.txt

    echo "System" | ${params.gmx_cmd} trjconv \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -o bound.xtc \\
        -b 0 -e ${bound_end} -tu ns \\
        2>&1 | tee -a phase_split_report.txt

    echo "System" | ${params.gmx_cmd} trjconv \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -o relocated.xtc \\
        -b ${reloc_start} -tu ns \\
        2>&1 | tee -a phase_split_report.txt

    N_BOUND=\$(${params.gmx_cmd} check -f bound.xtc 2>&1 | grep -E "^Last frame" | awk '{print \$3}' || echo "?")
    N_RELOC=\$(${params.gmx_cmd} check -f relocated.xtc 2>&1 | grep -E "^Last frame" | awk '{print \$3}' || echo "?")
    echo "" >> phase_split_report.txt
    echo "Frames bound.xtc     : \${N_BOUND}" | tee -a phase_split_report.txt
    echo "Frames relocated.xtc : \${N_RELOC}" | tee -a phase_split_report.txt
    """
}
