process LIGAND_TOPOLOGY {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/ligand_topo" }, mode: 'copy'

    input:
    tuple val(meta), path(ligand_mol2), path(ligand_str)

    output:
    tuple val(meta), path("unl.itp"), path("unl.prm"), path("posre_UNL.itp"),
                     path("unl_ini.pdb"), path("unl_ini.gro"), emit: topology

    script:
    """
    echo "=== LIGAND_TOPOLOGY: ${meta.id} (CGenFF -> GROMACS/CHARMM36) ===" >&2

    python3 ${projectDir}/bin/cgenff_charmm2gmx.py UNL \\
        ${ligand_mol2} ${ligand_str} \\
        ${projectDir}/ff/charmm36-mar2019.ff \\
        2>&1 | tee cgenff_charmm2gmx.log

    if [ ! -s unl.itp ]; then
        echo "ERRO: unl.itp nao foi gerado — ver cgenff_charmm2gmx.log"; exit 1
    fi

    # Converte a pose (todos os atomos, com H) para .gro, e gera a restricao
    # posicional do ligante (usada em NVT/NPT sob -DPOSRES_UNL)
    ${params.gmx_cmd} editconf -f unl_ini.pdb -o unl_ini.gro
    ${params.gmx_cmd} genrestr -f unl_ini.pdb -o posre_UNL.itp -fc 1000 1000 1000 <<< "0"

    cat >> unl.itp << 'EOF'

#ifdef POSRES_UNL
#include "posre_UNL.itp"
#endif
EOF
    echo "[OK] LIGAND_TOPOLOGY concluido para ${meta.id}" >&2
    """
}
